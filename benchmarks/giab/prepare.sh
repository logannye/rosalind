#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
DEST="${1:-$HERE/data}"
mkdir -p "$DEST/downloads" "$DEST/prepared"

for tool in curl samtools gzip awk shasum python3; do
  command -v "$tool" >/dev/null || { echo "missing required tool: $tool" >&2; exit 2; }
done

verify_sha256() {
  local expected="$1" file="$2" observed
  observed="$(shasum -a 256 "$file" | awk '{print $1}')"
  if [ "$observed" != "$expected" ]; then
    echo "SHA-256 mismatch for $file: expected $expected, got $observed" >&2
    exit 2
  fi
}

echo "Downloading and verifying pinned HG002 v5.0q resources..."
while IFS=$'\t' read -r id filename sha256 url; do
  case "$id" in ''|'#'*) continue ;; esac
  target="$DEST/downloads/$filename"
  if [ ! -f "$target" ]; then
    curl --fail --location --retry 3 --continue-at - --output "$target.part" "$url"
    mv "$target.part" "$target"
  fi
  verify_sha256 "$sha256" "$target"
  echo "  OK $id  $sha256"
done < "$HERE/resources.tsv"

downloads="$DEST/downloads"
prepared="$DEST/prepared"
reference="$downloads/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz"
reads="$downloads/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam"

samtools quickcheck "$reads"
samtools view -H "$reads" | grep -q $'^@SQ.*SN:chr20' || {
  echo "aligned reads do not declare chr20" >&2
  exit 2
}

# Re-extract through the index even though the pinned DeepVariant case-study BAM
# is already chr20-only. --no-PG avoids a tool-version-dependent header line.
samtools view --no-PG -b "$reads" chr20 -o "$prepared/HG002.chr20.bam"
samtools index "$prepared/HG002.chr20.bam"
samtools faidx "$reference" chr20 > "$prepared/GRCh38.chr20.fa"
gzip -dc "$downloads/HG002_GRCh38_v5.0q_smvar.vcf.gz" \
  | awk 'BEGIN{FS=OFS="\t"} /^#/ || $1=="chr20"' > "$prepared/HG002.v5.0q.chr20.vcf"
awk 'BEGIN{FS=OFS="\t"} $1=="chr20"' \
  "$downloads/HG002_GRCh38_v5.0q_smvar.benchmark.bed" \
  > "$prepared/HG002.v5.0q.chr20.bed"

for name in \
  GRCh38_AllTandemRepeatsandHomopolymers_slop5 \
  GRCh38_lowmappabilityall \
  GRCh38_segdups \
  GRCh38_alldifficultregions
do
  gzip -dc "$downloads/$name.bed.gz" \
    | awk 'BEGIN{FS=OFS="\t"} $1=="chr20"' \
    | gzip -n > "$prepared/$name.chr20.bed.gz"
done
cat > "$prepared/stratifications.tsv" <<'EOF'
low_complexity	/data/prepared/GRCh38_AllTandemRepeatsandHomopolymers_slop5.chr20.bed.gz
low_mappability	/data/prepared/GRCh38_lowmappabilityall.chr20.bed.gz
segmental_duplications	/data/prepared/GRCh38_segdups.chr20.bed.gz
all_difficult	/data/prepared/GRCh38_alldifficultregions.chr20.bed.gz
EOF

SAMTOOLS_VERSION="$(samtools --version | head -1)" \
DATA_ROOT="$DEST" RESOURCE_MANIFEST="$HERE/resources.tsv" python3 - <<'PY'
import hashlib, json, os
from pathlib import Path

root = Path(os.environ["DATA_ROOT"])
prepared = root / "prepared"
artifacts = {}
for path in sorted(prepared.iterdir()):
    if path.is_file():
        digest = hashlib.sha256()
        with path.open("rb") as stream:
            for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
                digest.update(chunk)
        artifacts[path.name] = {
            "bytes": path.stat().st_size,
            "sha256": digest.hexdigest(),
        }
sources = []
for line in Path(os.environ["RESOURCE_MANIFEST"]).read_text().splitlines():
    if not line or line.startswith("#"):
        continue
    identity, filename, sha256, url = line.split("\t")
    sources.append({"id": identity, "filename": filename, "sha256": sha256, "url": url})
manifest = {
    "schema": 1,
    "sample": "HG002/NA24385",
    "benchmark": "GIAB v5.0q small variants",
    "assembly": "GRCh38 GCA_000001405.15",
    "region": "chr20",
    "samtools": os.environ["SAMTOOLS_VERSION"],
    "sources": sources,
    "prepared_artifacts": artifacts,
}
(root / "data-manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
PY

echo "Prepared chr20 data and wrote $DEST/data-manifest.json"
