#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
DATA="${GIAB_DATA_DIR:-$HERE/../data}"
IMAGE="${GIAB_HAPPY_IMAGE:-}"
if [[ ! "$IMAGE" =~ ^[^[:space:]@]+@sha256:[0-9a-f]{64}$ ]]; then
  echo "GIAB_HAPPY_IMAGE must be an immutable NAME@sha256:DIGEST reference" >&2
  exit 2
fi
command -v docker >/dev/null || { echo "docker is required" >&2; exit 2; }

prepared="$DATA/prepared"
results="$DATA/results/happy"
rm -rf "$results"
mkdir -p "$results"
for file in GRCh38.chr20.fa HG002.v5.0q.chr20.vcf HG002.v5.0q.chr20.bed stratifications.tsv; do
  test -f "$prepared/$file" || { echo "missing prepared input: $prepared/$file" >&2; exit 2; }
done
test -f "$DATA/results/HG002.rosalind.chr20.vcf" || { echo "run the Rosalind GIAB workflow first" >&2; exit 2; }

awk 'BEGIN{FS=OFS="\t"} /^#/ || $7=="PASS"' \
  "$DATA/results/HG002.rosalind.chr20.vcf" > "$results/calls-pass.vcf"

run_happy() {
  local label="$1" calls="$2"
  local prefix="/data/results/happy/$label"
  docker run --rm --network none --platform linux/amd64 \
    -v "$prepared:/data/prepared:ro" \
    -v "$DATA/results:/data/results" \
    "$IMAGE" \
    /data/prepared/HG002.v5.0q.chr20.vcf "$calls" \
    -f /data/prepared/HG002.v5.0q.chr20.bed \
    -r /data/prepared/GRCh38.chr20.fa \
    --engine=vcfeval --engine-vcfeval-path=/opt/rtg-tools/rtg \
    --stratification /data/prepared/stratifications.tsv \
    -o "$prefix"
}

run_happy all /data/results/HG002.rosalind.chr20.vcf
run_happy pass /data/results/happy/calls-pass.vcf

IMAGE="$IMAGE" DATA_ROOT="$DATA" python3 - <<'PY'
import csv, json, os, subprocess
from pathlib import Path

root = Path(os.environ["DATA_ROOT"])
results = root / "results" / "happy"
def rows(label):
    with (results / f"{label}.summary.csv").open(newline="") as stream:
        return list(csv.DictReader(stream))
def snv(rows_):
    return next((row for row in rows_ if row.get("Type") in {"SNP", "SNV"} and row.get("Filter") in {None, "ALL"}), {})
report = {
    "schema": 1,
    "evaluator": {"name": "Illumina hap.py", "version": "0.3.15", "engine": "vcfeval", "rtg_tools": "3.12.1"},
    "container": os.environ["IMAGE"],
    "network": "disabled",
    "all_emitted": {"snv": snv(rows("all")), "stratified": rows("all")},
    "pass_only": {"snv": snv(rows("pass")), "stratified": rows("pass")},
    "argv": {
        "all": ["/opt/hap.py/bin/hap.py", "/data/prepared/HG002.v5.0q.chr20.vcf", "/data/results/HG002.rosalind.chr20.vcf", "-f", "/data/prepared/HG002.v5.0q.chr20.bed", "-r", "/data/prepared/GRCh38.chr20.fa", "--engine=vcfeval", "--engine-vcfeval-path=/opt/rtg-tools/rtg", "--stratification", "/data/prepared/stratifications.tsv", "-o", "/data/results/happy/all"],
        "pass": ["/opt/hap.py/bin/hap.py", "/data/prepared/HG002.v5.0q.chr20.vcf", "/data/results/happy/calls-pass.vcf", "-f", "/data/prepared/HG002.v5.0q.chr20.bed", "-r", "/data/prepared/GRCh38.chr20.fa", "--engine=vcfeval", "--engine-vcfeval-path=/opt/rtg-tools/rtg", "--stratification", "/data/prepared/stratifications.tsv", "-o", "/data/results/happy/pass"],
        "container": ["docker", "run", "--rm", "--network", "none", "--platform", "linux/amd64", "-v", f"{root / 'prepared'}:/data/prepared:ro", "-v", f"{root / 'results'}:/data/results", os.environ["IMAGE"]]
    },
}
(results / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
PY

echo "hap.py v0.3.15 / vcfeval report: $results/report.json"
