#!/usr/bin/env bash
set -euo pipefail
IMAGE="${1:?usage: smoke.sh LOCAL_IMAGE RESULTS_DIRECTORY}"
RESULTS="${2:?results directory required}"
mkdir -p "$RESULTS"
RESULTS="$(cd "$RESULTS" && pwd)"
python3 - "$RESULTS" <<'PY'
from pathlib import Path
import sys
root = Path(sys.argv[1])
(root / "reference.fa").write_text(">chr20\n" + "A" * 1000 + "\n")
vcf = ("##fileformat=VCFv4.2\n##contig=<ID=chr20,length=1000>\n"
       '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
       "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
       "chr20\t100\t.\tA\tC\t50\tPASS\t.\tGT\t0/1\n")
for name in ("truth.vcf", "query.vcf"):
    (root / name).write_text(vcf)
(root / "confident.bed").write_text("chr20\t0\t1000\n")
PY
docker run --rm --network none --platform linux/amd64 "$IMAGE" --version > "$RESULTS/happy-version.txt"
docker run --rm --network none --platform linux/amd64 --entrypoint /opt/rtg-tools/rtg "$IMAGE" version > "$RESULTS/rtg-version.txt"
# hap.py preprocessing requires an existing FAI. Build it with the evaluator's
# installed pysam, keeping both fixture preparation and evaluation offline.
docker run --rm --network none --platform linux/amd64 -v "$RESULTS:/data" \
  --entrypoint python "$IMAGE" -c "import pysam; pysam.faidx('/data/reference.fa')"
python3 - "$RESULTS/reference.fa.fai" <<'PY_INDEX'
from pathlib import Path
import sys
assert Path(sys.argv[1]).read_text().splitlines() == ["chr20\t1000\t7\t1000\t1001"]
PY_INDEX
if docker run --rm --network none --platform linux/amd64 -v "$RESULTS:/data" "$IMAGE" \
  /data/truth.vcf /data/query.vcf -r /data/reference.fa -f /data/confident.bed \
  --engine=vcfeval --engine-vcfeval-path=/opt/rtg-tools/rtg --threads 1 -o /data/evaluation \
  > "$RESULTS/evaluator.stdout.txt" 2> "$RESULTS/evaluator.stderr.txt"; then
  :
else
  status=$?
  cat "$RESULTS/evaluator.stderr.txt" >&2
  exit "$status"
fi
python3 - "$RESULTS" <<'PY'
import csv, json, sys
from pathlib import Path
root = Path(sys.argv[1])
rows = list(csv.DictReader((root / "evaluation.summary.csv").open()))
snps = [row for row in rows if row["Type"] == "SNP" and row["Filter"] == "ALL"]
assert len(snps) == 1, rows
row = snps[0]
assert int(row["TRUTH.TP"]) == 1 and int(row["TRUTH.FN"]) == 0 and int(row["QUERY.FP"]) == 0, row
(root / "smoke.json").write_text(json.dumps({"status": "passed", "scope": "one synthetic identical SNV; not GIAB validation", "snp": row}, indent=2) + "\n")
PY
echo "offline hap.py/vcfeval synthetic execution passed"
