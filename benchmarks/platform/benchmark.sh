#!/usr/bin/env bash
set -euo pipefail

REFERENCE="${REFERENCE:-/data/reference.fa}"
SOURCE_BAM="${BAM:-/data/input.bam}"
RESULTS="${RESULTS:-/results}"
BUDGET_MB="${BUDGET_MB:-4096}"
SHARDS="${SHARDS:-3}"
mkdir -p "$RESULTS"/{raw,outputs,inputs,shards}
cp "$REFERENCE" "$RESULTS/inputs/reference.fa"
REFERENCE="$RESULTS/inputs/reference.fa"
rosalind reference build --fasta "$REFERENCE" --output "$RESULTS/inputs/reference.rref"
# This proof point deliberately precedes every benchmark-side BAM copy, sort,
# index, or read. Planning has no alignment operand.
rosalind plan --reference-pack "$RESULTS/inputs/reference.rref" --budget-mb "$BUDGET_MB" --json > "$RESULTS/raw/rosalind-plan.json"

cp "$SOURCE_BAM" "$RESULTS/inputs/input.bam"
BAM="$RESULTS/inputs/input.bam"

samtools quickcheck "$BAM"
samtools sort -o "$RESULTS/inputs/sorted.bam" "$BAM"
BAM="$RESULTS/inputs/sorted.bam"

samtools faidx "$REFERENCE"
samtools index "$BAM"

measure() {
  local label="$1"; shift
  python3 - "$RESULTS/raw/$label.argv.json" "$@" <<'PY'
import json, sys
from pathlib import Path
Path(sys.argv[1]).write_text(json.dumps(sys.argv[2:]) + "\n")
PY
  set +e
  /usr/bin/time -v -o "$RESULTS/raw/$label.time.txt" "$@"
  local code=$?
  set -e
  printf '%s\n' "$code" > "$RESULTS/raw/$label.exit"
  test "$code" -eq 0
}

for run in 1 2 3; do
  measure "rosalind-$run" rosalind features --reference-pack "$RESULTS/inputs/reference.rref" \
    --alignments "$BAM" --max-depth 1000 --output "$RESULTS/outputs/rosalind-$run.tsv"
  measure "pysam-$run" python3 /opt/rosalind-platform/pysam_features.py \
    --reference "$REFERENCE" --bam "$BAM" --max-depth 1000 --output "$RESULTS/outputs/pysam-$run.tsv"
  measure "bcftools-$run" python3 /opt/rosalind-platform/bcftools_features.py \
    --reference "$REFERENCE" --bam "$BAM" --max-depth 1000 --output "$RESULTS/outputs/bcftools-$run.tsv"
done
rosalind verify --manifest "$RESULTS/outputs/rosalind-1.tsv.manifest.json" --json > "$RESULTS/raw/rosalind-verify.json"

rosalind features --reference-pack "$RESULTS/inputs/reference.rref" --alignments "$BAM" \
  --format arrow-ipc --output "$RESULTS/outputs/unsharded.arrow"
rosalind features --reference-pack "$RESULTS/inputs/reference.rref" --alignments "$BAM" \
  --format arrow-ipc --output "$RESULTS/outputs/unsharded-repeat.arrow"
for ((shard=0; shard<SHARDS; shard++)); do
  rosalind features --reference-pack "$RESULTS/inputs/reference.rref" --alignments "$BAM" \
    --format arrow-ipc --shard-count "$SHARDS" --shard-index "$shard" \
    --output "$RESULTS/shards/shard-$shard.arrow"
done
merge_args=()
for ((shard=0; shard<SHARDS; shard++)); do merge_args+=(--manifest "$RESULTS/shards/shard-$shard.arrow.manifest.json"); done
rosalind merge "${merge_args[@]}" --inputs "$RESULTS/shards" --output "$RESULTS/outputs/merged.arrow"

RESULTS="$RESULTS" python3 - <<'PY'
import hashlib, json, os
from pathlib import Path
import pyarrow.ipc as ipc
root = Path(os.environ["RESULTS"])
def digest(path):
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()
unsharded = root / "outputs/unsharded.arrow"
repeated = root / "outputs/unsharded-repeat.arrow"
merged = root / "outputs/merged.arrow"
with unsharded.open("rb") as stream:
    reader = ipc.open_stream(stream)
    batches = [batch.num_rows for batch in reader]
report = {
    "unsharded_sha256": digest(unsharded), "repeat_sha256": digest(repeated),
    "merged_sha256": digest(merged),
    "repeat_byte_identical": digest(unsharded) == digest(repeated),
    "merge_byte_identical": digest(unsharded) == digest(merged),
    "record_batch_rows": batches,
}
(root / "raw/arrow-checks.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
if not report["repeat_byte_identical"]: raise SystemExit("repeated Arrow output differs")
if not report["merge_byte_identical"]: raise SystemExit("sharded Arrow merge differs from unsharded output")
PY

# Separate declaration-refusal evidence from matched-budget container trials.
set +e
rosalind features --reference-pack "$RESULTS/inputs/reference.rref" --alignments "$BAM" \
  --memory-budget-mb 1 --enforce --output "$RESULTS/outputs/refused.tsv" \
  > "$RESULTS/raw/refusal.stdout" 2> "$RESULTS/raw/refusal.stderr"
refusal_code=$?
set -e
test "$refusal_code" -eq 3
test ! -e "$RESULTS/outputs/refused.tsv"
printf '{"declared_mb":1,"exit_code":%s,"output_created":false}\n' "$refusal_code" > "$RESULTS/raw/refusal.json"

dpkg-query -W -f='${Package}\t${Version}\n' | sort > "$RESULTS/environment-packages.tsv"
python3 -m pip freeze --all > "$RESULTS/environment-python.txt"
RESULTS="$RESULTS" python3 /opt/rosalind-platform/report.py
