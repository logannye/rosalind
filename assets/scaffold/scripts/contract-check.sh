#!/usr/bin/env bash
set -euo pipefail
ROSALIND_BIN="${ROSALIND_BIN:-rosalind}"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

"$ROSALIND_BIN" demo --output-dir "$WORK/demo" --json >/dev/null
cargo build --release
BIN="target/release/__PACKAGE_NAME__"

for n in 1 2; do
  "$BIN" run --index "$WORK/demo/ref.idx" --alignments "$WORK/demo/sorted.bam" \
    --memory-budget-mb 128 --enforce --output "$WORK/out$n.tsv"
done
cmp "$WORK/out1.tsv" "$WORK/out2.tsv"
"$ROSALIND_BIN" verify --manifest "$WORK/out1.tsv.manifest.json" --json >/dev/null
"$ROSALIND_BIN" receipt sanitize --manifest "$WORK/out1.tsv.manifest.json" \
  --output "$WORK/sanitized.json"
"$ROSALIND_BIN" reproduce --manifest "$WORK/out1.tsv.manifest.json" \
  --inputs "$WORK/demo" --binary "$BIN" --dry-run --json >/dev/null
"$ROSALIND_BIN" reproduce --manifest "$WORK/out1.tsv.manifest.json" \
  --inputs "$WORK/demo" --binary "$BIN" --json >/dev/null

"$BIN" run --index "$WORK/demo/ref.idx" --alignments "$WORK/demo/sorted.bam" \
  --scale 2 --output "$WORK/scaled.tsv"
set +e
diff_report=$("$ROSALIND_BIN" diff "$WORK/out1.tsv.manifest.json" \
  "$WORK/scaled.tsv.manifest.json")
diff_code=$?
set -e
test "$diff_code" -eq 1
printf '%s' "$diff_report" | grep -q 'analyzer.scale'

set +e
"$BIN" run --index "$WORK/demo/ref.idx" --alignments "$WORK/demo/sorted.bam" \
  --memory-budget-mb 1 --enforce --output "$WORK/refused.tsv"
code=$?
set -e
test "$code" -eq 3
test ! -e "$WORK/refused.tsv"

set +e
ROSALIND_FORCE_PEAK_RSS_BYTES=8589934592 "$BIN" run \
  --index "$WORK/demo/ref.idx" --alignments "$WORK/demo/sorted.bam" \
  --memory-budget-mb 128 --enforce --output "$WORK/breached.tsv"
code=$?
set -e
test "$code" -eq 4
test ! -e "$WORK/breached.tsv"
test -e "$WORK/breached.tsv.partial"

set +e
"$BIN" run --index "$WORK/demo/ref.idx" --alignments "$WORK/demo/sorted.bam" \
  --output "$WORK/out1.tsv"
code=$?
set -e
test "$code" -eq 2
