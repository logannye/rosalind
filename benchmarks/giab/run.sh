#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
DATA="${GIAB_DATA_DIR:-$HERE/data}"
BUDGET_MB="${GIAB_BUDGET_MB:-4096}"
UPDATE_BASELINE=false
if [ "${1:-}" = "--update-baseline" ]; then UPDATE_BASELINE=true; fi

for tool in python3; do
  command -v "$tool" >/dev/null || { echo "missing required tool: $tool" >&2; exit 2; }
done
test -f "$DATA/data-manifest.json" || {
  echo "prepare the opt-in data first: benchmarks/giab/prepare.sh '$DATA'" >&2
  exit 2
}

cd "$ROOT"
BIN="${ROSALIND_BIN:-target/release/rosalind}"
if [ ! -x "$BIN" ]; then cargo build --release --bin rosalind; fi
mkdir -p "$DATA/results"
REF="$DATA/prepared/GRCh38.chr20.fa"
BAM="$DATA/prepared/HG002.chr20.bam"
TRUTH="$DATA/prepared/HG002.v5.0q.chr20.vcf"
BED="$DATA/prepared/HG002.v5.0q.chr20.bed"
IDX="$DATA/results/GRCh38.chr20.idx"
CALLS="$DATA/results/HG002.rosalind.chr20.vcf"

"$BIN" index --reference "$REF" --output "$IDX"
"$BIN" variants --index "$IDX" --alignments "$BAM" \
  --memory-budget-mb "$BUDGET_MB" --enforce -o "$CALLS"
"$BIN" verify --manifest "$CALLS.manifest.json" --json > "$DATA/results/verify.json"
"$BIN" eval-germline --reference "$REF" --calls "$CALLS" --truth "$TRUTH" \
  --regions "$BED" --calls-filter all --json > "$DATA/results/metrics-all.json"
"$BIN" eval-germline --reference "$REF" --calls "$CALLS" --truth "$TRUTH" \
  --regions "$BED" --calls-filter pass --json > "$DATA/results/metrics-pass.json"

DATA_ROOT="$DATA" BUDGET_MB="$BUDGET_MB" python3 - <<'PY'
import json, os
from pathlib import Path

root = Path(os.environ["DATA_ROOT"])
results = root / "results"
receipt = json.loads((results / "HG002.rosalind.chr20.vcf.manifest.json").read_text())
report = {
    "schema": 1,
    "benchmark": "HG002 GIAB v5.0q small variants",
    "assembly": "GRCh38 GCA_000001405.15",
    "region": "chr20",
    "calls_filter_all": json.loads((results / "metrics-all.json").read_text()),
    "calls_filter_pass": json.loads((results / "metrics-pass.json").read_text()),
    "memory": {
        "budget_mb": int(os.environ["BUDGET_MB"]),
        "predicted_peak_rss_bytes": int(receipt["measurements"]["predicted_peak_rss_bytes"]),
        "peak_rss_bytes": int(receipt["measurements"]["peak_rss_bytes"]),
        "max_working_set_bytes": int(receipt["measurements"]["max_working_set_bytes"]),
    },
    "receipt_claim": receipt["params"]["manifest_blake3"],
    "command": receipt["params"].get("command"),
    "command_argv": json.loads(receipt["params"]["command_argv"]),
}
(results / "latest.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
PY

BASELINE="$HERE/baseline.json"
if grep -q '"status": "not-yet-established"' "$BASELINE"; then
  cp "$DATA/results/latest.json" "$BASELINE"
  echo "Stored the first honest GIAB baseline at $BASELINE"
elif ! python3 - "$BASELINE" "$DATA/results/latest.json" <<'PY'
import json, sys
a, b = (json.load(open(path)) for path in sys.argv[1:])
keys = ("calls_filter_all", "calls_filter_pass")
raise SystemExit(0 if all(a[k] == b[k] for k in keys) else 1)
PY
then
  if [ "$UPDATE_BASELINE" != true ]; then
    echo "GIAB metrics changed; rerun with --update-baseline after documenting the reason" >&2
    exit 1
  fi
  grep -q 'GIAB baseline update' CHANGELOG.md || {
    echo "baseline update requires a CHANGELOG entry containing 'GIAB baseline update'" >&2
    exit 1
  }
  cp "$DATA/results/latest.json" "$BASELINE"
  echo "Updated intentional GIAB baseline"
fi

echo "GIAB v5.0q chr20 benchmark reproduced exactly; latest report: $DATA/results/latest.json"
