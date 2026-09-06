#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
DATA="${GIAB_DATA_DIR:-$HERE/data}"
BUDGET_MB="${GIAB_BUDGET_MB:-4096}"
command -v python3 >/dev/null || { echo "missing required tool: python3" >&2; exit 2; }
LOCK_IMAGE="$(python3 - "$HERE/happy/lock.json" <<'PY'
import json, sys
image = json.load(open(sys.argv[1])).get("generated_image", {})
repository, digest = image.get("repository"), image.get("digest")
print(f"{repository}@{digest}" if repository and digest else "")
PY
)"
GIAB_HAPPY_IMAGE="${GIAB_HAPPY_IMAGE:-$LOCK_IMAGE}"
if [[ ! "$GIAB_HAPPY_IMAGE" =~ ^[^[:space:]@]+@sha256:[0-9a-f]{64}$ ]]; then
  echo "no committed immutable hap.py image; run cargo xtask giab image plan" >&2
  exit 2
fi
test -f "$DATA/data-manifest.json" || {
  echo "prepare the opt-in data first: benchmarks/giab/prepare.sh '$DATA'" >&2
  exit 2
}
python3 "$HERE/preflight.py" "$DATA"

cd "$ROOT"
BIN="${ROSALIND_BIN:-target/release/rosalind}"
if [ ! -x "$BIN" ]; then cargo build --release --bin rosalind; fi
mkdir -p "$DATA/results"
REF="$DATA/prepared/GRCh38.chr20.fa"
BAM="$DATA/prepared/HG002.chr20.bam"
TRUTH="$DATA/prepared/HG002.v5.0q.chr20.vcf"
BED="$DATA/prepared/HG002.v5.0q.chr20.bed"
RREF="$DATA/results/GRCh38.chr20.rref"
CALLS="$DATA/results/HG002.rosalind.chr20.vcf"

"$BIN" reference build --fasta "$REF" --output "$RREF" --force
"$BIN" verify --manifest "$RREF.manifest.json" --json > "$DATA/results/reference-verify.json"
"$BIN" variants --reference-pack "$RREF" --alignments "$BAM" \
  --memory-budget-mb "$BUDGET_MB" --enforce -o "$CALLS" --force
"$BIN" verify --manifest "$CALLS.manifest.json" --json > "$DATA/results/verify.json"
"$BIN" eval-germline --reference "$REF" --calls "$CALLS" --truth "$TRUTH" \
  --regions "$BED" --calls-filter all --json > "$DATA/results/metrics-all.json"
"$BIN" eval-germline --reference "$REF" --calls "$CALLS" --truth "$TRUTH" \
  --regions "$BED" --calls-filter pass --json > "$DATA/results/metrics-pass.json"
GIAB_HAPPY_IMAGE="$GIAB_HAPPY_IMAGE" "$HERE/happy/run.sh"

DATA_ROOT="$DATA" BUDGET_MB="$BUDGET_MB" python3 - <<'PY'
import json, os
from pathlib import Path

root = Path(os.environ["DATA_ROOT"])
results = root / "results"
receipt = json.loads((results / "HG002.rosalind.chr20.vcf.manifest.json").read_text())
reference_receipt = json.loads((results / "GRCh38.chr20.rref.manifest.json").read_text())
report = {
    "schema": 1,
    "benchmark": "HG002 GIAB v5.0q small variants",
    "assembly": "GRCh38 GCA_000001405.15",
    "region": "chr20",
    "calls_filter_all": json.loads((results / "metrics-all.json").read_text()),
    "calls_filter_pass": json.loads((results / "metrics-pass.json").read_text()),
    "external_happy_vcfeval": json.loads((results / "happy" / "report.json").read_text()),
    "data_manifest": json.loads((root / "data-manifest.json").read_text()),
    "reference_pack": {
        "receipt_claim": reference_receipt["params"]["manifest_blake3"],
        "content_blake3": reference_receipt["outputs"][0]["blake3"],
        "verify": json.loads((results / "reference-verify.json").read_text()),
        "command_argv": json.loads(reference_receipt["params"]["command_argv"]),
        "peak_rss_bytes": int(reference_receipt["measurements"]["peak_rss_bytes"]),
    },
    "memory": {
        "budget_mb": int(os.environ["BUDGET_MB"]),
        "predicted_peak_rss_bytes": int(receipt["measurements"]["predicted_peak_rss_bytes"]),
        "peak_rss_bytes": int(receipt["measurements"]["peak_rss_bytes"]),
        "max_working_set_bytes": int(receipt["measurements"]["max_working_set_bytes"]),
    },
    "receipt_claim": receipt["params"]["manifest_blake3"],
    "producer_identity": {key: value for key, value in receipt["params"].items() if key.startswith("producer.") or key in {"code_git_sha", "code_dirty", "rustc_version", "target_triple", "deps_lock_blake3"}},
    "command": receipt["params"].get("command"),
    "command_argv": json.loads(receipt["params"]["command_argv"]),
}
(results / "latest.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
PY

set +e
RESULT_STATUS="$(python3 "$HERE/compare.py" \
  --baseline "$HERE/baseline.json" \
  --latest "$DATA/results/latest.json" \
  --results "$DATA/results")"
COMPARE_CODE=$?
set -e
echo "GIAB v5.0q chr20 benchmark status: $RESULT_STATUS; latest report: $DATA/results/latest.json"
if [ "$COMPARE_CODE" -eq 1 ]; then
  echo "GIAB evidence diverged; review candidate and use cargo xtask giab baseline propose" >&2
fi
exit "$COMPARE_CODE"
