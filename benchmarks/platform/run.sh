#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
REFERENCE="${1:-$ROOT/examples/data/illumina_toy/reference.fa}"
BAM="${2:-$ROOT/examples/data/illumina_toy/alignments.bam}"
RESULTS="${3:-$HERE/results}"
IMAGE="${ROSALIND_PLATFORM_IMAGE:-rosalind-platform:local}"

docker build --platform linux/amd64 -f "$HERE/Dockerfile" -t "$IMAGE" "$ROOT"
rm -rf "$RESULTS"
mkdir -p "$RESULTS"
docker run --rm --network none --platform linux/amd64 \
  -v "$REFERENCE:/data/reference.fa:ro" -v "$BAM:/data/input.bam:ro" -v "$RESULTS:/results" \
  "$IMAGE"

# Retain explicit low-memory behavior for every implementation. These probes run
# separately so an OOM kill cannot erase the successful comparison artifacts.
mkdir -p "$RESULTS/memory-probes"
for implementation in rosalind pysam bcftools; do
  set +e
  docker run --rm --network none --platform linux/amd64 --memory 48m \
    -v "$REFERENCE:/data/reference.fa:ro" -v "$BAM:/data/input.bam:ro" \
    --entrypoint /opt/rosalind-platform/memory_probe.sh \
    "$IMAGE" "$implementation" > "$RESULTS/memory-probes/$implementation.stdout" 2> "$RESULTS/memory-probes/$implementation.stderr"
  code=$?
  set -e
  printf '{"exit_code":%s,"output_created":%s}\n' "$code" "$(test -s "$RESULTS/memory-probes/$implementation.stdout" && echo true || echo false)" \
    > "$RESULTS/memory-probes/$implementation.json"
done
RESULTS="$RESULTS" python3 - <<'PY'
import hashlib, json, os
from pathlib import Path
root = Path(os.environ["RESULTS"])
report_path = root / "report.json"
report = json.loads(report_path.read_text())
probes = {}
for implementation in ("rosalind", "pysam", "bcftools"):
    result = json.loads((root / "memory-probes" / f"{implementation}.json").read_text())
    for stream in ("stdout", "stderr"):
        path = root / "memory-probes" / f"{implementation}.{stream}"
        result[f"{stream}_bytes"] = path.stat().st_size
        result[f"{stream}_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    probes[implementation] = result
report["insufficient_48m_container"] = probes
report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
PY
echo "platform evidence: $RESULTS/report.json"
