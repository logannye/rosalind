#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
REFERENCE="${1:-$ROOT/examples/data/illumina_toy/reference.fa}"
BAM="${2:-$ROOT/examples/data/illumina_toy/alignments.bam}"
RESULTS="${3:-$HERE/results}"
IMAGE="${ROSALIND_PLATFORM_IMAGE:-rosalind-platform:local}"

# Evidence directories are create-new: reruns cannot erase a prior measurement.
if [ -e "$RESULTS" ]; then echo "results already exist: $RESULTS" >&2; exit 2; fi
mkdir -p "$RESULTS"
RESULTS="$(cd "$RESULTS" && pwd)"
docker build --platform linux/amd64 -f "$HERE/Dockerfile" -t "$IMAGE" "$ROOT"
docker image inspect "$IMAGE" > "$RESULTS/image.json"
docker run --rm --network none --platform linux/amd64 \
  -v "$REFERENCE:/data/reference.fa:ro" -v "$BAM:/data/input.bam:ro" -v "$RESULTS:/results" \
  "$IMAGE"

# Every implementation receives identical prepared inputs and actual cgroup
# memory/swap limits. The 1 MiB declaration-refusal check is recorded separately.
for budget in ${ROSALIND_PROBE_BUDGETS_MB:-48 128 512}; do
  [[ "$budget" =~ ^[1-9][0-9]*$ ]] || { echo "invalid probe budget: $budget" >&2; exit 2; }
  for implementation in rosalind pysam bcftools; do
    probe="$RESULTS/memory-probes/$budget/$implementation"
    mkdir -p "$probe"
    set +e
    docker run --rm --network none --platform linux/amd64 --memory "${budget}m" --memory-swap "${budget}m" \
      -v "$RESULTS/inputs:/data:ro" -v "$probe:/outputs" \
      --entrypoint /opt/rosalind-platform/memory_probe.sh \
      "$IMAGE" "$implementation" "$budget" > "$probe/stdout" 2> "$probe/stderr"
    code=$?
    set -e
    printf '%s\n' "$code" > "$probe/exit_code"
  done
done
RESULTS="$RESULTS" python3 - <<'PY'
import hashlib, json, os
from pathlib import Path
root = Path(os.environ["RESULTS"])
def digest(path):
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()
report_path = root / "report.json"
report = json.loads(report_path.read_text())
probes = []
for probe in sorted((root / "memory-probes").glob("*/*")):
    code = int((probe / "exit_code").read_text())
    output = probe / "features.tsv"
    result = {
        "implementation": probe.name, "budget_mb": int(probe.parent.name),
        "exit_code": code, "output_created": output.exists(),
        "status": {0: "completed", 3: "refused", 4: "breached", 137: "killed"}.get(code, "failed"),
        "stderr_sha256": digest(probe / "stderr"),
    }
    if code == 0:
        expected = root / "outputs" / f"{probe.name}-1.tsv"
        result["byte_identical_to_unconstrained"] = output.exists() and digest(output) == digest(expected)
    probes.append(result)
report["matched_cgroup_budgets"] = probes
report["declaration_refusal"] = json.loads((root / "raw/refusal.json").read_text())
report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
if any(p.get("byte_identical_to_unconstrained") is False for p in probes):
    raise SystemExit("memory budget changed successful output bytes")
if not any(p["implementation"] == "rosalind" and p["status"] == "completed" for p in probes):
    raise SystemExit("no matched-budget Rosalind trial completed; inspect retained probe evidence")
PY
echo "platform evidence: $RESULTS/report.json"
