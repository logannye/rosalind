#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PYTHON="${PYTHON:-python3}"
if [[ "$#" -lt 2 ]]; then
  echo "usage: $0 candidate.whl (--candidate-source PATH | --registry-sdk)" >&2
  exit 2
fi
case "$2" in
  --candidate-source)
    [[ "$#" -eq 3 ]] || { echo "--candidate-source needs one path" >&2; exit 2; }
    SOURCE="$("$PYTHON" -c 'import pathlib, sys; print(pathlib.Path(sys.argv[1]).resolve(strict=True))' "$3")"
    ORIGIN=(--candidate-source "$SOURCE")
    ;;
  --registry-sdk)
    [[ "$#" -eq 2 ]] || { echo "--registry-sdk takes no path" >&2; exit 2; }
    ORIGIN=(--registry-sdk)
    ;;
  *) echo "choose --candidate-source PATH or --registry-sdk" >&2; exit 2 ;;
esac
# Resolve before entering the isolated work directory; never install a different
# candidate merely because the caller supplied a relative path.
WHEEL="$("$PYTHON" -c 'import pathlib, sys; print(pathlib.Path(sys.argv[1]).resolve(strict=True))' "$1")"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
"$PYTHON" -m venv "$WORK/venv"
"$WORK/venv/bin/pip" install "$WHEEL"
"$WORK/venv/bin/pip" install 'pysam==0.23.3'
"$WORK/venv/bin/rosalind" --version
"$WORK/venv/bin/rosalind" sort --input "$ROOT/examples/data/illumina_toy/alignments.bam" --output "$WORK/sorted.bam"
"$WORK/venv/bin/rosalind" reference build --fasta "$ROOT/examples/data/illumina_toy/reference.fa" --output "$WORK/reference.rref"
# Run outside the checkout so import cannot accidentally use source Python files.
cd "$WORK"
"$WORK/venv/bin/python" - <<'PY'
import json
import rosalind
import pysam
import subprocess
from pathlib import Path
from rosalind.features import _normalized_version
binary = Path("venv/bin/rosalind").resolve()
native_version = subprocess.check_output([str(binary), "--version"], text=True).strip().split()[-1]
assert _normalized_version(native_version) == _normalized_version(rosalind.__version__), (
    native_version, rosalind.__version__)
run = rosalind.iter_features("reference.rref", "sorted.bam", workdir="features")
rows = sum(batch.num_rows for batch in run.batches())
assert rows > 0
assert run.result is not None and run.result.returncode == 0
receipt = rosalind.inspect_receipt(run.result.manifest_path)
assert receipt.params["feature.schema"] == "1"
print(f"installed rosalind {rosalind.__version__}: {rows} Arrow rows; receipt inspected")
pysam.index("sorted.bam")
with pysam.AlignmentFile("sorted.bam", "rb") as bam:
    Path("targets.bed").write_text(f"{bam.references[0]}\t0\t100\tfirst-100\n")
with rosalind.iter_evidence("reference.rref", "sorted.bam", regions="targets.bed", workdir="evidence") as evidence:
    assert sum(batch.num_rows for batch in evidence) == 100
    assert evidence.result is not None
materialized = rosalind.materialize_evidence("reference.rref", "sorted.bam", "evidence.arrow", regions="targets.bed")
subprocess.run([str(binary), "verify", "--manifest", str(materialized.manifest_path)], check=True)
with Path("evidence.replay.json").open("w") as report:
    subprocess.run([str(binary), "reproduce", "--manifest", str(materialized.manifest_path),
                    "--inputs", str(Path.cwd()), "--binary", str(binary), "--json"],
                   stdout=report, check=True)
summary = rosalind.panel_qc("sorted.bam", "targets.bed", "panel.tsv", reference="reference.rref")
assert summary.returncode == 0
print("installed version identity, evidence stream, materialization, offline replay and panel QC passed")

PY
# Stage the maintained SDK guide and validate its offline links, then execute the
# Python example from the installed wheel metadata and the SDK guide verbatim.
# Candidate source patches are an explicit option; registry mode has a fresh
# Cargo home and never substitutes an unpublished local SDK.
"$PYTHON" "$ROOT/scripts/onboarding.py" stage "$ROOT" "$WORK/bundle"
mkdir -p "$WORK/bundle/examples/data"
cp -R "$ROOT/examples/data/illumina_toy" "$WORK/bundle/examples/data/"
"$PYTHON" "$ROOT/scripts/onboarding.py" smoke --bundle "$WORK/bundle" \
  --binary "$WORK/venv/bin/rosalind" --python "$WORK/venv/bin/python" "${ORIGIN[@]}"
