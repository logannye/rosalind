#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
if [[ "$#" -ne 1 ]]; then
  echo "usage: $0 path/to/candidate.whl" >&2
  exit 2
fi
# Resolve before entering the isolated work directory; never install a different
# candidate merely because the caller supplied a relative path.
WHEEL="$(python -c 'import pathlib, sys; print(pathlib.Path(sys.argv[1]).resolve(strict=True))' "$1")"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
python -m venv "$WORK/venv"
"$WORK/venv/bin/pip" install "$WHEEL"
"$WORK/venv/bin/pip" install 'pysam==0.23.3'
"$WORK/venv/bin/rosalind" --version
"$WORK/venv/bin/rosalind" sort --input "$ROOT/examples/data/illumina_toy/alignments.bam" --output "$WORK/sorted.bam"
"$WORK/venv/bin/rosalind" reference build --fasta "$ROOT/examples/data/illumina_toy/reference.fa" --output "$WORK/reference.rref"
# Run outside the checkout so import cannot accidentally use source Python files.
cd "$WORK"
export ROSALIND_SMOKE_SOURCE_ROOT="$ROOT"
"$WORK/venv/bin/python" - <<'PY'
import json
import os
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

# The packaged executable supplies every scaffold file and runs conformance.
# Before registry publication, resolve its SDK dependencies to this candidate's
# source checkout. Fetch the generated SDK consumer's dependencies explicitly:
# on Linux maturin's container does not populate the host Cargo registry. The
# subsequent locked offline build proves compilation needs no further network.
source = Path(os.environ["ROSALIND_SMOKE_SOURCE_ROOT"])
scaffold = Path("wheel-smoke-analyzer").resolve()
subprocess.run([str(binary), "new", "analyzer", "wheel-smoke-analyzer", "--output", str(scaffold)], check=True)
assert f'version = "={native_version}"' in (scaffold / "Cargo.toml").read_text()
target = source / "target" / "wheel-sdk-smoke"
sdk_args = [
    "--manifest-path", str(scaffold / "Cargo.toml"),
    "--config", "patch.crates-io.rosalind-bio.path=" + json.dumps(str(source)),
    "--config", "patch.crates-io.rosalind-build-info.path=" + json.dumps(str(source / "crates/build-info")),
]
subprocess.run(["cargo", "fetch", *sdk_args], check=True)
assert (scaffold / "Cargo.lock").is_file(), "dependency preparation must lock the scaffold"
subprocess.run([
    "cargo", "build", "--locked", "--offline", *sdk_args,
    "--target-dir", str(target),
], check=True)
with Path("scaffold.conformance.json").open("w") as report:
    subprocess.run([str(binary), "conformance", "analyzer", "--binary",
                    str(target / "debug/wheel-smoke-analyzer"), "--json"], stdout=report, check=True)
assert json.loads(Path("scaffold.conformance.json").read_text())["passed"]
print("packaged scaffold built offline against candidate SDK and passed packaged conformance")
PY
