# Reuse evidence without the original alignments

This walkthrough extracts a small dataset once, moves its portable directory,
removes the generated input files, and answers candidate and panel questions from
the saved evidence. It checks that the saved-data answers equal fresh extraction.
No private inputs are used or removed.

Start with the [source CLI on PATH](installation.md) and run from the repository
root. The preparation step needs Python 3.9+, network access, and pinned pysam.
It uses the same [public NA18507 sources](../examples/research-filter/README.md)
as the researcher tutorial, with four independently generated candidate SNVs.

```sh
python3 -m venv /tmp/rosalind-reuse-env
/tmp/rosalind-reuse-env/bin/pip install pysam==0.23.3
```

The following block creates a new temporary directory for every run. It deletes
only the `inputs` directory that it just generated, after retaining the selection
files needed for offline queries. It leaves outputs in the printed directory.

<!-- smoke:reuse-quickstart -->
```sh
/tmp/rosalind-reuse-env/bin/python - <<'PY'
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

repository = Path.cwd()
binary = shutil.which("rosalind")
if binary is None:
    raise RuntimeError("Put the current source preview CLI on PATH first")
prepare = repository / "examples/research-filter/prepare.py"
if not prepare.is_file():
    raise RuntimeError("Run this example from the repository or native bundle root")

work = Path(tempfile.mkdtemp(prefix="rosalind-reuse-")).resolve()
inputs = work / "inputs"
subprocess.run([sys.executable, str(prepare), str(inputs)], check=True)

def run(*arguments):
    subprocess.run([binary, *map(str, arguments)], check=True)

reference = inputs / "ex1.fa"
alignment = inputs / "sample.bam"
targets = inputs / "targets.bed"
candidate_lines = (inputs / "candidates.vcf").read_text().splitlines(keepends=True)
headers = [line for line in candidate_lines if line.startswith("#")]
records = [line for line in candidate_lines if not line.startswith("#")]
assert len(records) == 4, "The pinned tutorial should produce four candidate records"
# Change the research shortlist to two of the supplied SNVs; invent no new calls.
candidates = work / "updated-candidates.vcf"
candidates.write_text("".join(headers + records[:2]))
common = ("--reference", reference, "--alignments", alignment,
          "--memory-budget-mb", "256")

# Keep fresh answers for comparison before the generated inputs disappear.
run("analyze", "evidence", *common, "--sites", candidates,
    "--fields", "depths,alleles", "--format", "arrow-ipc",
    "--output", work / "fresh-candidates.arrow")
run("analyze", "panel-qc", *common, "--regions", targets,
    "--min-callable-depth", "10", "--output", work / "fresh-panel.tsv")

# Store every target position and enough fields for both later questions.
run("analyze", "evidence", *common, "--regions", targets,
    "--fields", "depths,alleles,quality-sums", "--format", "arrow-ipc",
    "--cache-dir", work / "cache", "--output", work / "all-evidence.arrow")
receipt = json.loads((work / "all-evidence.arrow.manifest.json").read_text())
# "execution.evidence_dataset_manifest" is one literal key in measurements.
original_manifest = Path(receipt["measurements"]["execution.evidence_dataset_manifest"])
portable = work / "portable"
shutil.move(str(original_manifest.parent), str(portable))
shutil.copy2(targets, work / "targets.bed")
shutil.rmtree(inputs)
assert not inputs.exists()

dataset = portable / "evidence-dataset.manifest.json"
run("dataset", "verify", "--dataset", dataset)
run("dataset", "extract", "--dataset", dataset,
    "--sites", candidates, "--fields", "depths,alleles",
    "--memory-budget-mb", "256", "--format", "arrow-ipc",
    "--output", work / "reused-candidates.arrow")
run("dataset", "panel-qc", "--dataset", dataset,
    "--regions", work / "targets.bed", "--min-callable-depth", "10",
    "--memory-budget-mb", "256", "--output", work / "reused-panel.tsv")

for fresh, reused in (("fresh-candidates.arrow", "reused-candidates.arrow"),
                      ("fresh-panel.tsv", "reused-panel.tsv")):
    assert (work / fresh).read_bytes() == (work / reused).read_bytes(), (fresh, reused)
    run("verify", "--manifest", work / (reused + ".manifest.json"))
print("Verified identical candidate evidence and panel QC without original inputs.")
print(f"Inspect the dataset, outputs, and receipts in: {work}")
PY
```

## What this demonstrates

The original receipt stores the portable manifest path in `measurements` under
the literal key `execution.evidence_dataset_manifest`. Move or copy the **entire
directory containing that manifest**, not the manifest alone. The descriptor and
partition files are required.

The new queries use only the moved dataset and the retained VCF/BED selections.
The updated VCF narrows the original four candidates to a two-SNV shortlist;
the saved dataset covers all target positions, so that changed selection needs
no new alignment extraction.
They do not reopen the original BAM or FASTA. Equality of the result files checks
the candidate counts and target summaries for this pinned example; receipts
differ because the recorded recipes and dependencies differ.

The example stores `depths,alleles,quality-sums`. Candidate extraction projects
down to `depths,alleles`; panel QC requires `depths,quality-sums`. A later strand
analysis would need the `strands` group to have been stored too. Similarly, a
position outside the saved targets cannot be answered by this dataset alone.
Missing fields or positions cause refusal, not invented zero values.

This is an offline reuse check, not a performance or biological-accuracy
benchmark. Extraction, hashing, and verification have costs; measure repeated
workloads before claiming a speedup. Counts retain the original read-based
profile and cannot reconstruct observations excluded by its filters.

Continue with the [dataset reference](reusable-evidence.md) for filling missing
loci from verified alignments, or the
[Python, R, and SQL examples](../examples/persisted-evidence/README.md) for downstream
analysis. [Core concepts](concepts.md) explains profile compatibility and zero
versus missing evidence.
