# Synthetic cohort reanalysis

**Source preview on the isolated 0.6 development line.** This example shows a
three-sample cohort answering a second candidate question from saved evidence,
then explicitly filling two unmeasured loci into an immutable child snapshot.
The [CLI/Python preview guide](../../docs/cohort-cli-preview.md) explains the
interfaces; the [cohort contract](../../docs/cohort-contract.md) specifies their
scientific and storage semantics. These additions are separate from 0.5 release
closure and are not a claim of published availability, partner validation,
biological accuracy or representative performance.

All inputs are authored synthetic data, redistributable under this repository's
MIT OR Apache-2.0 license. There are no human sequences, participant identities,
network downloads or random seeds. The 64-base all-A reference and one-base `1M`
reads keep every observation directly inspectable. Do not treat this as a realistic
sequencing assay or a replacement for adversarial/native decoder fixtures.

## Files and independently specified expectations

- `reads.tsv` specifies every group of identical observations: specimen, one-based
  position, base, number of copies, MAPQ, BQ and SAM flag. Read names and order are
  deterministic; each BAM has its own named `@RG SM` sample.
- `members.tsv` maps asserted specimen IDs to named samples and stored BED selection.
  It is fixture metadata; `run_cohort.py` writes a separate native import table.
- `expected.json` contains hand-specified complete per-locus depths/counts, the
  partial-output rows with explicit nulls, and candidate-summary totals
  before and after filling missing loci. It is an oracle, not Rosalind output.
- `prepare.py` expands the reads into deterministic SAM, BAM/BAI and FASTA/FAI using
  pinned pysam. It records source/generated file hashes without paths or timestamps.
- `validate.py` independently parses these simple SAM reads without htslib, checks
  the authored arithmetic, and compares current Rosalind native/saved/reused output.
  Its tiny in-memory tables are fixture oracles, not a bounded cohort implementation.
- `run_cohort.py` exercises native cohort import, strict/partial planning, saved-only
  queries, relocation, explicit extension, verification and receipt replay. It
  compares the second query and completed child against fresh native extraction.

Every candidate has REF A. At the default technical depth threshold 10, a member is
ALT-supported only when depth-eligible and at least one requested ALT read exists.
This is read support, not a genotype frequency or calibrated confidence.

| Position / ALT | Specimen A depth / ALT | Specimen B depth / ALT | Specimen C depth / ALT |
|---|---|---|---|
| 10 / C | 10 / 4 | 9 / 1, below threshold | 10 / 0 |
| 20 / G | 0 / 0, stored | 12 / 4 | 3 / 0, below threshold |
| 30 / T | 10 / 1 | 10 / 5 | **Unmeasured**; source would yield 0 / 0 |
| 40 / C | 12 / 2 | **Unmeasured**; source would yield 10 / 3 | 0 / 0, stored |

Specimen B at position 10 also has one duplicate-flagged C, one low-MAPQ C,
one low-BQ C and one ambiguous N. Thus prefilter/aligned/callable depth is 13/11/9,
with one observation in each corresponding exclusive filter counter. None of these
four contributes to callable ALT support.

| Position | Requested | Observed | Depth-eligible | ALT-supported | Callable total | ALT total |
|---|---:|---:|---:|---:|---:|---:|
| 10 | 3 | 3 | 2 | 1 | 29 | 5 |
| 20 | 3 | 3 | 1 | 1 | 15 | 4 |
| 30 | 3 | 2 | 2 | 2 | 20 | 6 |
| 40 | 3 | 2 | 1 | 1 | 12 | 2 |

The last two totals include all observed members, including low-depth ones;
`expected.json` also specifies totals restricted to depth-eligible members.
Unmeasured counts/fractions/eligibility are null. Stored zero depth has exact zero
counts, false eligibility and an undefined fraction. A positive low-depth row has
a defined ALT/depth fraction but fails the technical screen. Filling C:30 changes
its status to observed zero; filling B:40 adds 10 callable and 3 ALT observations.

The [retained local validation note](VALIDATION.md) identifies the tested binary
and limits of the earlier single-sample primitive validation. A new cohort run
records its own exact binary identity and results in `cohort-validation/report.json`.

## Run the cohort preview

Build the binary from the isolated cohort development branch. Use a new fixture
directory; the script refuses existing output rather than replacing an earlier
success or failure:

```sh
cargo build --locked --bin rosalind
python3 -m venv /tmp/rosalind-cohort-fixture-env
/tmp/rosalind-cohort-fixture-env/bin/pip install pysam==0.23.3
/tmp/rosalind-cohort-fixture-env/bin/python examples/cohort-reanalysis/prepare.py \
  /tmp/rosalind-cohort-preview-fixture
/tmp/rosalind-cohort-fixture-env/bin/python examples/cohort-reanalysis/run_cohort.py \
  /tmp/rosalind-cohort-preview-fixture --binary "$PWD/target/debug/rosalind"
```

Keep that binary unchanged during the run; the script checks its identity so one
report cannot silently combine results from concurrent rebuilds.

The native processes use a cooperative 512 MiB budget by default; change it with
`--budget-mb`. This does not claim an OS memory cap. The Python demonstration loads
only the deliberately tiny fixture tables and is outside the native budget.

The demonstration:

1. Checks prepared-input hashes and the independent simple-SAM arithmetic; creates
   projected single-sample datasets and imports verified ordinary copies.
2. Plans the first four candidates, confirms strict refusal, then emits partial
   evidence and summaries matching the authored missingness and threshold oracle.
3. Relocates the cohort and makes both original input paths and imported dataset
   paths unavailable. A second list asks for G at positions 10 and 20, including a
   new ALT at a stored position; saved-only results verify and replay to identical
   TSV bytes.
4. Restores raw inputs and compares all six second-query sample/candidate rows
   against fresh native extraction, including every depth, allele and exclusive
   filter counter. The report records both costs for that matched question.
5. Supplies explicit local source mappings for B and C. Extension preserves A:4/0,
   B:3/1 and C:3/1 retained/computed loci, publishes a child last, and leaves the
   original snapshot's partial output byte-identical.
6. Compares all 12 child rows against fresh extraction and checks completed
   summary arithmetic. The child also queries, verifies and replays with original
   inputs unavailable. A `finally` block restores the prepared input directories.

Open `cohort-validation/REPORT.md` for the sample-by-candidate table.
`cohort-validation/report.json` retains every command and exit code, expected
refusals, binary/fixture/script hashes, native receipt measurements, import and
relocation costs, storage bytes, fresh-versus-saved elapsed times, and extension
hashing/record-visit counts. Failed runs retain their report and outputs for
inspection; rerun using a newly prepared directory.

These are single ordered measurements of tiny authored data. Process startup,
source hashing, verification, encoding and finalization are included; the script
does not establish representative speedups or an economic break-even point.
Independent team usage remains unmeasured.

With a matching preview Python installation, the resulting snapshot can be opened
without original alignments:

```python
import json
from pathlib import Path
from rosalind import open_cohort

work = Path("/tmp/rosalind-cohort-preview-fixture/cohort-validation")
report = json.loads((work / "report.json").read_text())
cohort = open_cohort(report["cohort"], report["child_snapshot"], binary=report["binary"])
print(cohort.plan(sites=work / "second.vcf"))
cohort.summarize(work / "python-summary.tsv", sites=work / "first.vcf", format="tsv")
```

## Prepare and validate current primitives

Use Python 3.9+ with `pysam==0.23.3` and an existing current-source Rosalind binary.
From the repository root, with a fresh output path:

```sh
python3 -m venv /tmp/rosalind-cohort-fixture-env
/tmp/rosalind-cohort-fixture-env/bin/pip install pysam==0.23.3
/tmp/rosalind-cohort-fixture-env/bin/python examples/cohort-reanalysis/prepare.py \
  /tmp/rosalind-cohort-fixture
/tmp/rosalind-cohort-fixture-env/bin/python examples/cohort-reanalysis/validate.py \
  /tmp/rosalind-cohort-fixture --binary /absolute/path/to/rosalind
```

Preparation and validation outputs are create-new. The validator:

1. Checks all 12 authored complete observation rows against an independent simple
   SAM parser, plus all authored partial rows and summary arithmetic.
2. Extracts each stored selection with `depths,alleles`, preserving true zeros;
   compares raw counts and specimen B's exclusive filter counters.
3. Copies intact portable leaf datasets to a new location and temporarily moves
   the generated original-input directory out of reach. Saved queries must match
   stored TSV bytes; queries asking B:40 or C:30 must refuse without output.
4. Checks absent-field refusal separately from missing loci, then restores inputs
   in a `finally` block.
5. Uses current serial `--reuse-dataset` to fill requested loci; compares complete
   output bytes and counts against fresh extraction, checks reused/computed locus
   counts (4/0, 3/1 and 3/1), and verifies the resulting artifacts.

`validation/report.json` records the actual binary hash/version, fixture hash,
executed commands, results and reuse counts. Generated BAMs/caches/reports stay in
the chosen output directory; do not commit them. Repeated preparation with the same
pinned tools should produce the same `preparation.json` hashes.

This older validator establishes single-sample primitives on a deliberately tiny
fixture. Use `run_cohort.py` above for the cohort demonstration. Broader resource,
cancellation, compatibility and adversarial cases are covered separately by the
native regression suite; neither fixture validator establishes independent usage.


## Explicit paired follow-up (paired development branch)

After a successful cohort demonstration, the separate
`codex/cohort-pairs-preview` branch can compare the supplied `pairs.tsv`: A→B,
B→A, and A→C. This ordering is explicit and does not infer biological pairing.

```sh
python3 examples/cohort-reanalysis/run_pairs.py \
  --demo-report /path/to/prepared/cohort-validation/report.json \
  --binary /absolute/path/to/rosalind \
  --output /new/path/to/pairs-validation
```

The runner checks every paired row against authored depths/ALT counts and Python
`fractions.Fraction`, including missingness, zero denominators, swapped direction
and low-depth eligibility, then verifies the native receipt. Reports retain full
commands and failures. The [paired contract](../../docs/cohort-pairs-contract.md)
explains the exact decimal-string difference columns. These synthetic examples
establish neither independent partner adoption nor biological interpretation.
