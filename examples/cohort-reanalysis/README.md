# Synthetic cohort contract fixture

**Design fixture, not a shipped cohort capability.** The proposed
[cohort contract](../../docs/cohort-contract.md) describes future work. This example
uses only existing single-sample `analyze evidence`, `dataset` and reuse commands.
There is no cohort runtime, public cohort schema or `open_cohort` implementation.
It does not establish partner validation, biological accuracy or representative
performance. Feature implementation may proceed on an isolated next-minor line
while 0.5 release closure and partner recruitment continue; neither publication
nor participant availability is a prerequisite for that development.

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
  It is fixture metadata, not the future cohort import format.
- `expected.json` contains hand-specified complete per-locus depths/counts, the
  proposed partial-output rows with explicit nulls, and candidate-summary totals
  before and after filling missing loci. It is an oracle, not Rosalind output.
- `prepare.py` expands the reads into deterministic SAM, BAM/BAI and FASTA/FAI using
  pinned pysam. It records source/generated file hashes without paths or timestamps.
- `validate.py` independently parses these simple SAM reads without htslib, checks
  the authored arithmetic, and compares current Rosalind native/saved/reused output.
  Its tiny in-memory tables are fixture oracles, not a bounded cohort implementation.

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
and exact limits of the observed result.

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

These checks establish current primitives on a deliberately tiny fixture. They do
not test future snapshot import, member compatibility, cross-member resource
planning, cancellation or cohort replay. Those are acceptance requirements in the
proposed contract and remain implementation work after the relevant gates.
