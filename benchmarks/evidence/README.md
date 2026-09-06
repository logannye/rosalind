# Exact-evidence budget curve

Prepare the content-locked [NA18507 tutorial](../../examples/research-filter/README.md).
Then run with the same Python environment containing pysam 0.23.3:

```sh
python benchmarks/evidence/run.py --binary target/debug/rosalind \
  --data /tmp/research-filter --output /tmp/evidence-curve --cache
```

Default: three repeats at 96/128/256MiB, three independent pysam extractions, and
optional cold/resumed cache runs. Results directories are create-new. Any differing
row, repeat, failed verification, or nonzero extraction exit makes the run fail;
the report and raw commands/timing/stderr remain available for diagnosis.

The baseline directly traverses CIGAR and applies the documented flag, MAPQ, BQ,
missing-quality, read-cycle, allele, strand, and histogram semantics. It streams
one locus at a time and retains only supplied candidate metadata. It does not use
the Rosalind implementation or sampled pileup defaults.

`report.json` retains input/binary hashes, source preparation provenance, versions,
exact argv, process time/RSS, raw timing I/O counters, artifact sizes/hashes, every
semantic comparison, receipt setup/hash/combined-analysis-encoding measurements,
and independently timed verification. The baseline does not implement provenance;
its extraction cost is labeled separately from Rosalind's end-to-end contract.

No cgroup is imposed by this portable harness. A declared budget is not a measured
hard limit for the baseline. The tiny four-locus example establishes repeatability
and agreement only; scale the supplied candidate set and aligned data, retaining
the same evidence, before making performance claims. Kernel-only encoding timing,
large-data curves, and independent-machine reproduction remain separate gates.

The synthetic projection pressure harness separately measures the storage change:

```sh
python benchmarks/evidence/projection.py --binary target/release/rosalind \
  --output /tmp/projection-pressure
```

It generates a deterministic 200-kb panel with approximately 100x read depth,
compares default panel fields with ALL at three budgets and three repetitions,
and retains exact commands, input/binary hashes, process RSS/time/I/O counters,
actual microtile widths, record visits, verification timing, and output equality.
This measures a controlled pressure case; it is not representative cohort or
whole-genome performance evidence.

## Representative HG002 workloads

`representative.py` measures an explicit case list against content-locked inputs.
It does not expand budgets, tiles, workers and formats into an implicit Cartesian
product. The maintained preparation script produces five workloads: up to 1,000 SNV
probes distributed across the full HG002 chromosome-20 BAM, and sparse SNV/dense
panel requests over an unfiltered 1-Mb window represented as both BAM and CRAM.
The source dictionary and matching full reference are preserved. Probe ALT bases
and panel intervals define research queries; they are not validated variant calls
or a target-enrichment assay. Preparation records public source URLs, checksums,
derivation commands, library versions and exact selection sizes.

```sh
# Use an environment containing pysam==0.23.3.
python benchmarks/evidence/prepare_representative.py \
  --download-cache /tmp/hg002-downloads --output /tmp/hg002-prepared

# Inspect the ordered invocation inventory before running the measurements.
python benchmarks/evidence/representative.py --binary target/release/rosalind \
  --manifest /tmp/hg002-prepared/workloads.json \
  --output /tmp/hg002-curve --plan-only

python benchmarks/evidence/representative.py --binary target/release/rosalind \
  --manifest /tmp/hg002-prepared/workloads.json --output /tmp/hg002-curve
```

The output directory must be new. `report.json` is updated after every invocation;
raw argv, stderr, OS timing/I/O counters, TSV outputs, receipts, partials and oracle
statistics remain on disk when a run fails or refuses its budget. Exact per-run
observations and medians/minima/maxima are both retained. Failed measurements are
excluded from successful timing summaries, and their requested/failed counts
remain visible. There are no checked-in performance conclusions until the
corresponding retained report has actually run.

The default prepared matrix uses three repeats and three serial budget/tile cases, with
additional 2-worker and 8-worker cases on the bounded window. Per-workload gates
check admitted worker counts and microtile widths in receipts, and require actual
microtile-count/record-visit differences when multiple widths are requested by a
workload gate. These counters establish a change in work; they do not measure
worker utilization or actual maximum fetch span. Merely passing different CLI
values does not establish that execution changed. A sparse
chromosome request can legitimately have one selected locus per fetch window;
the dense/window workloads provide the scheduling-pressure checks. Budgets and
tile widths change together in the default serial cases, so this is a configuration
sweep, not an isolated estimate of the effect of changing only RAM. A refusal or
failed verification stays a failed matrix result rather than disappearing from
the curve.

Three reuse modes are separately measured:

- `cold`: create a new application cache for this workload, case and repetition.
- `resumed`: run native extraction again with that cache. Inputs are still hashed;
  successful reuse must report zero alignment visits and zero computed partitions.
- `persisted`: run serial `dataset extract` using the portable dataset manifest
  published by the cold run. This reads and verifies stored evidence without
  reopening or rehashing the original alignment/reference. Its receipt must
  report zero alignment visits and `original_sources_rehashed=false`.

All three outputs must equal the independent oracle and the fresh native output.
“Cold” refers to Rosalind's application cache. The harness does not flush or make
claims about the operating system's filesystem cache. It deterministically
alternates independent case groups while keeping cold/resume/persisted order.

## Workload manifest contract

The JSON manifest is versioned (`schema: 1`), capped at 4 MiB, and supplies:

| Field | Meaning |
| --- | --- |
| `label`, `provenance` | Human-readable workload name and retained preparation evidence. |
| `repeats`, `seed` | At least three repeats; deterministic case-order seed. |
| `oracle` | `tile_bases` (1–16384), `max_read_len` (default250), and `max_record_bytes` (default1048576). |
| `cases` | Unique `id`, `budget_mib`, `tile_bases`, `workers`, and `cache`: `none`, `cold-resume`, or `cold-resume-dataset`. |
| `workloads` | Unique `id`; five file identities; selected case IDs; optional sample/pooling and expected row count; scheduling gates. |

Each workload contains `alignments`, `alignment_index`, `reference`,
`reference_fai`, and `selection` objects with `path` and lowercase `sha256`.
Paths are relative to the manifest unless absolute. `selection` additionally has
`kind: "sites"` or `kind: "regions"`. The oracle accepts coordinate-sorted SNV
VCF/VCF.gz/BCF or sorted BED; duplicate SNV ALTs are united and overlapping/adjacent
BED intervals are merged while streaming. Unsorted data fails explicitly instead
of being materialized and sorted in memory. File contents are checked before the
matrix starts. Native invocations perform their own source verification inside
the measured process.

The current harness checks captured input/binary/harness metadata between jobs
and rehashes those files and the workload manifest after completion. These
harness checks are outside timed invocations. A mutation fails the matrix;
unstarted jobs and a completed process whose metadata cannot be parsed remain
explicit records. Final hashes cannot prove that bytes never changed and reverted
between observations. Keep source files immutable throughout the run.

An optional `equivalence_group` joins scientifically equivalent BAM/CRAM
workloads for exact row comparisons. Source-bound receipt hashes legitimately
differ between distinct alignment encodings. `sample` names one declared sample;
`pool_samples: true` explicitly pools samples. Omit both to use the documented
automatic sample resolution. `expected_rows` checks complete denominators.

Workload `gates` specify `min_admitted_budgets` (default3),
`min_effective_tiles` (default1), and `workers` (default`[1]`). Gates apply to
native runs; a persisted query is always separately labeled serial execution.
`--smoke` permits a small diagnostic run with fewer repeats or incomplete
scheduling gates, and labels its report `diagnostic-smoke`/`smoke-passed`. Such a
run is not representative performance evidence.

## Oracle and measurement scope

`streaming_pysam.py` independently traverses each read's CIGAR. It performs one
indexed fetch per bounded selection window and accumulates fixed-width unsigned
integer arrays for the legacy field mask63: counts, strands, quality sums,
read-position sums, and both complete quality histograms. It retains no read map
or full output. Alignment flags, missing qualities, MAPQ20/BQ20, sample scope,
deletions/reference skips, overlapping mates as separate reads, zero-depth loci,
REF validation, and reverse-strand read positions match the documented evidence
profile. It uses no pileup sampling cap, BAQ adjustment, overlap collapsing or
hidden base-quality defaults. Native decoder/header/aux allocations remain
observed; post-decode envelope checks are not a hard memory barrier.
This oracle validates evidence for well-formed short-read inputs. Its malformed
CIGAR and invalid-quality rejection checks are less complete than the native
engine's; it is not a malformed-input conformance oracle.

The oracle's measured process performs selection parsing, indexed extraction and
TSV encoding. It does not implement input hashing, receipt sealing or artifact
verification. Native whole-process timing includes those production lifecycle
costs, and receipt phase timings expose hashing/setup/combined analysis-encoding/
finalization where available. Verification is timed independently for every
successful native or persisted artifact. This separation makes the cost of the
contract visible without presenting the oracle as a provenance-equivalent tool.
Setup includes hashing, so those two timing fields must not be added as disjoint
phases. Receipt budgets must match the case, both measured and receipt RSS must
fit the budget, and emitted-locus counts must match actual output rows. Prediction
underestimates are retained separately from budget breaches.

For retained-storage accounting, verification-time summaries, and a supplementary
identity/resource audit of a finished run:

```sh
python benchmarks/evidence/audit_representative.py \
  --report /tmp/hg002-curve/report.json \
  --manifest /tmp/hg002-prepared/workloads.json \
  --output /tmp/hg002-curve/audit.json
```

The audit never rewrites the original report. If harness sources were enhanced
after the run, pass `--harness-dir` pointing to an exact retained copy of all
original harness files; the original startup hashes must still match. Storage
figures describe retained files, not physical I/O or peak temporary disk use.

The report records binary/input/harness SHA-256 values, binary and dependency
versions, hardware/platform information, output sizes/hashes, process RSS, CPU and
wall time, and raw OS filesystem operation counters. Operation counts are not
converted to purported physical I/O bytes. Cooperative budget verification and
measured RSS do not establish a universal hard-RAM guarantee; the separate
[cgroup probe](cgroup-probe.md) exercises an actual container memory ceiling.

Focused, network-free checks:

```sh
ROSALIND_BENCH_BINARY="$PWD/target/debug/rosalind" \
  python -m unittest discover -s benchmarks/evidence -p test_representative.py -v
```

These exercise independent CIGAR/filter/zero-depth expectations, tile invariance,
BAM/CRAM parity, native field63 byte equality, canonical SNV/BED selections,
fingerprint failures, deterministic case dependencies, and observed-scheduling
gates. They use a synthetic fixture and do not substitute for the HG002 matrix.
