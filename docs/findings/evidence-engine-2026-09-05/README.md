# Evidence engine validation — 2026-09-05

These are local development results on macOS 26.2 arm64, from the edited
`codex/evidence-engine-roadmap` tree based on
`d6b40019a24912ebdd786f59638b4edc20c41a27`. They are not a published release,
independent-user validation, or a whole-genome benchmark. Exact executable hashes
and receipt producer metadata are retained; the benchmark and wheel executables
can have different hashes because wheel builds strip the binary.

## Functional and package checks

- **467 Rust tests passed**, no ignored tests, across 51 workspace test/doc-test
  suites. Includes adversarial CIGAR/quality/flag cases, historical receipts,
  1/2/8-worker execution, tile/budget invariance, corrupt/resumed caches, input
  mutation, failed-worker cancellation, output rollback, and dataset differences.
- Workspace/all-target Clippy with warnings denied; Rust 1.83 minimum-version
  check; formatting and documentation builds passed.
- **8 Python tests passed on 3.9.6 and 3.11.15**: lazy bounded streaming,
  finalization, early close, malformed producers, child errors, resource cleanup,
  and version matching.
- Actual local macOS arm64 wheel installed into fresh Python 3.9 and 3.11 virtual
  environments outside the checkout. Both completed legacy analysis, exact Arrow
  iteration, materialization, verification, offline replay, and panel QC.
  The installed CLI generated an analyzer and ran its conformance checks.
  Before publication, that scaffold used explicit local candidate SDK patches;
  this is not a fresh crates.io installation claim.
- All three publishable crates packaged and verified. Root package verification
  used local patches for unpublished receipt/build-info dependencies.
- Three fuzz targets ran for 15 seconds each (16 seconds elapsed): canonical
  receipts 7,171,419 runs; argv 4,782,132; legacy 4,609,338. No crashes. This is a
  smoke workload, not sustained fuzz assurance.
- Release-tooling tests, actionlint 1.7.12, static workflow policy, and shellcheck
  passed. The standalone evidence analyzer built. Nextflow 25.04.8 (3 tasks) and
  Snakemake 9.8.1 (4 rules) executed the real-origin panel example and verified
  outputs locally. Those executions were performed by the implementation team.

Raw local test/package/fuzz logs are in [validation](validation/).
[Package identity](package-identity.json) records the installable wheel hash.

## Independent evidence comparison and resource curve

The [retained report](curve/report.json) compares every one of the 34 TSV fields
against a direct CIGAR traversal implemented in pysam 0.23.3 / htslib 1.21.
The tutorial uses four candidate SNVs from the small NA18507 alignment slice in
samtools' pinned example data. It is not a variant truth benchmark.

Three repetitions at each budget produced identical rows; each Rosalind artifact
also passed byte verification. Median extraction-process measurements:

| Implementation | Declared budget | Wall seconds | Peak RSS MiB |
|---|---:|---:|---:|
| Rosalind | 96 MiB | 0.0256 | 8.38 |
| Rosalind | 128 MiB | 0.0251 | 8.39 |
| Rosalind | 256 MiB | 0.0260 | 8.38 |
| Streaming pysam | none | 0.0370 | 25.77 |

Startup dominates this tiny case. The baseline does not produce Rosalind receipts;
verification is separately timed. No cgroup was applied, no cache flush was
forced, and no large-data speedup follows from these measurements. The raw
commands, process time/I/O counters, package/input identities, receipts, outputs,
and comparisons are retained in [curve](curve/).

A cold local cache computed two partitions and visited 160 alignment records.
The resumed run reused both partitions, visited zero records, and preserved the
same evidence bytes. The measured native elapsed times were 31 ms and 15 ms;
these single cache trials demonstrate avoided work, not a stable speedup ratio.

## Specific memory regressions

The original legacy Arrow reproduction produces 70,000 rows, spanning a full
65,536-row batch and a 4,464-row tail. The corrected planner refuses 32 MiB before
creating output. At 256 MiB the run succeeds and verifies; observed peak RSS was
43 MiB in the retained manual reproduction. A generated equivalent is now a
permanent CLI regression. See [legacy-arrow](legacy-arrow/).

Packed-reference windowing was exercised over a generated 128 MiB FASTA, a
50,331,800-byte `.rref`, and 16,384 one-base requests spread across the full span.
At a declared 32 MiB, packed extraction completed at **14,434,304 bytes (13.77 MiB)**
peak RSS, versus 14,270,464 bytes for indexed FASTA. The plan predicted 33,553,608
bytes, below the 33,554,432-byte budget. Receipts and command logs are in
[reference-span](reference-span/); the large synthetic reference/output files
remain local rather than being checked into the repository. Unit tests separately
check `.rref`/`.idx`/FASTA windows across packed-word and contig boundaries.

## Remaining release evidence

Linux and Intel macOS candidate builds, actual Linux container/cgroup execution,
real filesystem exhaustion, large-record/data performance curves, a published
GIAB evaluator and reviewed caller baseline, seven-day RC soak for 0.5 onward,
independent replay, and external adoption have not been established here.
Protected publication settings and registry credentials remain prerequisites.
The [implementation status](../../implementation-status.md) tracks those gates.
