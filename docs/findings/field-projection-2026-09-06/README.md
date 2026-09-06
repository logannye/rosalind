# Physical evidence projection — local validation, 2026-09-06

Implemented in the source tree based on candidate preparation commit `b04dc7c`.
These measurements used local macOS arm64 release builds, with exact executable
hashes and producer identities retained in the reports and receipts. They are
source-development evidence, not published-package or independent-user evidence.

All 64 field masks preserve the values of present fields. Full schema-1 TSV and
Arrow hashes match a fixture captured before the storage refactor, including empty
streams and a canonical Arrow batch boundary. Projected schema-2 reads validate
mask/version, available counter equations, frame layout, and expected fields before
record allocation. Missing groups remain absent. Dataset diff refuses differing
masks and compares matching projected Arrow/TSV through bounded streams.

Python extraction, materialization, cache reuse, byte replay, and locus comparison
pass with 1/2/8 workers at 512/768/1024 MiB and different execution tile widths.
This sparse fixture crosses eight ownership boundaries and includes uncovered loci.
The public preview batch API and standalone analyzer example use borrowed optional
groups; full owned rows require explicit materialization.

## Pressure reproduction and correction

The deterministic 200-kb synthetic panel has approximately 100x read depth. Each
case runs three times; the harness alternates projection order and records process
wall time, RSS, I/O counters, actual admitted tile widths, input/output/executable
hashes, receipts, and independently timed verification. See [before](before/report.json)
and [after](after/report.json), including raw files alongside them.

Before buffer reuse, the full-evidence 64 MiB case exceeded its admitted budget in
all three repetitions when a smaller boundary tile allocated alongside retained
allocator pages. Each run failed with resource status and a partial artifact.
The engine now determines maximum selected local pressure from normalized
intervals, reserves its group buffers once, and reuses two bounded identity vectors.
The planner explicitly includes both identity vectors. The same 18-case matrix
then passed with identical panel TSV bytes and successful byte verification.

| Fields | Budget MiB | Median peak RSS MiB | Median extraction wall seconds |
|---|---:|---:|---:|
| Full | 32 | 21.52 | 0.386 |
| Full | 64 | 52.62 | 0.390 |
| Full | 96 | 56.72 | 0.381 |
| Depths + quality sums | 32 | 11.31 | 0.231 |
| Depths + quality sums | 64 | 11.27 | 0.235 |
| Depths + quality sums | 96 | 11.31 | 0.235 |

The permanent CLI regression runs a long zero-heavy panel through full and projected
storage at 64/96/128 MiB and checks process peaks and exact output equality.
Sparse-window tests separately show fewer indexed record visits while preserving
full output bytes and zero loci.

This controlled synthetic test isolates storage pressure. It does not establish
representative biological performance, an OS-enforced bound, cross-platform RSS,
or cohort-scale speed. Real BAM/CRAM workload curves and external use remain the
later adoption milestone. The original failing data is retained without relabeling
it as a passing run.

Validation: 495 workspace tests passed, followed by the new permanent
boundary-buffer regression. Thirteen Python tests also passed on fresh Python 3.9 and 3.11 environments
after buffer reuse, alongside workspace Clippy, Rust 1.83 all-target checks, and
the standalone analyzer crate test.
