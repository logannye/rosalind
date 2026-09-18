# First decoder correction: retained failed attempt

Source `9862692685ee865aa8678b185d60a45380c51ac0` corrected the CRAM allocation
model but exposed a native index-lifetime bug. **This candidate is not verified
for release.** Its complete results remain here, including failures.

The same [HG002 workload](../../representative-evidence-2026-09-06/README.md)
requested 108 runs. [report.json](report.json) records 106 completed runs, one
eight-worker CRAM panel extraction terminated by SIGSEGV, and one dependent saved
query that could not run. Every completed output matches the earlier baseline's
bytes ([comparison](baseline-comparison.json)). The [audit](audit.json) reports
zero memory-plan underestimates among completed runs and unchanged final input,
binary and harness identities, but correctly fails the overall study.

The [sanitized crash](native-crash.json) was localized to
`cram_index_free_recurse`, called during index destruction. In rust-htslib 0.44.1,
`IndexedReader::drop` closes the CRAM file before dropping its index. The pinned
HTSlib index destructor then accesses the already freed CRAM handle. This matches
the [upstream fix](https://github.com/rust-bio/rust-htslib/pull/518). Its GitHub tag
was not available from crates.io when checked; no dependency migration is claimed.

[matrix-raw.tar.xz](matrix-raw.tar.xz) retains all outputs, partial cache contents,
receipts, logs, timings, failed invocations and the executed harness snapshot.
The original genomic inputs are referenced by their verified hashes, not embedded.
[candidate-identities.json](candidate-identities.json) identifies the local wheel,
installed binary, source archive and unchanged workload manifest.

Independent checks on this same source passed:

- 577 workspace tests, full Clippy, Rust 1.83 and 38 installed-binary harness tests.
- 23 installed-package Python tests on each of Python 3.9 and 3.11; real-origin
  BAM/CRAM byte equality, verification, replay and panel/position verification.
  [Package report](package-smoke.json) and [raw evidence](package-smoke-raw.tar.xz).
- 29 Linux focused tests and seven real cgroup scenarios for each encoding,
  including verified 128 MiB completions, typed failures and established native
  kernel OOMs. [Linux provenance](linux-validation.json), [BAM report](cgroup-bam.json),
  [CRAM report](cgroup-cram.json), [BAM raw evidence](cgroup-bam-raw.tar.xz), and
  [CRAM raw evidence](cgroup-cram-raw.tar.xz).

These passing checks do not override the parallel crash. The Linux VM uses the
same physical development workstation. Cgroup charges include file cache and
are not process RSS. Timings include concurrent lightweight validation activity;
they are not an isolated performance ranking. Public release, independent users,
high-depth assays and 30-day return usage remain separate requirements.
