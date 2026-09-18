# Local fixture validation — 2026-09-18

These results validate existing single-sample primitives on the authored synthetic
fixture. They do not validate a cohort implementation, performance, independent
users or a release candidate.

- CLI: `rosalind 0.5.0` from the local source build.
- CLI SHA256: `f6c18ea49e3b5b924bb95b6bf6bd8654b11957e7c0b03ecca3bf61783a501643`.
- Expected-values SHA256: `c69e822c613b279636ced4b805889835ac05ccdadd030612d63357b70fd129fe`.
- Preparation: Python 3.11.15, pysam 0.23.3; no downloads or original biological data.
- Two separate preparations produced byte-identical `preparation.json` values,
  including every generated SAM/BAM/index/reference hash.
- All 12 full observation rows match the independent simple-SAM parser and current
  Rosalind extraction. Authored partial rows and both summary tables pass their
  independent arithmetic checks.
- Specimen B's four exclusive filter counters match their authored expectations.
- All three relocated saved queries match original TSV bytes while original input
  paths are unavailable. Missing-locus and missing-field requests refuse without
  successful output.
- Fresh and serial reuse output bytes match for every member. Reused/computed
  locus counts are A:4/0, B:3/1, C:3/1. All three reuse receipts verify.
- Local Markdown destinations and Python compilation checks pass.

The [sanitized validation record](../../docs/findings/roadmap-foundation-2026-09-18/cohort-fixture-validation.json)
retains fixture/preparation hashes, every executed command with portable path
placeholders, refusal outcomes and observed reuse counts. It contains no generated
BAMs or caches.

Run the [documented preparation and validator](README.md#prepare-and-validate-current-primitives)
to retain a fresh `validation/report.json` with the actual tested binary identity,
commands and results. Machine-local run directories and generated binary data are
not committed. Future snapshot, comparison, lifecycle and replay requirements
remain in the [proposed contract](../../docs/cohort-contract.md).
