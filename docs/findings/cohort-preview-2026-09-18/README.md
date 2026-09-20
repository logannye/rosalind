# Cohort preview: retained engineering validation

The unpublished **0.6.0** preview is developed separately from the 0.5 release in
[PR #145](https://github.com/logannye/rosalind/pull/145). This report records authored
synthetic engineering checks, not independent adoption, assay validation, or an
economic performance claim.

## Identified source and full workflow

The final local demonstration ran a clean native build from
`ba6cd43fe4583eb8c148b7f66782e16ff1ad2943` on macOS arm64. Its result receipts record
that exact `code_git_sha` and `code_dirty=false`. The executable SHA256 is
`5c51678758c3da47693f0dabe9df61703f2e00677460a47dbc467e24abc266e1`.
The [retained report](validated-ba6cd43.json) includes the fixture/script hashes,
all 28 commands, intentional strict-coverage refusal, outputs, timings and source
hashing/extension measurements. Local path prefixes are redacted; the original
report SHA256 is retained. Original receipts were not rewritten.

The [executable tutorial](../../../examples/cohort-reanalysis/README.md) establishes:

- Three named synthetic samples agree with a separate single-base SAM arithmetic
  oracle, including exclusive filter counters, observed zeros and low depth.
- Imported snapshots work after relocation while original alignments, reference
  and original dataset locations are unavailable.
- Six cells for a different candidate/ALT question agree with fresh native
  extraction across every shared depth/allele/filter field.
- Explicit extension computes only B:40 and C:30; the unaffected sample performs
  no native record visits. Twelve resulting cells agree with fresh extraction.
- The original snapshot's output bytes remain unchanged. Second-question and
  extended-query receipts replay to identical TSV bytes after relocation.
- The matching 0.6 Python example reads the same snapshot and produces the same
  native summary. Separate installed-wheel validation passed all 30 Python tests,
  including seven cohort tests, outside the checkout without a version override.

The measured one-time import and verification was about 84 ms, and the imported
store occupied 30,467 bytes (51,503 after extension). The second saved query took
about 44 ms; three serial fresh native queries took about 79 ms. These are one
ordered run over a 64-base synthetic reference, dominated by startup and
verification overhead. They do not establish a speedup or a break-even point for
real cohorts. Preparation, first extraction, copies, hashing, verification and
storage are all separate costs in the report.

## Correctness and execution checks

The source tests cover exact uint64 counts/overflow, reference/profile/sample
compatibility, required-field projections, strict and partial coverage, true zero
versus null, immutable imports/parents, late inventory mutation, publication
failure and cancellation. Managed extraction/summary tests compare both Arrow and
TSV bytes across 384/512/768 MiB admitted budgets and widths 1/3/16,384.

CLI integration tests cover VCF/VCF.gz/BCF selection, explicit extension, all four
extract/summary × Arrow/TSV relocated replay combinations, altered consumed
partitions, and a normalized 2,048-candidate query larger than the inline replay
envelope. Large replay operands use a bounded request file rather than relying on
OS argument limits; unrelated commands and mismatched enforcement are rejected.

Local validation included the workspace tests, Clippy for all targets, Rust 1.83
checks, doctests and guide packaging checks. The first workspace run exposed a
release-tool test hardcoded to 0.5; it now takes the policy's version and all 38
release-tool tests pass. Two test-only Clippy clone warnings were corrected.
The [platform CI run](https://github.com/logannye/rosalind/actions/runs/35384722832)
passed all ten jobs, including the authored workflow and saved-only Linux cgroup
probe. GitHub tested PR head `ba6cd43fe4583eb8c148b7f66782e16ff1ad2943` through
its synthetic merge commit `769d5f071508abbfe1181c38a9f713baad13e24e`; the Linux
receipts identify the latter commit and record `code_dirty=false`.

The [Linux cgroup report](linux-cgroup.json) passed three cases: extraction and
summary under an observed 512 MiB cgroup-v2 limit with swap disabled, and refusal
of an inadmissible 1 MiB declared budget without successful or partial artifacts.
Both successful outputs exactly match the host baseline, have verified receipts,
and record zero original-alignment decodes. Only saved objects, the candidate
file and binary were mounted; original alignments and the reference were absent.
There were no OOM events. This demonstrates these operations on the small fixture,
not a universal memory bound or large-cohort scalability.

The [retained raw CI archive](linux-ci-evidence.tar.gz) contains controller logs,
receipts, verification results, cgroup observations and the full workflow report.
The corresponding [GitHub artifact](https://github.com/logannye/rosalind/actions/runs/35384722832/artifacts/10563732042)
was retained before its expiration. Archive SHA256:
`abeb45ae05c5be94cf9b872173cbe75e5e43083dcb722d02bde6d50a7c579966`.

A second clean local run staged the offline guides outside the checkout and ran
the entire tutorial from that bundle using the installed 0.6 wheel and the clean
`ba6cd43` binary. All 28 commands passed. This exposed and corrected missing
cohort fixture/script assets in guide packaging; CI now executes the staged
example as well.

## Retained failures and corrections

- [Initial tutorial attempt](attempt-demo-invocation.json): evidence and offline
  queries passed, but the demonstration invoked `reproduce` with a positional
  receipt instead of `--manifest`. The script was corrected and rerun from a new
  fixture directory; the failed report remains separate.
- [Cross-build failure](cross-build-before-fix.json): an intact older snapshot was
  rejected because the reader compared recorded producer fields with its own
  build identity. The fix in `64349f7` verifies the recorded receipt self-hash,
  semantic keys and exact lineage while preserving original producer provenance.
  The [corrected check](cross-build-after-fix.json) opens and verifies that same old
  snapshot with a later binary, without modifying its snapshot or receipt bytes.
  Regression tests also reject unsealed tampering and incorrectly resealed lineage.

Independent researchers, recurring cohort teams, clinical interpretation, and
representative large-cohort economics remain unestablished. Their absence is not
a software shipping gate, and is not replaced by these authored checks.
