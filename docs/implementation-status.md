# Evidence engine implementation status

Updated 2026-09-06. The evidence-engine foundation was merged in PR #97. The current
development package is **0.5.0**, preparing the first feature-bearing RC; no 0.4
release was published. M1–M3 releases and independent adoption remain separate
from implementation. Current follow-up work is recorded in
[ADOPTION_ROADMAP.md](ADOPTION_ROADMAP.md).

Sample-scope and panel-default corrections now pass focused Rust, CLI, Python,
cache/replay, and existing engine conformance checks. No new candidate has been
published by this follow-up work.

| Phase | Code state | Verification state | Released | External adoption |
|---|---|---|---|---|
| M0 stabilization | Implemented: exact-or-fail pileup, streaming merge, Arrow/resource corrections, replay and release infrastructure | Local regression/package/fuzz checks pass; platform/container/GIAB release evidence pending | No new release from this work | Not established |
| M1 exact evidence | Implemented: indexed tiled engine, bounded FASTA/pack windows, semantics, CLI, Arrow/TSV, Python, tutorial | Independent field comparison, adversarial cases, budget/tile invariance and fresh local wheel checks pass | No | Not established |
| M2 panel/SDK | Implemented: target summaries, fused output, declared SDK requirements, standalone crate, workflows | Oracle/fused-output tests, separate crate and both workflow examples pass locally | No | Not established |
| M3 datasets | Implemented: verified cache/resume, bounded workers, streaming dataset comparison | Cache integrity/resume/input-mutation/1–8 worker tests and small resource curve pass; larger/external evidence pending | No | Not established |

## Current code and compatibility

Original `features`, `analyze coverage`, experimental `variants`, `.rref`, legacy
`.idx`, schema-5 receipts, analyzer scaffold, and first-party shard merge remain.
New legacy-pileup runs use exact-or-fail capacity semantics. Historical receipt
verification remains supported; replay of old sampling semantics requires its
matching producer.

`analyze evidence` and `analyze panel-qc` form the indexed path. Their scientific
profile is defined in [SEMANTICS.md](SEMANTICS.md), separately from legacy feature
and caller defaults. Schema v1 emits full rows; field requirements do not yet imply
smaller physically projected rows.

## Validation results

These checks apply to the locally edited tree, not a published candidate.
Full results and retained raw logs are in the
[validation report](findings/evidence-engine-2026-09-05/README.md).

- 467 Rust workspace tests pass, including historical schema fixtures; no ignored tests.
- Workspace/all-target Clippy, Rust 1.83 minimum-version check, formatting, and docs pass.
- Eight Python lifecycle/version tests pass on Python 3.9 and 3.11.
- Fresh local macOS arm64 wheel installs on both Python versions complete analysis,
  bounded evidence iteration, panel QC, verification, offline replay, and packaged
  scaffold/conformance. Unpublished SDK dependencies use explicit candidate source
  patches; registry-install validation remains a release gate.
- All three publishable crates package and verify with the candidate dependency patches.
- Three 15-second fuzz workloads completed without crashes; raw logs are retained.
- M0 release tooling: 24 xtask unit and 5 publication integration tests pass.
- Six Python release-helper and four streaming-baseline tests pass.
- Legacy Python API lifecycle/version tests pass, including distinct RC versions.
- Workflow policy, actionlint 1.7.12, and shellcheck pass for edited release and
  benchmark scripts.
- Standalone analyzer build and Nextflow/Snakemake real-origin examples pass locally.
- Independent pysam comparison matches all 34 fields across three budgets and
  three repetitions. Cached reuse preserves output with zero alignment-record visits.
- The original 70,000-row Arrow memory case now refuses an insufficient budget and
  completes an admitted run. Packed-reference sparse queries across a 128 MiB
  reference complete at 13.77 MiB peak RSS under a declared 32 MiB budget.

Historical July checks tested the preceding overhaul. They are not fresh evidence
for the new semantics, encoders, or scheduling changes.

### Subsequent GitHub validation (2026-09-06)

- [Candidate wheel CI](https://github.com/logannye/rosalind/actions/runs/34009965756)
  passed on Linux x86_64, macOS arm64, and macOS x86_64, including fresh Python
  3.9/3.11 installations, native/Python version agreement, analysis, verification,
  offline replay, and scaffold/conformance. These are retained CI artifacts;
  they have not been published to a package registry.
- [Application CI](https://github.com/logannye/rosalind/actions/runs/34010203922)
  passed the full Linux Rust suite, independent analyzer, package verification,
  Python/CLI workflows, advisory audit, Clippy, MSRV, contract harness, release
  automation checks, and regenerated browser verifier check.
- [Workflow integration CI](https://github.com/logannye/rosalind/actions/runs/34010203990)
  passed Nextflow and Snakemake execution on Linux. These checks are performed by
  repository automation, not independent external users.
- [Evaluator candidate CI](https://github.com/logannye/rosalind/actions/runs/34010511992)
  built the pinned Linux image and completed an offline hap.py/vcfeval comparison
  of one synthetic SNV. This verifies evaluator setup and execution; the image
  remains unpublished and does not establish the GIAB caller baseline.

## External prerequisites and unmeasured claims

No registry publication, tag, container push, or GitHub release was performed by
this work. Release tooling checks registry/trusted-publisher configuration and
protected environment approval before publication. Wheel CI results above
establish candidate installation on the three supported platforms. Standalone
native release assets and registry installations retain their publication gates.

The local Docker daemon is unavailable. Container builds, cgroup probes, and the
pinned platform benchmark have not run here. Static harness validation is not a
benchmark result.

The inspected repository has no `CARGO_REGISTRY_TOKEN`. The protected `rc` and
`release` environments retain required review, but their custom deployment
policies currently allow zero refs. PyPI trusted-publisher configuration is
unverified. These settings were inspected without changing them.

The GIAB baseline remains `not-yet-established`; the evaluator lock needs an
immutable published image. Scheduled readiness records the blocked prerequisite
without claiming evaluation. Changed caller defaults, implementation, or watched
shared processing require reviewed evidence from tested source before release.
Unrelated evidence-only changes do not automatically require a caller score.

No whole-genome speedup, universal hard RAM bound, competitive accuracy,
second-machine replay, external adoption, or 30-day repeat usage is established by
implementation. [ROADMAP.md](ROADMAP.md) separates executable and adoption gates.
