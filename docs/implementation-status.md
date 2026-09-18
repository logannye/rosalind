# Capability and release status

Updated **2026-09-18**. This is the authoritative availability inventory. The
[roadmap](ROADMAP.md) tracks delivery tasks; retained reports establish only the
source, workload and platform they actually tested.

## Which version can I use?

| Channel | Availability | Appropriate use |
|---|---|---|
| Public stable | [v0.1.0, June 2](https://github.com/logannye/rosalind/releases/tag/v0.1.0) | Historical legacy feature/caller commands; not the current evidence/dataset SDK tutorials |
| Current development | Package version **0.5.0** | Exact-evidence, dataset and SDK research previews, built from an identified source checkout |
| Feature-bearing RC | Not published at this review | The older unpublished RC workflow was canceled as superseded; no soak has started |
| Cohort APIs | Implemented in separate, unpublished 0.6 preview [PR #145](https://github.com/logannye/rosalind/pull/145) | Not in 0.5/main; review preview-specific validation before use |

No 0.4 release was published. Native/wheel build artifacts, author-run examples
and green CI are not registry publication or independent adoption. Follow the
[installation guide](installation.md) for the channel actually available.

## Capability inventory

“Source” means implemented in development, not necessarily published. **No
independent adoption is established for the current evidence platform.**

| Capability | Implementation | Validation evidence and limits | Published availability | Independent use |
|---|---|---|---|---|
| Exact read-level SNV evidence | Source; indexed BAM/BAI/CSI and checked CRAM preview | Independent integer observation oracle; CIGAR/filter/sample fixtures; admitted budget/tile/worker invariance. See foundation and representative reports below | Not in current public stable | Not established |
| Panel QC and fused position evidence | Source; full target denominators and uncovered positions | Independent target summaries and fused/separate equality; configurable technical thresholds have no clinical meaning | Not published | Not established |
| VCF/VCF.gz/BCF selection and annotation | Source; record, genotype and allele-order preservation | All 128 field masks, allele-sum oracle, output/annotation/replay tests in [PR101](https://github.com/logannye/rosalind/pull/101) | Not published | Not established |
| Physical field projection and sparse fetch coalescing | Source; omitted groups do not allocate their accumulators/outputs | Corrected synthetic pressure matrix in [PR99](https://github.com/logannye/rosalind/pull/99); no universal speed claim | Not published | Not established |
| Verified cache/resume and workers | Source; canonical first-party partition extraction | Integrity, mutation, resume and worker tests. Native CRAM resume still performs complete source validation | Not published | Not established |
| Portable evidence and subset queries | Source; query stored evidence without original alignments | Relocated offline extraction/panel QC, integer preservation and projection checks in [PR102](https://github.com/logannye/rosalind/pull/102) | Not published | Not established |
| Partial-overlap reuse | Source; serial missing-locus computation | Fresh/reused output equality; explicit compatible source identity. Expanded Arrow/TSV output does not publish an expanded cache | Not published | Not established |
| Arrow/TSV and Parquet interchange | Source; bounded export and Python/R/SQL examples | uint64 preservation, error cleanup and independently executed adapters. Native directory replay for Parquet is unsupported | Evidence APIs not published; historical feature TSV is in v0.1.0 | Not established |
| Python interface | Source; matching CLI bundled with wheel | Fresh Python 3.9/3.11 installs on Linux x86_64 and macOS arm64/x86_64; native/Python identity and lifecycle tests | Current mixed package not published | Not established |
| Managed evidence analyzer SDK | Source; factory, canonical batches, cancellation, transactional artifacts, receipts/replay | 19 separately compiled conformance checks and fresh candidate source-patched builds in [PR103](https://github.com/logannye/rosalind/pull/103). These are not non-author users | Not published | Not established |
| Nextflow/Snakemake | Source; exact-evidence examples | Author-run local and Linux CI integration; individual site executors/limits need their own evidence | Guides available in source | Not established |
| GitHub Action | Existing legacy variants/features/coverage integration | Does not expose the new evidence or dataset command family; defaults preserved | Historical Action available | Current adoption not established |
| Receipts, verification and replay | Historical verification plus additive evidence recipes | Schema 1–5 compatibility, tamper/artifact checks, relocated replay. Unsigned receipts do not authenticate an author or prove biology | Historical receipts in v0.1.0; new recipes source-only | Not established for new recipes |
| Legacy features, pileup and experimental caller | Available; new legacy runs exact-or-fail | Capacity/Arrow/merge regressions; historical downsampling replay needs its original producer. GIAB caller baseline is unestablished | v0.1.0 contains historical behavior, not all current corrections | Not established here |
| Local cohort snapshots and cross-sample reanalysis | Separate unpublished 0.6 preview, [PR #145](https://github.com/logannye/rosalind/pull/145) | [Source-bound validation](https://github.com/logannye/rosalind/blob/codex/cohort-preview/docs/findings/cohort-preview-2026-09-18/README.md); synthetic fixtures, exact oracles, relocation, extension and replay | Not published; not in 0.5/main | Not established |

## Validation records and their identities

- **Foundation, September 5:** [retained report](https://github.com/logannye/rosalind/blob/main/docs/findings/evidence-engine-2026-09-05/README.md)
  records 467 Rust tests, historical receipts, Python 3.9/3.11 and fresh wheel
  checks, plus a small NA18507 34-field oracle study. Those are historical counts,
  not a claim about every later commit.
- **Projection, annotation, portable datasets and SDK, September 6:** merged
  PRs 99/101/102/103 retain their own source-specific evidence. The
  [portable report](https://github.com/logannye/rosalind/blob/main/docs/findings/portable-evidence-2026-09-06/README.md)
  includes R and DuckDB checks; the
  [SDK report](https://github.com/logannye/rosalind/blob/main/docs/findings/evidence-sdk-2026-09-06/README.md)
  includes the standalone analyzer and installed candidate smoke.
- **Representative baseline:** [PR #104](https://github.com/logannye/rosalind/pull/104)
  retains 108 requested HG002 runs over chromosome-sparse and matched BAM/CRAM
  windows. All observation comparisons passed and RSS remained below requested
  budgets, but **eight CRAM plan underestimates** were found. WGS-derived panel
  positions are not a high-depth capture-assay benchmark. Linux probes use a VM
  on the same physical workstation, not an independent physical machine.
- **First decoder correction, source 9862692:** the retained failed attempt has
  106 completed runs, one eight-worker CRAM SIGSEGV and one dependent query that
  could not run. Completed outputs match the original baseline and had no plan
  underestimates; that does not make the failed candidate release-ready.
- **Native lifetime correction, source 9b8e12f:** destroys the native CRAM index
  before its borrowed file handle; preserves the public engine's Send contract.
  The [verified report](https://github.com/logannye/rosalind/blob/main/docs/findings/cram-decoder-validation-2026-09-06/verified-9b8e12f-2026-09-18/README.md)
  records all 108 packaged representative runs passing with no baseline output
  differences or prediction underestimates, 14 real Linux cgroup cases, and
  installed Python 3.9/3.11 plus analyzer onboarding checks. The new regression
  reproduces the old SIGSEGV as a negative control. The failed attempt remains
  separately retained. Merge/publication gates remain distinct from these results.
- **Onboarding and cohort preparation:** the
  [September 18 report](https://github.com/logannye/rosalind/blob/main/docs/findings/roadmap-foundation-2026-09-18/README.md)
  records executable researcher/reuse/builder checks and a small independently
  specified multi-sample fixture. The [cohort contract](https://github.com/logannye/rosalind/blob/main/docs/cohort-contract.md) is a
  0.5 planning document. The separate [0.6 preview PR](https://github.com/logannye/rosalind/pull/145) implements cohort commands; independent partner usage remains unestablished.

Findings are linked remotely so onboarding bundles do not accidentally include
large raw benchmark archives. Retained raw failures must not be overwritten by
later successful runs.

## Platforms and assurance

Supported candidate packaging targets are Linux x86_64 and macOS arm64/x86_64.
Linux ARM and Windows remain deferred. Fresh candidate wheel CI is evidence of
installation on its tested platform/version, not evidence of public publication.
See [candidate wheel CI](https://github.com/logannye/rosalind/actions/runs/34009965756)
and each PR's current checks.

Memory accounting is cooperative and sampled RSS is not a native-allocation
sandbox. In particular, native CRAM validation may allocate before a checkpoint.
An existing Linux cgroup-v2 limit is a separate, explicitly checked assurance.
Retained Python/R/SQL consumer memory is outside the native process's budget.

“Exact” describes the declared short-read counting profile. Overlapping mates
count separately; there is no UMI/fragment consensus, indel calling, long-read
model, methylation model or calibrated clinical inference. Saved aggregates
cannot undo their original filtering. See [semantics](SEMANTICS.md).

## Release and adoption gates

The September 18 configuration update allows `main` in both protected RC/release
environments and preserves their required reviewers. The maintainer confirms
TestPyPI setup; a successful protected upload still needs to verify that mapping.
The Cargo credential for stable crates is absent from the secret-name inventory;
PyPI stable publishing remains to be verified. Anonymous GHCR retrieval of the
published hap.py evaluator is verified; see the [evaluator record](findings/evaluator-publication-2026-09-18/README.md). R02 records
configuration work separately from code completion; never place credentials here.

A real immutable prerelease starts the seven-day soak, and changes to watched
caller or shared source require reviewed evidence. Independent feedback is advisory:
submitted records are validated and their absence remains visible, but missing
participants do not block publication or isolated next-minor development. Explicit
unresolved release-blocking defects in submitted feedback still block technical
eligibility. Candidate
SDK source patches do not establish registry installation. Stable promotion must
pass fresh registry-only smoke.

The [adoption kit](adoption-validation.md) distinguishes researcher/builder/workflow
tasks from the advisory release-persona scenarios. Three non-authors, one
independent analyzer, one integrated workflow and two teams returning after 30 days are targets, not
completed achievements. The 30-day return milestone does not block the first
stable release. No whole-genome speedup, universal hard RAM bound, competitive
caller accuracy or commercial purchase is claimed here.
