# Engineering checkpoint — 2026-09-18

Work is paused at a reviewable checkpoint at the maintainer's request to conserve
compute. This is not a claim that stable publication or the entire adoption
roadmap is complete.

## What is ready

- Main's README, landing page, researcher/builder guides, executable onboarding,
  capability inventory, and linked issue roadmap are updated. Main is the 0.5
  source preview; the latest published stable release remains v0.1.0.
- [PR #145](https://github.com/logannye/rosalind/pull/145) implements C01–C11's
  engineering workflow on the separate 0.6 branch: immutable local cohorts,
  comparison checks, strict/partial missingness, exact summaries, bounded managed
  execution, explicit extension, CLI/Python, and relocated replay. The
  [retained validation](findings/cohort-preview-2026-09-18/README.md) identifies
  tested source commits, authored inputs, successes, failures and limitations.
- [PR #146](https://github.com/logannye/rosalind/pull/146) fixes the missing CMake
  container prerequisite. Its actual Docker build and offline non-root smoke
  passed, as did its primary CI suite.
- [PR #147](https://github.com/logannye/rosalind/pull/147) fixes stale build identity
  in linked worktrees. Focused tests exercise real incremental Cargo rebuilds;
  local Clippy and Rust 1.83 checks passed. Remaining remote checks are visible
  on the PR.
- The evaluator image was successfully published, smoked and attested in
  [run 35378796613](https://github.com/logannye/rosalind/actions/runs/35378796613).
  Only automatic PR creation failed. Recover its generated lock through an
  ordinary PR; do not rebuild or republish the image merely to recover that PR.

The cohort platform CI at `ba6cd43` passed its ten jobs, including actual Linux
cgroup extraction, summary and refusal checks. The subsequent `47baf5f` change
adds retained evidence and packaged tutorial assets. Its redundant full CI and
wheel runs (35385938344 and 35385938360) were intentionally canceled to conserve
compute; local link, packaging and staged-example checks passed. Cancellation is
not recorded as a successful CI result.

## Resume release closure without repeating finished work

1. Inspect the existing checks on #146/#147; merge when their applicable checks
   are satisfied. Keep 0.6 cohort code outside the 0.5 candidate.
2. Review/import the already-published evaluator digest and retained attestation.
   Resolve or document the manual-PR handoff while preserving environment rules.
3. Reuse the maintainer control plane to prepare a fresh immutable **0.5.0-rc.3**
   from the resulting main commit. RC2 was canceled after its container failure;
   do not assign its number to different source. Run the existing publication
   workflow once, review its protected jobs, and verify the resulting public
   bundles, wheels and container outside the checkout.
4. Run the existing pinned HG002 chr20 GIAB benchmark using the reviewed image;
   retain all results and review the baseline required by watched caller/shared
   source changes. The synthetic evaluator smoke is not this benchmark.
5. Complete stable account configuration: `CARGO_REGISTRY_TOKEN` in GitHub and
   the PyPI trusted publisher for `release.yml`, environment `release`. The
   maintainer confirmed TestPyPI setup; actual OIDC upload remains the proof.
6. Observe the existing seven-day soak **from actual GitHub prerelease
   publication**, then use the stable planner and ordered registry publication.
   Only after success update the default installation/availability labels and
   record the release-bound demonstration.

No RC or seven-day soak is established by a successful local build. Do not
weaken correctness, artifact identity or technical release checks to label work
complete. Independent-user feedback is advisory and does not block shipment.

## Preview and later value work

Keep the cohort PR reviewable independently of release closure. Paired/longitudinal
comparison work is checkpointed separately on `codex/cohort-pairs-preview`;
its exact arithmetic and bounded core tests passed, but its adapter/demo/platform
validation must be reported separately before promotion.

Independent use, two recurring cohort teams, 30-day return usage, measured
real-cohort economic value, the next partner-selected reducer, and a molecule
feasibility study remain unestablished. Preserve these as adoption/demand
milestones rather than fabricating evidence or inventing shipping prerequisites.
