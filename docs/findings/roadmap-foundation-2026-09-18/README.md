# Roadmap and onboarding foundation — 2026-09-18

This is maintainer/engineering validation, not independent adoption or public
package publication. The runtime correction is source `9b8e12f3b6e802d2103c5ebeb5ff89a08c8d85ce`;
the separate release-contract inventory change is `48c63f3`. Documentation and
onboarding-harness changes are identified by the commit containing this report.

## Implemented

- One implementation/published/validated/adopted inventory and two README paths.
- Version-aware installation, researcher, reuse, builder, concepts, troubleshooting,
  workflow and historical-context guides. Static landing updated; Receipt Studio
  stays at `/verify/`. The Action keeps its legacy behavior.
- Thirty-nine linked delivery tasks across six milestones. Existing publication
  issues are reused, superseded infrastructure issues closed, and deferred
  platform work labeled. `docs/roadmap.json` is the editable source;
  `scripts/roadmap.py --render --check` validates and renders the backlog.
- Dataset command leaves and both analyzer API modes added to the frozen CLI
  inventory; 49 discovered public entrypoints match release policy.
- Both primary scientific tutorials now run from staged package assets outside
  the source checkout. Selected interpreters, native binaries and temporary
  directories are explicit; native smoke uses a private preparation environment.

## Checks and scope

- All 33 Python tests under `scripts/` passed, including roadmap mapping-conflict
  refusal before GitHub writes, preservation of existing labels/assignees, and
  isolated tutorial launchers with shell-sensitive paths.
- Onboarding staging/link closure passed for 45 guides. Large benchmark archives
  remain remote links and are excluded from bundles.
- Research tutorial: four supplied candidates, record-preserving annotation,
  receipt verification and byte-identical replay.
- Reuse tutorial: a different two-candidate query plus panel QC over 3,159 target
  loci, comparing saved evidence with fresh outputs. The portable artifact is
  relocated and generated tutorial source inputs removed before saved queries.
- Both supplied-Python and native private-environment modes passed against a
  local debug CLI, SHA-256
  `f6c18ea49e3b5b924bb95b6bf6bd8654b11957e7c0b03ecca3bf61783a501643`.
  This checks the changed guides/harness, not a released package or final clean
  candidate identity. Clean installed-runtime checks are retained separately.
- The builder guide's customized `zero_callable_loci` reducer compiled in an
  isolated copy, passed four scientific/lifecycle tests and all 19 evidence
  conformance checks. These include native/saved agreement, budget/tile byte
  invariance and relocated replay. No independent developer is asserted.
- Release inventory: 29 xtask unit tests and five integration tests passed; a
  real help walk matched all 49 declared entrypoints.
- `cargo fmt --all -- --check`, `git diff --check`, shell syntax check for
  `install.sh`, and the non-mutating release-workflow policy checker passed.
- The landing page was opened in a browser and its main layout, audience links,
  availability notice and retained verifier link were inspected.

## Outstanding gates

Additional retained checks bind the updated harness to the clean installed
`9b8e12f` wheel/runtime: [workflow log](packaged-research-workflows.log) and
[identities and command record](packaged-research-workflows.metadata.json).
Both scientific workflows and the primary README example pass; an independent
pysam check confirms candidate counts and full target denominators. The
[cohort fixture check](cohort-fixture-validation.json) separately records the
authored synthetic oracle, strict refusal, saved-only and missing-locus reuse
checks. It establishes no cohort runtime. See the
[fixture guide](../../../examples/cohort-reanalysis/README.md).

The supported Linux/macOS CI matrix must validate the final changed package
smoke. The public stable release remains v0.1.0. A source build or passing local
wheel does not satisfy R04/R06 publication.

On September 18, protected `rc` and `release` environments retained required
reviewers but their custom deployment policies had no allowed refs. No
`CARGO_REGISTRY_TOKEN` name appeared in repository or environment inventories.
TestPyPI/PyPI private publisher state was not established. Account configuration
and environment approvals are maintainer-owned; no settings or credentials were
changed by these inspections.

No independent user sessions, cohort design partners or 30-day return usage have
been supplied. A proposed cohort contract and synthetic fixture prepare C01;
they do not ship cohort commands or satisfy partner validation. Stable release,
cohort APIs and subsequent biological expansion retain their roadmap gates.
