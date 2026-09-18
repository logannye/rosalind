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

The public v0.1.0 macOS/arm64 installer was also executed from a fresh directory.
Its binary predates `--version`: the initial new completion probe failed after
a successful checksum-verified download. The probe now handles that historical
interface, and the fresh retry exits successfully with explicit legacy/source
guidance. Both the [failed attempt](public-stable-installer-before.log) and
[corrected result](public-stable-installer.log) are retained. A regression check
exercises historical and current CLI interfaces without network downloads.

Additional retained checks bind the updated harness to the clean installed
`9b8e12f` wheel/runtime: [workflow log](packaged-research-workflows.log) and
[identities and command record](packaged-research-workflows.metadata.json).
Both scientific workflows and the primary README example pass; an independent
pysam check confirms candidate counts and full target denominators. The
[cohort fixture check](cohort-fixture-validation.json) separately records the
authored synthetic oracle, strict refusal, saved-only and missing-locus reuse
checks. It establishes no cohort runtime. See the
[fixture guide](../../../examples/cohort-reanalysis/README.md).

The supported Linux/macOS wheel matrix must validate the final changed package
smoke. Actual native tarball staging and outside-checkout smoke on all supported
targets remain a separate D07 acceptance gate. Those checks already run in the
RC build stage before publication; D07 does not prevent dispatching that validation
stage once its account and source prerequisites are ready. The public stable
release remains v0.1.0. A source build or passing local wheel does not satisfy
R04/R06 publication.

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

## Merged delivery and acceptance review

- PRs [#104](https://github.com/logannye/rosalind/pull/104) and
  [#105](https://github.com/logannye/rosalind/pull/105) are merged. The corrected
  runtime head `9b8e12f` passed all 15 checks, including the
  [three-platform wheel matrix](https://github.com/logannye/rosalind/actions/runs/35369668546).
  The #105 merge is `6dc920b9706348e319273fa2db1532360a994f3f`.
- The documentation/onboarding implementation
  [#143](https://github.com/logannye/rosalind/pull/143), head `bc6e9e6`, passed all
  15 checks, including Python 3.9 and 3.11 wheel smoke on Linux x86_64 and both
  supported macOS architectures. Its merge is
  `f6d0cbb6d075e96f1a4bae81eb0600aef969bc21`.
- The installer and Python contract follow-up
  [#144](https://github.com/logannye/rosalind/pull/144) is merged as
  `ad7675d83d68849adbcbe854bf5d1df3023351b4`. Its source, workflow and Linux/macOS
  ARM wheel checks passed before merge; its separate Intel wheel check was still
  running at merge time. The scientific runtime and Python package source are
  unchanged from the fully checked onboarding head. The
  [complete contract comparison](../python-contract-2026-09-18/README.md)
  records the added Python fingerprint and unchanged existing scientific components.
- [Pages deployment](https://github.com/logannye/rosalind/actions/runs/35375394358)
  succeeded from `f6d0cbb6`. The public landing page was opened in a browser:
  audience routes, source installation and unpublished-preview notice were present.
  The preserved `/verify/` route loaded its sample chain and reported intact
  receipt hashes and lineage. No user files were uploaded.
- Acceptance review closes G01–G07, D01–D03, D05–D06, R01 and R03. The May target
  architecture now has an explicit archive banner. D04 remains open for independent
  developer completion; D07 remains open for the actual native-bundle matrix.
  Maintainer-authored fixtures, successful CI and these engineering checks do not
  count as independent adoption.

The [linked roadmap](../../ROADMAP.md) records all remaining work and evidence.
The next publication action is maintainer-owned R02 setup, followed by a fresh RC
plan and its build/validation stage. This is not a declaration that 0.5 is released
or that cohort APIs are implemented.
