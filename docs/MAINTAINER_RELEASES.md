# Programmatic release completion

Rosalind's public binary contains no publishing credentials or release commands.
Maintainers use the source-only `cargo xtask` control plane. Read-only commands
generate deterministic plans; confirmed dispatch commands start workflows whose
mutations remain protected by GitHub environments.

## Bootstrap

```sh
cargo xtask doctor --json
cargo install cargo-public-api --version '=0.50.1' --locked
rustup toolchain install 1.95.0 --profile minimal
rustup toolchain install nightly-2026-04-15 --profile minimal
gh secret set CARGO_REGISTRY_TOKEN --repo logannye/rosalind
```

Create `rc` and `release` repository environments with required reviewers and
deployment rules allowing the intended dispatch branch, normally `main`, before
dispatching. The control plane dispatches workflows from the default branch while
their checkout is pinned to the planned candidate SHA. A custom deployment policy
with no allowed refs blocks every protected job.

Keep the hap.py GHCR package public after its first successful build. The `release`
environment protects crate uploads, image pushes, both Python indexes, and
automated baseline PRs; `rc` protects prerelease tag creation. `doctor` checks names
only and never reads or prints secret values. Its local Docker/tool diagnostics
are distinct from the prerequisites for a build on a GitHub-hosted runner.

## Python trusted publishers

Wheel builds and fresh-install checks run in reusable `wheels.yml` without OIDC
or publication privileges. Final uploads run directly in `rc.yml` and
`release.yml`, because PyPI currently does not support a reusable workflow as a
[trusted publisher](https://docs.pypi.org/trusted-publishers/troubleshooting/#reusable-workflows-on-github).
Keep the protected `release` environment on both upload jobs.

Configure these exact publisher fields in each index's authenticated project
publishing settings, or as a pending publisher if `rosalind-bio` does not yet exist:

| Field | TestPyPI candidate publisher | PyPI stable publisher |
|---|---|---|
| Project name | `rosalind-bio` | `rosalind-bio` |
| Repository owner | `logannye` | `logannye` |
| Repository name | `rosalind` | `rosalind` |
| Workflow filename | `rc.yml` | `release.yml` |
| Environment | `release` | `release` |

The workflow filename is the top-level caller, not `wheels.yml`, and does not
include `.github/workflows/`. Account settings are separate on
[TestPyPI](https://test.pypi.org/manage/account/publishing/) and
[PyPI](https://pypi.org/manage/account/publishing/). GitHub secret-name inspection
and a public project API response cannot establish whether a private pending
publisher is configured. No long-lived PyPI token is used by these workflows.

Before upload, `wheel-artifacts.py` requires exactly the three supported target
reports, verifies the requested commit/native/Python versions, checks platform and
distribution metadata, and verifies each tested wheel's SHA-256 and size. It
rejects missing, extra, changed, or mismatched wheels before staging upload files.
`verify-wheel-upload.py` then permits retrying an existing filename only when the
index reports identical bytes. GitHub release assets retain the candidate reports;
PyPI's upload action retains its trusted-publishing attestation behavior.

Every `plan` is create-new JSON with one common report schema. Its BLAKE3 `plan_id`
authenticates the commit, gate outcomes, contract fingerprint, metadata, and intended
actions. A dispatch rejects a changed body, a blocked plan, or a confirmation typo;
the workflow then recomputes the plan from a clean checkout before doing anything.

Exit codes are stable: `0` ready/success, `2` invalid configuration, `3` an unmet
external gate, `4` a failed remote operation, and `5` an integrity or contract
fingerprint failure.

The v0.4 foundation was checkpointed with 336 existing workspace tests. That count
is recorded in `release/policy.toml`; new maintainer and benchmark tests are additive,
and the original suite remains a non-regression floor.

## Image and benchmark

```sh
cargo xtask giab image plan --output image-plan.json --json
cargo xtask giab image dispatch --plan image-plan.json --confirm PLAN_ID
cargo xtask giab image status --json

cargo xtask giab benchmark plan --output benchmark-plan.json --json
cargo xtask giab benchmark dispatch --plan benchmark-plan.json --confirm PLAN_ID
cargo xtask giab benchmark status --json
```

Routine benchmark runs emit an attested candidate and never edit the committed
baseline. After review, `giab baseline propose` dispatches a PR-producing workflow.

```sh
# Download and verify the newest attested evidence without overwriting a directory.
cargo xtask giab benchmark download --output evidence

# First add a specific "GIAB baseline update" changelog entry, then:
cargo xtask giab baseline propose --run-id RUN_ID --reason "WHY" --json
cargo xtask giab baseline propose --run-id RUN_ID --reason "WHY" --confirm PLAN_ID
```

The committed evaluator lock is the default. A manual benchmark image override must
be `NAME@sha256:DIGEST`; both the OCI platform and GitHub attestation are checked
before evaluation. Downloads/prepared inputs are cached independently from results,
and every restored input is re-hashed by `prepare.sh`.

## RC and stable release

```sh
cargo xtask rc plan --version 0.5.0 --number 1 --ref FULL_CANDIDATE_SHA \
  --output /tmp/rosalind-rc-plan.json --json
cargo xtask rc dispatch --plan /tmp/rosalind-rc-plan.json --confirm PLAN_ID
cargo xtask rc status --tag v0.5.0-rc.1 --json

cargo xtask release plan --version 0.5.0 --rc-tag v0.5.0-rc.1 \
  --ref FULL_CANDIDATE_SHA --output /tmp/rosalind-release-plan.json --json
cargo xtask release dispatch --plan /tmp/rosalind-release-plan.json --confirm PLAN_ID
```

The stable workflow enforces an identical public-contract fingerprint, idempotent
publication, and fresh-cache downstream installation before creating the stable
tag. The requested version must match the committed root package and release
policy, and the clean candidate must be reachable from the pushed default branch.
Each `PLAN_ID` is the exact ID returned by the corresponding plan. Plans are kept
outside the checkout so they do not make the release tree dirty.

The feature-bearing evidence-engine candidate uses the 0.5.0 release sequence.
The server-timestamped seven-day soak and three partner personas apply at 0.5.0
and later. The soak begins when the GitHub prerelease is published, after its
TestPyPI wheels and OCI image pass their protected jobs. The historical 0.4.0
stabilization exemption does not justify relabeling the evidence-engine release.
Changed caller/shared-processing source additionally requires the reviewed GIAB
baseline described below before stable promotion. Missing crate credentials or
that baseline do not prevent preparing an RC that is clearly labeled experimental.

Publication can be rerun after interruption. Each package is recreated locally and
the downloaded registry `.crate` must have identical SHA-256 bytes before it is
skipped. A mismatch is permanent integrity failure. The stable tag and GitHub
Release are created only after Linux and macOS install the published crate into an
empty `CARGO_HOME`, run the demo, build a scaffold, and pass conformance.

## Design-partner gate

Generate private interview packets locally:

```sh
cargo xtask partners init --persona analyzer-builder --output release/private-design-partners/analyzer
cargo xtask partners init --persona workflow-hpc --output release/private-design-partners/hpc
cargo xtask partners init --persona constrained-offline --output release/private-design-partners/offline
```

Only copy the finalized anonymized `feedback.json` records into
`release/design-partners/`. Validation requires exact persona scenarios, a full
tested commit, the frozen contract fingerprint, consent to publish the sanitized
record, and no unresolved release-blocking finding. At promotion time the tested
commit must be an ancestor of the candidate and its fingerprint must equal the RC.
Names, contact details, organizations, credentials, genomic data, and raw notes are
rejected or kept outside the repository.

## What remains intentionally human

Automation narrows but does not counterfeit external evidence. A maintainer must
provide the crates.io token, approve protected environments, make the GHCR package
public, provide benchmark compute, wait the full server-timestamped seven days, and
conduct three real partner engagements. The CLI reports those as blockers rather
than silently weakening them.

## Candidate publication and evaluator readiness

RC/stable workflows invoke wheel and OCI publication explicitly with the exact
candidate SHA; they do not depend on a downstream release event from GITHUB_TOKEN.
Linux x86_64 and macOS arm64/x86_64 wheels run outside-checkout install checks on
Python 3.11 and 3.9, including exact evidence, panel QC, and artifact verification.
The installed binary version must equal the Python distribution version after RC
normalization. Persisted evidence is replayed through that installed binary with
local inputs. The packaged scaffold is generated, built offline against the
candidate SDK source, and checked by the packaged conformance runner. This source
patch is explicit pre-publication validation; fresh registry SDK installation is
still checked separately after crate publication. The candidate smoke explicitly
prepares its dependency cache before the locked offline scaffold build. RC build
checkouts derive numbered native and Python versions, for example `0.5.0-rc.N` and
`0.5.0rcN`, from the committed stable base version. Derivation records source and
manifest hashes; it does not change the committed candidate. Generated build
reports live outside the source checkout before compilation, then wheel hashes
are captured after testing. Reusable builds never publish by themselves.

Caller evidence is change-based: the selected variants dispatch/defaults,
transitive local helpers, watched caller/shared-input files, and dependency lock
are compared to the tested source. A root package version bump or unrelated CLI
addition alone does not require GIAB. Shared pileup/filter changes do. Missing,
stale, dirty, or incomplete baseline evidence fails closed. An unchanged bootstrap
caller remains explicitly experimental with no accuracy claim.

`happy-candidate.yml` builds the evaluator without registry credentials and runs a
synthetic identical-SNV case through the actual offline hap.py/vcfeval path. It
runs on evaluator PRs and when a scheduled GIAB run lacks a published evaluator.
The protected image workflow requires that candidate gate, then publishes and
attests a digest and proposes its lock. A candidate build is not a GIAB result;
HG002 execution still requires the reviewed published image and prepared data.
Python 2.7 is required by pinned hap.py 0.3.15; compatible bx-python/six wheels are
hash-locked rather than resolved from floating dependencies.

Read-only prerequisite inspection on 2026-09-06 found no repository, rc, or release
secrets: CARGO_REGISTRY_TOKEN is absent. Both protected environments have the owner
as required reviewer, but custom deployment policies contain zero allowed
branch/tag rules. Configure the intended release workflow refs before dispatch.
PyPI/TestPyPI trusted publishers must be configured with the top-level caller and
protected environment listed above; GitHub secret-name inspection cannot establish
that external state. The public PyPI/TestPyPI project APIs returned 404, which does
not rule out a pending publisher. GitHub's setting for Actions to create/approve
pull requests was disabled, so automatic evaluator-lock/baseline PR creation also
needs explicit maintainer handling. No secret values were read and no settings
were changed by that inspection.
