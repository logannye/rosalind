# Programmatic release completion

Rosalind's public binary contains no publishing credentials or release commands.
Maintainers use the source-only `cargo xtask` control plane. Read-only commands
generate deterministic plans; confirmed dispatch commands start workflows whose
mutations remain protected by GitHub environments.

## Bootstrap

```sh
cargo xtask doctor --json
cargo install cargo-public-api --version 0.50.1
rustup toolchain install 1.95.0 --profile minimal
rustup toolchain install nightly-2026-04-15 --profile minimal
gh secret set CARGO_REGISTRY_TOKEN --repo logannye/rosalind
```

Create `rc` and `release` repository environments with required reviewers before
dispatching. Keep the hap.py GHCR package public after its first successful build.
The `release` environment protects crate uploads, image pushes, and automated
baseline PRs; `rc` protects prerelease tag creation. `doctor` checks names only and
never reads or prints secret values.

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
cargo xtask rc plan --version 0.4.0 --number 1 --output rc-plan.json --json
cargo xtask rc dispatch --plan rc-plan.json --confirm PLAN_ID
cargo xtask rc status --tag v0.4.0-rc.1 --json

cargo xtask release plan --version 0.4.0 --rc-tag v0.4.0-rc.1 \
  --output release-plan.json --json
cargo xtask release dispatch --plan release-plan.json --confirm PLAN_ID
```

The stable workflow enforces the server-timestamped seven-day soak, an identical
public-contract fingerprint, three anonymized partner personas, idempotent crates.io
publication, and fresh-cache downstream installation before creating the stable tag.

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
