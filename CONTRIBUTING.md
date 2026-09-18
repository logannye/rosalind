# Contributing to Rosalind

Start with [architecture](ARCHITECTURE.md), [semantics](docs/SEMANTICS.md) and
[the resource contract](CONTRACT.md). For a custom statistic, use the
[builder quickstart](docs/builder-quickstart.md). Current evidence source/candidate
APIs are unpublished; public stable v0.1.0 predates them.

## Build and validate

```sh
cargo build --locked --bin rosalind
cargo test --workspace
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
python3 scripts/onboarding.py check .
```

Use Rust 1.83 or newer and the [native build prerequisites](docs/analyzer-sdk.md#build-prerequisites).
Run focused checks during development, then affected integration/platform gates.
Recorded results must identify the tested source/build. A source smoke does not
establish registry installation, representative performance or external adoption.

## Preserve the contracts

1. **Scientific semantics.** Budgets, microtiles and workers cannot alter successful
   exact-evidence values/bytes. Preserve zero-depth denominators, exclusive filters,
   sample scope and absent-field behavior. Define new metrics and their identity.
2. **Determinism.** Use checked integer reduction and canonical ordering. Preserve
   relevant repeat-run, budget/tile/worker and saved/native equality. Historical
   fixtures retain their documented semantics. See [determinism](docs/determinism.md).
3. **Resources and publication.** Account for decoder, metadata, consumer, encoder
   and finalization state. Models can fail; cooperative checks do not prove every
   allocation fits. CRAM validation can allocate before a checkpoint. Failure and
   cancellation must never publish success.

## Extend the evidence SDK

Implement `EvidenceArtifactFactory` and a bounded `EvidenceAnalyzer`, then call
`run_evidence_artifact`. Native and persisted sources share its reducer interface;
the runner owns the artifact lifecycle. Declare fields and a conservative retained
bound, and record scientific parameters in both receipt metadata and replay.
See [the SDK guide](docs/analyzer-sdk.md) and [standalone example](examples/evidence-analyzer/).

`ColumnAnalyzer` remains a legacy API. Its lower-level `run_bounded_whole_genome`
helper does not create receipts or govern arbitrary consumer allocations. A
streaming iterator alone does not confer those guarantees on new code.

## Review and bug reports

For non-trivial changes, record behavior, scientific assumptions, failure modes and
acceptance checks before implementation. Existing specs/plans live under
`docs/superpowers/`. Keep changes reviewable and require relevant CI before merge.

Bug reports should include the command, `rosalind --version`, platform, claim hash
(`manifest_blake3`), error output and a shareable minimal fixture. Identify native
versus persisted input, budget, CRAM validation, cache reuse and worker settings.
Do not attach sensitive sequencing data to a public issue.
