# Contributing to Rosalind

Thanks for building on Rosalind. This guide gets you productive fast and explains the two
guarantees every change must preserve.

## Build, test, format

```bash
cargo build --release          # the rosalind binary + library
cargo test                     # full unit + integration suite
cargo fmt --check              # formatting (run `cargo fmt` to fix)
```

## Two invariants you must not break

1. **Determinism.** Primary artifacts (VCF, gVCF, feature TSV) are byte-identical across repeated runs
   on identical inputs. See [`docs/determinism.md`](docs/determinism.md). If you touch the calling or
   egress path, add or keep a repeat-run byte-equality test.
2. **The memory contract.** `rosalind plan` predicts a peak that must upper-bound the realized peak;
   `--enforce` honors it (refuse exit 3, or the runtime governor exits 4); the receipt is self-hashing.
   See [`CONTRACT.md`](CONTRACT.md). A streaming stage's working set must stay bounded by coverage, not
   input size.

## How we work: brainstorm → spec → plan → TDD

Non-trivial changes go through a short design loop. Specs live in `docs/superpowers/specs/`, plans in
`docs/superpowers/plans/`. Write a failing test first; keep commits small; one branch per increment;
open a PR for review (CI must be green before merge).

## Extending the engine

The fastest way to add your own per-locus analytic is the **ColumnKit SDK**: implement one
`ColumnAnalyzer` method and run it through `run_bounded_whole_genome` to inherit bounded memory,
determinism, and a verifiable receipt — no contract code to re-derive. See
[`examples/columnkit_coverage.rs`](examples/columnkit_coverage.rs) and [`CONTRACT.md`](CONTRACT.md).

## Reporting bugs

Use the bug-report template. Include the exact command, the inputs, and the receipt hash
(`manifest_blake3`) + `rosalind --version` — Rosalind's determinism makes most bugs exactly
reproducible.
