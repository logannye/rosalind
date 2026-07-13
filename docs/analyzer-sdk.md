# Analyzer SDK

Implement `ColumnAnalyzer` when the output can be derived independently from each
bounded pileup column. The contract runner owns reference/BAM validation, prediction,
refusal, governance, transactional output, receipt sealing, and replay.

Start with:

```sh
rosalind new analyzer locus-qc --output ./locus-qc
cd locus-qc
cargo test
```

An analyzer provides an optional deterministic header, claim-protected parameters,
and `on_column`. State must be O(1) or have a conservative declared bound. Avoid
timestamps, randomized identifiers, host metadata, unordered-map iteration, locale-
dependent formatting, and unbounded genome-wide accumulation.

Use `AnalyzerMemoryModel::Fixed` only for a real conservative maximum. Unknown
bounds are valid in record-only mode and refused under enforcement. Treat
`ContractRunError::Refused` and `::Breached` as distinct outcomes and never publish
a partial artifact under the successful name.

`run_column_analysis(analyzer, spec)` remains the whole-genome, sequential-BAM
entry point. Use `run_column_analysis_selected(analyzer, spec, selection)` with an
`IntervalSet` or deterministic `AnalysisSelection::shard` for indexed BAM+BAI
execution. Selection normalization, partition claims, replay arguments, and memory
prediction are inherited from the same transactional runner.

Conformance is available without a Rosalind source checkout:

```sh
rosalind conformance analyzer --binary ./target/release/locus-qc --json
```

The scaffold includes repeat-byte, receipt-integrity, relocated-path, refusal, and
external-replay tests. See `src/call/columnkit.rs` and `src/contract.rs` for the
complete API.
