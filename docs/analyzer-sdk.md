# Analyzer SDK

Use `rosalind::evidence` for new exact short-read evidence consumers. Open an
`EvidenceRequest`, resolve its VCF/BED selection against the engine's contigs, and
pass an `EvidenceAnalyzer` to `run`. The independent
[example crate](../examples/evidence-analyzer/) demonstrates the full public path.

An analyzer declares `EvidenceRequirements`: field capabilities, whether real
reference bases are required, flanking context, and a conservative peak retained
byte bound. Call `plan_for_analyzer` before exposing the final plan or creating an
output. Unknown memory is allowed without a declared budget and rejected for
budgeted runs. Schema v1 supports zero flanking context and physically emits all
fields; a smaller consumer capability set is not physical projection.

`on_batch` borrows ordered exact rows. Copies retained after a callback belong to
the analyzer's memory model. Use checked integer reducers and deterministic order;
do not retain a genome table, timestamps, random identifiers, or unordered output.
`PanelQcAnalyzer` keeps state proportional to original target count and declares
that memory. `FusedAnalyzers` combines two consumers over one traversal.

The engine supplies exact indexed extraction and resource planning. It does not
automatically supply every standalone binary with transactional files, receipts,
or OS enforcement. The CLI and Python `materialize_evidence` provide the first-party
artifact lifecycle. Custom binaries must explicitly add the lifecycle they claim.
See [SEMANTICS.md](SEMANTICS.md) for depth/filter/coordinate rules.

## Legacy ColumnAnalyzer compatibility

`ColumnAnalyzer` remains for existing pileup-column consumers and the scaffold:

```sh
rosalind new analyzer locus-qc --output ./locus-qc
cd locus-qc
cargo test
rosalind conformance analyzer --binary ./target/release/locus-qc --json
```

`run_column_analysis(analyzer, spec)` owns reference/BAM validation, prediction,
refusal, governance, atomic output, receipt sealing, and replay.
`run_column_analysis_selected` adds the legacy normalized interval/shard path.
The lower-level `run_bounded_whole_genome` only drives the column stream and writes
to its supplied sink; it does not itself create receipts or govern total process
memory. An analyzer cannot inherit a bound for arbitrary additional retained state.

Use `AnalyzerMemoryModel::Fixed` only for a real maximum including encoder flush
transients. Unknown bounds are valid for observation but refused under enforcement.
Distinguish `ContractRunError::Refused` from `::Breached`; neither publishes a
successful artifact. New pileup runs fail at capacity rather than sampling reads.
The adapter exposes base/quality/strand/read-position observations; it does not
provide UMI groups, modification tags, methylation calls, phased haplotypes, or
other unsupported biological capabilities.

The scaffold includes repeat bytes, receipt integrity, relocated inputs, refusal,
and explicit external replay tests. It remains a compatibility option while the
new evidence API and its release/adoption gates mature.
