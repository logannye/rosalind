# Rosalind architecture

Rosalind separates analysis references, the bounded computation kernel, analyzer
logic, and the execution contract. Search indexing is an adjacent workload, not a
dependency of per-locus analysis.

```text
FASTA ── reference build ──> .rref ──┐
legacy .idx ── compatibility ────────┤
                                     ├─> bounded pileup walk ─> ColumnAnalyzer
coordinate-sorted BAM ───────────────┘              │
                                                    v
plan/refuse ── governor ── transactional output ── receipt/replay/diff

BAM + BAI ── normalized region/BED or reference-span shard ──┐
                                                              ├─> same bounded kernel
complete shard receipts ── strict canonical merge ───────────> merged artifact + receipt
```

## Layers

| Layer | Code | Responsibility |
|---|---|---|
| Coordinates and budgets | `src/core/` | Contigs, loci, reads, memory budgets, process-wide governor, typed errors |
| Analysis reference | `src/genomics/reference_pack.rs` | `.rref`, streaming builder, mmap reader, `ReferenceSequence` and `ReferenceProvider` |
| Legacy search index | `src/genomics/index/` | `.idx`, FM-index lookup, compatible embedded reference provider |
| Selection | `src/selection.rs` | Region/BED normalization and `reference-span-v1` ownership |
| Alignment input | `src/io/bam.rs` | Sequential BAM stream or bounded indexed interval fetch |
| Bounded kernel | `src/pileup/` | Filtered/depth-capped CIGAR-aware pileup columns |
| Analyzer SDK | `src/call/columnkit.rs` | `ColumnAnalyzer` and maintained feature/coverage analyzers |
| Execution contract | `src/contract.rs` | Prediction, refusal, enforcement, output lifetime, measurements, receipt sealing |
| Canonical merge | `src/merge.rs` | Parent evidence validation and first-party TSV/Arrow/VCF/gVCF merge codecs |
| Trust and replay | `crates/receipt/`, `src/reproduce.rs`, `src/receipt_tools.rs` | Canonical claims, verification, replay, certificates, diff and inspection |
| CLI adapter | `src/main.rs` | User-facing argument validation and stable exit-code mapping |

## Invariants

1. An analyzer cannot require an FM-index merely to read reference bases.
2. Reference and BAM contig dictionaries must agree before output creation.
3. Bounded analyzers stream columns and declare any additional memory bound.
4. Enforced refusal happens before the destination exists.
5. Successful artifacts publish atomically; breaches never take the successful name.
6. Identical inputs, parameters, analyzer, and Rosalind release produce identical bytes.
7. Output-affecting choices are claim-protected; host measurements remain separate.
8. Replay is tokenized and allowlisted. A receipt is tamper-evident, not proof of authorship.

## Compatibility

`AnalysisReference` opens `.rref` or legacy `.idx` by magic. `ReferenceSequence`
keeps existing borrowed `ReferenceView` callers source-compatible, while
`ReferenceProvider` adds contig metadata and source identity for complete analysis
artifacts. Receipt schemas 1–5 and legacy `.idx` replay remain supported.

The RSS governor is process-wide; enforced runs must be serialized within a process.
Deterministic sharding, rather than nondeterministic internal threading, is the
parallel execution model. Indexed reads may contribute on both sides of a shard
boundary, but each output locus is owned exactly once by its reference span.
