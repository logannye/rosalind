# Rosalind architecture

Current source turns indexed short-read alignments into exact per-locus summaries,
persists those summaries, and supplies the same evidence to bounded reducers.
Scientific selection and filtering are separate from execution tiling. These APIs
are an unpublished candidate; public stable v0.1.0 predates them.
[Implementation status](docs/implementation-status.md) separates implementation,
validation, publication and adoption.

```text
indexed BAM / supported CRAM + local reference + SNV sites / BED
                  |
       sample scope + scientific profile + fields
                  |
       admitted aggregate microtiles --> canonical evidence batches
                  |                              |
       verified partition datasets --------------+
                                                 |
                              EvidenceAnalyzer / panel QC / encoders
                                                 |
                       managed admission + cancellation + publication
                                                 |
                             artifact + receipt + verification/replay
```

## Evidence layers

| Layer | Code | Responsibility |
|---|---|---|
| Scientific request | `src/evidence/mod.rs`, `sample.rs`, `selection.rs` | Loci, sample scope, filters, fields and consumer requirements |
| Reference access | `src/evidence/reference.rs`, `src/genomics/reference_pack.rs` | Indexed FASTA, `.rref` and compatible `.idx` windows |
| Decoder admission | `src/evidence/cram.rs`, `engine.rs` | Indexed input, CRAM metadata envelopes and whole-file record validation |
| Exact aggregation | `src/evidence/engine.rs`, `batch.rs` | Budget-selected microtiles, checked summaries and physical field projection |
| Encoding and consumers | `src/evidence/encoding.rs`, `panel.rs` | Canonical Arrow/TSV, full target denominators and fused consumers |
| Managed SDK | `src/evidence/artifact.rs` | Native/persisted sources, admission, cancellation, atomic artifacts, receipts and replay |
| Saved evidence | `src/dataset.rs`, `src/dataset/` | Partition workers, cache/resume, portable datasets, offline queries, partial reuse and Parquet |
| Trust and replay | `crates/receipt/`, `src/reproduce.rs`, `src/receipt_tools.rs` | Claims, byte verification, tokenized replay and inspection |
| Adapters | `src/evidence_cli.rs`, `src/dataset_cli.rs`, `python/rosalind/` | CLI/Python interfaces over the same scientific semantics |

BAM accepts BAI or CSI indexes; the supported CRAM profile requires CRAI and an
explicit local FASTA with adjacent FAI. CRAM performs complete cooperative record
validation before indexed extraction.

Canonical ownership tiles span 16,384 bases. Budget-selected microtiles can be
smaller and revisit indexed reads. Managed/encoded batches contain at most 1,024
rows independently of computation width. Workers produce first-party partitions;
reducers consume them serially in canonical order. Eligible observations are never
discarded to fit a successful exact-evidence run.

Saved datasets retain profile, sample scope, selection, fields and source identity.
Queries can project available fields, select stored loci and request new SNV ALT
annotations from A/C/G/T counts. They cannot reconstruct filtered observations or
invent absent loci. Partial-overlap reuse checks unchanged original inputs and
computes missing loci. Its current one-worker path emits a materialized artifact,
not an expanded portable cache.

## Invariants and boundaries

1. Loci, filters, sample scope and counting unit determine scientific meaning.
   Admitted budgets, microtiles and workers do not change successful evidence bytes.
2. Alignment/reference dictionaries agree before successful output creation.
   Analysis does not require building an FM-index.
3. Consumers declare capabilities and retained memory. Unknown bounds cannot
   support enforcement. Absent groups do not mean zero.
4. Native/persisted sources feed the same reducer contract. Query lineage represents
   only consumed partitions as locally verified.
5. Refusal, failure and cancellation do not publish a successful artifact pair.
6. Scientific identity, full receipt claim and measurements are distinct. Receipts
   are tamper-evident records, not signatures or biological validation.
7. Inputs remain immutable. Content hashes establish identity; metadata/inode
   guards detect ordinary subsequent changes.
8. Checks are cooperative. Native decoding, including CRAM validation, can allocate
   before a checkpoint. A Linux cgroup supplies separately recorded OS assurance.

See [semantics](docs/SEMANTICS.md), [the contract](CONTRACT.md),
[reusable evidence](docs/reusable-evidence.md) and [the SDK](docs/analyzer-sdk.md).

## Legacy analysis and search

```text
.rref / compatible .idx + coordinate-sorted BAM
                  |
         exact-or-fail pileup columns --> ColumnAnalyzer
                  |
  whole-genome / selected interval / reference-span shard
                  |
       contract runner --> artifact + receipt
       complete shards --> strict canonical merge
```

`src/pileup/`, `src/call/columnkit.rs`, `src/contract.rs`, `src/selection.rs` and
`src/merge.rs` maintain this compatibility path. Current runs fail at declared
capacities rather than downsampling; their defaults differ from the evidence
profile. `src/genomics/index/` retains FM-index search as an adjacent workload.

Historical receipt and legacy `.idx` verification remain supported; replaying old
behavior requires its matching producer. Serialize managed runs within a process
or use separate processes because the governor/cancellation scope is process-wide.
Evidence partition workers and legacy reference-span shards are distinct parallel
execution paths.
