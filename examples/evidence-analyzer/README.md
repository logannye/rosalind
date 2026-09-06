# Standalone evidence analyzer

This is an independent Cargo crate with its own workspace boundary. It imports
only public Rosalind APIs, declares the fields/reference it needs and its 24-byte
retained state, and consumes exact batches with checked integer sums. No read or
genome-wide table is retained by the analyzer. The engine also reserves its own
input, tile, and output memory;24bytes is not the whole process budget.

After [preparing the real-origin tutorial](../research-filter/README.md):

```sh
cargo run --manifest-path examples/evidence-analyzer/Cargo.toml -- \
  /tmp/research-filter/sample.bam /tmp/research-filter/ex1.fa \
  /tmp/research-filter/candidates.vcf
cargo test --manifest-path examples/evidence-analyzer/Cargo.toml
```

The local dependency path selects this unreleased source tree. Once a matching
release exists, an external repository can replace it with that registry version
and retain its own Cargo.lock. The current output is a three-integer research
summary, not a variant caller. On the pinned tutorial it reports 4 loci, 141 callable
read observations, and 67 observations of supplied ALT alleles.

`EvidenceEngine` supplies exact extraction and admission planning. This minimal
binary does not add transactional file receipts or OS enforcement automatically.
For persistent verifiable datasets use `rosalind analyze evidence` or Python
`materialize_evidence`; for the complete legacy analyzer receipt/replay lifecycle
use the `rosalind new analyzer` scaffold. See [SDK contracts](../../docs/analyzer-sdk.md).
