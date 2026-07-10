# Current architecture

Rosalind separates the bounded computation kernel from the products built on it:

1. `IndexReader` memory-maps a deterministic multi-contig reference index.
2. `StreamingBamSource` admits a coordinate-sorted BAM one record at a time.
3. The pileup engine emits deterministic, depth-capped `PileupColumn` values.
4. A `ColumnAnalyzer` converts each column into a streaming artifact.
5. `contract::run_column_analysis` owns prediction, upfront refusal, live RSS
   enforcement, output lifetime, measurement capture, and receipt sealing.
6. The CLI is an adapter: `features` and `analyze` delegate to the public runner and
   map typed refusal/breach outcomes to exits 3/4. The library never exits a process.

Producer and analyzer identities live in `ContractRunSpec`, keeping the existing
`ColumnAnalyzer` trait source-compatible. Replay is tokenized argv, and receipts
bind both identities into the deterministic claim. Downstream binaries can therefore
inherit the contract without forking Rosalind's CLI.

The RSS governor is process-wide. Only one enforced contract runner may be active in
a process; a concurrent attempt returns a typed configuration error. Record-only
runs do not make a concurrency guarantee and callers should serialize them when
process-wide RSS attribution matters.

The browser-side Receipt Studio consumes the leaf `rosalind-receipt` crate compiled
to WASM. It performs receipt inspection, streaming artifact hashing, causal diff,
and chain traversal locally, with no upload or third-party network request.
