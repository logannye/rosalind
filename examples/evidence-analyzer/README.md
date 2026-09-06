# Standalone evidence artifact analyzer

This independent Cargo crate imports only public Rosalind APIs. Its
`EvidenceArtifactFactory` creates a reducer borrowing the runner's output sink;
`run_evidence_artifact` owns resource admission, cooperative enforcement,
cancellation, verified inputs, atomic output/receipt publication, and replay.
The reducer requires depth and allele groups and real reference bases, retains
three checked integer counters, and declares a conservative 128-byte bound for
its state and allocation. The runner separately reserves its input, tile, output
and receipt memory.

After [preparing the real-origin tutorial](../research-filter/README.md), run from
the repository root:

```sh
cargo test --manifest-path examples/evidence-analyzer/Cargo.toml
cargo build --manifest-path examples/evidence-analyzer/Cargo.toml --locked
examples/evidence-analyzer/target/debug/rosalind-example-evidence-analyzer run \
  --alignments /tmp/research-filter/sample.bam --reference /tmp/research-filter/ex1.fa \
  --sites /tmp/research-filter/candidates.vcf --memory-budget-mb 128 --enforce \
  --output /tmp/research-filter/summary.tsv
rosalind verify --manifest /tmp/research-filter/summary.tsv.manifest.json
rosalind reproduce --manifest /tmp/research-filter/summary.tsv.manifest.json \
  --inputs /tmp/research-filter \
  --binary examples/evidence-analyzer/target/debug/rosalind-example-evidence-analyzer --dry-run
```

The pinned tutorial reports 4 loci, 141 callable read observations, and 67
observations of supplied ALT alleles. Zero-depth selected loci remain in the locus
count. Multiallelic SNVs count every requested ALT. The summary is a research
reducer; it does not make variant calls. `--scale 1` is an example analyzer parameter
recorded in the receipt and replay recipe; scaling uses checked arithmetic.

The same binary can read a [portable evidence dataset](../persisted-evidence/README.md)
without opening its original BAM or reference:

```sh
examples/evidence-analyzer/target/debug/rosalind-example-evidence-analyzer run \
  --dataset /path/to/bundle/evidence-dataset.manifest.json \
  --sites /tmp/research-filter/candidates.vcf --memory-budget-mb 128 --enforce \
  --output /tmp/research-filter/offline-summary.tsv
```

The dataset must cover every requested locus and contain depth and allele fields
with a real stored reference. Its filter and sample identity is preserved; native
filter/sample flags are rejected with `--dataset`. Omit the selection or use
`--whole-dataset` to summarize all stored loci. `--regions` accepts a BED union.

Receipts from external analyzers always require an explicitly selected `--binary`
for reproduction. `--force` opts into atomic replacement. SIGINT/SIGTERM triggers
cooperative cancellation; failed execution/finalization publishes no successful
artifact/receipt pair. `--require-os-limit` additionally requires an existing
matching Linux cgroup-v2 memory limit.

The source manifest selects this candidate via path dependencies for both Rosalind
and build identity support. After publication, replace those dependencies with
`rosalind-bio = { version = "=0.5.0" }` and `rosalind-build-info = "=0.1.0"` and
retain your own Cargo.lock. Native bundles render those replacements explicitly
and record them in `ONBOARDING-BUNDLE.json`. See the [SDK guide](../../docs/analyzer-sdk.md)
for generated projects and candidate-versus-registry validation.
