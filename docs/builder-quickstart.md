# Build a useful evidence analyzer

Build a reducer that reports selected loci, callable read observations and supplied
ALT support, then extend it to expose zero-coverage loci. The managed runner handles
resources, cancellation, publication and replay. These APIs require the current
source/unpublished evidence candidate; public stable v0.1.0 predates them.

## Build the source example

Use Rust 1.83+ and the [SDK prerequisites](analyzer-sdk.md#build-prerequisites).
From the Rosalind repository root:

```sh
cargo build --locked --bin rosalind
cargo test --manifest-path examples/evidence-analyzer/Cargo.toml
cargo build --manifest-path examples/evidence-analyzer/Cargo.toml --locked
```

The [independent example](../examples/evidence-analyzer/) already uses local path
dependencies. For a separate project, generate `rosalind new analyzer candidate-qc
--api evidence --output ./candidate-qc` and apply the SDK's
[candidate dependency patches](analyzer-sdk.md#candidate-dependency-setup) before
running Cargo. Generated projects use registry pins, which require those patches
until the matching crates are published. Retain the generated lockfile.

Prepare the [real-origin tutorial](../examples/research-filter/README.md), which
produces `/tmp/research-filter/` with reference, indexed BAM and candidate sites.
Use a fresh tutorial/output directory for repeated runs. From the repository root:

```sh
examples/evidence-analyzer/target/debug/rosalind-example-evidence-analyzer run \
  --alignments /tmp/research-filter/sample.bam \
  --reference /tmp/research-filter/ex1.fa \
  --sites /tmp/research-filter/candidates.vcf --memory-budget-mb 128 --enforce \
  --output /tmp/research-filter/summary.tsv
target/debug/rosalind verify --manifest /tmp/research-filter/summary.tsv.manifest.json
```

With `--scale 1` (the default), the unmodified pinned example emits 4 selected loci,
141 callable read observations and 67 supplied-ALT observations. These are a small
research reducer's counts, not a caller-accuracy result.

## Add a zero-coverage statistic

In your copy of `examples/evidence-analyzer/src/main.rs` (or the generated
`src/main.rs`), extend the existing `CandidateSummary`; keep `run_evidence_artifact`
as lifecycle owner. This statistic reports selected loci with no callable read
observations, which makes coverage gaps visible alongside total read support.

1. Change `counts: [u64; 3]` to `counts: [u64; 4]`, its factory initializer to
   `[0; 4]`, and the overflow test initializer to `[u64::MAX, 0, 0, 0]`.
2. In `SummaryFactory::params()`, change the metric identity from
   `candidate-summary-v1` to `candidate-summary-zero-depth-v1`. Bump your analyzer
   package version when distributing the changed statistic.
3. Immediately after adding `depths.callable_depth` in `on_batch`, add:

```rust
add(&mut self.counts[3], u64::from(depths.callable_depth == 0))?;
```

4. Replace `CandidateSummary::finish` with:

```rust
fn finish(&mut self) -> Result<(), EvidenceError> {
    let values = self.counts.map(|value| {
        value.checked_mul(self.scale).ok_or(EvidenceError::CounterOverflow)
    });
    let [loci, callable, alternates, zero_depth] = values;
    let (loci, callable, alternates, zero_depth) =
        (loci?, callable?, alternates?, zero_depth?);
    writeln!(
        self.out,
        "selected_loci\tcallable_read_observations\tcandidate_alt_read_observations\tzero_callable_loci"
    )?;
    writeln!(self.out, "{loci}\t{callable}\t{alternates}\t{zero_depth}")?;
    Ok(())
}
```

5. Update the existing mixed zero/nonzero fixture's expected output to:

```rust
b"selected_loci\tcallable_read_observations\tcandidate_alt_read_observations\tzero_callable_loci\n2\t10\t4\t1\n"
```

The existing fixture already includes a zero-depth row and multiple ALT bases,
so it checks the new statistic independently of extraction. Keep the missing-field
and overflow assertions. The field requirements remain `depths,alleles`; no new
observations are needed. The extra checked `u64` still fits this example's declared
128-byte reducer allowance; reassess the declaration if you retain more state.
`--scale` remains the example's checked multiplier; use 1 for literal counts.

Rerun the example tests/build above, then check the executable contract:

```sh
target/debug/rosalind conformance analyzer --api evidence \
  --binary examples/evidence-analyzer/target/debug/rosalind-example-evidence-analyzer --json
```

Conformance checks runner behavior and native/saved equality; your fixture checks
what the statistic means. Run the changed binary with a new output filename and
verify that receipt. A stored zero-depth locus is measurable; a locus absent from
a dataset is refused, not counted as zero.

## Run the same reducer on saved evidence

Persist the required evidence once:

```sh
target/debug/rosalind analyze evidence --reference /tmp/research-filter/ex1.fa \
  --alignments /tmp/research-filter/sample.bam --sites /tmp/research-filter/candidates.vcf \
  --fields depths,alleles --cache-dir /tmp/research-filter/evidence-cache \
  --memory-budget-mb 128 --format arrow-ipc --output /tmp/research-filter/saved.arrow
```

The artifact receipt's `measurements.execution.evidence_dataset_manifest` records
the portable manifest location. Use that exact path below; the cache root itself
is not necessarily the dataset directory. Move/copy the whole portable directory
together when relocating it.

```sh
examples/evidence-analyzer/target/debug/rosalind-example-evidence-analyzer run \
  --dataset /path/to/portable/evidence-dataset.manifest.json \
  --sites /tmp/research-filter/candidates.vcf --memory-budget-mb 128 --enforce \
  --output /tmp/research-filter/offline-summary.tsv
```

This path does not open original alignments/reference. It requires complete locus
coverage, depth/allele fields and a real stored reference. Stored filters and sample
scope cannot change during the query. Compare output bytes with a native run of the
same customized binary and selection.

## Verify and replay the custom artifact

```sh
target/debug/rosalind verify --manifest /tmp/research-filter/offline-summary.tsv.manifest.json
target/debug/rosalind reproduce --manifest /tmp/research-filter/offline-summary.tsv.manifest.json \
  --inputs /tmp/research-filter \
  --binary examples/evidence-analyzer/target/debug/rosalind-example-evidence-analyzer --dry-run
```

Keep the portable dataset and selection under the replay input tree; if it lives
elsewhere, choose a tree containing those dependencies. Inspect the dry-run recipe,
then omit `--dry-run` to execute byte replay. An external analyzer always requires
an explicitly selected binary. Receipts bind bytes and provenance; they do not
prove biological accuracy or authorship.

Continue with [SDK lifecycle details](analyzer-sdk.md), [reusable evidence](reusable-evidence.md)
and [workflow integration](workflow-integration.md).
