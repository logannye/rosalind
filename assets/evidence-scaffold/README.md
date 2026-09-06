# __PACKAGE_NAME__

A standalone batch evidence analyzer that inherits Rosalind's resource planning,
cooperative enforcement, cancellation, verified inputs, atomic artifact/receipt
publication, and explicit replay. Customize `SummaryFactory` and `CandidateSummary`
in `src/main.rs`; retain `run_evidence_artifact` as the lifecycle owner.

The example counts selected loci (including zero coverage), callable reads and
reads supporting requested SNV ALT alleles. It requires depth and allele groups
and a real reference. The default `--scale 1` is recorded as an analyzer parameter;
all counters and scaling use checked integer arithmetic. The analyzer declares
128 retained bytes, including its borrowed writer and allocation allowance. The
runner reserves its own encoders, input buffers and receipt state separately.

```sh
cargo test
cargo build --release --locked
# Native input: BAM/CRAM plus index; reference FASTA+.fai, .rref or .idx.
target/release/__PACKAGE_NAME__ run \
  --alignments sample.bam --reference genome.fa --sites candidates.vcf \
  --memory-budget-mb 128 --enforce --output summary.tsv
rosalind verify --manifest summary.tsv.manifest.json
# Persisted input: the portable bundle supplies reference/profile/sample identity.
target/release/__PACKAGE_NAME__ run \
  --dataset cache/evidence-dataset.manifest.json --sites candidates.vcf \
  --memory-budget-mb 128 --enforce --output offline-summary.tsv
# External receipts always require an explicitly chosen executable for replay.
rosalind reproduce --manifest summary.tsv.manifest.json --inputs input-directory \
  --binary target/release/__PACKAGE_NAME__ --dry-run
rosalind conformance analyzer --api evidence --binary target/release/__PACKAGE_NAME__ --json
```

`--whole-dataset` selects all stored loci; `--regions` selects a BED union. Requested
loci absent from a dataset are refused, and stored filters cannot be changed during
an offline query. `--force` explicitly opts into atomic replacement. `--require-os-limit`
requires an existing matching Linux cgroup-v2 memory limit. SIGINT/SIGTERM cancellation
and failed finalization leave no successfully published artifact/receipt pair.

The manifest pins Rosalind `=__ROSALIND_VERSION__`. Before that version is published,
use explicit `[patch.crates-io]` entries for both the candidate `rosalind-bio` checkout
and `rosalind-build-info` (`crates/build-info`). Commit the generated Cargo.lock;
`cargo test` resolves it once and subsequent builds can use `--locked --offline`.
The CI registry installation is intended for the published version.
