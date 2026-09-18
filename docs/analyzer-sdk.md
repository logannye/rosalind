# Analyzer SDK

## Build prerequisites

An installed native binary or wheel needs no Rust compiler. Building Rosalind or
an analyzer needs Rust/Cargo, a C and C++ compiler, CMake, and pkg-config; htslib's
compression dependencies also need development headers. Use Rust 1.83 or newer.
The supported dependency constraints are tested with both the repository lockfile
and a newly resolved generated analyzer on Rust 1.83. Retain the generated lockfile
to keep subsequent builds reproducible.

On Debian/Ubuntu, install the native prerequisites before running Cargo:

```sh
sudo apt-get update
sudo apt-get install -y build-essential clang cmake pkg-config libbz2-dev liblzma-dev zlib1g-dev
```

On macOS, install Xcode Command Line Tools (`xcode-select --install`) and, with
Homebrew, `brew install cmake pkg-config xz`. Python source wheel builds additionally
need Python 3.9+ and Maturin; see [Python installation](../python/README.md).

## Exact evidence consumers

Use `rosalind::evidence::run_evidence_artifact` for new standalone exact short-read
analyzers. Implement `EvidenceArtifactFactory` to declare requirements and scientific
parameters before allocation; its `create` method returns an `EvidenceAnalyzer`
borrowing the runner-owned output sink. The runner validates the created analyzer's
requirements, admits the complete working set, checks cancellation and immutable
inputs, finalizes the analyzer, and publishes the output and receipt together.

`EvidenceArtifactSource::Native` reads indexed BAM/CRAM plus a local reference.
CRAM uses the [checked decoder profile](SEMANTICS.md#cram-decoder-admission),
including one complete cooperative validation pass before indexed extraction.
Its cost and decoder envelope are recorded separately from queried record visits.
`EvidenceArtifactSource::Dataset` reads a verified portable dataset without reopening
its original BAM/reference. Both feed canonical batches to the same analyzer. Dataset
queries require complete locus coverage and available fields; stored profile/sample
identity cannot be changed. Selection files are hashed and guarded before parsing.
The independent [example crate](../examples/evidence-analyzer/) exercises this public
artifact path with a small checked candidate-summary reducer.

<!-- smoke:evidence-scaffold -->
```sh
rosalind new analyzer candidate-qc --api evidence --output ./candidate-qc
cd candidate-qc
cargo test
cargo build --release --locked --offline
rosalind conformance analyzer --api evidence --binary ./target/release/candidate-qc --json > conformance.json
```

The source example uses local Cargo dependencies. Native bundles include its
manifest, lockfile, build script, and Rust source, with dependencies rendered to the
bundle's exact registry SDK versions. `ONBOARDING-BUNDLE.json` records those changes.
An unpublished candidate needs the explicit source patches described below.

An analyzer declares `EvidenceRequirements`: field capabilities, whether real
reference bases are required, flanking context, and a conservative peak retained
byte bound covering `create`, `on_batch`, `finish`, and encoder flush transients.
The artifact runner plans before creating the analyzer or output. Unknown memory
is refused under enforcement. `RecordOnly` can record a budget with an unknown
bound, but does not enforce that budget. Flanking context remains zero.
Set `request.fields` explicitly to
match the desired groups: omitted groups have no accumulator or encoder buffer.
The engine validates that the request contains every analyzer requirement; it does
not silently expand scientific metrics. Full output retains schema 1; physical
projections use schema 2 with a versioned field mask.

`on_batch` borrows ordered exact rows through `batch.rows()` and `batch.row(index)`.
Managed runs emit canonical batches of at most 1,024 rows, with boundaries that
do not depend on execution microtile width or persisted partition layout.
Groups such as `row.depths` and `row.alleles` are optional borrowed values. Check
capabilities before reading metrics; absence never means zero. Copies retained
after a callback belong to
the analyzer's memory model. Use checked integer reducers and deterministic order;
do not retain a genome table, timestamps, random identifiers, or unordered output.
`PanelQcAnalyzer` keeps state proportional to original target count and declares
that memory. `FusedAnalyzers` combines two consumers over one traversal.

The unpublished batch API changed from a public `Vec<EvidenceRow>` to private
contiguous group storage. Migrate `batch.rows.len()` to `batch.len()` and
`batch.rows.iter()` to `batch.rows()`. `row.try_to_full_row()` copies a full owned
row only when every group is available. `EvidenceBatch::from_full_rows(...)`
provides an explicit fixture/compatibility constructor. Existing `ColumnAnalyzer`
APIs are unchanged.

Use `EvidenceArrowWriter::with_fields(output, fields)` or the TSV equivalent;
`new(output)` remains ALL. Configure fields before writing, including empty
streams. `read_evidence_batches_expected_fields` validates the expected schema
before allocating record buffers; use its field-specific reader memory bound for
budgeted persisted consumption. The metadata-returning reader reports fields even
for an empty stream. Custom consumers without a known input mask should reserve
the full reader bound.

Set `EvidenceArtifactSpec::handle_signals = true` for executable entry points;
SIGINT/SIGTERM then requests cooperative cancellation through the shared scope.
Library hosts may supply a cancellation token. The runner governs callbacks and
finalization, but a callback must return or check cancellation to cooperate promptly.
An arbitrary native allocation can exceed an envelope before a cooperative check;
`RequireOsLimit` additionally verifies an existing Linux cgroup-v2 limit.
The managed resource/cancellation scope is process-wide: only one managed runner
may be active per process. Serialize independent runs or use separate processes.

Factories record analyzer parameters through `params()` and a tokenized
`ReplayInvocation`. Every parameter changing bytes or scientific interpretation
must be represented in both. External receipts never choose an executable implicitly:

```sh
rosalind reproduce --manifest summary.tsv.manifest.json --inputs input-directory \
  --binary /absolute/path/to/candidate-qc --dry-run
```

`EvidenceArtifactError::exit_code()` distinguishes invalid input/output (2), refusal
(3), resource breach (4), integrity failure (5), and cancellation (130). Failed runs
publish no successful artifact/receipt pair. Use `OutputPolicy::ReplaceAtomic` only
when the user has selected replacement (`--force` in the scaffold).

`EvidenceEngine` remains a lower-level extraction/planning API for hosts that own
their own lifecycle. It does not itself publish an atomic receipted artifact.
See [SEMANTICS.md](SEMANTICS.md) for depth/filter/coordinate rules.

## Legacy ColumnAnalyzer compatibility

`ColumnAnalyzer` remains for existing pileup-column consumers and the scaffold:

<!-- smoke:legacy-scaffold -->
```sh
rosalind new analyzer locus-qc --output ./locus-qc
cd locus-qc
cargo test
cargo build --release --locked --offline
rosalind conformance analyzer --binary ./target/release/locus-qc --json > conformance.json
```

`cargo test` creates debug test artifacts; the explicit release build creates the
binary used by conformance, using the dependencies already fetched by the tests.
Run this sequence in a fresh directory with the desired `rosalind` on `PATH`,
and leave `CARGO_TARGET_DIR` unset so the documented
binary path is correct. The generated manifest pins the SDK version to the
generating executable. It resolves from crates.io only after that exact version
is published.

For an unpublished checkout or candidate, add explicit local patches to the
generated analyzer `Cargo.toml` before running Cargo. Replace the example paths
with absolute paths to that same candidate checkout:

```toml
[patch.crates-io]
rosalind-bio = { path = "/absolute/path/to/rosalind" }
rosalind-build-info = { path = "/absolute/path/to/rosalind/crates/build-info" }
```

This validates the candidate SDK, not a registry installation. Remove those
patches and regenerate the lockfile to validate a subsequently published SDK.
The maintained package smokes distinguish `--candidate-source PATH` from
`--registry-sdk`; only the latter tests registry-only resolution.

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

The legacy scaffold preserves its existing API and conformance suite. New batch
consumers should choose `--api evidence` for the exact artifact runner.
