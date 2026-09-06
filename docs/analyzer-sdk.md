# Analyzer SDK

## Build prerequisites

An installed native binary or wheel needs no Rust compiler. Building Rosalind or
an analyzer needs Rust/Cargo, a C and C++ compiler, CMake, and pkg-config; htslib's
compression dependencies also need development headers. Use a current stable
Rust toolchain for a newly resolved analyzer dependency graph. The source checkout
declares its minimum supported Rust version in `Cargo.toml` and tests it with the
repository lockfile.

On Debian/Ubuntu, install the native prerequisites before running Cargo:

```sh
sudo apt-get update
sudo apt-get install -y build-essential clang cmake pkg-config libbz2-dev liblzma-dev zlib1g-dev
```

On macOS, install Xcode Command Line Tools (`xcode-select --install`) and, with
Homebrew, `brew install cmake pkg-config xz`. Python source wheel builds additionally
need Python 3.9+ and Maturin; see [Python installation](../python/README.md).

## Exact evidence consumers

Use `rosalind::evidence` for new exact short-read evidence consumers. Open an
`EvidenceRequest`, resolve its VCF/BED selection against the engine's contigs, and
pass an `EvidenceAnalyzer` to `run`. The independent
[example crate](../examples/evidence-analyzer/) demonstrates the full public path.
The source checkout uses a local Cargo dependency. Native bundles include its
manifest, lockfile, and Rust source, with the manifest rendered to the bundle's
exact registry SDK version; `ONBOARDING-BUNDLE.json` records that transformation.
An unpublished candidate still needs the explicit source patches described below.

An analyzer declares `EvidenceRequirements`: field capabilities, whether real
reference bases are required, flanking context, and a conservative peak retained
byte bound. Call `plan_for_analyzer` before exposing the final plan or creating an
output. Unknown memory is allowed without a declared budget and rejected for
budgeted runs. Flanking context remains zero. Set `request.fields` explicitly to
match the desired groups: omitted groups have no accumulator or encoder buffer.
The engine validates that the request contains every analyzer requirement; it does
not silently expand scientific metrics. Full output retains schema 1; physical
projections use schema 2 with a versioned field mask.

`on_batch` borrows ordered exact rows through `batch.rows()` and `batch.row(index)`.
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

The engine supplies exact indexed extraction and resource planning. It does not
automatically supply every standalone binary with transactional files, receipts,
or OS enforcement. The CLI and Python `materialize_evidence` provide the first-party
artifact lifecycle. Custom binaries must explicitly add the lifecycle they claim.
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
generated `locus-qc/Cargo.toml` before running Cargo. Replace the example paths
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

The scaffold includes repeat bytes, receipt integrity, relocated inputs, refusal,
and explicit external replay tests. It remains a compatibility option while the
new evidence API and its release/adoption gates mature.
