mod evidence_cli;

use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

use anyhow::{anyhow, bail, Context, Result};
use clap::{Args, CommandFactory, FromArgMatches, Parser, Subcommand, ValueEnum};
use rosalind::core::MemoryBudget;
use rosalind::genomics::{
    compare_callsets, create_bam_writer, read_vcf_variants, render_plan_line,
    sort_bam_deterministic, AlignedRead, BWTAligner, BedIndex, BuildMemoryModel, CigarOp,
    CigarOpKind, GenomeIndex, IndexBuildReport, IndexReader, IndexWriter,
};
use rosalind::io::decompress::open_input;
use rosalind::io::fasta::{FastaReader, FastaRecord};
use rosalind::io::fastq::{FastqReader, FastqRecord};
use rosalind::util::rss::peak_rss_bytes;
use rust_htslib::bam::Read as BamRead;
use rust_htslib::bam::{
    self, record::Aux, record::Cigar as BamCigar, record::CigarString, record::Record,
};

#[derive(Parser, Debug)]
#[command(
    name = "rosalind",
    version,
    about = "Deterministic, resource-governed per-locus genomics analyses with portable receipts"
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Args, Debug, Clone, Default)]
struct SelectionArgs {
    /// One samtools-style 1-based inclusive interval (`chr:start-end`).
    #[arg(long, conflicts_with_all = ["regions", "shard_count", "shard_index"])]
    region: Option<String>,
    /// BED intervals (zero-based half-open); requires an alignment index.
    #[arg(long, conflicts_with_all = ["region", "shard_count", "shard_index"])]
    regions: Option<PathBuf>,
    /// Total deterministic reference-span shards.
    #[arg(long, requires = "shard_index", conflicts_with_all = ["region", "regions"])]
    shard_count: Option<u32>,
    /// Zero-based deterministic shard index.
    #[arg(long, requires = "shard_count", conflicts_with_all = ["region", "regions"])]
    shard_index: Option<u32>,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Build, inspect, or convert lightweight analysis reference packs.
    Reference {
        #[command(subcommand)]
        action: ReferenceAction,
    },
    /// Canonically merge a complete compatible first-party shard set.
    Merge {
        /// Shard receipt; repeat once per shard.
        #[arg(long, required = true)]
        manifest: Vec<PathBuf>,
        /// Root used to content-locate relocated shard artifacts.
        #[arg(long)]
        inputs: Vec<PathBuf>,
        /// Canonical merged artifact destination.
        #[arg(short, long)]
        output: PathBuf,
        /// Merge receipt destination (default: `<output>.manifest.json`).
        #[arg(long)]
        output_manifest: Option<PathBuf>,
        /// Atomically replace existing output and receipt.
        #[arg(long)]
        force: bool,
    },
    /// Align reads against a reference genome and emit SAM records.
    Align {
        /// Reference genome in FASTA format (only the first record is used).
        #[arg(long)]
        reference: PathBuf,
        /// Reads file in FASTQ format (single-end).
        #[arg(long)]
        reads: Option<PathBuf>,
        /// Reads R1 file in FASTQ format (paired-end).
        #[arg(long)]
        reads_r1: Option<PathBuf>,
        /// Reads R2 file in FASTQ format (paired-end).
        #[arg(long)]
        reads_r2: Option<PathBuf>,
        /// Maximum mismatches permitted when seeding alignments.
        #[arg(long, default_value_t = 2)]
        max_mismatches: usize,
        /// Offset applied to reported reference positions (1-based in SAM).
        #[arg(long, default_value_t = 0)]
        reference_offset: u32,
        /// Output format for the alignment.
        #[arg(long, value_enum, default_value_t = OutputFormat::Sam)]
        format: OutputFormat,
        /// Optional path to write the output (stdout if omitted for SAM).
        #[arg(short, long)]
        output: Option<PathBuf>,
        /// Atomically replace existing output and receipt destinations.
        #[arg(long)]
        force: bool,
        /// Receipt path (default: `<output>.manifest.json` for file output).
        #[arg(long)]
        manifest: Option<PathBuf>,
    },
    /// Call germline variants from aligned reads (streaming pileup engine +
    /// abstention-aware genotype-likelihood caller).
    Variants {
        /// Persisted index (`rosalind index`); calls all contigs, reference from
        /// the index. Mutually exclusive with `--reference`.
        #[arg(long, conflicts_with_all = ["reference", "reference_pack"], required_unless_present_any = ["reference", "reference_pack"])]
        index: Option<PathBuf>,
        /// Reference genome (FASTA) — single-contig path. Mutually exclusive with `--index`.
        #[arg(long, conflicts_with = "reference_pack", required_unless_present_any = ["index", "reference_pack"])]
        reference: Option<PathBuf>,
        /// Lightweight analysis reference built by `rosalind reference build`.
        #[arg(long, required_unless_present_any = ["index", "reference"])]
        reference_pack: Option<PathBuf>,
        /// Alignments in SAM or BAM format (coordinate-sorted for `--index`).
        #[arg(long)]
        alignments: PathBuf,
        /// Chromosome name (single-contig `--reference` path only; defaults to the
        /// first FASTA record). Not allowed with `--index`.
        #[arg(long)]
        chrom: Option<String>,
        /// Starting offset (0-based) for the reference region (`--reference` only).
        #[arg(long, default_value_t = 0)]
        region_start: u32,
        /// Minimum MAPQ required for a read to be considered.
        #[arg(long, default_value_t = 0)]
        mapq_threshold: u8,
        /// Optional VCF output path (stdout if omitted).
        #[arg(short, long)]
        output: Option<PathBuf>,
        /// Deprecated and ignored.
        #[arg(long, default_value_t = 1024, hide = true)]
        block_size: usize,
        /// Minimum quality threshold for reporting variants.
        #[arg(long, default_value_t = 30.0)]
        quality_threshold: f32,
        /// Declared memory budget (MiB) for the run — records a plan/peak line.
        /// With `--enforce` it is honored (exit 3 refuse / exit 4 breach). (`--index` path.)
        #[arg(long)]
        memory_budget_mb: Option<u64>,
        /// Declare the active read capacity; exceeding it fails with partial output;
        /// the bound `plan`/`--enforce` rely on. `0` = uncapped.
        #[arg(long, default_value_t = 1000)]
        max_depth: u32,
        /// Max read length assumed by the pre-run `--enforce` estimate AND enforced
        /// at ingest under `--enforce` (a longer read aborts the run).
        #[arg(long, default_value_t = 250)]
        max_read_len: u32,
        /// Honor the budget: refuse up front if predicted peak exceeds it (exit 3),
        /// or fail after the run if the realized peak does (exit 4). Requires
        /// `--memory-budget-mb` and `--max-depth > 0`.
        #[arg(long, default_value_t = false)]
        enforce: bool,
        /// Require an existing Linux cgroup-v2 memory.max at or below the budget.
        #[arg(long, requires = "enforce")]
        require_os_limit: bool,
        /// Atomically replace an existing output and receipt.
        #[arg(long)]
        force: bool,
        /// Where to write the reproducibility receipt. Default: `<output>.manifest.json`
        /// for file output, or `./rosalind.variants.manifest.json` for stdout output.
        #[arg(long)]
        manifest: Option<PathBuf>,
        /// Emit a banded gVCF (every callable locus → a variant or a `<NON_REF>`
        /// reference block) instead of a sites-only VCF. Output is bounded and
        /// byte-reproducible; cohort integration is not yet qualified.
        #[arg(long)]
        gvcf: bool,
        #[command(flatten)]
        selection: SelectionArgs,
    },
    /// Stream a bounded, deterministic per-locus FEATURE table (TSV) over a
    /// persisted index — the same memory contract as `variants`, but every
    /// callable locus is emitted as ML-ready features. Byte-identical run-to-run.
    Features {
        /// Legacy persisted search index. Prefer `--reference-pack` for analysis.
        #[arg(
            long,
            conflicts_with = "reference_pack",
            required_unless_present = "reference_pack"
        )]
        index: Option<PathBuf>,
        /// Lightweight analysis reference built by `rosalind reference build`.
        #[arg(long, required_unless_present = "index")]
        reference_pack: Option<PathBuf>,
        /// Coordinate-sorted alignments (BAM).
        #[arg(long)]
        alignments: PathBuf,
        /// Minimum MAPQ required for a read to be considered.
        #[arg(long, default_value_t = 0)]
        mapq_threshold: u8,
        /// Declared memory budget (MiB). With `--enforce` it is honored (exit 3/4).
        #[arg(long)]
        memory_budget_mb: Option<u64>,
        /// Active-read capacity (exact-or-fail). `0` = uncapped.
        #[arg(long, default_value_t = 1000)]
        max_depth: u32,
        /// Max read length assumed by the `--enforce` estimate and enforced at ingest.
        #[arg(long, default_value_t = 250)]
        max_read_len: u32,
        /// Honor the budget: refuse up front (exit 3) / fail after (exit 4).
        #[arg(long, default_value_t = false)]
        enforce: bool,
        /// Require an existing Linux cgroup-v2 memory.max at or below the budget.
        #[arg(long, requires = "enforce")]
        require_os_limit: bool,
        /// Atomically replace an existing output and receipt.
        #[arg(long)]
        force: bool,
        /// Output TSV path (stdout if omitted).
        #[arg(short, long)]
        output: Option<PathBuf>,
        /// Where to write the reproducibility receipt (default: `<output>.manifest.json`).
        #[arg(long)]
        manifest: Option<PathBuf>,
        #[command(flatten)]
        selection: SelectionArgs,
        /// Feature artifact encoding.
        #[arg(long, value_enum, default_value_t = FeatureFormat::Tsv)]
        format: FeatureFormat,
    },
    /// Extract indexed exact evidence or panel QC, or run a legacy column analyzer,
    /// with a verifiable receipt recording the analysis parameters.
    Analyze {
        /// Which analyzer to run.
        #[arg(value_enum)]
        kind: AnalyzerKind,
        /// Legacy persisted search index. Prefer `--reference-pack` for analysis.
        #[arg(long, conflicts_with_all = ["reference_pack", "reference"])]
        index: Option<PathBuf>,
        /// Lightweight analysis reference built by `rosalind reference build`.
        #[arg(long, conflicts_with = "reference")]
        reference_pack: Option<PathBuf>,
        /// Indexed local FASTA for exact evidence; also accepted: .rref or .idx.
        #[arg(long)]
        reference: Option<PathBuf>,
        /// Coordinate-sorted alignments (BAM; CRAM for evidence/panel-qc).
        #[arg(long)]
        alignments: PathBuf,
        /// Minimum MAPQ required for a read to be considered.
        #[arg(long)]
        mapq_threshold: Option<u8>,
        /// Memory budget (MiB). Evidence/panel QC enforce it automatically;
        /// legacy analyzers require `--enforce` (refusal 3, runtime limit 4).
        #[arg(long)]
        memory_budget_mb: Option<u64>,
        /// Active-read capacity (exact-or-fail). `0` = uncapped.
        #[arg(long, default_value_t = 1000)]
        max_depth: u32,
        /// Max read length assumed by the `--enforce` estimate and enforced at ingest.
        #[arg(long, default_value_t = 250)]
        max_read_len: u32,
        /// Honor the budget: refuse up front (exit 3) / fail after (exit 4).
        #[arg(long, default_value_t = false)]
        enforce: bool,
        /// Require an existing Linux cgroup-v2 memory.max at or below the budget.
        #[arg(long, requires = "enforce")]
        require_os_limit: bool,
        /// Atomically replace an existing output and receipt.
        #[arg(long)]
        force: bool,
        /// Output path (stdout if omitted).
        #[arg(short, long)]
        output: Option<PathBuf>,
        /// Where to write the reproducibility receipt (default: `<output>.manifest.json`).
        #[arg(long)]
        manifest: Option<PathBuf>,
        #[command(flatten)]
        selection: SelectionArgs,
        #[command(flatten)]
        evidence: evidence_cli::EvidenceOptions,
    },
    /// Deterministically coordinate-sort a BAM file using bounded memory.
    Sort {
        /// Input BAM path.
        #[arg(long)]
        input: PathBuf,
        /// Output BAM path.
        #[arg(short, long)]
        output: PathBuf,
        /// Memory budget (MiB) for in-memory sorting chunks.
        #[arg(long, default_value_t = 1024)]
        memory_mb: usize,
        /// Atomically replace existing output and receipt destinations.
        #[arg(long)]
        force: bool,
        /// Receipt path (default: `<output>.manifest.json`).
        #[arg(long)]
        manifest: Option<PathBuf>,
    },
    /// End-to-end tumor/normal somatic calling (align + sort + call).
    Somatic {
        /// Reference FASTA (single-contig for now).
        #[arg(long)]
        reference: PathBuf,
        /// Tumor FASTQ reads (single-end).
        #[arg(long)]
        tumor: Option<PathBuf>,
        /// Tumor FASTQ reads R1 (paired-end).
        #[arg(long)]
        tumor_r1: Option<PathBuf>,
        /// Tumor FASTQ reads R2 (paired-end).
        #[arg(long)]
        tumor_r2: Option<PathBuf>,
        /// Normal FASTQ reads (single-end).
        #[arg(long)]
        normal: Option<PathBuf>,
        /// Normal FASTQ reads R1 (paired-end).
        #[arg(long)]
        normal_r1: Option<PathBuf>,
        /// Normal FASTQ reads R2 (paired-end).
        #[arg(long)]
        normal_r2: Option<PathBuf>,
        /// Output VCF path.
        #[arg(short, long)]
        output: PathBuf,
        /// Working directory for intermediate BAMs.
        #[arg(long)]
        workdir: Option<PathBuf>,
        /// Sorting memory budget (MiB).
        #[arg(long, default_value_t = 1024)]
        memory_mb: usize,
        /// Atomically replace existing outputs, receipts, and work artifacts.
        #[arg(long)]
        force: bool,
    },
    /// Compare a called somatic VCF against a truth VCF (optionally masked by BED).
    EvalSomatic {
        /// Reference FASTA (single contig slice used by the VCFs).
        #[arg(long)]
        reference: PathBuf,
        /// Called VCF path.
        #[arg(long)]
        calls: PathBuf,
        /// Truth VCF path.
        #[arg(long)]
        truth: PathBuf,
        /// Optional BED mask (0-based half-open).
        #[arg(long)]
        regions: Option<PathBuf>,
    },
    /// Compare a called germline VCF against a truth VCF (e.g. a GIAB benchmark),
    /// optionally masked by a high-confidence BED. Reports precision/recall/F1.
    EvalGermline {
        /// Reference FASTA (the contigs the VCFs use).
        #[arg(long)]
        reference: PathBuf,
        /// Called VCF path.
        #[arg(long)]
        calls: PathBuf,
        /// Truth VCF path.
        #[arg(long)]
        truth: PathBuf,
        /// Optional high-confidence BED mask (0-based half-open).
        #[arg(long)]
        regions: Option<PathBuf>,
        /// Which emitted calls participate in the comparison.
        #[arg(long, value_enum, default_value_t = CallsFilter::All)]
        calls_filter: CallsFilter,
        /// Emit a machine-readable metrics object.
        #[arg(long)]
        json: bool,
    },
    /// Build a reference index once into a portable, memory-mappable artifact.
    Index {
        /// Reference genome in FASTA (plain or gzip; `-` for stdin). All contigs.
        #[arg(long)]
        reference: PathBuf,
        /// Output path for the index artifact.
        #[arg(short, long)]
        output: PathBuf,
        /// Declared memory budget (MiB) for the build. Records a plan line; does
        /// not enforce (enforcement is a later phase).
        #[arg(long)]
        memory_budget_mb: Option<u64>,
        /// Atomically replace existing index and receipt destinations.
        #[arg(long)]
        force: bool,
    },
    /// Scaffold a downstream project that inherits Rosalind's contract.
    New {
        #[command(subcommand)]
        action: NewAction,
    },
    /// Run the complete contract story offline on embedded toy data.
    Demo {
        /// Directory to create. It must be absent or empty.
        #[arg(long, default_value = "rosalind-demo")]
        output_dir: PathBuf,
        /// Declared budget for the enforced demo run.
        #[arg(long, default_value_t = 128)]
        budget_mb: u64,
        /// Emit one JSON summary after the demo completes.
        #[arg(long)]
        json: bool,
    },
    /// Preflight index/BAM compatibility, output safety, and memory feasibility.
    Doctor {
        /// Legacy persisted search index. Prefer `--reference-pack`.
        #[arg(
            long,
            conflicts_with = "reference_pack",
            required_unless_present = "reference_pack"
        )]
        index: Option<PathBuf>,
        /// Lightweight analysis reference.
        #[arg(long, required_unless_present = "index")]
        reference_pack: Option<PathBuf>,
        /// Coordinate-sorted BAM to validate against the index.
        #[arg(long)]
        alignments: PathBuf,
        /// Planned output path whose collision safety should be checked.
        #[arg(long)]
        output: Option<PathBuf>,
        /// Intended memory budget in MiB.
        #[arg(long)]
        budget_mb: Option<u64>,
        /// Scan every mapped read to prove coordinate order and maximum read length.
        #[arg(long)]
        deep: bool,
        /// Emit a stable machine-readable report.
        #[arg(long)]
        json: bool,
        #[command(flatten)]
        selection: SelectionArgs,
    },
    /// Open the embedded, loopback-only Receipt Studio with optional preloaded receipts.
    Studio {
        /// Run receipts and reproduction certificates to preload.
        receipts: Vec<PathBuf>,
        /// Print the URL without opening the system browser.
        #[arg(long)]
        no_open: bool,
        /// Loopback port; zero chooses an unused port.
        #[arg(long, default_value_t = 0)]
        port: u16,
        /// Emit one startup JSON object before serving.
        #[arg(long)]
        json: bool,
    },
    /// Walk a directory of receipts as a provenance DAG and verify it offline.
    Chain {
        #[command(subcommand)]
        action: ChainAction,
    },
    /// Inspect, sanitize, or export reproducibility receipts.
    Receipt {
        #[command(subcommand)]
        action: ReceiptAction,
    },
    /// Validate an external analyzer against Rosalind's complete builder contract.
    Conformance {
        #[command(subcommand)]
        action: ConformanceAction,
    },
    /// Locate exact occurrences of a pattern in a prebuilt index (load + query).
    Locate {
        /// Index artifact built by `rosalind index`.
        #[arg(long)]
        index: PathBuf,
        /// Pattern to locate (ASCII A/C/G/T/N; case-insensitive).
        #[arg(long)]
        pattern: String,
        /// Maximum number of candidate hits to locate.
        #[arg(long, default_value_t = 1024)]
        max_hits: usize,
    },
    /// Predict whether a job fits a declared memory budget, before committing.
    Plan {
        /// Persisted index (`rosalind index`): predict the bounded whole-genome
        /// `variants` peak. Mutually exclusive with `--reference`.
        #[arg(long, conflicts_with_all = ["reference", "reference_pack"], required_unless_present_any = ["reference", "reference_pack"])]
        index: Option<PathBuf>,
        /// Reference FASTA: predict the index BUILD peak (advisory — build is
        /// O(reference); Phase D enforces). Mutually exclusive with `--index`.
        #[arg(long, conflicts_with = "reference_pack", required_unless_present_any = ["index", "reference_pack"])]
        reference: Option<PathBuf>,
        /// Analysis reference pack: predict a bounded analyzer run.
        #[arg(long, required_unless_present_any = ["index", "reference"])]
        reference_pack: Option<PathBuf>,
        /// Max active depth assumed for the `variants` working-set bound.
        #[arg(long, default_value_t = 1000)]
        max_depth: u32,
        /// Max read length assumed for the `variants` working-set bound.
        #[arg(long, default_value_t = 250)]
        max_read_len: u32,
        /// Declared memory budget (MiB) to check feasibility against.
        #[arg(long)]
        budget_mb: Option<u64>,
        /// Include the canonical Arrow feature encoder's additional working set.
        #[arg(long, value_enum, default_value_t = FeatureFormat::Tsv)]
        format: FeatureFormat,
        /// Emit the predicted-peak breakdown as one-line JSON (for a scheduler/CI
        /// to read), instead of the human-readable table. `--index` only.
        #[arg(long)]
        json: bool,
        #[command(flatten)]
        selection: SelectionArgs,
    },
    /// Pack many bounded `variants` jobs onto fixed-size nodes by their PREDICTED
    /// peaks — show a co-location fits within budget before launching a byte. Each
    /// job's peak is read from its index header (no run); peaks are additive, so the
    /// sum is a conservative bound a scheduler can refuse on. Exit 3 if no packing fits.
    Pack {
        /// A jobs file: one job per line, `<index_path>[\t<max_depth>[\t<max_read_len>]]`
        /// (TSV or whitespace; blank lines and `#` comments ignored).
        #[arg(long)]
        jobs: PathBuf,
        /// Per-node memory capacity, in MiB (the RAM each node can give a co-located batch).
        #[arg(long)]
        node_mb: u64,
        /// Cap on the number of nodes; refuse (exit 3) if the jobs need more. Unset = as many as needed.
        #[arg(long)]
        nodes: Option<usize>,
        /// Default max active depth for jobs that do not specify one.
        #[arg(long, default_value_t = 1000)]
        max_depth: u32,
        /// Default max read length for jobs that do not specify one.
        #[arg(long, default_value_t = 250)]
        max_read_len: u32,
        /// Emit the schedule as JSON instead of the human-readable plan.
        #[arg(long)]
        json: bool,
    },
    /// Re-check a reproducibility receipt without re-running: re-hash its inputs
    /// and outputs and confirm the realized peak landed within the budget.
    Verify {
        /// Path to a `*.manifest.json` written by a previous run.
        #[arg(long)]
        manifest: PathBuf,
        /// Budget (MiB) to check the recorded peak against (overrides the
        /// `memory_budget_mb` recorded in the manifest, if any).
        #[arg(long)]
        budget_mb: Option<u64>,
        /// Assert the receipt was built from exactly this commit SHA (prefix ok).
        /// Fails verify on a mismatch, or on a clean match from a dirty build.
        #[arg(long)]
        expect_code: Option<String>,
        /// Emit a stable JSON report instead of human-readable lines.
        #[arg(long)]
        json: bool,
    },
    /// Localize how two run receipts' claims differ, bucketed by causal role: inputs /
    /// code-identity / params (causes), outputs (effect), measurements (noise). Exit
    /// 0 = identical claims, 1 = claims differ, 2 = read/parse error.
    Diff {
        /// First receipt (`*.manifest.json`).
        a: PathBuf,
        /// Second receipt (`*.manifest.json`).
        b: PathBuf,
        /// Emit a compact JSON summary instead of the human report.
        #[arg(long)]
        json: bool,
        /// Create a streaming per-locus metric delta TSV from verified exact evidence.
        #[arg(long)]
        loci_output: Option<PathBuf>,
    },
    /// Re-derive a recorded result from its receipt and content-located inputs, and
    /// write a chainable reproduction certificate. The verdict is over output bytes
    /// (exit 0 REPRODUCED / 6 DIVERGED / 7 INCONCLUSIVE / 5 tampered receipt).
    Reproduce {
        /// Path to a `*.manifest.json` from a previous run.
        #[arg(long)]
        manifest: PathBuf,
        /// Directory holding the recorded inputs (located by content hash).
        #[arg(long)]
        inputs: PathBuf,
        /// Do not write a `.repro.json` reproduction certificate.
        #[arg(long, default_value_t = false)]
        no_attest: bool,
        /// Where to write the certificate (default: `<manifest>.repro.json`).
        #[arg(short, long)]
        output: Option<PathBuf>,
        /// Executable to use for replay. Required for third-party analyzer receipts;
        /// a path recorded inside a receipt is never executed automatically.
        #[arg(long)]
        binary: Option<PathBuf>,
        /// Validate and print the isolated execution plan without running it.
        #[arg(long)]
        dry_run: bool,
        /// Atomically replace an existing reproduction certificate.
        #[arg(long)]
        force: bool,
        /// Emit a stable JSON report instead of human-readable lines.
        #[arg(long)]
        json: bool,
    },
    /// Emit a self-hosted status badge — a shields.io endpoint JSON or a static SVG.
    /// An intact receipt is blue; only a valid linked reproduction certificate is green.
    Badge {
        /// The run's `*.manifest.json`.
        #[arg(long)]
        manifest: PathBuf,
        /// An optional `*.repro.json` whose verdict backs the "reproducible" claim.
        #[arg(long)]
        repro: Option<PathBuf>,
        /// Output path; a `.svg` extension emits the static SVG, else a shields JSON.
        #[arg(short, long)]
        output: PathBuf,
        /// Atomically replace an existing badge.
        #[arg(long)]
        force: bool,
    },
}

#[derive(Subcommand, Debug)]
enum ChainAction {
    /// Verify every node self-hashes and every internal edge resolves by content hash.
    Verify {
        /// Directory of `*.manifest.json` receipts to walk.
        dir: PathBuf,
        /// Emit a compact JSON report instead of human-readable lines.
        #[arg(long)]
        json: bool,
    },
}

#[derive(Subcommand, Debug)]
enum ReceiptAction {
    /// Explain independent trust dimensions and exactly which evidence is missing.
    Inspect {
        /// Run receipt to inspect.
        #[arg(long)]
        manifest: PathBuf,
        /// Artifact file to content-match; repeat for multiple inputs and outputs.
        #[arg(long)]
        artifact: Vec<PathBuf>,
        /// Optional reproduction certificate to validate and link.
        #[arg(long)]
        certificate: Option<PathBuf>,
        /// Emit a stable machine-readable trust report.
        #[arg(long)]
        json: bool,
    },
    /// Replace recorded paths with stable role labels without changing the claim ID.
    Sanitize {
        /// Schema-3 through schema-5 run receipt.
        #[arg(long)]
        manifest: PathBuf,
        /// New sanitized receipt path; existing files are never overwritten.
        #[arg(short, long)]
        output: PathBuf,
    },
    /// Export an unsigned in-toto Statement v1 using the Rosalind predicate.
    ExportIntoto {
        /// Intact native schema-5 run receipt.
        #[arg(long)]
        manifest: PathBuf,
        /// New statement path; existing files are never overwritten.
        #[arg(short, long)]
        output: PathBuf,
    },
}

#[derive(Subcommand, Debug)]
enum ConformanceAction {
    /// Run embedded determinism, contract, receipt, replay, diff, and collision checks.
    Analyzer {
        /// External analyzer executable to test explicitly.
        #[arg(long)]
        binary: PathBuf,
        /// Emit a stable badge-ready result.
        #[arg(long)]
        json: bool,
    },
}

#[derive(Subcommand, Debug)]
enum NewAction {
    /// Create a standalone Rust `ColumnAnalyzer` binary and contract tests.
    Analyzer {
        /// Lowercase kebab-case package and analyzer name.
        name: String,
        /// Destination directory (must be absent or empty).
        #[arg(short, long)]
        output: PathBuf,
    },
}

#[derive(Copy, Clone, Debug, ValueEnum, Eq, PartialEq)]
enum OutputFormat {
    Sam,
    Bam,
}

#[derive(Copy, Clone, Debug, ValueEnum, Eq, PartialEq)]
enum CallsFilter {
    All,
    Pass,
}

/// Registered per-locus analyzers for `rosalind analyze <kind>`. A compile-time
/// registry — adding a kind is one variant + one dispatch arm.
#[derive(Copy, Clone, Debug, ValueEnum, Eq, PartialEq)]
enum AnalyzerKind {
    Features,
    Coverage,
    Evidence,
    PanelQc,
}

#[derive(Copy, Clone, Debug, ValueEnum, Eq, PartialEq)]
enum FeatureFormat {
    Tsv,
    ArrowIpc,
}

#[derive(Subcommand, Debug)]
enum ReferenceAction {
    /// Stream a FASTA into a deterministic, mmap-friendly `.rref`.
    Build {
        /// Input FASTA path (plain or gzip). A path is required for the two-pass build.
        #[arg(long)]
        fasta: PathBuf,
        /// Destination reference pack.
        #[arg(short, long)]
        output: PathBuf,
        /// Receipt destination (default: `<output>.manifest.json`).
        #[arg(long)]
        manifest: Option<PathBuf>,
        /// Atomically replace an existing destination.
        #[arg(long)]
        force: bool,
    },
    /// Validate and describe an analysis reference pack.
    Inspect {
        /// Reference pack to inspect.
        #[arg(long)]
        reference_pack: PathBuf,
        /// Emit stable machine-readable JSON.
        #[arg(long)]
        json: bool,
    },
    /// Extract the reference section of a legacy `.idx` into `.rref`.
    Convert {
        /// Legacy search index.
        #[arg(long)]
        index: PathBuf,
        /// Destination reference pack.
        #[arg(short, long)]
        output: PathBuf,
        /// Receipt destination (default: `<output>.manifest.json`).
        #[arg(long)]
        manifest: Option<PathBuf>,
        /// Atomically replace an existing destination.
        #[arg(long)]
        force: bool,
    },
}

#[derive(Debug, Clone)]
struct FastqPair {
    name: String,
    r1: FastqRecord,
    r2: FastqRecord,
}

#[derive(Debug)]
enum ResolvedReads {
    Single(Vec<FastqRecord>),
    Paired(Vec<FastqPair>),
}

struct AlignmentCandidate {
    position: usize,
    mismatches: usize,
    mapq: u8,
    is_reverse: bool,
    cigar: Vec<CigarOp>,
    as_score: i32,
    md: String,
}

fn main() -> Result<()> {
    // Install the build-identity (baked at compile time by build.rs) into the receipt
    // crate, so every receipt records exactly which code/toolchain/deps produced the run.
    rosalind::provenance::set_build_identity(rosalind::provenance::BuildIdentity {
        code_git_sha: env!("ROSALIND_GIT_SHA").to_string(),
        code_dirty: env!("ROSALIND_GIT_DIRTY").to_string(),
        rustc_version: env!("ROSALIND_RUSTC_VERSION").to_string(),
        target_triple: env!("ROSALIND_TARGET").to_string(),
        deps_lock_blake3: env!("ROSALIND_DEPS_LOCK_BLAKE3").to_string(),
    });
    let matches = Cli::command().get_matches();
    if let Some(("analyze", args)) = matches.subcommand() {
        let exact = matches!(
            args.get_one::<AnalyzerKind>("kind"),
            Some(AnalyzerKind::Evidence | AnalyzerKind::PanelQc)
        );
        let evidence_only = [
            "reference",
            "sites",
            "base_quality_threshold",
            "format",
            "tile_bases",
            "max_record_bytes",
            "workers",
            "cache_dir",
            "resume",
            "cram_reference",
            "alignment_index",
            "reference_fai",
            "cram_reference_fai",
            "min_callable_depth",
            "position_output",
            "plan",
        ];
        for name in evidence_only {
            if !exact && args.value_source(name) == Some(clap::parser::ValueSource::CommandLine) {
                clap::Error::raw(
                    clap::error::ErrorKind::ArgumentConflict,
                    format!(
                        "--{} requires analyze evidence or analyze panel-qc",
                        name.replace('_', "-")
                    ),
                )
                .exit();
            }
        }
        if exact && args.value_source("max_depth") == Some(clap::parser::ValueSource::CommandLine) {
            clap::Error::raw(clap::error::ErrorKind::ArgumentConflict,
                "--max-depth applies to legacy analyzers; exact evidence accumulates all eligible reads").exit();
        }
        if args.get_one::<AnalyzerKind>("kind") == Some(&AnalyzerKind::Evidence)
            && args.value_source("min_callable_depth")
                == Some(clap::parser::ValueSource::CommandLine)
        {
            clap::Error::raw(
                clap::error::ErrorKind::ArgumentConflict,
                "--min-callable-depth requires analyze panel-qc",
            )
            .exit();
        }
    }
    let cli = Cli::from_arg_matches(&matches)?;

    match cli.command {
        Commands::Reference { action } => run_reference(action)?,
        Commands::Merge {
            manifest,
            inputs,
            output,
            output_manifest,
            force,
        } => match rosalind::merge_shards(
            &manifest,
            &inputs,
            &output,
            output_manifest.as_deref(),
            force,
        ) {
            Ok(outcome) => eprintln!(
                "merged {} shards: {}; receipt: {}",
                outcome.shard_count,
                outcome.output.display(),
                outcome.manifest.display()
            ),
            Err(error) => {
                eprintln!("merge: {error}");
                std::process::exit(error.exit_code());
            }
        },
        Commands::Align {
            reference,
            reads,
            reads_r1,
            reads_r2,
            max_mismatches,
            reference_offset,
            format,
            output,
            force,
            manifest,
        } => run_align(
            reference,
            reads,
            reads_r1,
            reads_r2,
            max_mismatches,
            reference_offset,
            format,
            output,
            force,
            manifest,
        )?,
        Commands::Variants {
            index,
            reference,
            reference_pack,
            alignments,
            chrom,
            region_start,
            mapq_threshold,
            output,
            block_size: _,
            quality_threshold,
            memory_budget_mb,
            max_depth,
            max_read_len,
            enforce,
            require_os_limit,
            force,
            manifest,
            gvcf,
            selection,
        } => {
            if index.is_some() || reference_pack.is_some() {
                if chrom.is_some() || region_start != 0 {
                    bail!("--chrom/--region-start are not valid with --index/--reference-pack (the whole reference is called)");
                }
                run_variants_index(
                    select_analysis_reference(index, reference_pack),
                    alignments,
                    mapq_threshold,
                    output,
                    quality_threshold,
                    memory_budget_mb,
                    max_depth,
                    max_read_len,
                    enforce,
                    require_os_limit,
                    force,
                    manifest,
                    gvcf,
                    selection,
                )?
            } else {
                if selection_requested(&selection) {
                    bail!("--region/--regions/--shard-* require --index or --reference-pack");
                }
                if gvcf {
                    bail!("--gvcf requires --index (the bounded whole-genome path)");
                }
                let reference = reference.expect("clap guarantees one of --index/--reference");
                run_variants(
                    reference,
                    alignments,
                    chrom,
                    region_start,
                    mapq_threshold,
                    output,
                    1024,
                    quality_threshold,
                    force,
                    manifest,
                )?
            }
        }
        Commands::Features {
            index,
            reference_pack,
            alignments,
            mapq_threshold,
            memory_budget_mb,
            max_depth,
            max_read_len,
            enforce,
            require_os_limit,
            force,
            output,
            manifest,
            selection,
            format,
        } => run_features(
            select_analysis_reference(index, reference_pack),
            alignments,
            mapq_threshold,
            memory_budget_mb,
            max_depth,
            max_read_len,
            enforce,
            require_os_limit,
            force,
            output,
            manifest,
            selection,
            format,
        )?,
        Commands::Analyze {
            kind,
            index,
            reference_pack,
            reference,
            alignments,
            mapq_threshold,
            memory_budget_mb,
            max_depth,
            max_read_len,
            enforce,
            require_os_limit,
            force,
            output,
            manifest,
            selection,
            evidence,
        } => {
            if matches!(kind, AnalyzerKind::Evidence | AnalyzerKind::PanelQc) {
                evidence_cli::run(evidence_cli::EvidenceCommand {
                    panel: kind == AnalyzerKind::PanelQc,
                    reference: reference.or(reference_pack).or(index),
                    alignments,
                    mapq_threshold,
                    memory_budget_mb,
                    max_read_len,
                    enforce,
                    require_os_limit,
                    force,
                    output,
                    manifest,
                    selection,
                    options: evidence,
                })?;
            } else {
                if reference.is_some() || (index.is_none() && reference_pack.is_none()) {
                    bail!("legacy analyzers require --index or --reference-pack");
                }
                evidence.reject_legacy_options()?;
                let mapq_threshold = mapq_threshold.unwrap_or(0);
                let label = match kind {
                    AnalyzerKind::Features => "analyze features",
                    AnalyzerKind::Coverage => "analyze coverage",
                    AnalyzerKind::Evidence | AnalyzerKind::PanelQc => unreachable!(),
                };
                let analysis_reference = select_analysis_reference(index, reference_pack);
                match kind {
                    AnalyzerKind::Features => {
                        let mut a = rosalind::call::FeatureAnalyzer::default();
                        run_bounded_analysis(
                            label,
                            "analyzer.",
                            &mut a,
                            analysis_reference,
                            alignments,
                            mapq_threshold,
                            memory_budget_mb,
                            max_depth,
                            max_read_len,
                            enforce,
                            require_os_limit,
                            force,
                            output,
                            manifest,
                            selection,
                        )?
                    }
                    AnalyzerKind::Coverage => {
                        let mut a = rosalind::call::CoverageTrack;
                        run_bounded_analysis(
                            label,
                            "analyzer.",
                            &mut a,
                            analysis_reference,
                            alignments,
                            mapq_threshold,
                            memory_budget_mb,
                            max_depth,
                            max_read_len,
                            enforce,
                            require_os_limit,
                            force,
                            output,
                            manifest,
                            selection,
                        )?
                    }
                    AnalyzerKind::Evidence | AnalyzerKind::PanelQc => unreachable!(),
                }
            }
        }
        Commands::Sort {
            input,
            output,
            memory_mb,
            force,
            manifest,
        } => run_sort(input, output, memory_mb, force, manifest)?,
        Commands::Somatic {
            reference,
            tumor,
            tumor_r1,
            tumor_r2,
            normal,
            normal_r1,
            normal_r2,
            output,
            workdir,
            memory_mb,
            force,
        } => {
            run_somatic(
                reference, tumor, tumor_r1, tumor_r2, normal, normal_r1, normal_r2, output,
                workdir, memory_mb, force,
            )?;
        }
        Commands::EvalSomatic {
            reference,
            calls,
            truth,
            regions,
        } => {
            run_eval(reference, calls, truth, regions, CallsFilter::All, false)?;
        }
        Commands::EvalGermline {
            reference,
            calls,
            truth,
            regions,
            calls_filter,
            json,
        } => {
            run_eval(reference, calls, truth, regions, calls_filter, json)?;
        }
        Commands::Index {
            reference,
            output,
            memory_budget_mb,
            force,
        } => run_index(reference, output, memory_budget_mb, force)?,
        Commands::New { action } => match action {
            NewAction::Analyzer { name, output } => run_new_analyzer(&name, &output)?,
        },
        Commands::Demo {
            output_dir,
            budget_mb,
            json,
        } => run_demo(output_dir, budget_mb, json)?,
        Commands::Doctor {
            index,
            reference_pack,
            alignments,
            output,
            budget_mb,
            deep,
            json,
            selection,
        } => run_doctor_command(
            select_analysis_reference(index, reference_pack),
            alignments,
            output,
            budget_mb,
            deep,
            json,
            selection,
        )?,
        Commands::Studio {
            receipts,
            no_open,
            port,
            json,
        } => rosalind::serve_studio(&rosalind::StudioSpec {
            receipts,
            no_open,
            port,
            json,
        })?,
        Commands::Chain { action } => match action {
            ChainAction::Verify { dir, json } => run_chain_verify(dir, json)?,
        },
        Commands::Receipt { action } => match action {
            ReceiptAction::Inspect {
                manifest,
                artifact,
                certificate,
                json,
            } => run_receipt_inspect(manifest, artifact, certificate, json)?,
            ReceiptAction::Sanitize { manifest, output } => run_receipt_sanitize(manifest, output)?,
            ReceiptAction::ExportIntoto { manifest, output } => {
                run_receipt_export_intoto(manifest, output)?
            }
        },
        Commands::Conformance { action } => match action {
            ConformanceAction::Analyzer { binary, json } => run_conformance_analyzer(binary, json)?,
        },
        Commands::Locate {
            index,
            pattern,
            max_hits,
        } => run_locate(index, pattern, max_hits)?,
        Commands::Plan {
            index,
            reference,
            reference_pack,
            max_depth,
            max_read_len,
            budget_mb,
            format,
            json,
            selection,
        } => run_plan(
            index.or(reference_pack),
            reference,
            max_depth,
            max_read_len,
            budget_mb,
            json,
            selection,
            format,
        )?,
        Commands::Pack {
            jobs,
            node_mb,
            nodes,
            max_depth,
            max_read_len,
            json,
        } => run_pack(jobs, node_mb, nodes, max_depth, max_read_len, json)?,
        Commands::Verify {
            manifest,
            budget_mb,
            expect_code,
            json,
        } => run_verify(manifest, budget_mb, expect_code, json)?,
        Commands::Diff {
            a,
            b,
            json,
            loci_output,
        } => run_diff(a, b, json, loci_output)?,
        Commands::Reproduce {
            manifest,
            inputs,
            no_attest,
            output,
            binary,
            dry_run,
            force,
            json,
        } => run_reproduce(
            manifest, inputs, no_attest, output, binary, dry_run, force, json,
        )?,
        Commands::Badge {
            manifest,
            repro,
            output,
            force,
        } => run_badge(manifest, repro, output, force)?,
    }

    Ok(())
}

fn run_reference(action: ReferenceAction) -> Result<()> {
    use rosalind::genomics::{ReferencePackBuilder, ReferencePackReader, ReferenceProvider};

    let render_hash =
        |hash: &[u8; 32]| -> String { hash.iter().map(|byte| format!("{byte:02x}")).collect() };
    match action {
        ReferenceAction::Build {
            fasta,
            output,
            manifest,
            force,
        } => {
            use rosalind::provenance::CommandCapture;
            use rosalind::util::atomic::write_atomic;

            let receipt_path = manifest.unwrap_or_else(|| sidecar_path(&output, ".manifest.json"));
            require_safe_cli_destination(&receipt_path, force, "reference-pack receipt");
            let metadata = ReferencePackBuilder::build(&fasta, &output, force)
                .with_context(|| format!("failed to build {}", output.display()))?;
            let mut receipt = new_run_manifest("reference build");
            let mut command = CommandCapture::from_argv_prefix(["reference", "build"]);
            command.input("--fasta", &fasta)?;
            command.flag_if(force, "--force");
            command.output("--output", &output)?;
            command.record_into(&mut receipt);
            receipt
                .params
                .insert("artifact.input.0.role".into(), "reference-fasta".into());
            receipt.params.insert(
                "artifact.output.0.role".into(),
                "analysis-reference-pack".into(),
            );
            receipt
                .params
                .insert("artifact.output.0.format".into(), "rref-v1".into());
            receipt.params.insert(
                "reference.source_blake3".into(),
                render_hash(&metadata.source_reference_blake3),
            );
            receipt.params.insert(
                "reference.total_bases".into(),
                metadata.total_bases.to_string(),
            );
            receipt.record_measurement("peak_rss_bytes", peak_rss_bytes().to_string());
            receipt.finalize();
            write_atomic(&receipt_path, receipt.to_canonical_json().as_bytes(), force)?;
            println!(
                "built {}: {} contigs, {} bases, source {}",
                output.display(),
                metadata.contig_count,
                metadata.total_bases,
                render_hash(&metadata.source_reference_blake3)
            );
        }
        ReferenceAction::Inspect {
            reference_pack,
            json,
        } => {
            let reader = ReferencePackReader::open(&reference_pack)
                .with_context(|| format!("failed to inspect {}", reference_pack.display()))?;
            let metadata = reader.metadata();
            if json {
                println!(
                    "{{\"path\":\"{}\",\"format\":\"rref\",\"format_version\":{},\"contig_count\":{},\"total_bases\":{},\"source_reference_blake3\":\"{}\",\"content_blake3\":\"{}\"}}",
                    reference_pack.display(),
                    metadata.format_version,
                    metadata.contig_count,
                    metadata.total_bases,
                    render_hash(&metadata.source_reference_blake3),
                    render_hash(&metadata.content_blake3)
                );
            } else {
                println!("reference pack: {}", reference_pack.display());
                println!("  format version : {}", metadata.format_version);
                println!(
                    "  contigs / bases: {} / {}",
                    metadata.contig_count, metadata.total_bases
                );
                println!(
                    "  source BLAKE3  : {}",
                    render_hash(&metadata.source_reference_blake3)
                );
                for contig in reader.contigs().iter() {
                    println!("  {}\t{}", contig.name, contig.length);
                }
            }
        }
        ReferenceAction::Convert {
            index,
            output,
            manifest,
            force,
        } => {
            use rosalind::provenance::CommandCapture;
            use rosalind::util::atomic::write_atomic;

            let receipt_path = manifest.unwrap_or_else(|| sidecar_path(&output, ".manifest.json"));
            require_safe_cli_destination(&receipt_path, force, "reference-pack receipt");
            let metadata = ReferencePackBuilder::convert(&index, &output, force)
                .with_context(|| format!("failed to convert {}", index.display()))?;
            let mut receipt = new_run_manifest("reference convert");
            let mut command = CommandCapture::from_argv_prefix(["reference", "convert"]);
            command.input("--index", &index)?;
            command.flag_if(force, "--force");
            command.output("--output", &output)?;
            command.record_into(&mut receipt);
            receipt.params.insert(
                "artifact.input.0.role".into(),
                "legacy-reference-index".into(),
            );
            receipt.params.insert(
                "artifact.output.0.role".into(),
                "analysis-reference-pack".into(),
            );
            receipt
                .params
                .insert("artifact.output.0.format".into(), "rref-v1".into());
            receipt.params.insert(
                "reference.source_blake3".into(),
                render_hash(&metadata.source_reference_blake3),
            );
            receipt.params.insert(
                "reference.total_bases".into(),
                metadata.total_bases.to_string(),
            );
            receipt.record_measurement("peak_rss_bytes", peak_rss_bytes().to_string());
            receipt.finalize();
            write_atomic(&receipt_path, receipt.to_canonical_json().as_bytes(), force)?;
            println!(
                "converted {} to {}: {} contigs, {} bases",
                index.display(),
                output.display(),
                metadata.contig_count,
                metadata.total_bases
            );
        }
    }
    Ok(())
}

fn run_new_analyzer(name: &str, output: &std::path::Path) -> Result<()> {
    let report = rosalind::scaffold::create_analyzer_project(name, output)
        .with_context(|| format!("failed to scaffold analyzer {name:?}"))?;
    println!("created analyzer project: {}", report.root.display());
    for path in report.files {
        println!("  {}", path.display());
    }
    println!("next: cd {} && cargo test", report.root.display());
    Ok(())
}

fn run_conformance_analyzer(binary: PathBuf, json: bool) -> Result<()> {
    let report = rosalind::conform_analyzer(&binary)?;
    if json {
        println!("{}", report.to_json());
    } else {
        println!(
            "analyzer conformance: {}",
            if report.passed { "PASS" } else { "FAIL" }
        );
        for (name, passed) in &report.checks {
            println!("  {} {name}", if *passed { "✓" } else { "✗" });
        }
        for failure in &report.failures {
            eprintln!("  {failure}");
        }
    }
    if report.passed {
        Ok(())
    } else {
        std::process::exit(1);
    }
}

fn run_doctor_command(
    index: PathBuf,
    alignments: PathBuf,
    output: Option<PathBuf>,
    budget_mb: Option<u64>,
    deep: bool,
    json: bool,
    selection_args: SelectionArgs,
) -> Result<()> {
    let selection = resolve_selection(&index, &selection_args)?;
    let report = rosalind::doctor::run_doctor_selected(
        &rosalind::DoctorSpec {
            index,
            alignments,
            output,
            budget_mb,
            deep,
        },
        &selection,
    );
    if json {
        println!("{}", report.to_json());
    } else {
        println!(
            "doctor: {}",
            if report.ok {
                "READY"
            } else {
                "ACTION REQUIRED"
            }
        );
        println!(
            "  index / alignments : {} / {}",
            status_word(report.index_readable),
            status_word(report.alignments_readable)
        );
        println!(
            "  declared sort      : {}",
            report.declared_sort_order.as_deref().unwrap_or("unknown")
        );
        if let Some(proven) = report.coordinate_order_proven {
            println!("  deep order proof   : {}", status_word(proven));
        }
        println!(
            "  predicted memory   : {} MiB (minimum budget {} MiB)",
            report.predicted_peak_rss_bytes.div_ceil(1 << 20),
            report.required_budget_mb
        );
        println!("  assurance available: {}", report.available_assurance);
        for issue in &report.issues {
            eprintln!("  issue: {issue}");
        }
        for action in &report.remediation {
            println!("  next: {action}");
        }
    }
    if report.ok {
        Ok(())
    } else {
        std::process::exit(2);
    }
}

fn status_word(value: bool) -> &'static str {
    if value {
        "ok"
    } else {
        "failed"
    }
}

fn sidecar_path(path: &std::path::Path, suffix: &str) -> PathBuf {
    let mut value = path.as_os_str().to_os_string();
    value.push(suffix);
    PathBuf::from(value)
}

fn new_run_manifest(subcommand: &str) -> rosalind::provenance::RunManifest {
    let mut manifest = rosalind::provenance::RunManifest::new(subcommand);
    manifest.tool_version = env!("CARGO_PKG_VERSION").to_string();
    manifest
}

fn require_safe_cli_destination(path: &std::path::Path, force: bool, kind: &str) {
    if !force && path.exists() {
        eprintln!(
            "{kind} already exists: {} (choose a new path or pass --force for atomic replacement)",
            path.display()
        );
        std::process::exit(2);
    }
}

fn run_demo(output_dir: PathBuf, budget_mb: u64, json: bool) -> Result<()> {
    use rosalind::provenance::RunManifest;

    if output_dir.exists() && std::fs::read_dir(&output_dir)?.next().is_some() {
        bail!(
            "demo output directory is not empty: {} (choose another --output-dir)",
            output_dir.display()
        );
    }
    std::fs::create_dir_all(&output_dir)?;
    let reference = output_dir.join("reference.fa");
    let reads = output_dir.join("reads.fastq");
    let raw_bam = output_dir.join("raw.bam");
    let sorted_bam = output_dir.join("sorted.bam");
    let index = output_dir.join("ref.idx");
    let calls = output_dir.join("calls.vcf");
    let manifest = output_dir.join("calls.vcf.manifest.json");
    rosalind::util::atomic::write_atomic(
        &reference,
        include_str!("../assets/demo/reference.fa").as_bytes(),
        false,
    )?;
    rosalind::util::atomic::write_atomic(
        &reads,
        include_str!("../assets/demo/reads.fastq").as_bytes(),
        false,
    )?;

    let exe = std::env::current_exe().context("locating the rosalind executable")?;
    let run = |step: &str, args: &[String]| -> Result<std::process::Output> {
        if !json {
            println!("demo: {step}");
        }
        let output = std::process::Command::new(&exe)
            .args(args)
            .output()
            .with_context(|| format!("demo step {step} could not start"))?;
        if !json {
            print!("{}", String::from_utf8_lossy(&output.stdout));
            eprint!("{}", String::from_utf8_lossy(&output.stderr));
        }
        if !output.status.success() {
            bail!(
                "demo step {step} failed with exit {:?}: {}",
                output.status.code(),
                String::from_utf8_lossy(&output.stderr).trim()
            );
        }
        Ok(output)
    };
    let p = |path: &std::path::Path| path.display().to_string();

    run(
        "build the portable index",
        &[
            "index".into(),
            "--reference".into(),
            p(&reference),
            "--output".into(),
            p(&index),
        ],
    )?;
    run(
        "align the embedded reads",
        &[
            "align".into(),
            "--reference".into(),
            p(&reference),
            "--reads".into(),
            p(&reads),
            "--format".into(),
            "bam".into(),
            "--output".into(),
            p(&raw_bam),
        ],
    )?;
    run(
        "sort deterministically",
        &[
            "sort".into(),
            "--input".into(),
            p(&raw_bam),
            "--output".into(),
            p(&sorted_bam),
        ],
    )?;
    run(
        "predict before running",
        &[
            "plan".into(),
            "--index".into(),
            p(&index),
            "--budget-mb".into(),
            budget_mb.to_string(),
        ],
    )?;
    run(
        "honor the declared budget",
        &[
            "variants".into(),
            "--index".into(),
            p(&index),
            "--alignments".into(),
            p(&sorted_bam),
            "--memory-budget-mb".into(),
            budget_mb.to_string(),
            "--enforce".into(),
            "-o".into(),
            p(&calls),
        ],
    )?;
    run(
        "verify the receipt and artifacts",
        &[
            "verify".into(),
            "--manifest".into(),
            p(&manifest),
            "--json".into(),
        ],
    )?;
    run(
        "reproduce the result byte-for-byte",
        &[
            "reproduce".into(),
            "--manifest".into(),
            p(&manifest),
            "--inputs".into(),
            p(&output_dir),
            "--json".into(),
        ],
    )?;
    run(
        "verify the local provenance chain",
        &[
            "chain".into(),
            "verify".into(),
            p(&output_dir),
            "--json".into(),
        ],
    )?;

    let receipt = RunManifest::from_canonical_json(&std::fs::read_to_string(&manifest)?)
        .map_err(|error| anyhow!("demo receipt could not be parsed: {error}"))?;
    let claim = receipt.content_hash();
    let reproduction = sidecar_path(&manifest, ".repro.json");
    let studio_command = format!(
        "rosalind studio {} {} {} {} {}",
        p(&output_dir.join("ref.idx.manifest.json")),
        p(&output_dir.join("raw.bam.manifest.json")),
        p(&output_dir.join("sorted.bam.manifest.json")),
        p(&manifest),
        p(&reproduction),
    );
    if json {
        let escape = |s: &str| s.replace('\\', "\\\\").replace('"', "\\\"");
        println!(
            "{{\"schema\":2,\"ok\":true,\"trust\":\"reproduced\",\"claim\":\"{}\",\"output_dir\":\"{}\",\"index\":\"{}\",\"alignments\":\"{}\",\"output\":\"{}\",\"manifest\":\"{}\",\"reproduction_certificate\":\"{}\",\"studio_command\":\"{}\"}}",
            claim,
            escape(&p(&output_dir)),
            escape(&p(&index)),
            escape(&p(&sorted_bam)),
            escape(&p(&calls)),
            escape(&p(&manifest)),
            escape(&p(&reproduction)),
            escape(&studio_command),
        );
    } else {
        println!(
            "demo: COMPLETE — reproduced · fits {budget_mb} MiB · claim {}",
            &claim[..claim.len().min(10)]
        );
        println!("demo: artifacts {}", output_dir.display());
        println!("demo: {studio_command}");
    }
    Ok(())
}

/// Compare a called VCF against a truth VCF over a reference, optionally masked by
/// a BED. VCF-agnostic — used by both `eval-somatic` and `eval-germline` (and the
/// drop-in interface for a real GIAB germline benchmark).
fn run_eval(
    reference_path: PathBuf,
    calls_path: PathBuf,
    truth_path: PathBuf,
    regions_path: Option<PathBuf>,
    calls_filter: CallsFilter,
    json: bool,
) -> Result<()> {
    let references = read_fasta_map(&reference_path)
        .with_context(|| format!("failed to read reference from {}", reference_path.display()))?;

    let calls_txt = std::fs::read_to_string(&calls_path)
        .with_context(|| format!("failed to read calls VCF {}", calls_path.display()))?;
    let truth_txt = std::fs::read_to_string(&truth_path)
        .with_context(|| format!("failed to read truth VCF {}", truth_path.display()))?;

    let mut calls = read_vcf_variants(&calls_txt)
        .with_context(|| format!("failed to parse calls VCF {}", calls_path.display()))?;
    if calls_filter == CallsFilter::Pass {
        calls.retain(|call| call.filter == "PASS");
    }
    let truth = read_vcf_variants(&truth_txt)
        .with_context(|| format!("failed to parse truth VCF {}", truth_path.display()))?;

    let bed = if let Some(path) = regions_path {
        let bed_txt = std::fs::read_to_string(&path)
            .with_context(|| format!("failed to read BED {}", path.display()))?;
        Some(BedIndex::from_str(&bed_txt).with_context(|| "failed to parse BED")?)
    } else {
        None
    };

    let report = compare_callsets(&references, &calls, &truth, bed.as_ref())?;
    let (p, r) = (report.precision(), report.recall());
    let f1 = if p + r == 0.0 {
        0.0
    } else {
        2.0 * p * r / (p + r)
    };
    if json {
        println!(
            "{{\"schema\":1,\"calls_filter\":\"{}\",\"truth_total\":{},\"calls_total\":{},\"tp\":{},\"fp\":{},\"fn\":{},\"precision\":{p:.9},\"recall\":{r:.9},\"f1\":{f1:.9},\"genotype_concordant\":{},\"genotype_discordant\":{},\"genotype_unknown\":{},\"genotype_concordance\":{:.9}}}",
            match calls_filter {
                CallsFilter::All => "all",
                CallsFilter::Pass => "pass",
            },
            report.total_truth,
            report.total_calls,
            report.true_positive,
            report.false_positive,
            report.false_negative,
            report.genotype_concordant,
            report.genotype_discordant,
            report.genotype_unknown,
            report.genotype_concordance(),
        );
    } else {
        println!("truth_total={}", report.total_truth);
        println!("calls_total={}", report.total_calls);
        println!("tp={}", report.true_positive);
        println!("fp={}", report.false_positive);
        println!("fn={}", report.false_negative);
        println!("precision={p:.6}");
        println!("recall={r:.6}");
        println!("f1={f1:.6}");
        println!("genotype_concordant={}", report.genotype_concordant);
        println!("genotype_discordant={}", report.genotype_discordant);
        println!("genotype_unknown={}", report.genotype_unknown);
        println!("genotype_concordance={:.6}", report.genotype_concordance());
        for (ty, (tp, fp, fn_)) in report.by_type.iter() {
            println!("type={:?} tp={} fp={} fn={}", ty, tp, fp, fn_);
        }
    }
    Ok(())
}

/// Build a multi-contig index from a FASTA and persist it (B3c).
fn run_index(
    reference: PathBuf,
    output: PathBuf,
    memory_budget_mb: Option<u64>,
    force: bool,
) -> Result<()> {
    use rosalind::util::atomic::{write_atomic, AtomicFile};

    let receipt_path = sidecar_path(&output, ".manifest.json");
    require_safe_cli_destination(&output, force, "index output");
    require_safe_cli_destination(&receipt_path, force, "index receipt");
    let mut atomic_output = AtomicFile::create(&output)
        .with_context(|| format!("failed to reserve index output {}", output.display()))?;
    atomic_output.close_for_path_writer();
    // Read every FASTA record (all contigs) into (name, sequence) pairs.
    let fasta_reader = open_input(&reference)
        .with_context(|| format!("failed to open reference {}", reference.display()))?;
    let records: Vec<FastaRecord> = FastaReader::new(fasta_reader)
        .collect::<std::result::Result<Vec<_>, _>>()
        .with_context(|| format!("failed to parse FASTA {}", reference.display()))?;
    if records.is_empty() {
        bail!(
            "reference {} contains no FASTA records",
            reference.display()
        );
    }
    let total_bp: u64 = records.iter().map(|r| r.sequence.len() as u64).sum();
    // The code-grounded build-memory model the plan line AND the build receipt share
    // (anti-drift): the up-front prediction and the realized accounting use one model.
    let model = BuildMemoryModel::from_reference_len(total_bp);

    // Record-only budget plan line, printed BEFORE the build. Never refuses.
    if let Some(mb) = memory_budget_mb {
        eprintln!(
            "{}",
            render_plan_line(model.working_set(), MemoryBudget::from_mb(mb))
        );
    }

    let named: Vec<(String, Vec<u8>)> = records.into_iter().map(|r| (r.name, r.sequence)).collect();
    let index = GenomeIndex::from_named_sequences(&named)
        .with_context(|| format!("failed to build index from {}", reference.display()))?;

    IndexWriter::create(atomic_output.temporary_path())
        .with_context(|| format!("failed to create index file {}", output.display()))?
        .write_genome_index(&index)
        .with_context(|| format!("failed to write index to {}", output.display()))?;
    atomic_output
        .commit(force)
        .with_context(|| format!("failed to commit index {}", output.display()))?;

    // Deterministic build receipt → stdout.
    let index_bytes = std::fs::metadata(&output)
        .with_context(|| format!("failed to stat index file {}", output.display()))?
        .len();
    let reference_blake3 = *blake3::hash(index.reference()).as_bytes();
    let report = IndexBuildReport {
        index_path: output.display().to_string(),
        contigs: index
            .contigs()
            .iter()
            .map(|c| (c.name.to_string(), c.length))
            .collect(),
        total_bp,
        reference_blake3,
        index_bytes,
    };
    print!("{}", report.render());

    // Build receipt: realized peak RSS vs the modeled n-scale SA-IS build memory —
    // the D0 measure-first probe. The realized peak is machine-dependent; the
    // breakdown + attribution ratio are the analysis payload.
    let peak = peak_rss_bytes();
    let denom = total_bp.max(1);
    eprintln!(
        "build: realized peak RSS {} MiB ({} B/base) over {} bp",
        peak / (1 << 20),
        peak / denom,
        total_bp
    );
    eprint!("{}", model.render(total_bp));
    let ratio = if peak > 0 {
        model.total_bytes as f64 / peak as f64
    } else {
        0.0
    };
    eprintln!(
        "build: model/realized attribution = {:.2} [{}]",
        ratio,
        if ratio >= 0.70 {
            "CONFIRM ≥0.70"
        } else {
            "below 0.70"
        }
    );

    // Content-addressed receipt: makes the index a chainable, verifiable node — the
    // root of every downstream provenance chain. Mirrors the `variants`/`somatic` path.
    // The `--output` operand records blake3(.idx file) into outputs[]; that digest is
    // bit-identical to what `variants` records as its `--index` input, so the chain edge
    // resolves by construction. `--reference` records blake3(FASTA file) as the root.
    {
        use rosalind::provenance::{blake3_hex, CommandCapture};

        let mut manifest = new_run_manifest("index");
        let mut cmd = CommandCapture::new("index");
        cmd.input("--reference", &reference)?;
        if let Some(mb) = memory_budget_mb {
            cmd.opt("--memory-budget-mb", mb);
        }
        cmd.output("--output", &output)?;
        cmd.record_into(&mut manifest);
        manifest
            .params
            .insert("artifact.input.0.role".to_string(), "reference".to_string());
        manifest.params.insert(
            "artifact.output.0.role".to_string(),
            "reference-index".to_string(),
        );
        // `reference_blake3` is the in-memory NORMALIZED sequence hash (a "what genome"
        // id, stable across FASTA reformatting) — informational, NOT the chain edge.
        manifest.params.insert(
            "reference_blake3".to_string(),
            blake3_hex(&reference_blake3),
        );
        manifest
            .params
            .insert("total_bp".to_string(), total_bp.to_string());
        // Realized build peak is machine-dependent → a MEASUREMENT (relocated out of the
        // claim by finalize). The build is still O(reference) RAM; this records the cost,
        // it does NOT claim the budget was honored.
        manifest.record_measurement("peak_rss_bytes", peak.to_string());
        manifest.finalize();
        write_atomic(
            &receipt_path,
            manifest.to_canonical_json().as_bytes(),
            force,
        )
        .with_context(|| format!("failed to write index receipt for {}", output.display()))?;
        eprintln!("wrote reproducibility receipt: {}", receipt_path.display());
    }
    Ok(())
}

fn run_sort(
    input: PathBuf,
    output: PathBuf,
    memory_mb: usize,
    force: bool,
    manifest_out: Option<PathBuf>,
) -> Result<()> {
    use rosalind::provenance::CommandCapture;
    use rosalind::util::atomic::{write_atomic, AtomicFile};

    let receipt_path = manifest_out.unwrap_or_else(|| sidecar_path(&output, ".manifest.json"));
    require_safe_cli_destination(&output, force, "sorted BAM output");
    require_safe_cli_destination(&receipt_path, force, "sort receipt");
    let mut atomic_output = AtomicFile::create(&output)
        .with_context(|| format!("failed to reserve sorted BAM {}", output.display()))?;
    atomic_output.close_for_path_writer();
    let bytes = memory_mb.saturating_mul(1024 * 1024).max(1024 * 1024);
    sort_bam_deterministic(&input, atomic_output.temporary_path(), bytes)
        .context("sorting BAM failed")?;
    atomic_output
        .commit(force)
        .with_context(|| format!("failed to commit sorted BAM {}", output.display()))?;

    let mut receipt = new_run_manifest("sort");
    let mut command = CommandCapture::new("sort");
    command.input("--input", &input)?;
    command.opt("--memory-mb", memory_mb);
    command.flag_if(force, "--force");
    command.output("--output", &output)?;
    command.record_into(&mut receipt);
    receipt.params.insert(
        "artifact.input.0.role".to_string(),
        "raw-alignments".to_string(),
    );
    receipt.params.insert(
        "artifact.output.0.role".to_string(),
        "sorted-alignments".to_string(),
    );
    receipt.record_measurement("peak_rss_bytes", peak_rss_bytes().to_string());
    receipt.finalize();
    write_atomic(&receipt_path, receipt.to_canonical_json().as_bytes(), force)
        .with_context(|| format!("failed to write sort receipt {}", receipt_path.display()))?;
    eprintln!("wrote reproducibility receipt: {}", receipt_path.display());
    Ok(())
}

/// Predict whether a job fits a declared budget, before committing. `--index`
/// predicts the bounded whole-genome `variants` peak (largest contig + active set
/// @ the declared cap, atop the measured process baseline). `--reference`
/// predicts the index build peak (advisory; build is O(reference)).
#[allow(clippy::too_many_arguments)] // Each argument maps to a distinct public planning option.
fn run_plan(
    index: Option<PathBuf>,
    reference: Option<PathBuf>,
    max_depth: u32,
    max_read_len: u32,
    budget_mb: Option<u64>,
    json: bool,
    selection_args: SelectionArgs,
    format: FeatureFormat,
) -> Result<()> {
    use rosalind::call::plan::{predicted_peak_rss_bytes, render_variants_plan};
    use rosalind::genomics::{AnalysisReference, ReferenceProvider};

    if let Some(index_path) = index {
        let loaded = AnalysisReference::open(&index_path).with_context(|| {
            format!("failed to open analysis reference {}", index_path.display())
        })?;
        let selection = selection_from_provider(&loaded, &selection_args)?;
        let largest = selection.largest_reference_span(loaded.contigs());
        // Measure the process baseline now (binary + libs + index mmap header);
        // the per-contig reference decode + active set are modeled on top.
        let baseline = peak_rss_bytes();
        let encoder_bytes = if format == FeatureFormat::ArrowIpc {
            rosalind::call::arrow::feature_arrow_memory_bytes(
                loaded
                    .contigs()
                    .iter()
                    .map(|c| c.name.len())
                    .max()
                    .unwrap_or(0),
            )
        } else {
            0
        };
        if json {
            // Machine-readable: the predicted-peak fields a scheduler/CI reads.
            let predicted = predicted_peak_rss_bytes(largest, max_depth, max_read_len, baseline)
                .saturating_add(encoder_bytes);
            let verdict = match budget_mb {
                Some(mb) => {
                    if MemoryBudget::from_mb(mb).admits(predicted) {
                        "fits"
                    } else {
                        "refuse"
                    }
                }
                None => "none",
            };
            let budget_field = budget_mb
                .map(|mb| mb.to_string())
                .unwrap_or_else(|| "null".to_string());
            println!(
                "{{\"index\":\"{}\",\"predicted_peak_rss_bytes\":{},\"baseline_rss_bytes\":{},\
                 \"largest_contig_len\":{},\"max_depth\":{},\"max_read_len\":{},\
                 \"budget_mb\":{},\"verdict\":\"{}\",\"encoder_additional_bytes\":{}}}",
                index_path.display(),
                predicted,
                baseline,
                largest,
                max_depth,
                max_read_len,
                budget_field,
                verdict,
                encoder_bytes,
            );
        } else {
            print!(
                "{}",
                render_variants_plan(
                    largest,
                    max_depth,
                    max_read_len,
                    baseline.saturating_add(encoder_bytes),
                    budget_mb
                )
            );
            if encoder_bytes != 0 {
                println!(
                    "Arrow encoder additional working set: {encoder_bytes} bytes (included above)"
                );
            }
        }
    } else {
        if selection_requested(&selection_args) {
            bail!("analysis selection flags require --index or --reference-pack");
        }
        let reference = reference.expect("clap guarantees one of --index/--reference");
        let fasta_reader = open_input(&reference)
            .with_context(|| format!("failed to open reference {}", reference.display()))?;
        let total_bp: u64 = FastaReader::new(fasta_reader)
            .collect::<std::result::Result<Vec<_>, _>>()
            .with_context(|| format!("failed to parse FASTA {}", reference.display()))?
            .iter()
            .map(|r| r.sequence.len() as u64)
            .sum();
        // Predict the build peak from the code-grounded BuildMemoryModel (~45 B/base,
        // envelopes the D0-measured 41 B/base realized) — NOT the old 12 B/base estimate,
        // which under-predicted ~3.4x. Advisory: the build is still O(reference) until D1a.
        let model = BuildMemoryModel::from_reference_len(total_bp);
        print!("{}", model.render(total_bp));
        match budget_mb {
            Some(mb) => println!(
                "{}",
                render_plan_line(model.working_set(), MemoryBudget::from_mb(mb))
            ),
            None => println!(
                "plan: est. build peak ~{} MiB (advisory; build is O(reference)) [no budget]",
                model.total_bytes / (1 << 20)
            ),
        }
    }
    Ok(())
}

/// Pack many bounded `variants` jobs onto fixed-size nodes by their PREDICTED
/// peaks — each read from the job's index header (no run, no read I/O). Peaks are
/// conservative upper bounds and additive, so the printed schedule shows every
/// node within capacity by predicted peak before a single job launches — the
/// contract turned into a placement decision. Exits 3 when no safe packing exists.
fn run_pack(
    jobs_path: PathBuf,
    node_mb: u64,
    nodes: Option<usize>,
    default_max_depth: u32,
    default_max_read_len: u32,
    json: bool,
) -> Result<()> {
    use rosalind::call::plan::predicted_peak_rss_bytes;
    use rosalind::call::{first_fit_decreasing, PackJob, PackOutcome};
    use rosalind::genomics::IndexReader;

    let text = std::fs::read_to_string(&jobs_path)
        .with_context(|| format!("failed to read jobs file {}", jobs_path.display()))?;

    // A nominal per-process baseline (binary + libs + index header): each
    // co-located job runs in its own process with roughly this floor. Measured
    // once and applied per job, so the per-job predicted peaks are comparable.
    let baseline = peak_rss_bytes();

    let mut pack_jobs: Vec<PackJob> = Vec::new();
    for (lineno, raw) in text.lines().enumerate() {
        let line = raw.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut fields = line.split_whitespace();
        let index_path = fields.next().expect("a non-empty line has a first field");
        let max_depth = match fields.next() {
            Some(s) => s.parse::<u32>().with_context(|| {
                format!(
                    "jobs file {}:{}: invalid max_depth '{s}'",
                    jobs_path.display(),
                    lineno + 1
                )
            })?,
            None => default_max_depth,
        };
        let max_read_len = match fields.next() {
            Some(s) => s.parse::<u32>().with_context(|| {
                format!(
                    "jobs file {}:{}: invalid max_read_len '{s}'",
                    jobs_path.display(),
                    lineno + 1
                )
            })?,
            None => default_max_read_len,
        };
        let loaded = IndexReader::open(std::path::Path::new(index_path)).with_context(|| {
            format!(
                "jobs file {}:{}: failed to open index {index_path}",
                jobs_path.display(),
                lineno + 1
            )
        })?;
        let largest = loaded
            .contigs()
            .iter()
            .map(|c| c.length as u64)
            .max()
            .unwrap_or(0);
        let predicted = predicted_peak_rss_bytes(largest, max_depth, max_read_len, baseline);
        pack_jobs.push(PackJob {
            label: index_path.to_string(),
            predicted_peak_bytes: predicted,
        });
    }

    if pack_jobs.is_empty() {
        bail!("jobs file {} has no jobs", jobs_path.display());
    }

    let node_cap = node_mb.saturating_mul(1 << 20);
    const MIB: u64 = 1 << 20;
    match first_fit_decreasing(&pack_jobs, node_cap, nodes) {
        PackOutcome::Packed(assignments) => {
            if json {
                let mut s = format!("{{\"node_mb\":{node_mb},\"nodes\":[");
                for (i, n) in assignments.iter().enumerate() {
                    if i > 0 {
                        s.push(',');
                    }
                    let labels = n
                        .job_labels
                        .iter()
                        .map(|l| format!("\"{l}\""))
                        .collect::<Vec<_>>()
                        .join(",");
                    s.push_str(&format!(
                        "{{\"node\":{},\"used_bytes\":{},\"jobs\":[{}]}}",
                        n.node, n.used_bytes, labels
                    ));
                }
                s.push_str("]}");
                println!("{s}");
            } else {
                println!(
                    "pack: {} job(s) → {} node(s) of {} MiB — every node within capacity by predicted peak",
                    pack_jobs.len(),
                    assignments.len(),
                    node_mb
                );
                for n in &assignments {
                    println!(
                        "  node {}: {} / {} MiB  [{}]",
                        n.node,
                        n.used_bytes / MIB,
                        node_mb,
                        n.job_labels.join(", ")
                    );
                }
            }
            Ok(())
        }
        PackOutcome::NoFit { reason } => {
            eprintln!("pack: REFUSE — {reason}");
            std::process::exit(3);
        }
    }
}

/// Re-check a reproducibility receipt without re-running: parse it, re-hash each
/// listed input/output and confirm the digests match, and confirm the recorded
/// realized peak RSS landed within the budget (supplied, or recorded in the
/// manifest). Exits non-zero with a per-check report on any mismatch.
fn run_verify(
    manifest_path: PathBuf,
    budget_mb: Option<u64>,
    expect_code: Option<String>,
    json: bool,
) -> Result<()> {
    use rosalind::provenance::{verify_receipt, VerifyOpts};

    let text = std::fs::read_to_string(&manifest_path)
        .with_context(|| format!("failed to read manifest {}", manifest_path.display()))?;
    let report = verify_receipt(
        &text,
        &VerifyOpts {
            budget_mb,
            expect_code,
            rehash_files: true,
        },
    );
    if json {
        println!("{}", report.to_json());
    } else {
        for n in &report.notes {
            println!("verify: {n}");
        }
    }
    if report.ok {
        if !json {
            let m = report
                .manifest
                .as_ref()
                .expect("a passing report always carries the parsed manifest");
            println!(
                "verify: OK — {} input(s), {} output(s) match",
                m.inputs.len(),
                m.outputs.len()
            );
        }
        Ok(())
    } else {
        if !json {
            for p in &report.problems {
                eprintln!("verify: FAIL — {p}");
            }
        }
        std::process::exit(5);
    }
}

fn run_receipt_inspect(
    manifest: PathBuf,
    artifacts: Vec<PathBuf>,
    certificate: Option<PathBuf>,
    json: bool,
) -> Result<()> {
    use rosalind::provenance::TrustState;

    let report = rosalind::inspect_receipt(&manifest, &artifacts, certificate.as_deref())?;
    if json {
        println!("{}", report.to_json());
    } else {
        let facets = [
            ("receipt integrity", &report.trust.receipt_integrity),
            ("artifact completeness", &report.trust.artifact_completeness),
            ("resource contract", &report.trust.resource_contract),
            ("reproduction evidence", &report.trust.reproduction_evidence),
            ("signature", &report.trust.signature),
        ];
        println!(
            "receipt: {}",
            report
                .trust
                .claim_id
                .as_deref()
                .map(|claim| &claim[..claim.len().min(12)])
                .unwrap_or("unavailable")
        );
        for (name, facet) in facets {
            println!("  {name}: {}", facet.status);
            println!("    why: {}", facet.explanation);
        }
        if !report.missing_artifacts.is_empty() {
            println!("  missing content hashes:");
            for hash in &report.missing_artifacts {
                println!("    {hash}");
            }
        }
    }
    let failed = [
        &report.trust.receipt_integrity,
        &report.trust.artifact_completeness,
        &report.trust.resource_contract,
        &report.trust.reproduction_evidence,
    ]
    .iter()
    .any(|facet| facet.state == TrustState::Failed);
    if failed {
        std::process::exit(5);
    }
    Ok(())
}

fn run_receipt_sanitize(manifest: PathBuf, output: PathBuf) -> Result<()> {
    let claim = rosalind::sanitize_receipt(&manifest, &output)?;
    eprintln!("wrote sanitized receipt: {}", output.display());
    eprintln!("claim unchanged: {}", &claim[..claim.len().min(12)]);
    eprintln!(
        "warning: claim-bearing extension parameters were retained and may still contain sensitive text"
    );
    Ok(())
}

fn run_receipt_export_intoto(manifest: PathBuf, output: PathBuf) -> Result<()> {
    let claim = rosalind::export_intoto(&manifest, &output)?;
    eprintln!("wrote unsigned in-toto Statement v1: {}", output.display());
    eprintln!("native claim: {}", &claim[..claim.len().min(12)]);
    eprintln!("note: this export is neither signed nor a claim of SLSA compliance");
    Ok(())
}

/// Walk a directory of receipts as a provenance DAG: every node self-hashes and every
/// expected-internal edge (`--index`) resolves by content hash. External inputs
/// (reads/alignments/reference FASTA) are integrity-verified, not byte-reproduced.
/// Exit 0 = CHAIN INTACT, 5 = CHAIN BROKEN (matching `verify`).
fn run_chain_verify(dir: PathBuf, json: bool) -> Result<()> {
    use rosalind::provenance::{walk_chain, EdgeStatus, RunManifest};
    use std::collections::HashMap;

    let mut receipts: Vec<RunManifest> = Vec::new();
    for entry in std::fs::read_dir(&dir)
        .with_context(|| format!("failed to read directory {}", dir.display()))?
    {
        let path = entry?.path();
        let is_manifest = path
            .file_name()
            .and_then(|n| n.to_str())
            .map(|n| n.ends_with(".manifest.json"))
            .unwrap_or(false);
        if !is_manifest {
            continue;
        }
        let text = std::fs::read_to_string(&path)
            .with_context(|| format!("failed to read {}", path.display()))?;
        // Skip anything that is not a parseable run manifest (e.g. a stray file).
        if let Ok(m) = RunManifest::from_canonical_json(&text) {
            receipts.push(m);
        }
    }
    if receipts.is_empty() {
        bail!("no receipts (*.manifest.json) found in {}", dir.display());
    }
    // Stable order (independent of read_dir) so the report is deterministic.
    receipts.sort_by_key(|m| m.content_hash());

    let report = walk_chain(&receipts);

    if json {
        println!("{}", report.to_json());
    } else {
        let sub: HashMap<&str, &str> = report
            .nodes
            .iter()
            .map(|n| (n.id.as_str(), n.subcommand.as_str()))
            .collect();
        let tampered = report
            .nodes
            .iter()
            .filter(|n| n.self_hash == Some(false))
            .count();
        if tampered == 0 {
            println!("chain: {} nodes, all self-hash OK", report.nodes.len());
        } else {
            println!(
                "chain: {} nodes, {tampered} TAMPERED (self-hash mismatch)",
                report.nodes.len()
            );
        }
        for e in &report.edges {
            match &e.status {
                EdgeStatus::Resolved { parent_id } => {
                    let p = sub.get(parent_id.as_str()).copied().unwrap_or("?");
                    println!(
                        "edge: {}  {}-->  {p}  [resolved]",
                        e.child_subcommand, e.flag
                    );
                }
                EdgeStatus::External => {
                    println!(
                        "edge: {}  {}-->  (external)  [integrity-only]",
                        e.child_subcommand, e.flag
                    );
                }
                EdgeStatus::Broken => {
                    println!(
                        "edge: {}  {}-->  (unresolved)  [BROKEN]",
                        e.child_subcommand, e.flag
                    );
                }
            }
        }
        let resolved = report
            .edges
            .iter()
            .filter(|e| matches!(e.status, EdgeStatus::Resolved { .. }))
            .count();
        if report.intact {
            println!(
                "VERDICT: CHAIN INTACT ({} nodes, {resolved} internal edges resolve)",
                report.nodes.len()
            );
        } else {
            println!("VERDICT: CHAIN BROKEN");
        }
    }

    if report.intact {
        Ok(())
    } else {
        std::process::exit(5);
    }
}

/// Localize how two receipts' claims differ (claim-level only; it does not re-derive or
/// diff output bytes). Exit 0 = identical, 1 = differ, 2 = read/parse error.
fn run_diff(a: PathBuf, b: PathBuf, json: bool, loci_output: Option<PathBuf>) -> Result<()> {
    use rosalind::provenance::{diff_receipts, RunManifest};

    let parse = |p: &PathBuf| -> RunManifest {
        let text = match std::fs::read_to_string(p) {
            Ok(t) => t,
            Err(e) => {
                eprintln!("diff: cannot read {}: {e}", p.display());
                std::process::exit(2);
            }
        };
        match RunManifest::from_canonical_json(&text) {
            Ok(m) => m,
            Err(e) => {
                eprintln!("diff: cannot parse {}: {e}", p.display());
                std::process::exit(2);
            }
        }
    };
    let ma = parse(&a);
    let mb = parse(&b);
    if let Some(output) = &loci_output {
        match rosalind::dataset_diff::diff_evidence_datasets(&a, &b, output) {
            Ok(summary) => eprintln!(
                "locus deltas: {} changed, {} added, {} removed; {} metric changes written to {}",
                summary.changed_loci,
                summary.added_loci,
                summary.removed_loci,
                summary.metric_changes,
                output.display()
            ),
            Err(error) => {
                eprintln!("diff: {error}");
                std::process::exit(error.exit_code());
            }
        }
    }
    let report = diff_receipts(&ma, &mb);

    let short = |o: &Option<String>| -> String {
        match o {
            Some(s) if s.len() > 12 => format!("{}…", &s[..12]),
            Some(s) => s.clone(),
            None => "(absent)".to_string(),
        }
    };

    if json {
        println!("{}", report.to_json());
    } else {
        println!("diff: {}  vs  {}", ma.subcommand, mb.subcommand);
        for c in &report.code_identity {
            println!(
                "CAUSE  — code-identity {}  {} → {}",
                c.key,
                short(&c.a),
                short(&c.b)
            );
        }
        for c in &report.inputs {
            println!(
                "CAUSE  — input  {}  {} → {}",
                c.flag,
                short(&c.a),
                short(&c.b)
            );
        }
        for c in &report.science_params {
            println!(
                "CAUSE  — param  {}  {} → {}",
                c.key,
                short(&c.a),
                short(&c.b)
            );
        }
        for c in &report.execution_params {
            println!(
                "EXECUTION — setting {}  {} → {}",
                c.key,
                short(&c.a),
                short(&c.b)
            );
        }
        for c in &report.outputs {
            println!(
                "EFFECT — output {}  {} → {}",
                c.flag,
                short(&c.a),
                short(&c.b)
            );
        }
        for c in &report.measurements {
            println!(
                "noise  — measurement {}  {} → {}",
                c.key,
                short(&c.a),
                short(&c.b)
            );
        }
        println!("VERDICT: {}", report.verdict());
    }
    std::process::exit(report.exit_code());
}

/// Re-derive a recorded result and report REPRODUCED / DIVERGED / INCONCLUSIVE (and
/// TAMPERED for a modified receipt), exiting 0 / 6 / 7 / 5 respectively. Writes a
/// `.repro.json` reproduction certificate next to the receipt unless `--no-attest`.
#[allow(clippy::too_many_arguments)] // CLI entry point: each replay safety switch is independent.
fn run_reproduce(
    manifest: PathBuf,
    inputs: PathBuf,
    no_attest: bool,
    output: Option<PathBuf>,
    binary: Option<PathBuf>,
    dry_run: bool,
    force: bool,
    json: bool,
) -> Result<()> {
    use rosalind::provenance::{ReproOutput, ReproReceipt};

    if dry_run {
        match rosalind::reproduce::plan_reproduction(&manifest, &inputs, binary.as_deref()) {
            Ok(plan) => {
                if json {
                    println!("{}", plan.to_json());
                } else {
                    println!("reproduction plan: validated");
                    println!("  binary: {}", plan.binary.display());
                    println!("  argv: {}", plan.argv.join(" "));
                    println!(
                        "  inputs: {}; outputs: {}",
                        plan.inputs.len(),
                        plan.outputs.len()
                    );
                }
                std::fs::remove_dir_all(&plan.work_dir).ok();
                return Ok(());
            }
            Err(error) => {
                eprintln!("reproduction plan rejected: {error}");
                std::process::exit(7);
            }
        }
    }

    let report = rosalind::reproduce::reproduce_with_binary(&manifest, &inputs, binary.as_deref())?;
    if json {
        println!("{}", report.to_json());
    } else {
        for line in &report.lines {
            println!("{line}");
        }
    }

    // Mint a chainable reproduction certificate when a real comparison happened
    // (REPRODUCED / DIVERGED) — unless the caller opted out.
    if report.compared && !no_attest {
        let outputs: Vec<ReproOutput> = report
            .outputs
            .iter()
            .map(|c| ReproOutput {
                role: c.role.clone(),
                recorded_blake3: c.recorded.clone(),
                observed_blake3: c.observed.clone(),
                matched: c.matched,
            })
            .collect();
        let mut cert = ReproReceipt::build_with_identities(
            &report.parent_claim,
            &report.parent_subcommand,
            &report.verdict_label,
            1,
            &outputs,
            report.resource_here.peak_rss_bytes,
            report.resource_here.declared_budget_mb,
            &report.original_code,
            &report.reproducer_code,
        );
        cert.set_tool_version(env!("CARGO_PKG_VERSION"));
        let dest = match &output {
            Some(p) => p.clone(),
            None => {
                let mut s = manifest.as_os_str().to_os_string();
                s.push(".repro.json");
                PathBuf::from(s)
            }
        };
        require_safe_cli_destination(&dest, force, "reproduction certificate");
        rosalind::util::atomic::write_atomic(&dest, cert.to_canonical_json().as_bytes(), force)
            .with_context(|| {
                format!(
                    "failed to write reproduction certificate {}",
                    dest.display()
                )
            })?;
        if !json {
            println!(
                "  -> wrote reproduction certificate: {} (chains to {})",
                dest.display(),
                &report.parent_claim[..report.parent_claim.len().min(10)]
            );
        }
    }

    if report.exit_code != 0 {
        std::process::exit(report.exit_code);
    }
    Ok(())
}

/// Emit a self-hosted status badge for a run without conflating an intact receipt with
/// evidence that another execution reproduced its output bytes.
fn run_badge(
    manifest: PathBuf,
    repro: Option<PathBuf>,
    output: PathBuf,
    force: bool,
) -> Result<()> {
    use rosalind::provenance::{
        badge_json_for, badge_svg_for, ArtifactEvidence, BadgeStatus, CertificateEvidence,
        ReproReceipt, RunManifest, TrustReport, TrustState,
    };

    let text = std::fs::read_to_string(&manifest)
        .with_context(|| format!("failed to read receipt {}", manifest.display()))?;
    let m = RunManifest::from_canonical_json(&text)
        .map_err(|e| anyhow!("failed to parse receipt {}: {e}", manifest.display()))?;

    let certificate = repro.as_ref().map(|path| {
        std::fs::read_to_string(path)
            .map_err(|error| format!("cannot read {}: {error}", path.display()))
            .and_then(|text| {
                ReproReceipt::from_canonical_json(&text)
                    .map_err(|error| format!("cannot parse {}: {error}", path.display()))
            })
    });
    let certificate_evidence = match &certificate {
        None => CertificateEvidence::NotSupplied,
        Some(Ok(certificate)) => CertificateEvidence::Parsed(certificate),
        Some(Err(error)) => CertificateEvidence::Invalid(error.clone()),
    };
    let trust = TrustReport::evaluate(&m, ArtifactEvidence::NotChecked, certificate_evidence);
    let status = BadgeStatus::from_trust(&trust);
    let fits_mb = (trust.resource_contract.state == TrustState::Satisfied)
        .then_some(trust.budget_mb)
        .flatten();

    let is_svg = output
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.eq_ignore_ascii_case("svg"))
        .unwrap_or(false);
    let body = if is_svg {
        badge_svg_for(status, fits_mb)
    } else {
        badge_json_for(status, fits_mb)
    };
    require_safe_cli_destination(&output, force, "badge output");
    rosalind::util::atomic::write_atomic(&output, body.as_bytes(), force)
        .with_context(|| format!("failed to write badge {}", output.display()))?;
    eprintln!("wrote badge: {}", output.display());
    Ok(())
}

/// Load a prebuilt index and print exact-match loci for `pattern` (B3c). This is
/// a memory-mapped load + exact match — it never rebuilds the index.
fn run_locate(index: PathBuf, pattern: String, max_hits: usize) -> Result<()> {
    let loaded = IndexReader::open(&index)
        .with_context(|| format!("failed to open index {}", index.display()))?;
    let view = loaded
        .genome_view()
        .with_context(|| format!("failed to view index {}", index.display()))?;

    let loci = view.locate_exact(pattern.as_bytes(), max_hits);
    if loci.is_empty() {
        eprintln!("no hits");
        return Ok(());
    }

    let mut stdout = io::BufWriter::new(io::stdout().lock());
    for locus in loci {
        let name = view
            .contigs()
            .by_id(locus.contig)
            .map(|c| c.name.to_string())
            .unwrap_or_else(|| locus.contig.to_string());
        writeln!(stdout, "{name}\t{}", locus.pos.0)?;
    }
    stdout.flush()?;
    Ok(())
}

#[allow(clippy::too_many_arguments)] // a CLI entry point: each flag is a parameter
fn run_somatic(
    reference_path: PathBuf,
    tumor_fastq: Option<PathBuf>,
    tumor_r1: Option<PathBuf>,
    tumor_r2: Option<PathBuf>,
    normal_fastq: Option<PathBuf>,
    normal_r1: Option<PathBuf>,
    normal_r2: Option<PathBuf>,
    output_vcf: PathBuf,
    workdir: Option<PathBuf>,
    memory_mb: usize,
    force: bool,
) -> Result<()> {
    use rosalind::util::atomic::{write_atomic, AtomicFile};

    let start_total = Instant::now();
    let fasta = read_fasta(&reference_path)
        .with_context(|| format!("failed to read reference from {}", reference_path.display()))?;
    let reference: Arc<[u8]> = Arc::from(fasta.sequence.clone().into_boxed_slice());

    let workdir = workdir.unwrap_or_else(|| {
        output_vcf
            .parent()
            .map(|p| p.to_path_buf())
            .unwrap_or_else(|| PathBuf::from("."))
    });
    std::fs::create_dir_all(&workdir)
        .with_context(|| format!("failed to create workdir {}", workdir.display()))?;

    let tumor_bam = workdir.join("tumor.bam");
    let tumor_sorted = workdir.join("tumor.sorted.bam");
    let normal_bam = workdir.join("normal.bam");
    let normal_sorted = workdir.join("normal.sorted.bam");
    let manifest_path = sidecar_path(&output_vcf, ".manifest.json");
    let perf_path = workdir.join("somatic.perf.txt");
    for (path, kind) in [
        (&output_vcf, "somatic VCF"),
        (&manifest_path, "somatic receipt"),
        (&tumor_bam, "tumor BAM"),
        (&tumor_sorted, "sorted tumor BAM"),
        (&normal_bam, "normal BAM"),
        (&normal_sorted, "sorted normal BAM"),
        (&perf_path, "performance report"),
    ] {
        require_safe_cli_destination(path, force, kind);
    }

    // Align tumor reads.
    let start_align_tumor = Instant::now();
    {
        let mut artifact = AtomicFile::create(&tumor_bam)?;
        artifact.close_for_path_writer();
        let mut writer = create_bam_writer(
            artifact.temporary_path(),
            fasta.name.as_str(),
            fasta.sequence.len(),
        )
        .with_context(|| format!("failed to create tumor BAM {}", tumor_bam.display()))?;
        let mut aligner = BWTAligner::new(&reference)?;

        match resolve_reads(
            tumor_fastq.clone(),
            tumor_r1.clone(),
            tumor_r2.clone(),
            "tumor",
        )? {
            ResolvedReads::Single(reads) => {
                let alignments = align_reads(&mut aligner, &reads, 10)?;
                write_bam_alignments(&mut writer, 0, &reads, &alignments)?;
            }
            ResolvedReads::Paired(pairs) => {
                let alignments = align_pairs(&mut aligner, &pairs, 10)?;
                write_bam_alignments_paired(&mut writer, 0, &pairs, &alignments)?;
            }
        }
        drop(writer);
        artifact.commit(force)?;
    }
    let dur_align_tumor = start_align_tumor.elapsed();

    // Align normal reads.
    let start_align_normal = Instant::now();
    {
        let mut artifact = AtomicFile::create(&normal_bam)?;
        artifact.close_for_path_writer();
        let mut writer = create_bam_writer(
            artifact.temporary_path(),
            fasta.name.as_str(),
            fasta.sequence.len(),
        )
        .with_context(|| format!("failed to create normal BAM {}", normal_bam.display()))?;
        let mut aligner = BWTAligner::new(&reference)?;

        match resolve_reads(
            normal_fastq.clone(),
            normal_r1.clone(),
            normal_r2.clone(),
            "normal",
        )? {
            ResolvedReads::Single(reads) => {
                let alignments = align_reads(&mut aligner, &reads, 10)?;
                write_bam_alignments(&mut writer, 0, &reads, &alignments)?;
            }
            ResolvedReads::Paired(pairs) => {
                let alignments = align_pairs(&mut aligner, &pairs, 10)?;
                write_bam_alignments_paired(&mut writer, 0, &pairs, &alignments)?;
            }
        }
        drop(writer);
        artifact.commit(force)?;
    }
    let dur_align_normal = start_align_normal.elapsed();

    // Sort BAMs deterministically.
    let start_sort = Instant::now();
    let bytes = memory_mb.saturating_mul(1024 * 1024).max(1024 * 1024);
    let mut tumor_sorted_artifact = AtomicFile::create(&tumor_sorted)?;
    tumor_sorted_artifact.close_for_path_writer();
    sort_bam_deterministic(&tumor_bam, tumor_sorted_artifact.temporary_path(), bytes)?;
    tumor_sorted_artifact.commit(force)?;
    let mut normal_sorted_artifact = AtomicFile::create(&normal_sorted)?;
    normal_sorted_artifact.close_for_path_writer();
    sort_bam_deterministic(&normal_bam, normal_sorted_artifact.temporary_path(), bytes)?;
    normal_sorted_artifact.commit(force)?;
    let dur_sort = start_sort.elapsed();

    // Call somatic SNVs on the new engine (indels deferred to Phase C).
    use rosalind::call::{call_somatic_region, SomaticParams};
    use rosalind::core::ContigSet;
    use rosalind::io::bam::BamSource;
    use rosalind::io::vcf::write_somatic_vcf;
    use rosalind::pileup::PileupParams;
    use rosalind::provenance::CommandCapture;

    let start_call = Instant::now();
    let mut contigs = ContigSet::new();
    let contig_id = contigs.push(fasta.name.clone(), reference.len() as u32);
    let region = 0..(reference.len() as u32);

    let tumor_src = BamSource::new(&tumor_sorted, &contigs)
        .map_err(|e| anyhow!("failed to read tumor BAM {}: {e}", tumor_sorted.display()))?;
    let normal_src = BamSource::new(&normal_sorted, &contigs)
        .map_err(|e| anyhow!("failed to read normal BAM {}: {e}", normal_sorted.display()))?;
    let calls = call_somatic_region(
        tumor_src,
        normal_src,
        Arc::clone(&reference),
        contig_id,
        region,
        PileupParams::default(),
        &SomaticParams::default(),
    )
    .map_err(|e| anyhow!("somatic calling failed: {e}"))?;
    let dur_call = start_call.elapsed();

    // Write spec-valid somatic VCF (TUMOR/NORMAL).
    {
        let mut artifact = AtomicFile::create(&output_vcf)?;
        let mut writer = io::BufWriter::new(artifact.file_mut());
        write_somatic_vcf(&mut writer, &contigs, &calls)?;
        writer.flush()?;
        drop(writer);
        artifact.commit(force)?;
    }

    // Reproducibility receipt (BLAKE3, canonical JSON).
    let tumor_inputs: Vec<PathBuf> = match (&tumor_fastq, &tumor_r1, &tumor_r2) {
        (Some(p), None, None) => vec![p.clone()],
        (None, Some(r1), Some(r2)) => vec![r1.clone(), r2.clone()],
        _ => Vec::new(),
    };
    let normal_inputs: Vec<PathBuf> = match (&normal_fastq, &normal_r1, &normal_r2) {
        (Some(p), None, None) => vec![p.clone()],
        (None, Some(r1), Some(r2)) => vec![r1.clone(), r2.clone()],
        _ => Vec::new(),
    };
    let mut manifest = new_run_manifest("somatic");
    let mut cmd = CommandCapture::new("somatic");
    cmd.input("--reference", &reference_path)?;
    for p in tumor_inputs.iter() {
        if p.exists() {
            cmd.input("--tumor", p)?;
        }
    }
    for p in normal_inputs.iter() {
        if p.exists() {
            cmd.input("--normal", p)?;
        }
    }
    cmd.flag_if(force, "--force");
    cmd.output("-o", &output_vcf)?;
    cmd.record_into(&mut manifest);
    manifest
        .params
        .insert("artifact.input.0.role".to_string(), "reference".to_string());
    for index in 1..manifest.inputs.len() {
        manifest.params.insert(
            format!("artifact.input.{index}.role"),
            if index <= tumor_inputs.len() {
                "tumor-reads"
            } else {
                "normal-reads"
            }
            .to_string(),
        );
    }
    manifest.params.insert(
        "artifact.output.0.role".to_string(),
        "somatic-calls".to_string(),
    );
    manifest
        .params
        .insert("somatic_snv_only".to_string(), "true".to_string());
    manifest.finalize();
    write_atomic(
        &manifest_path,
        manifest.to_canonical_json().as_bytes(),
        force,
    )?;
    eprintln!("wrote reproducibility receipt: {}", manifest_path.display());

    // Performance/RSS report (kept separate from determinism checks).
    let mut perf_artifact = AtomicFile::create(&perf_path)?;
    let mut f = io::BufWriter::new(perf_artifact.file_mut());
    writeln!(f, "align_tumor_ms={}", dur_align_tumor.as_millis())?;
    writeln!(f, "align_normal_ms={}", dur_align_normal.as_millis())?;
    writeln!(f, "sort_ms={}", dur_sort.as_millis())?;
    writeln!(f, "call_ms={}", dur_call.as_millis())?;
    writeln!(f, "total_ms={}", start_total.elapsed().as_millis())?;
    writeln!(f, "peak_rss_bytes={}", peak_rss_bytes())?;
    f.flush()?;
    drop(f);
    perf_artifact.commit(force)?;

    Ok(())
}

#[allow(clippy::too_many_arguments)] // a CLI entry point: each flag is a parameter
fn run_align(
    reference_path: PathBuf,
    reads_path: Option<PathBuf>,
    reads_r1: Option<PathBuf>,
    reads_r2: Option<PathBuf>,
    max_mismatches: usize,
    reference_offset: u32,
    format: OutputFormat,
    output: Option<PathBuf>,
    force: bool,
    manifest_out: Option<PathBuf>,
) -> Result<()> {
    use rosalind::util::atomic::{write_atomic, AtomicFile};

    if format == OutputFormat::Bam && output.is_none() {
        bail!("--output <FILE> must be provided when writing BAM output");
    }
    let receipt_path = manifest_out.clone().or_else(|| {
        output
            .as_ref()
            .map(|path| sidecar_path(path, ".manifest.json"))
    });
    if let Some(path) = &output {
        require_safe_cli_destination(path, force, "alignment output");
    }
    if let Some(path) = &receipt_path {
        require_safe_cli_destination(path, force, "alignment receipt");
    }
    let mut atomic_output = output
        .as_ref()
        .map(|path| AtomicFile::create(path))
        .transpose()
        .with_context(|| "failed to reserve atomic alignment output")?;
    if format == OutputFormat::Bam {
        if let Some(file) = &mut atomic_output {
            file.close_for_path_writer();
        }
    }
    let receipt_reads = [
        ("--reads", reads_path.clone()),
        ("--reads-r1", reads_r1.clone()),
        ("--reads-r2", reads_r2.clone()),
    ];
    let fasta = read_fasta(&reference_path)
        .with_context(|| format!("failed to read reference from {}", reference_path.display()))?;
    let reads = resolve_reads(reads_path, reads_r1, reads_r2, "reads")
        .context("failed to resolve reads inputs")?;
    let mut aligner = BWTAligner::new(&fasta.sequence)
        .context("failed to initialize the FM-index aligner for the reference")?;

    match reads {
        ResolvedReads::Single(records) => {
            let alignments = align_reads(&mut aligner, &records, max_mismatches)
                .context("aligning reads failed")?;
            match format {
                OutputFormat::Sam => {
                    if let Some(file) = &mut atomic_output {
                        let mut writer = io::BufWriter::new(file.file_mut());
                        write_sam_alignments(
                            &mut writer,
                            &fasta.name,
                            fasta.sequence.len(),
                            reference_offset,
                            &records,
                            &alignments,
                        )?;
                        writer.flush()?;
                    } else {
                        let stdout = io::stdout();
                        let mut handle = stdout.lock();
                        write_sam_alignments(
                            &mut handle,
                            &fasta.name,
                            fasta.sequence.len(),
                            reference_offset,
                            &records,
                            &alignments,
                        )?;
                    }
                }
                OutputFormat::Bam => {
                    let temporary = atomic_output
                        .as_ref()
                        .expect("BAM output was validated")
                        .temporary_path();
                    let mut writer =
                        create_bam_writer(temporary, &fasta.name, fasta.sequence.len())
                            .with_context(|| {
                                format!("failed to create BAM writer for {}", temporary.display())
                            })?;
                    write_bam_alignments(&mut writer, reference_offset, &records, &alignments)?;
                }
            }
        }
        ResolvedReads::Paired(pairs) => {
            let alignments = align_pairs(&mut aligner, &pairs, max_mismatches)
                .context("aligning reads failed")?;
            match format {
                OutputFormat::Sam => {
                    if let Some(file) = &mut atomic_output {
                        let mut writer = io::BufWriter::new(file.file_mut());
                        write_sam_alignments_paired(
                            &mut writer,
                            &fasta.name,
                            fasta.sequence.len(),
                            reference_offset,
                            &pairs,
                            &alignments,
                        )?;
                        writer.flush()?;
                    } else {
                        let stdout = io::stdout();
                        let mut handle = stdout.lock();
                        write_sam_alignments_paired(
                            &mut handle,
                            &fasta.name,
                            fasta.sequence.len(),
                            reference_offset,
                            &pairs,
                            &alignments,
                        )?;
                    }
                }
                OutputFormat::Bam => {
                    let temporary = atomic_output
                        .as_ref()
                        .expect("BAM output was validated")
                        .temporary_path();
                    let mut writer =
                        create_bam_writer(temporary, &fasta.name, fasta.sequence.len())
                            .with_context(|| {
                                format!("failed to create BAM writer for {}", temporary.display())
                            })?;
                    write_bam_alignments_paired(
                        &mut writer,
                        reference_offset,
                        &pairs,
                        &alignments,
                    )?;
                }
            }
        }
    }

    if let Some(file) = atomic_output {
        file.commit(force)
            .with_context(|| "failed to commit atomic alignment output")?;
    }
    if let Some(receipt_path) = receipt_path {
        use rosalind::provenance::CommandCapture;

        let mut receipt = new_run_manifest("align");
        let mut command = CommandCapture::new("align");
        command.input("--reference", &reference_path)?;
        for (flag, path) in &receipt_reads {
            if let Some(path) = path {
                command.input(flag, path)?;
            }
        }
        command.opt("--max-mismatches", max_mismatches);
        command.opt("--reference-offset", reference_offset);
        command.opt(
            "--format",
            match format {
                OutputFormat::Sam => "sam",
                OutputFormat::Bam => "bam",
            },
        );
        command.flag_if(force, "--force");
        if let Some(output) = &output {
            command.output("--output", output)?;
        }
        command.record_into(&mut receipt);
        receipt
            .params
            .insert("artifact.input.0.role".to_string(), "reference".to_string());
        for index in 1..receipt.inputs.len() {
            receipt.params.insert(
                format!("artifact.input.{index}.role"),
                if receipt.inputs.len() == 2 {
                    "reads".to_string()
                } else if index == 1 {
                    "reads-r1".to_string()
                } else {
                    "reads-r2".to_string()
                },
            );
        }
        if !receipt.outputs.is_empty() {
            receipt.params.insert(
                "artifact.output.0.role".to_string(),
                "raw-alignments".to_string(),
            );
        }
        receipt.record_measurement("peak_rss_bytes", peak_rss_bytes().to_string());
        receipt.finalize();
        write_atomic(&receipt_path, receipt.to_canonical_json().as_bytes(), force)
            .with_context(|| format!("failed to write align receipt {}", receipt_path.display()))?;
        eprintln!("wrote reproducibility receipt: {}", receipt_path.display());
    }
    Ok(())
}

fn align_reads(
    aligner: &mut BWTAligner,
    records: &[FastqRecord],
    max_mismatches: usize,
) -> Result<Vec<Option<AlignmentCandidate>>> {
    let mut results = Vec::with_capacity(records.len());
    for record in records {
        let alignment = aligner
            .align_read(&record.sequence)
            .with_context(|| format!("failed to align read {}", record.name))?;

        if alignment.has_candidates() && alignment.mismatches <= max_mismatches {
            let position = alignment
                .primary_position
                .map(|p| p as usize)
                .or_else(|| {
                    if alignment.interval.is_empty() {
                        None
                    } else {
                        Some(aligner.fm_index().sa_at(alignment.interval.lower as usize))
                    }
                })
                .ok_or_else(|| {
                    anyhow!("alignment reported candidates but no position could be located")
                })?;
            results.push(Some(AlignmentCandidate {
                position,
                mismatches: alignment.mismatches,
                mapq: alignment.mapq,
                is_reverse: alignment.is_reverse,
                cigar: alignment.cigar.clone(),
                as_score: alignment.as_score,
                md: alignment.md.clone(),
            }));
        } else {
            results.push(None);
        }
    }

    Ok(results)
}

type PairAlignments = (Option<AlignmentCandidate>, Option<AlignmentCandidate>);

fn align_pairs(
    aligner: &mut BWTAligner,
    pairs: &[FastqPair],
    max_mismatches: usize,
) -> Result<Vec<PairAlignments>> {
    let mut out = Vec::with_capacity(pairs.len());
    for pair in pairs {
        let a1 = aligner
            .align_read(&pair.r1.sequence)
            .with_context(|| format!("failed to align read {}", pair.name))?;
        let a2 = aligner
            .align_read(&pair.r2.sequence)
            .with_context(|| format!("failed to align read {}", pair.name))?;

        let c1 = if a1.has_candidates() && a1.mismatches <= max_mismatches {
            Some(AlignmentCandidate {
                position: a1
                    .primary_position
                    .map(|p| p as usize)
                    .or_else(|| {
                        if a1.interval.is_empty() {
                            None
                        } else {
                            Some(aligner.fm_index().sa_at(a1.interval.lower as usize))
                        }
                    })
                    .ok_or_else(|| {
                        anyhow!("alignment reported candidates but no position could be located")
                    })?,
                mismatches: a1.mismatches,
                mapq: a1.mapq,
                is_reverse: a1.is_reverse,
                cigar: a1.cigar.clone(),
                as_score: a1.as_score,
                md: a1.md.clone(),
            })
        } else {
            None
        };

        let c2 = if a2.has_candidates() && a2.mismatches <= max_mismatches {
            Some(AlignmentCandidate {
                position: a2
                    .primary_position
                    .map(|p| p as usize)
                    .or_else(|| {
                        if a2.interval.is_empty() {
                            None
                        } else {
                            Some(aligner.fm_index().sa_at(a2.interval.lower as usize))
                        }
                    })
                    .ok_or_else(|| {
                        anyhow!("alignment reported candidates but no position could be located")
                    })?,
                mismatches: a2.mismatches,
                mapq: a2.mapq,
                is_reverse: a2.is_reverse,
                cigar: a2.cigar.clone(),
                as_score: a2.as_score,
                md: a2.md.clone(),
            })
        } else {
            None
        };

        out.push((c1, c2));
    }
    Ok(out)
}

#[allow(clippy::too_many_arguments)] // a CLI entry point: each flag is a parameter
fn run_variants(
    reference_path: PathBuf,
    alignments_path: PathBuf,
    chrom: Option<String>,
    region_start: u32,
    mapq_threshold: u8,
    output: Option<PathBuf>,
    _block_size: usize,
    quality_threshold: f32,
    force: bool,
    manifest_out: Option<PathBuf>,
) -> Result<()> {
    use rosalind::call::{call_germline_region, GermlineParams};
    use rosalind::core::ContigSet;
    use rosalind::io::bam::BamSource;
    use rosalind::io::vcf::{write_germline_vcf, GermlineRow};
    use rosalind::pileup::{PileupParams, SliceSource};
    use rosalind::provenance::CommandCapture;
    use rosalind::util::atomic::{write_atomic, AtomicFile};

    let receipt_dest = manifest_out.clone().or_else(|| {
        output
            .as_ref()
            .map(|path| sidecar_path(path, ".manifest.json"))
    });
    if let Some(path) = &output {
        require_safe_cli_destination(path, force, "VCF output");
    }
    if let Some(path) = &receipt_dest {
        require_safe_cli_destination(path, force, "VCF receipt");
    }

    let fasta = read_fasta(&reference_path)
        .with_context(|| format!("failed to read reference from {}", reference_path.display()))?;
    let chrom_name = chrom.unwrap_or_else(|| fasta.name.clone());
    let chrom_arc: Arc<str> = chrom_name.clone().into();
    let reference: Arc<[u8]> = Arc::from(fasta.sequence.into_boxed_slice());

    // Single-contig run: one contig spanning the reference window.
    let mut contigs = ContigSet::new();
    let contig_id = contigs.push(chrom_name.clone(), region_start + reference.len() as u32);

    let region = region_start..(region_start + reference.len() as u32);
    let pileup_params = PileupParams {
        min_mapq: mapq_threshold, // now actually honored by the engine
        ..PileupParams::default()
    };
    let germline_params = GermlineParams {
        min_qual: quality_threshold as f64,
        ..GermlineParams::default()
    };

    let is_bam = alignments_path
        .extension()
        .and_then(|ext| ext.to_str())
        .map(|ext| ext.eq_ignore_ascii_case("bam"))
        .unwrap_or(false);

    let sites = if is_bam {
        let source = BamSource::new(&alignments_path, &contigs)
            .map_err(|e| anyhow!("failed to read BAM {}: {e}", alignments_path.display()))?;
        call_germline_region(
            source,
            Arc::clone(&reference),
            contig_id,
            region,
            pileup_params,
            &germline_params,
        )
        .map_err(|e| anyhow!("variant calling failed (BAM): {e}"))?
    } else {
        let legacy =
            read_alignment_file(&alignments_path, Some(&chrom_arc)).with_context(|| {
                format!(
                    "failed to read SAM alignments from {}",
                    alignments_path.display()
                )
            })?;
        let core_reads: Vec<rosalind::core::AlignedRead> = legacy
            .iter()
            .map(|r| legacy_read_to_core(r, contig_id))
            .collect();
        let source = SliceSource::new(core_reads);
        call_germline_region(
            source,
            Arc::clone(&reference),
            contig_id,
            region,
            pileup_params,
            &germline_params,
        )
        .map_err(|e| anyhow!("variant calling failed (SAM): {e}"))?
    };

    let rows: Vec<GermlineRow> = sites
        .into_iter()
        .map(|(locus, ref_base, call)| GermlineRow {
            locus,
            ref_base,
            call,
        })
        .collect();

    match &output {
        Some(path) => {
            let mut file = AtomicFile::create(path)
                .with_context(|| format!("failed to reserve VCF file {}", path.display()))?;
            let mut writer = io::BufWriter::new(file.file_mut());
            write_germline_vcf(&mut writer, &contigs, &chrom_name, &rows)?;
            writer.flush()?;
            drop(writer);
            file.commit(force)
                .with_context(|| format!("failed to commit VCF file {}", path.display()))?;
        }
        None => {
            let stdout = io::stdout();
            let mut handle = stdout.lock();
            write_germline_vcf(&mut handle, &contigs, &chrom_name, &rows)?;
        }
    }

    if let Some(receipt_dest) = receipt_dest {
        let mut manifest = new_run_manifest("variants");
        manifest
            .params
            .insert("pileup.semantics".into(), "exact-or-fail-v1".into());
        let mut cmd = CommandCapture::new("variants");
        cmd.input("--reference", &reference_path)?;
        cmd.input("--alignments", &alignments_path)?;
        cmd.opt("--chrom", &chrom_name);
        cmd.opt("--region-start", region_start);
        cmd.opt("--mapq-threshold", mapq_threshold);
        cmd.opt("--quality-threshold", quality_threshold as f64);
        cmd.flag_if(force, "--force");
        if let Some(path) = &output {
            cmd.output("-o", path)?;
        }
        cmd.record_into(&mut manifest);
        manifest
            .params
            .insert("artifact.input.0.role".to_string(), "reference".to_string());
        manifest.params.insert(
            "artifact.input.1.role".to_string(),
            "alignments".to_string(),
        );
        if !manifest.outputs.is_empty() {
            manifest
                .params
                .insert("artifact.output.0.role".to_string(), "calls".to_string());
        }
        manifest
            .params
            .insert("model.germline".to_string(), "baseq-mapq-v1".to_string());
        manifest.finalize();
        write_atomic(
            &receipt_dest,
            manifest.to_canonical_json().as_bytes(),
            force,
        )?;
        eprintln!("wrote reproducibility receipt: {}", receipt_dest.display());
    }

    Ok(())
}

/// Call germline variants across all contigs of a persisted index (B4), reading
/// the reference from the index and streaming the (coordinate-sorted) BAM.
/// Stream a bounded, deterministic per-locus FEATURE table (TSV) over a persisted
/// index. Same engine + memory contract as `run_variants_index`, but every callable
/// locus is emitted as ML-ready features. Byte-identical run-to-run.
#[allow(clippy::too_many_arguments)]
fn run_bounded_analysis(
    subcommand: &str,
    param_prefix: &str,
    analyzer: &mut dyn rosalind::call::ColumnAnalyzer,
    analysis_reference: PathBuf,
    alignments_path: PathBuf,
    mapq_threshold: u8,
    memory_budget_mb: Option<u64>,
    max_depth: u32,
    max_read_len: u32,
    enforce: bool,
    require_os_limit: bool,
    force: bool,
    output: Option<PathBuf>,
    manifest_out: Option<PathBuf>,
    selection_args: SelectionArgs,
) -> Result<()> {
    use rosalind::contract::{
        AnalyzerIdentity, AnalyzerMemoryModel, ContractRunError, ContractRunSpec, ContractVerdict,
        EnforcementMode, OutputPolicy, OutputTarget, ProducerIdentity, ReplayInvocation,
    };

    let output_target = output
        .clone()
        .map(OutputTarget::File)
        .unwrap_or(OutputTarget::Stdout);
    let spec = ContractRunSpec {
        producer: ProducerIdentity::rosalind(),
        analyzer: AnalyzerIdentity::new(
            subcommand.strip_prefix("analyze ").unwrap_or(subcommand),
            env!("CARGO_PKG_VERSION"),
        )
        .with_param_prefix(param_prefix),
        analyzer_memory: AnalyzerMemoryModel::Fixed {
            model_id: "fixed-additional-v1".to_string(),
            max_additional_bytes: 0,
        },
        invocation: ReplayInvocation::new(subcommand.split_whitespace()),
        index: analysis_reference,
        alignments: alignments_path,
        output: output_target,
        output_policy: if force {
            OutputPolicy::ReplaceAtomic
        } else {
            OutputPolicy::CreateNewAtomic
        },
        manifest: manifest_out,
        mapq_threshold,
        max_depth,
        max_read_len,
        memory_budget_mb,
        enforcement: if require_os_limit {
            EnforcementMode::RequireOsLimit
        } else if enforce {
            EnforcementMode::Cooperative
        } else {
            EnforcementMode::RecordOnly
        },
    };

    let render_outcome = |outcome: &rosalind::contract::ContractRunOutcome| {
        if let Some(path) = &outcome.manifest_path {
            eprintln!("wrote reproducibility receipt: {}", path.display());
        } else {
            eprintln!(
                "no receipt written (stdout output) — pass --manifest <path> or -o <tsv> to persist one"
            );
        }
        eprintln!(
            "{subcommand}: peak RSS {} MiB; max pileup working set {} KiB",
            outcome.peak_rss_bytes / (1 << 20),
            outcome.max_working_set_bytes / 1024
        );
    };

    match rosalind::contract::run_column_analysis_resolving(analyzer, spec, |reference| {
        selection_from_provider(reference, &selection_args)
            .map_err(|error| rosalind::core::CoreError::MalformedRecord(error.to_string()))
    }) {
        Ok(outcome) => {
            render_outcome(&outcome);
            if enforce && outcome.verdict == ContractVerdict::Within {
                eprintln!(
                    "contract: OK — realized peak {} MiB within declared {} MiB",
                    outcome.peak_rss_bytes / (1 << 20),
                    memory_budget_mb.expect("enforced run has a budget")
                );
            }
            Ok(())
        }
        Err(ContractRunError::Refused(report)) => {
            eprintln!(
                "contract: REFUSE — declared {} MiB, predicted peak ~{} MiB (largest contig {} MiB \
                 + active @ max-depth {} / max-read-len {} atop a {} MiB baseline). Raise \
                 --memory-budget-mb, lower --max-depth, or drop --enforce.",
                report.budget_mb,
                report.predicted_peak_rss_bytes / (1 << 20),
                report.largest_contig_bytes / (1 << 20),
                report.max_depth,
                report.max_read_len,
                report.baseline_rss_bytes / (1 << 20),
            );
            std::process::exit(3);
        }
        Err(ContractRunError::Breached(outcome)) => {
            render_outcome(&outcome);
            if let Some(failure) = &outcome.capacity_exceeded {
                eprintln!("pileup capacity exceeded at contig {}, zero-based position {}: {} active reads exceed --max-depth {}; exact partial output and receipt preserved; raise --max-depth and rerun", failure.contig, failure.position, failure.required, failure.capacity);
            } else {
                eprintln!(
                    "contract: VIOLATED — realized peak {} MiB (partial output + receipt written)",
                    outcome.peak_rss_bytes / (1 << 20)
                );
            }
            std::process::exit(4);
        }
        Err(ContractRunError::OutputExists(path)) => {
            eprintln!(
                "output collision: {} (choose a new path or pass --force for atomic replacement)",
                path.display()
            );
            std::process::exit(2);
        }
        Err(ContractRunError::UnknownAnalyzerBound) => {
            eprintln!("contract: REFUSE — the analyzer has no declared memory bound");
            std::process::exit(3);
        }
        Err(ContractRunError::OsEnforcementUnavailable(message)) => {
            eprintln!("contract: REFUSE — OS enforcement unavailable: {message}");
            std::process::exit(3);
        }
        Err(error) => Err(anyhow!(error)),
    }
}

/// The shipped `features` egress: the FeatureAnalyzer through the one bounded-analysis
/// path. `param_prefix=""` reproduces the historic `feature_rows` claim key exactly, so
/// the receipt stays byte-identical.
#[allow(clippy::too_many_arguments)]
fn run_features(
    analysis_reference: PathBuf,
    alignments_path: PathBuf,
    mapq_threshold: u8,
    memory_budget_mb: Option<u64>,
    max_depth: u32,
    max_read_len: u32,
    enforce: bool,
    require_os_limit: bool,
    force: bool,
    output: Option<PathBuf>,
    manifest_out: Option<PathBuf>,
    selection_args: SelectionArgs,
    format: FeatureFormat,
) -> Result<()> {
    match format {
        FeatureFormat::Tsv => {
            let mut analyzer = rosalind::call::FeatureAnalyzer::default();
            run_bounded_analysis(
                "features",
                "",
                &mut analyzer,
                analysis_reference,
                alignments_path,
                mapq_threshold,
                memory_budget_mb,
                max_depth,
                max_read_len,
                enforce,
                require_os_limit,
                force,
                output,
                manifest_out,
                selection_args,
            )
        }
        FeatureFormat::ArrowIpc => {
            let mut analyzer = rosalind::call::FeatureArrowAnalyzer::new()?;
            run_bounded_analysis(
                "features",
                "",
                &mut analyzer,
                analysis_reference,
                alignments_path,
                mapq_threshold,
                memory_budget_mb,
                max_depth,
                max_read_len,
                enforce,
                require_os_limit,
                force,
                output,
                manifest_out,
                selection_args,
            )
        }
    }
}

fn select_analysis_reference(index: Option<PathBuf>, reference_pack: Option<PathBuf>) -> PathBuf {
    if let Some(reference_pack) = reference_pack {
        reference_pack
    } else {
        let index = index.expect("clap requires --index or --reference-pack");
        eprintln!(
            "migration: --index remains supported; build a smaller analysis reference with `rosalind reference convert --index {} --output reference.rref`",
            index.display()
        );
        index
    }
}

fn selection_requested(args: &SelectionArgs) -> bool {
    args.region.is_some()
        || args.regions.is_some()
        || args.shard_count.is_some()
        || args.shard_index.is_some()
}

fn resolve_selection(path: &Path, args: &SelectionArgs) -> Result<rosalind::AnalysisSelection> {
    use rosalind::genomics::AnalysisReference;
    let reference = AnalysisReference::open(path)
        .with_context(|| format!("failed to open analysis reference {}", path.display()))?;
    selection_from_provider(&reference, args)
}

fn selection_from_provider(
    reference: &dyn rosalind::ReferenceProvider,
    args: &SelectionArgs,
) -> Result<rosalind::AnalysisSelection> {
    if let Some(region) = &args.region {
        return Ok(rosalind::AnalysisSelection::Intervals(
            rosalind::IntervalSet::parse_region(region, reference.contigs())?,
        ));
    }
    if let Some(path) = &args.regions {
        return Ok(rosalind::AnalysisSelection::Intervals(
            rosalind::IntervalSet::from_bed(path, reference.contigs())?,
        ));
    }
    if let (Some(count), Some(index)) = (args.shard_count, args.shard_index) {
        return Ok(rosalind::AnalysisSelection::shard(
            count,
            index,
            reference.contigs(),
        )?);
    }
    Ok(rosalind::AnalysisSelection::WholeGenome)
}

fn record_selection_claims(
    params: &mut std::collections::BTreeMap<String, String>,
    selection: &rosalind::AnalysisSelection,
) {
    params.insert("partition.kind".into(), selection.kind().into());
    if let Some(intervals) = selection.intervals() {
        params.insert(
            "partition.interval_count".into(),
            intervals.intervals().len().to_string(),
        );
        params.insert(
            "partition.total_bases".into(),
            intervals.total_bases().to_string(),
        );
        params.insert("partition.intervals_blake3".into(), intervals.blake3());
    }
    if let rosalind::AnalysisSelection::Shard { count, index, .. } = selection {
        params.insert("partition.algorithm".into(), "reference-span-v1".into());
        params.insert("partition.shard_count".into(), count.to_string());
        params.insert("partition.shard_index".into(), index.to_string());
    }
}

#[allow(clippy::too_many_arguments)] // a CLI entry point: each flag is a parameter
fn run_variants_index(
    analysis_reference: PathBuf,
    alignments_path: PathBuf,
    mapq_threshold: u8,
    output: Option<PathBuf>,
    quality_threshold: f32,
    memory_budget_mb: Option<u64>,
    max_depth: u32,
    max_read_len: u32,
    enforce: bool,
    require_os_limit: bool,
    force: bool,
    manifest_out: Option<PathBuf>,
    gvcf: bool,
    selection_args: SelectionArgs,
) -> Result<()> {
    use rosalind::util::atomic::AtomicFile;

    let receipt_dest = manifest_out.clone().or_else(|| {
        output
            .as_ref()
            .map(|path| sidecar_path(path, ".manifest.json"))
    });
    if !force {
        if let Some(path) = &output {
            require_safe_cli_destination(path, false, "VCF output");
            require_safe_cli_destination(
                &sidecar_path(path, ".partial"),
                false,
                "partial VCF output",
            );
        }
        if let Some(path) = &receipt_dest {
            require_safe_cli_destination(path, false, "VCF receipt");
        }
    }
    if enforce && memory_budget_mb.is_none() {
        bail!("--enforce requires --memory-budget-mb");
    }
    if enforce && max_depth == 0 {
        bail!("--enforce requires --max-depth > 0");
    }
    if let (Some(artifact), Some(receipt)) = (&output, &receipt_dest) {
        if artifact == receipt || sidecar_path(artifact, ".partial") == *receipt {
            bail!("artifact, partial artifact, and receipt destinations must differ");
        }
    }
    let os_limit_bytes = if require_os_limit {
        let budget = memory_budget_mb.expect("clap requires --enforce; validation requires budget");
        let required = budget.saturating_mul(1 << 20);
        match rosalind::contract::detected_os_memory_limit_bytes() {
            Some(limit) if limit <= required => Some(limit),
            Some(limit) => {
                eprintln!(
                    "contract: REFUSE — OS enforcement unavailable: cgroup v2 memory.max {} MiB exceeds the declared {budget} MiB budget. Lower the container/systemd memory limit to at most {budget} MiB.",
                    limit / (1 << 20)
                );
                std::process::exit(3);
            }
            None => {
                eprintln!(
                    "contract: REFUSE — OS enforcement unavailable: run inside a cgroup-v2 container or systemd scope with memory.max at most {budget} MiB. Rosalind does not create privileged cgroups."
                );
                std::process::exit(3);
            }
        }
    } else {
        None
    };
    use rosalind::call::{
        call_germline_selected_bam, call_germline_whole_genome, stream_gvcf_selected_bam,
        stream_gvcf_whole_genome, write_gvcf_header, GermlineParams,
    };
    use rosalind::core::governor::MemoryGovernor;
    use rosalind::core::PILEUP_IO_RSS_OVERHEAD;
    use rosalind::genomics::{AnalysisReference, ReferenceProvider};
    use rosalind::io::bam::{find_bai, StreamingBamSource};
    use rosalind::io::vcf::{write_germline_header, write_germline_row, GermlineRow};
    use rosalind::pileup::PileupParams;
    use rosalind::provenance::CommandCapture;

    let startup_baseline = peak_rss_bytes();
    // Live RSS source for the governor: a test seam (ROSALIND_FORCE_LIVE_RSS_BYTES)
    // standing in for live RSS, else the real getrusage high-water mark. Distinct
    // from ROSALIND_FORCE_PEAK_RSS_BYTES, which overrides only the POST-run peak.
    let live_rss = || {
        std::env::var("ROSALIND_FORCE_LIVE_RSS_BYTES")
            .ok()
            .and_then(|v| v.parse::<u64>().ok())
            .unwrap_or_else(peak_rss_bytes)
    };
    let poll_ms = std::env::var("ROSALIND_GOVERNOR_POLL_MS")
        .ok()
        .and_then(|v| v.parse::<u64>().ok())
        .unwrap_or(100);
    // Under --enforce, the governor fails the run LOUD the moment live RSS crosses
    // the budget (exit 4 with partial output + receipt). Held
    // for the duration of the calling pass; dropped (thread stopped) at scope end.
    let _governor_guard = if enforce {
        let mb = memory_budget_mb.expect("--enforce requires --memory-budget-mb (checked above)");
        Some(
            MemoryGovernor::start(
                MemoryBudget::from_mb(mb).bytes,
                std::time::Duration::from_millis(poll_ms),
                live_rss,
            )
            .map_err(|e| anyhow!("failed to start memory governor: {e}"))?,
        )
    } else {
        None
    };

    let loaded = AnalysisReference::open(&analysis_reference).with_context(|| {
        format!(
            "failed to open analysis reference {}",
            analysis_reference.display()
        )
    })?;
    let contigs = loaded.contigs();
    let selection = selection_from_provider(&loaded, &selection_args)?;
    let baseline = peak_rss_bytes().max(startup_baseline);

    let pileup_params = PileupParams {
        min_mapq: mapq_threshold,
        // `--max-depth 0` opts out of the cap (then the working set is unbounded
        // and `--enforce` is rejected below).
        max_depth: if max_depth == 0 {
            None
        } else {
            Some(max_depth)
        },
        // Under `--enforce` the predicted envelope assumes reads <= max_read_len;
        // check it at ingest so a longer read aborts loudly instead of voiding it.
        max_read_len: if enforce { Some(max_read_len) } else { None },
        ..PileupParams::default()
    };
    let germline_params = GermlineParams {
        min_qual: quality_threshold as f64,
        ..GermlineParams::default()
    };

    // `--index` streams a coordinate-sorted BAM (the bounded, multi-contig WGS
    // path). SAM under `--index` is unsupported — the legacy SAM reader is
    // single-contig; use a sorted BAM (`rosalind sort`), or `--reference` for
    // single-contig SAM.
    let is_bam = alignments_path
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.eq_ignore_ascii_case("bam"))
        .unwrap_or(false);
    if !is_bam {
        bail!(
            "--index/--reference-pack requires a coordinate-sorted BAM (use `rosalind sort`); \
             for single-contig SAM use --reference"
        );
    }
    if !matches!(selection, rosalind::AnalysisSelection::WholeGenome) {
        find_bai(&alignments_path).ok_or_else(|| {
            anyhow!(
                "sparse analysis requires {}.bai or {}",
                alignments_path.display(),
                alignments_path.with_extension("bai").display()
            )
        })?;
    }

    // Predict the peak RSS up front: a measured baseline (binary + libs + index/
    // BAM open) plus the depth-capped working set plus an RSS margin. Computed
    // unconditionally and recorded in the receipt — it is the contract's up-front
    // claim, which the post-run check and `verify` assert the realized peak honors.
    // Under `--enforce` it also gates the run: refuse cleanly before any work,
    // with the cooperative assurance recorded in the receipt.
    let largest = selection.largest_reference_span(contigs);
    let predicted_peak =
        rosalind::call::plan::predicted_peak_rss_bytes(largest, max_depth, max_read_len, baseline);
    if enforce {
        if memory_budget_mb.is_none() {
            bail!("--enforce requires --memory-budget-mb");
        }
        if max_depth == 0 {
            bail!(
                "--enforce requires --max-depth > 0 (an uncapped active set has no a-priori bound)"
            );
        }
        let mb = memory_budget_mb.unwrap();
        if !MemoryBudget::from_mb(mb).admits(predicted_peak) {
            eprintln!(
                "contract: REFUSE — declared {} MiB, predicted peak ~{} MiB \
                 (largest contig {} MiB + active @ max-depth {} / max-read-len {} \
                 atop a {} MiB baseline). Raise --memory-budget-mb, lower --max-depth, \
                 or drop --enforce.",
                mb,
                predicted_peak / (1 << 20),
                largest / (1 << 20),
                max_depth,
                max_read_len,
                baseline / (1 << 20),
            );
            std::process::exit(3);
        }
    }

    // Stream straight to the writer (header once, then one record per callable
    // locus in gVCF mode, or per variant otherwise) so no genome-wide row buffer
    // accumulates. The returned WorkingSet is the high-water (reference + active
    // set), captured per contig — the gVCF banding state is O(1), so the bound
    // holds for both modes. A small macro keeps the file/stdout arms DRY while
    // each passes a concrete (sized) writer.
    macro_rules! drive {
        ($writer:expr) => {{
            let w = $writer;
            let r = if gvcf {
                write_gvcf_header(&mut *w, contigs, "SAMPLE")?;
                match &selection {
                    rosalind::AnalysisSelection::WholeGenome => {
                        let source = StreamingBamSource::new(&alignments_path, contigs)?;
                        stream_gvcf_whole_genome(
                            source,
                            &loaded,
                            contigs,
                            pileup_params,
                            &germline_params,
                            &mut *w,
                        )
                    }
                    _ => stream_gvcf_selected_bam(
                        &alignments_path,
                        &loaded,
                        &selection,
                        pileup_params,
                        &germline_params,
                        &mut *w,
                    ),
                }
            } else {
                write_germline_header(&mut *w, contigs, "SAMPLE")?;
                let mut emit = |(locus, ref_base, call)| {
                    write_germline_row(
                        &mut *w,
                        contigs,
                        &GermlineRow {
                            locus,
                            ref_base,
                            call,
                        },
                    )
                    .map_err(rosalind::core::CoreError::from)
                };
                match &selection {
                    rosalind::AnalysisSelection::WholeGenome => {
                        let source = StreamingBamSource::new(&alignments_path, contigs)?;
                        call_germline_whole_genome(
                            source,
                            &loaded,
                            contigs,
                            pileup_params,
                            &germline_params,
                            &mut emit,
                        )
                    }
                    _ => call_germline_selected_bam(
                        &alignments_path,
                        &loaded,
                        &selection,
                        pileup_params,
                        &germline_params,
                        &mut emit,
                    ),
                }
            };
            // Flush even on a governed abort so partial output survives; surface a
            // genuine flush failure only when the calling pass itself succeeded.
            if r.is_ok() {
                w.flush()?;
            } else {
                let _ = w.flush();
            }
            r
        }};
    }
    let mut atomic_output = output
        .as_ref()
        .map(|path| AtomicFile::create(path))
        .transpose()
        .with_context(|| "failed to reserve transactional VCF output")?;
    let drive_result: Result<
        (rosalind::core::WorkingSet, rosalind::pileup::SkipCounts),
        rosalind::core::CoreError,
    > = match &mut atomic_output {
        Some(file) => {
            let mut writer = io::BufWriter::new(file.file_mut());
            drive!(&mut writer)
        }
        None => {
            let stdout = io::stdout();
            let mut handle = stdout.lock();
            drive!(&mut handle)
        }
    };
    // A governor trip is the one error we do NOT bail on: we still write the proof
    // receipt (verdict=over, governor=tripped) and exit 4 via the existing post-run
    // check. Any other error is a genuine failure.
    let (max_ws, skips, mut breached, breach_peak, capacity_failure) = match drive_result {
        Ok((ws, sk)) => (ws, sk, false, 0u64, None),
        Err(rosalind::core::CoreError::BudgetExceeded { needed, .. }) => (
            rosalind::core::WorkingSet { bytes: 0 },
            rosalind::pileup::SkipCounts::default(),
            true,
            needed,
            None,
        ),
        Err(rosalind::core::CoreError::CapacityExceeded {
            contig,
            position,
            capacity,
            required,
        }) => (
            rosalind::core::WorkingSet { bytes: 0 },
            rosalind::pileup::SkipCounts::default(),
            false,
            0,
            Some(rosalind::contract::CapacityFailure {
                contig,
                position,
                capacity,
                required,
            }),
        ),
        Err(e) => return Err(anyhow!("variant calling failed: {e}")),
    };
    // Hash under the same governor, before final measurements and publication.
    let output_hash = atomic_output
        .as_ref()
        .map(|file| rosalind::provenance::blake3_file(file.temporary_path()))
        .transpose()?;
    let input_hashes = if receipt_dest.is_some() {
        let mut hashes = vec![
            rosalind::provenance::blake3_file(&analysis_reference)?,
            rosalind::provenance::blake3_file(&alignments_path)?,
        ];
        if let Some(path) = selection
            .intervals()
            .and_then(|intervals| intervals.bed_origin())
        {
            hashes.push(rosalind::provenance::blake3_file(path)?);
        }
        hashes
    } else {
        Vec::new()
    };
    // Realized peak: the governor's tripping peak on a live breach, else the
    // post-run high-water (monotonic). Test-only seam: ROSALIND_FORCE_PEAK_RSS_BYTES
    // overrides ONLY this post-run realized peak (never the pre-run baseline at the
    // --enforce gate), so the exit-4 backstop can be exercised without allocating.
    let mut peak_rss = if breached {
        breach_peak
    } else {
        std::env::var("ROSALIND_FORCE_PEAK_RSS_BYTES")
            .ok()
            .and_then(|v| v.parse::<u64>().ok())
            .unwrap_or_else(peak_rss_bytes)
    };

    if let Err(rosalind::core::CoreError::BudgetExceeded { needed, .. }) =
        rosalind::core::governor::checkpoint()
    {
        breached = true;
        peak_rss = peak_rss.max(needed);
    }

    // Compute the contract verdict before writing the receipt (so it records it).
    let mut verdict = match memory_budget_mb.map(|mb| MemoryBudget::from_mb(mb).admits(peak_rss)) {
        None => "unset",
        Some(true) => "within",
        Some(false) => "over",
    };
    let mut is_breach = capacity_failure.is_some() || breached || (enforce && verdict == "over");
    let receipt_output = output.as_ref().map(|path| {
        if is_breach {
            sidecar_path(path, ".partial")
        } else {
            path.clone()
        }
    });
    // Distinguish a live-governed abort from a post-run-detected overrun in the receipt.
    let mut governor_state = if breached {
        "tripped"
    } else if enforce {
        "enforced"
    } else {
        "record-only"
    };
    // Measured RSS residual telemetry: the real I/O + allocator slack this run
    // incurred above the modeled working set, recorded so the fixed 8 MiB
    // prediction margin can later be re-tuned with evidence (Sprint 1.1).
    let baseline_rss_bytes = baseline;
    let rss_residual_bytes = peak_rss
        .saturating_sub(max_ws.bytes)
        .saturating_sub(baseline_rss_bytes);

    // Reproducibility + memory receipt. Written when there is a destination — an
    // explicit --manifest path, or a sidecar next to a `-o` VCF. A stdout run
    // without --manifest writes NO file (no surprise cwd write, no race on a fixed
    // filename) but says how to persist one.
    if let Some(dest) = receipt_dest {
        let mut manifest = new_run_manifest("variants");
        let mut cmd = CommandCapture::new("variants");
        let reference_flag = if loaded.is_legacy_index() {
            "--index"
        } else {
            "--reference-pack"
        };
        cmd.input_hashed(
            reference_flag,
            &analysis_reference.display().to_string(),
            &input_hashes[0],
        );
        cmd.input_hashed(
            "--alignments",
            &alignments_path.display().to_string(),
            &input_hashes[1],
        );
        match &selection {
            rosalind::AnalysisSelection::WholeGenome => {}
            rosalind::AnalysisSelection::Intervals(intervals) => {
                if let Some(region) = intervals.region_origin() {
                    cmd.opt("--region", region);
                } else if let Some(path) = intervals.bed_origin() {
                    cmd.input_hashed("--regions", &path.display().to_string(), &input_hashes[2]);
                }
            }
            rosalind::AnalysisSelection::Shard { count, index, .. } => {
                cmd.opt("--shard-count", *count);
                cmd.opt("--shard-index", *index);
            }
        }
        cmd.opt("--mapq-threshold", mapq_threshold);
        cmd.opt("--quality-threshold", quality_threshold as f64);
        cmd.opt("--max-depth", max_depth);
        cmd.opt("--max-read-len", max_read_len);
        cmd.flag_if(enforce, "--enforce");
        cmd.flag_if(gvcf, "--gvcf");
        if let Some(mb) = memory_budget_mb {
            cmd.opt("--memory-budget-mb", mb);
        }
        cmd.flag_if(require_os_limit, "--require-os-limit");
        cmd.flag_if(force, "--force");
        if let (Some(path), Some(hash)) = (&receipt_output, &output_hash) {
            cmd.output_hashed("-o", &path.display().to_string(), hash);
        }
        cmd.record_into(&mut manifest);
        manifest
            .params
            .insert("pileup.semantics".into(), "exact-or-fail-v1".into());
        manifest.params.insert(
            "artifact.input.0.role".to_string(),
            if loaded.is_legacy_index() {
                "reference-index"
            } else {
                "analysis-reference-pack"
            }
            .to_string(),
        );
        manifest.params.insert(
            "artifact.input.1.role".to_string(),
            "sorted-alignments".to_string(),
        );
        manifest.params.insert(
            "artifact.output.0.format".to_string(),
            if gvcf { "gvcf" } else { "vcf-sites" }.to_string(),
        );
        if selection
            .intervals()
            .and_then(|intervals| intervals.bed_origin())
            .is_some()
        {
            manifest.params.insert(
                "artifact.input.2.role".to_string(),
                "regions-bed".to_string(),
            );
        }
        record_selection_claims(&mut manifest.params, &selection);
        if !manifest.outputs.is_empty() {
            manifest.params.insert(
                "artifact.output.0.role".to_string(),
                if is_breach { "partial-calls" } else { "calls" }.to_string(),
            );
        }
        manifest
            .params
            .insert("model.germline".to_string(), "baseq-mapq-v1".to_string());
        manifest
            .params
            .insert("producer.name".to_string(), "rosalind".to_string());
        manifest.params.insert(
            "producer.version".to_string(),
            env!("CARGO_PKG_VERSION").to_string(),
        );
        manifest
            .params
            .insert("producer.binary".to_string(), "rosalind".to_string());
        manifest.params.insert(
            "analyzer.id".to_string(),
            if gvcf {
                "germline-gvcf"
            } else {
                "germline-sites"
            }
            .to_string(),
        );
        manifest.params.insert(
            "analyzer.version".to_string(),
            env!("CARGO_PKG_VERSION").to_string(),
        );
        manifest.params.insert(
            "analyzer.memory_model".to_string(),
            "fixed-additional-v1".to_string(),
        );
        manifest
            .params
            .insert("analyzer.max_additional_bytes".to_string(), "0".to_string());
        manifest.params.insert(
            "contract.assurance".to_string(),
            if require_os_limit {
                "cgroup-v2"
            } else if enforce {
                "declared-bound-cooperative"
            } else {
                "observed-only"
            }
            .to_string(),
        );
        if let Some(limit) = os_limit_bytes {
            manifest
                .params
                .insert("os.memory_limit_bytes".to_string(), limit.to_string());
        }
        manifest.params.insert(
            "run_status".to_string(),
            if is_breach { "breached" } else { "completed" }.to_string(),
        );
        manifest
            .params
            .insert("peak_rss_bytes".to_string(), peak_rss.to_string());
        manifest.params.insert(
            "predicted_peak_rss_bytes".to_string(),
            predicted_peak.to_string(),
        );
        manifest.params.insert(
            "max_working_set_bytes".to_string(),
            max_ws.bytes.to_string(),
        );
        manifest
            .params
            .insert("governor".to_string(), governor_state.to_string());
        manifest.params.insert(
            "baseline_rss_bytes".to_string(),
            baseline_rss_bytes.to_string(),
        );
        manifest.params.insert(
            "rss_residual_bytes".to_string(),
            rss_residual_bytes.to_string(),
        );
        manifest.params.insert(
            "io_rss_overhead_assumed_bytes".to_string(),
            PILEUP_IO_RSS_OVERHEAD.to_string(),
        );
        // The deterministic working-set PREDICTION (index header + declared caps, no
        // baseline) — a CLAIM field (not in MEASUREMENT_KEYS), so it is cross-machine
        // stable and `verify` can re-check `predicted >= realized` offline.
        manifest.params.insert(
            "predicted_working_set_bytes".to_string(),
            rosalind::call::plan::estimate_variants_working_set(largest, max_depth, max_read_len)
                .bytes
                .to_string(),
        );
        manifest.params.insert(
            "over_max_depth".to_string(),
            skips.over_max_depth.to_string(),
        );
        manifest
            .params
            .insert("reads_skipped_total".to_string(), skips.total().to_string());
        manifest
            .params
            .insert("contract_verdict".to_string(), verdict.to_string());
        if let Some(failure) = &capacity_failure {
            manifest
                .params
                .insert("failure.kind".into(), "capacity-exceeded".into());
            manifest
                .params
                .insert("failure.contig".into(), failure.contig.to_string());
            manifest
                .params
                .insert("failure.position0".into(), failure.position.to_string());
            manifest.params.insert(
                "failure.capacity_reads".into(),
                failure.capacity.to_string(),
            );
            manifest.params.insert(
                "failure.required_reads".into(),
                failure.required.to_string(),
            );
        }
        let mut staged = AtomicFile::create(&dest)?;
        for pass in 0..2 {
            peak_rss = peak_rss.max(
                std::env::var("ROSALIND_FORCE_PEAK_RSS_BYTES")
                    .ok()
                    .and_then(|value| value.parse::<u64>().ok())
                    .unwrap_or_else(peak_rss_bytes),
            );
            if let Err(rosalind::core::CoreError::BudgetExceeded { needed, .. }) =
                rosalind::core::governor::checkpoint()
            {
                breached = true;
                peak_rss = peak_rss.max(needed);
                governor_state = "tripped";
            }
            verdict = match memory_budget_mb.map(|mb| MemoryBudget::from_mb(mb).admits(peak_rss)) {
                None => "unset",
                Some(true) => "within",
                Some(false) => "over",
            };
            is_breach = capacity_failure.is_some() || breached || (enforce && verdict == "over");
            manifest
                .params
                .insert("peak_rss_bytes".into(), peak_rss.to_string());
            manifest
                .params
                .insert("contract_verdict".into(), verdict.into());
            manifest
                .params
                .insert("governor".into(), governor_state.into());
            manifest.params.insert(
                "rss_residual_bytes".into(),
                peak_rss
                    .saturating_sub(max_ws.bytes)
                    .saturating_sub(baseline)
                    .to_string(),
            );
            manifest.params.insert(
                "run_status".into(),
                if capacity_failure.is_some() {
                    "capacity-exceeded"
                } else if is_breach {
                    "breached"
                } else {
                    "completed"
                }
                .into(),
            );
            if let (Some(entry), Some(path)) = (manifest.outputs.first_mut(), &output) {
                entry.path = if is_breach {
                    sidecar_path(path, ".partial")
                } else {
                    path.clone()
                }
                .display()
                .to_string();
                manifest.params.insert(
                    "artifact.output.0.role".into(),
                    if is_breach { "partial-calls" } else { "calls" }.into(),
                );
            }
            manifest.finalize();
            if pass != 0 {
                use std::io::{Seek, SeekFrom};
                staged.file_mut().set_len(0)?;
                staged.file_mut().seek(SeekFrom::Start(0))?;
            }
            staged
                .file_mut()
                .write_all(manifest.to_canonical_json().as_bytes())?;
            staged.file_mut().flush()?;
        }
        staged
            .commit(force)
            .with_context(|| format!("failed to write manifest {}", dest.display()))?;
        eprintln!("wrote reproducibility receipt: {}", dest.display());
    } else {
        eprintln!(
            "no receipt written (stdout output) — pass --manifest <path> or -o <vcf> to persist one"
        );
    }
    if let Some(file) = atomic_output {
        if is_breach {
            file.commit_as(
                &sidecar_path(output.as_ref().expect("file output"), ".partial"),
                force,
            )?;
        } else {
            file.commit(force)?;
        }
    }
    if let Some(failure) = capacity_failure {
        eprintln!("pileup capacity exceeded at contig {}, zero-based position {}: {} active reads exceed --max-depth {}; exact partial output and receipt preserved; raise --max-depth and rerun", failure.contig, failure.position, failure.required, failure.capacity);
        std::process::exit(4);
    }
    // Memory receipt: the bounded contract, made visible + verifiable.
    eprintln!(
        "memory: peak RSS {} MiB; max pileup working set {} KiB",
        peak_rss / (1 << 20),
        max_ws.bytes / 1024
    );
    let other_skipped = skips.total() - skips.over_max_depth;
    if other_skipped > 0 {
        eprintln!(
            "pileup: skipped {other_skipped} reads by filter \
             (unmapped/wrong-contig/secondary/supplementary/duplicate/low-mapq)"
        );
    }
    if let Some(mb) = memory_budget_mb {
        let budget = MemoryBudget::from_mb(mb);
        let within = budget.admits(peak_rss);
        if enforce {
            if within {
                eprintln!(
                    "contract: OK — realized peak {} MiB within declared {mb} MiB",
                    peak_rss / (1 << 20)
                );
            } else {
                eprintln!(
                    "contract: VIOLATED — realized peak {} MiB exceeded declared {mb} MiB (partial output + receipt written)",
                    peak_rss / (1 << 20)
                );
                std::process::exit(4);
            }
        } else if within {
            eprintln!("memory: within budget ({mb} MiB)");
        } else {
            eprintln!(
                "memory: EXCEEDED budget {mb} MiB (realized peak {} MiB) — record-only, run completed",
                peak_rss / (1 << 20)
            );
        }
    }
    Ok(())
}

/// Read a reference FASTA (plain or gzip; `-` = stdin). Phase B1 keeps the
/// single-contig CLI policy: only the first record is used; additional records
/// are warned about (multi-contig consumption is a later phase). The streaming
/// parser itself lives in `io::fasta`.
fn read_fasta(path: &PathBuf) -> Result<FastaRecord> {
    let reader = open_input(path).with_context(|| format!("failed to open {}", path.display()))?;
    let mut records = FastaReader::new(reader);
    let first = records
        .next()
        .ok_or_else(|| anyhow!("FASTA file {} is missing a record", path.display()))?
        .with_context(|| format!("failed to parse FASTA {}", path.display()))?;
    if records.next().is_some() {
        eprintln!(
            "warning: {}: only the first FASTA record is currently used; ignoring the rest \
             (multi-contig lands in a later phase)",
            path.display()
        );
    }
    Ok(first)
}

/// Read ALL FASTA records into a contig-name → sequence map. Used by `eval-*`,
/// which must normalize each variant against its OWN contig — loading only the
/// first record (the single-contig `read_fasta` policy) silently miscompares or
/// crashes on any multi-contig benchmark.
fn read_fasta_map(path: &PathBuf) -> Result<std::collections::BTreeMap<String, Vec<u8>>> {
    let reader = open_input(path).with_context(|| format!("failed to open {}", path.display()))?;
    let mut map = std::collections::BTreeMap::new();
    for rec in FastaReader::new(reader) {
        let rec = rec.with_context(|| format!("failed to parse FASTA {}", path.display()))?;
        if map.insert(rec.name.clone(), rec.sequence).is_some() {
            bail!(
                "FASTA {} has a duplicate contig name '{}'",
                path.display(),
                rec.name
            );
        }
    }
    if map.is_empty() {
        bail!("FASTA file {} is missing a record", path.display());
    }
    Ok(map)
}

/// Read a FASTQ file (plain or gzip; `-` = stdin) into a vector of records.
/// The streaming parser lives in `io::fastq`.
fn read_fastq(path: &PathBuf) -> Result<Vec<FastqRecord>> {
    let reader = open_input(path)
        .with_context(|| format!("failed to open FASTQ file {}", path.display()))?;
    FastqReader::new(reader)
        .collect::<std::result::Result<Vec<_>, _>>()
        .with_context(|| format!("failed to parse FASTQ {}", path.display()))
}

fn normalize_read_name(raw: &str) -> String {
    let name = raw.trim();
    let name = name.split_whitespace().next().unwrap_or(name);
    if name.ends_with("/1") || name.ends_with("/2") {
        name[..name.len() - 2].to_string()
    } else {
        name.to_string()
    }
}

fn read_fastq_pairs(r1: &PathBuf, r2: &PathBuf) -> Result<Vec<FastqPair>> {
    let mut r1_records =
        read_fastq(r1).with_context(|| format!("failed to read FASTQ R1 from {}", r1.display()))?;
    let mut r2_records =
        read_fastq(r2).with_context(|| format!("failed to read FASTQ R2 from {}", r2.display()))?;

    if r1_records.len() != r2_records.len() {
        bail!(
            "paired FASTQ length mismatch: {} has {} records, {} has {} records",
            r1.display(),
            r1_records.len(),
            r2.display(),
            r2_records.len()
        );
    }

    let mut out = Vec::with_capacity(r1_records.len());
    for (mut rec1, mut rec2) in r1_records.drain(..).zip(r2_records.drain(..)) {
        let name1 = normalize_read_name(&rec1.name);
        let name2 = normalize_read_name(&rec2.name);
        if name1 != name2 {
            bail!("mate name mismatch: '{}' vs '{}'", name1, name2);
        }
        rec1.name = name1.clone();
        rec2.name = name1.clone();
        out.push(FastqPair {
            name: name1,
            r1: rec1,
            r2: rec2,
        });
    }
    Ok(out)
}

fn resolve_reads(
    single: Option<PathBuf>,
    r1: Option<PathBuf>,
    r2: Option<PathBuf>,
    label: &str,
) -> Result<ResolvedReads> {
    match (single, r1, r2) {
        (Some(path), None, None) => Ok(ResolvedReads::Single(read_fastq(&path)?)),
        (None, Some(r1), Some(r2)) => Ok(ResolvedReads::Paired(read_fastq_pairs(&r1, &r2)?)),
        (None, None, None) => bail!(
            "missing reads for {}: provide --{} <FASTQ> or --{}-r1/--{}-r2",
            label,
            label,
            label,
            label
        ),
        _ => bail!(
            "invalid reads args for {}: provide either single-end or paired-end (not both/incomplete)",
            label
        ),
    }
}

fn write_sam_alignments<W: Write>(
    writer: &mut W,
    reference_name: &str,
    reference_len: usize,
    reference_offset: u32,
    reads: &[FastqRecord],
    alignments: &[Option<AlignmentCandidate>],
) -> Result<()> {
    writeln!(writer, "@HD\tVN:1.6\tSO:unknown")?;
    writeln!(writer, "@SQ\tSN:{reference_name}\tLN:{reference_len}")?;
    writeln!(
        writer,
        "@PG\tID:rosalind\tPN:rosalind\tVN:{}",
        env!("CARGO_PKG_VERSION")
    )?;

    if reads.len() != alignments.len() {
        bail!("internal error: read-alignment count mismatch");
    }

    for (record, alignment) in reads.iter().zip(alignments.iter()) {
        let seq_str = String::from_utf8(record.sequence.clone())
            .map_err(|_| anyhow!("FASTQ sequence for {} is not valid ASCII", record.name))?;
        let qual_str = String::from_utf8(record.qualities.clone())
            .map_err(|_| anyhow!("FASTQ qualities for {} are not valid ASCII", record.name))?;

        if let Some(hit) = alignment {
            let mapq = hit.mapq;
            let pos = hit.position + reference_offset as usize + 1;
            let flag = if hit.is_reverse { 0x10 } else { 0 };
            let cigar = if hit.cigar.is_empty() {
                format!("{}M", record.sequence.len())
            } else {
                format_cigar(&hit.cigar)
            };
            writeln!(
                writer,
                "{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{seq}\t{qual}\tNM:i:{nm}\tAS:i:{as_score}\tMD:Z:{md}",
                qname = record.name,
                flag = flag,
                rname = reference_name,
                pos = pos,
                cigar = cigar,
                seq = seq_str,
                qual = qual_str,
                nm = hit.mismatches,
                as_score = hit.as_score,
                md = hit.md
            )?;
        } else {
            writeln!(
                writer,
                "{qname}\t{flag}\t*\t0\t0\t*\t*\t0\t0\t{seq}\t{qual}",
                qname = record.name,
                flag = 4,
                seq = seq_str,
                qual = qual_str
            )?;
        }
    }

    writer.flush()?;
    Ok(())
}

fn write_sam_alignments_paired<W: Write>(
    writer: &mut W,
    reference_name: &str,
    reference_len: usize,
    reference_offset: u32,
    reads: &[FastqPair],
    alignments: &[PairAlignments],
) -> Result<()> {
    writeln!(writer, "@HD\tVN:1.6\tSO:unknown")?;
    writeln!(writer, "@SQ\tSN:{reference_name}\tLN:{reference_len}")?;
    writeln!(
        writer,
        "@PG\tID:rosalind\tPN:rosalind\tVN:{}",
        env!("CARGO_PKG_VERSION")
    )?;

    if reads.len() != alignments.len() {
        bail!("internal error: read-alignment count mismatch");
    }

    for (pair, (a1, a2)) in reads.iter().zip(alignments.iter()) {
        write_one_sam_mate(
            writer,
            reference_name,
            reference_offset,
            &pair.r1,
            a1,
            a2,
            true,
        )?;
        write_one_sam_mate(
            writer,
            reference_name,
            reference_offset,
            &pair.r2,
            a2,
            a1,
            false,
        )?;
    }

    writer.flush()?;
    Ok(())
}

fn write_one_sam_mate<W: Write>(
    writer: &mut W,
    reference_name: &str,
    reference_offset: u32,
    record: &FastqRecord,
    this: &Option<AlignmentCandidate>,
    mate: &Option<AlignmentCandidate>,
    is_first: bool,
) -> Result<()> {
    let seq_str = String::from_utf8(record.sequence.clone())
        .map_err(|_| anyhow!("FASTQ sequence for {} is not valid ASCII", record.name))?;
    let qual_str = String::from_utf8(record.qualities.clone())
        .map_err(|_| anyhow!("FASTQ qualities for {} are not valid ASCII", record.name))?;

    let mut flag: u16 = 0x1; // paired
    flag |= if is_first { 0x40 } else { 0x80 };

    let (rname, pos, mapq, cigar) = if let Some(hit) = this {
        if hit.is_reverse {
            flag |= 0x10;
        }
        let pos = hit.position + reference_offset as usize + 1;
        let cigar = if hit.cigar.is_empty() {
            format!("{}M", record.sequence.len())
        } else {
            format_cigar(&hit.cigar)
        };
        (reference_name.to_string(), pos, hit.mapq, cigar)
    } else {
        flag |= 0x4;
        ("*".to_string(), 0usize, 0u8, "*".to_string())
    };

    let (rnext, pnext, tlen, mate_reverse_flag) = if let Some(m) = mate {
        let pnext = m.position + reference_offset as usize + 1;
        let tlen = compute_tlen(this.as_ref(), mate.as_ref(), record.sequence.len());
        let mate_is_rev = m.is_reverse;
        ("=".to_string(), pnext, tlen, mate_is_rev)
    } else {
        flag |= 0x8;
        ("*".to_string(), 0usize, 0i64, false)
    };
    if mate_reverse_flag {
        flag |= 0x20;
    }

    let nm = this.as_ref().map(|h| h.mismatches).unwrap_or(0);
    let as_score = this.as_ref().map(|h| h.as_score).unwrap_or(0);
    let md = this.as_ref().map(|h| h.md.as_str()).unwrap_or("");

    writeln!(
        writer,
        "{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}\tNM:i:{nm}\tAS:i:{as_score}\tMD:Z:{md}",
        qname = record.name,
        flag = flag,
        rname = rname,
        pos = pos,
        mapq = mapq,
        cigar = cigar,
        rnext = rnext,
        pnext = pnext,
        tlen = tlen,
        seq = seq_str,
        qual = qual_str,
        nm = nm,
        as_score = as_score,
        md = md
    )?;
    Ok(())
}

fn compute_tlen(
    this: Option<&AlignmentCandidate>,
    mate: Option<&AlignmentCandidate>,
    read_len: usize,
) -> i64 {
    let (Some(a), Some(b)) = (this, mate) else {
        return 0;
    };
    let a_start = a.position as i64;
    let b_start = b.position as i64;
    let a_end = a_start + read_len as i64;
    let b_end = b_start + read_len as i64;
    let left = a_start.min(b_start);
    let right = a_end.max(b_end);
    let tlen = right - left;
    if a_start <= b_start {
        tlen
    } else {
        -tlen
    }
}

/// Convert FASTQ-ASCII quality bytes (Phred+33, e.g. `b'I'` = 73 for Q40) to
/// raw Phred values (e.g. 40) as required by the BAM QUAL field.
///
/// The FASTQ parser stores quality bytes verbatim from the quality line, so
/// in-memory `FastqRecord.qualities` are ASCII-encoded (Phred+33). BAM/htslib
/// expects raw Phred (0–93) in `Record::set`; this function decodes them.
fn fastq_quals_to_phred(ascii: &[u8]) -> Vec<u8> {
    ascii.iter().map(|q| q.saturating_sub(33)).collect()
}

fn write_bam_alignments(
    writer: &mut bam::Writer,
    reference_offset: u32,
    reads: &[FastqRecord],
    alignments: &[Option<AlignmentCandidate>],
) -> Result<()> {
    if reads.len() != alignments.len() {
        bail!("internal error: read-alignment count mismatch");
    }

    for (record, alignment) in reads.iter().zip(alignments.iter()) {
        let phred_quals = fastq_quals_to_phred(&record.qualities);
        if let Some(hit) = alignment {
            let mut bam_record = Record::new();
            let cigar_ops = if hit.cigar.is_empty() {
                vec![BamCigar::Match(record.sequence.len() as u32)]
            } else {
                hit.cigar
                    .iter()
                    .map(|op| match op.kind {
                        CigarOpKind::Match => BamCigar::Match(op.len),
                        CigarOpKind::Insertion => BamCigar::Ins(op.len),
                        CigarOpKind::Deletion => BamCigar::Del(op.len),
                        CigarOpKind::SoftClip => BamCigar::SoftClip(op.len),
                        CigarOpKind::HardClip => BamCigar::HardClip(op.len),
                    })
                    .collect()
            };
            let cigar = CigarString::from(cigar_ops);
            bam_record.set(
                record.name.as_bytes(),
                Some(&cigar),
                &record.sequence,
                &phred_quals,
            );
            bam_record.set_tid(0);
            bam_record.set_pos((hit.position + reference_offset as usize) as i64);
            let flags = if hit.is_reverse { 0x10 } else { 0 };
            bam_record.set_flags(flags);
            bam_record.set_mapq(hit.mapq);
            bam_record.set_mtid(-1);
            bam_record.set_mpos(-1);
            bam_record.set_insert_size(0);
            bam_record.push_aux(b"NM", Aux::I32(hit.mismatches as i32))?;
            bam_record.push_aux(b"AS", Aux::I32(hit.as_score))?;
            if !hit.md.is_empty() {
                bam_record.push_aux(b"MD", Aux::String(hit.md.as_str()))?;
            }
            writer.write(&bam_record)?;
        } else {
            let mut bam_record = Record::new();
            bam_record.set(record.name.as_bytes(), None, &record.sequence, &phred_quals);
            bam_record.set_tid(-1);
            bam_record.set_pos(-1);
            bam_record.set_flags(0x4);
            bam_record.set_mapq(0);
            bam_record.set_mtid(-1);
            bam_record.set_mpos(-1);
            bam_record.set_insert_size(0);
            writer.write(&bam_record)?;
        }
    }

    Ok(())
}

fn write_bam_alignments_paired(
    writer: &mut bam::Writer,
    reference_offset: u32,
    reads: &[FastqPair],
    alignments: &[PairAlignments],
) -> Result<()> {
    if reads.len() != alignments.len() {
        bail!("internal error: read-alignment count mismatch");
    }

    for (pair, (a1, a2)) in reads.iter().zip(alignments.iter()) {
        write_one_bam_mate(writer, reference_offset, &pair.r1, a1, a2, true)?;
        write_one_bam_mate(writer, reference_offset, &pair.r2, a2, a1, false)?;
    }
    Ok(())
}

fn write_one_bam_mate(
    writer: &mut bam::Writer,
    reference_offset: u32,
    record: &FastqRecord,
    this: &Option<AlignmentCandidate>,
    mate: &Option<AlignmentCandidate>,
    is_first: bool,
) -> Result<()> {
    let mut bam_record = Record::new();

    let mut flags: u16 = 0x1; // paired
    flags |= if is_first { 0x40 } else { 0x80 };

    let (tid, pos) = if let Some(hit) = this {
        if hit.is_reverse {
            flags |= 0x10;
        }
        (0i32, (hit.position + reference_offset as usize) as i64)
    } else {
        flags |= 0x4;
        (-1i32, -1i64)
    };

    if mate.is_none() {
        flags |= 0x8;
    } else if mate.as_ref().is_some_and(|m| m.is_reverse) {
        flags |= 0x20;
    }

    let cigar_ops = if let Some(hit) = this {
        if hit.cigar.is_empty() {
            vec![BamCigar::Match(record.sequence.len() as u32)]
        } else {
            hit.cigar
                .iter()
                .map(|op| match op.kind {
                    CigarOpKind::Match => BamCigar::Match(op.len),
                    CigarOpKind::Insertion => BamCigar::Ins(op.len),
                    CigarOpKind::Deletion => BamCigar::Del(op.len),
                    CigarOpKind::SoftClip => BamCigar::SoftClip(op.len),
                    CigarOpKind::HardClip => BamCigar::HardClip(op.len),
                })
                .collect()
        }
    } else {
        Vec::new()
    };
    let cigar = if cigar_ops.is_empty() {
        None
    } else {
        Some(CigarString::from(cigar_ops))
    };

    let phred_quals = fastq_quals_to_phred(&record.qualities);
    bam_record.set(
        record.name.as_bytes(),
        cigar.as_ref(),
        &record.sequence,
        &phred_quals,
    );
    bam_record.set_flags(flags);
    bam_record.set_tid(tid);
    bam_record.set_pos(pos);
    bam_record.set_mapq(this.as_ref().map(|h| h.mapq).unwrap_or(0));

    if let Some(m) = mate {
        bam_record.set_mtid(0);
        bam_record.set_mpos((m.position + reference_offset as usize) as i64);
        bam_record.set_insert_size(compute_tlen(
            this.as_ref(),
            mate.as_ref(),
            record.sequence.len(),
        ));
    } else {
        bam_record.set_mtid(-1);
        bam_record.set_mpos(-1);
        bam_record.set_insert_size(0);
    }

    if let Some(hit) = this {
        bam_record.push_aux(b"NM", Aux::I32(hit.mismatches as i32))?;
        bam_record.push_aux(b"AS", Aux::I32(hit.as_score))?;
        if !hit.md.is_empty() {
            bam_record.push_aux(b"MD", Aux::String(hit.md.as_str()))?;
        }
    }

    writer.write(&bam_record)?;
    Ok(())
}

fn format_cigar(cigar: &[CigarOp]) -> String {
    let mut out = String::new();
    for op in cigar {
        out.push_str(&op.len.to_string());
        let ch = match op.kind {
            CigarOpKind::Match => 'M',
            CigarOpKind::Insertion => 'I',
            CigarOpKind::Deletion => 'D',
            CigarOpKind::SoftClip => 'S',
            CigarOpKind::HardClip => 'H',
        };
        out.push(ch);
    }
    out
}

fn read_alignment_file(
    path: &PathBuf,
    target_chrom: Option<&Arc<str>>,
) -> Result<Vec<AlignedRead>> {
    if path
        .extension()
        .and_then(|ext| ext.to_str())
        .map(|ext| ext.eq_ignore_ascii_case("bam"))
        .unwrap_or(false)
    {
        return read_bam_alignment_file(path, target_chrom);
    }

    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut reads = Vec::new();

    for (line_no, line) in reader.lines().enumerate() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('@') {
            continue;
        }

        let fields: Vec<&str> = trimmed.split('\t').collect();
        if fields.len() < 11 {
            bail!(
                "SAM record on line {} has {} fields (expected ≥ 11)",
                line_no + 1,
                fields.len()
            );
        }

        let flag: u16 = fields[1]
            .parse()
            .with_context(|| format!("invalid FLAG '{}' on line {}", fields[1], line_no + 1))?;
        if flag & 0x4 != 0 {
            continue;
        }

        let rname = fields[2];
        if rname == "*" {
            continue;
        }
        if let Some(target) = target_chrom {
            if rname != target.as_ref() {
                continue;
            }
        }

        let pos: u32 = fields[3]
            .parse()
            .with_context(|| format!("invalid POS '{}' on line {}", fields[3], line_no + 1))?;
        let mapq: u8 = fields[4]
            .parse()
            .with_context(|| format!("invalid MAPQ '{}' on line {}", fields[4], line_no + 1))?;
        let cigar = parse_cigar(fields[5], fields[9].len())
            .with_context(|| format!("invalid CIGAR '{}' on line {}", fields[5], line_no + 1))?;

        let sequence = fields[9].to_ascii_uppercase().into_bytes();
        let qual_field = fields[10].as_bytes();
        if qual_field.len() != sequence.len() {
            bail!(
                "sequence/quality length mismatch on line {} ({} vs {})",
                line_no + 1,
                sequence.len(),
                qual_field.len()
            );
        }
        let qualities: Vec<u8> = qual_field.iter().map(|q| q.saturating_sub(33)).collect();

        let is_reverse = flag & 0x10 != 0;

        reads.push(AlignedRead::new(
            rname.to_string(),
            pos.saturating_sub(1),
            mapq,
            cigar,
            sequence,
            qualities,
            is_reverse,
        ));
    }

    Ok(reads)
}

/// Convert a legacy `genomics::AlignedRead` into a canonical `core::AlignedRead`
/// for the new calling vertical. The legacy type has no RefSkip/Pad CIGAR ops.
fn legacy_read_to_core(r: &AlignedRead, contig: u32) -> rosalind::core::AlignedRead {
    use rosalind::core::{CigarOp as CoreOp, CigarOpKind as CoreKind};
    let cigar = r
        .cigar
        .iter()
        .map(|op| {
            let kind = match op.kind {
                CigarOpKind::Match => CoreKind::Match,
                CigarOpKind::Insertion => CoreKind::Insertion,
                CigarOpKind::Deletion => CoreKind::Deletion,
                CigarOpKind::SoftClip => CoreKind::SoftClip,
                CigarOpKind::HardClip => CoreKind::HardClip,
            };
            CoreOp::new(kind, op.len)
        })
        .collect();
    let flags = if r.is_reverse {
        rosalind::core::SamFlags(rosalind::core::SamFlags::REVERSE)
    } else {
        rosalind::core::SamFlags::default()
    };
    rosalind::core::AlignedRead {
        contig,
        pos: rosalind::core::Position(r.pos),
        mapq: r.mapq,
        flags,
        cigar,
        seq: std::sync::Arc::clone(&r.sequence),
        qual: std::sync::Arc::clone(&r.qualities),
    }
}

fn read_bam_alignment_file(
    path: &PathBuf,
    target_chrom: Option<&Arc<str>>,
) -> Result<Vec<AlignedRead>> {
    let mut reader = bam::Reader::from_path(path)
        .with_context(|| format!("failed to open BAM {}", path.display()))?;
    let header = reader.header().to_owned();

    let mut reads = Vec::new();
    for record_result in reader.records() {
        let record = record_result?;

        if record.is_unmapped() {
            continue;
        }

        let tid = record.tid();
        if tid < 0 {
            continue;
        }
        let rname_bytes = header.tid2name(tid as u32);
        let rname = std::str::from_utf8(rname_bytes)
            .map_err(|_| anyhow!("BAM reference name is not valid UTF-8"))?;
        if let Some(target) = target_chrom {
            if rname != target.as_ref() {
                continue;
            }
        }

        let pos0 = record.pos();
        if pos0 < 0 {
            continue;
        }

        let mut cigar = Vec::new();
        for c in record.cigar().iter() {
            let (kind, len) = match *c {
                BamCigar::Match(l) | BamCigar::Equal(l) | BamCigar::Diff(l) => {
                    (CigarOpKind::Match, l)
                }
                BamCigar::Ins(l) => (CigarOpKind::Insertion, l),
                BamCigar::Del(l) => (CigarOpKind::Deletion, l),
                BamCigar::SoftClip(l) => (CigarOpKind::SoftClip, l),
                BamCigar::HardClip(l) => (CigarOpKind::HardClip, l),
                _ => continue,
            };
            cigar.push(CigarOp::new(kind, len));
        }

        let seq = record.seq().as_bytes();
        let sequence: Vec<u8> = seq.iter().map(|b| b.to_ascii_uppercase()).collect();
        let qualities: Vec<u8> = record.qual().to_vec();

        reads.push(AlignedRead::new(
            rname.to_string(),
            pos0 as u32,
            record.mapq(),
            cigar,
            sequence,
            qualities,
            record.is_reverse(),
        ));
    }

    Ok(reads)
}

fn parse_cigar(cigar: &str, read_len: usize) -> Result<Vec<CigarOp>> {
    if cigar == "*" {
        bail!("variant calling requires mapped reads (CIGAR cannot be '*')");
    }

    let mut ops = Vec::new();
    let mut run = String::new();
    let mut read_consuming: usize = 0;

    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            run.push(ch);
            continue;
        }
        let len: u32 = run
            .parse()
            .with_context(|| format!("invalid CIGAR run length in '{cigar}'"))?;
        run.clear();
        // Map operations the same way `read_bam_alignment_file` does; skip ops
        // that `CigarOpKind` does not represent (e.g. N/P), mirroring that reader.
        let kind = match ch {
            'M' | '=' | 'X' => CigarOpKind::Match,
            'I' => CigarOpKind::Insertion,
            'D' => CigarOpKind::Deletion,
            'S' => CigarOpKind::SoftClip,
            'H' => CigarOpKind::HardClip,
            _ => continue,
        };
        // M/=/X, I, and S consume read bases; D and H do not.
        if matches!(
            kind,
            CigarOpKind::Match | CigarOpKind::Insertion | CigarOpKind::SoftClip
        ) {
            read_consuming += len as usize;
        }
        ops.push(CigarOp::new(kind, len));
    }

    if !run.is_empty() {
        bail!("CIGAR '{cigar}' ends with a run length but no operation");
    }
    if ops.is_empty() {
        bail!("CIGAR '{cigar}' contains no supported operations");
    }
    if read_consuming != read_len {
        bail!("CIGAR '{cigar}' consumes {read_consuming} read bases but the read has {read_len}");
    }

    Ok(ops)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::env;
    use std::fs;
    use std::path::PathBuf;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn temp_file_path(suffix: &str) -> PathBuf {
        let timestamp = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("time went backwards")
            .as_nanos();
        let mut path = env::temp_dir();
        path.push(format!("rosalind-test-{suffix}-{timestamp}.tmp"));
        path
    }

    fn write_temp_file(contents: &str, suffix: &str) -> PathBuf {
        let path = temp_file_path(suffix);
        fs::write(&path, contents).expect("failed to write temp file");
        path
    }

    #[test]
    fn parse_cigar_handles_multi_op_indels() {
        // The exact CIGAR the SAM reader previously choked on (regression).
        // 146M + 2I + 2M consumes 150 read bases.
        let ops = parse_cigar("146M2I2M", 150).expect("multi-op CIGAR should parse");
        assert_eq!(
            ops,
            vec![
                CigarOp::new(CigarOpKind::Match, 146),
                CigarOp::new(CigarOpKind::Insertion, 2),
                CigarOp::new(CigarOpKind::Match, 2),
            ]
        );
    }

    #[test]
    fn parse_cigar_handles_match_softclip_and_deletion() {
        assert_eq!(
            parse_cigar("150M", 150).unwrap(),
            vec![CigarOp::new(CigarOpKind::Match, 150)]
        );
        // Soft-clips consume read bases (5 + 140 + 5 = 150).
        let sc = parse_cigar("5S140M5S", 150).unwrap();
        assert_eq!(sc[0], CigarOp::new(CigarOpKind::SoftClip, 5));
        assert_eq!(sc[2], CigarOp::new(CigarOpKind::SoftClip, 5));
        // Deletions consume reference only, so read-consuming length is 20.
        let del = parse_cigar("10M2D10M", 20).unwrap();
        assert_eq!(del[1], CigarOp::new(CigarOpKind::Deletion, 2));
    }

    #[test]
    fn parse_cigar_rejects_length_mismatch_and_star() {
        assert!(parse_cigar("100M", 150).is_err());
        assert!(parse_cigar("*", 0).is_err());
    }

    #[test]
    fn fasta_parser_extracts_primary_name() {
        let path = write_temp_file(">chr1 CP068277.2 description\nACGT\n", "fasta");
        let record = read_fasta(&path).expect("FASTA read should succeed");
        assert_eq!(record.name, "chr1");
        let _ = fs::remove_file(path);
    }

    #[test]
    fn fastq_parser_trims_after_space() {
        let contents = "@read1 1:N:0:CG\nAC\n+\n!!\n";
        let path = write_temp_file(contents, "fastq");
        let records = read_fastq(&path).expect("FASTQ read should succeed");
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].name, "read1");
        let _ = fs::remove_file(path);
    }

    #[test]
    fn align_reads_with_fm_index() {
        let reference = b"ACGTACGT";
        let mut aligner =
            BWTAligner::new(reference).expect("should build FM-index aligner for test");
        let records = vec![FastqRecord {
            name: "read1".to_string(),
            sequence: b"ACGT".to_vec(),
            qualities: b"IIII".to_vec(),
        }];
        let alignments = align_reads(&mut aligner, &records, 2).expect("alignment should run");
        assert!(alignments[0].is_some());
    }

    // ── fastq_quals_to_phred unit tests ────────────────────────────────────────
    // BAM/htslib expect raw Phred (0–93) in Record::set; the FASTQ parser keeps
    // ASCII-encoded (Phred+33) bytes. These tests guard the conversion helper.

    #[test]
    fn fastq_quals_to_phred_decodes_single_byte() {
        // ASCII b'I' = 73; Phred = 73 - 33 = 40 (Q40, a typical high-quality base).
        assert_eq!(fastq_quals_to_phred(b"I"), vec![40]);
    }

    #[test]
    fn fastq_quals_to_phred_decodes_lowest_quality() {
        // ASCII b'!' = 33; Phred = 33 - 33 = 0.
        assert_eq!(fastq_quals_to_phred(b"!"), vec![0]);
    }

    #[test]
    fn fastq_quals_to_phred_decodes_run() {
        // b"!I~" → [0, 40, 93]  (b'~' = 126; 126 - 33 = 93, the maximum valid Phred score).
        assert_eq!(fastq_quals_to_phred(b"!I~"), vec![0, 40, 93]);
    }

    #[test]
    fn fastq_quals_to_phred_empty_input_returns_empty() {
        assert_eq!(fastq_quals_to_phred(b""), Vec::<u8>::new());
    }

    #[test]
    fn fastq_quals_to_phred_saturating_does_not_underflow() {
        // A byte below 33 should saturate to 0 rather than wrapping (defensive).
        assert_eq!(fastq_quals_to_phred(&[0u8]), vec![0]);
        assert_eq!(fastq_quals_to_phred(&[32u8]), vec![0]);
    }
}
