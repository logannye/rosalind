//! CLI lifecycle for exact evidence. Scientific identity excludes execution knobs.

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{self, Read, Write};
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};

use anyhow::{bail, Result};
use clap::Args;
use rosalind::core::governor::{checkpoint, MemoryGovernor};
use rosalind::core::CoreError;
use rosalind::dataset::InputSnapshot;
use rosalind::evidence::*;
use rosalind::provenance::{CommandCapture, RunManifest};
use rosalind::util::atomic::{ensure_destination, AtomicFile};
use rosalind::util::rss::peak_rss_bytes;

use crate::{FeatureFormat, SelectionArgs};
mod output;
use output::execute_outputs;

#[derive(Args, Debug, Clone)]
pub(crate) struct EvidenceOptions {
    /// Supplied plain-text SNV VCF; mutually exclusive with BED selection.
    #[arg(long, conflicts_with_all = ["regions", "region", "shard_count", "shard_index"])]
    sites: Option<PathBuf>,
    /// Base quality threshold for exact evidence (unavailable qualities excluded).
    #[arg(long, default_value_t = 20)]
    base_quality_threshold: u8,
    /// Exact evidence output format (panel summaries use TSV).
    #[arg(long, value_enum, default_value_t = FeatureFormat::Tsv)]
    format: FeatureFormat,
    /// Maximum execution microtile width; does not change scientific results.
    #[arg(long, default_value_t = 16_384)]
    tile_bases: u32,
    /// Maximum decoded record envelope; larger records cause a resource failure.
    #[arg(long, default_value_t = 1_048_576)]
    max_record_bytes: usize,
    /// Bounded workers for independent canonical partitions.
    #[arg(long, default_value_t = 1)]
    workers: usize,
    /// Opt-in verified local evidence cache.
    #[arg(long)]
    cache_dir: Option<PathBuf>,
    /// Reuse compatible verified completed partitions; requires --cache-dir.
    #[arg(long, requires = "cache_dir")]
    resume: bool,
    /// Explicit local FASTA used to decode CRAM when the analysis reference is a pack.
    #[arg(long)]
    cram_reference: Option<PathBuf>,
    /// Explicit BAI/CSI/CRAI path, useful when replay inputs were relocated.
    #[arg(long)]
    alignment_index: Option<PathBuf>,
    /// Explicit FAI for an indexed FASTA analysis reference.
    #[arg(long)]
    reference_fai: Option<PathBuf>,
    /// Explicit FAI for the CRAM decoder FASTA.
    #[arg(long)]
    cram_reference_fai: Option<PathBuf>,
    /// Minimum callable A/C/G/T read depth for a callable panel position.
    #[arg(long, default_value_t = 10)]
    min_callable_depth: u64,
    /// Optional exact per-position evidence from the same panel traversal.
    #[arg(long)]
    position_output: Option<PathBuf>,
    /// Print the complete admitted evidence plan as JSON without creating output.
    #[arg(long)]
    plan: bool,
}

impl EvidenceOptions {
    pub(crate) fn reject_legacy_options(&self) -> Result<()> {
        if self.sites.is_some()
            || self.cram_reference.is_some()
            || self.alignment_index.is_some()
            || self.reference_fai.is_some()
            || self.cram_reference_fai.is_some()
            || self.position_output.is_some()
            || self.cache_dir.is_some()
            || self.resume
            || self.workers != 1
            || self.plan
            || self.format != FeatureFormat::Tsv
            || self.base_quality_threshold != 20
            || self.tile_bases != 16_384
            || self.max_record_bytes != 1_048_576
            || self.min_callable_depth != 10
        {
            bail!("evidence-specific options require analyze evidence or analyze panel-qc");
        }
        Ok(())
    }
}

pub(crate) struct EvidenceCommand {
    pub panel: bool,
    pub reference: Option<PathBuf>,
    pub alignments: PathBuf,
    pub mapq_threshold: Option<u8>,
    pub memory_budget_mb: Option<u64>,
    pub max_read_len: u32,
    pub enforce: bool,
    pub require_os_limit: bool,
    pub force: bool,
    pub output: Option<PathBuf>,
    pub manifest: Option<PathBuf>,
    pub selection: SelectionArgs,
    pub options: EvidenceOptions,
}

pub(crate) fn run(command: EvidenceCommand) -> Result<()> {
    match run_inner(command) {
        Ok(()) => Ok(()),
        Err(error) => {
            if let Some(dataset_error) = error.downcast_ref::<rosalind::dataset::DatasetError>() {
                eprintln!("{dataset_error}");
                std::process::exit(dataset_error.exit_code());
            }
            let code = match error.downcast_ref::<EvidenceError>() {
                Some(EvidenceError::Refused { .. }) => 3,
                Some(EvidenceError::RecordLimit(_))
                | Some(EvidenceError::Core(CoreError::BudgetExceeded { .. })) => 4,
                Some(EvidenceError::Analyzer(_)) | Some(EvidenceError::CounterOverflow) => 1,
                _ => 2,
            };
            eprintln!("{error:#}");
            std::process::exit(code);
        }
    }
}

fn hash_file(path: &Path) -> Result<String> {
    let mut input = File::open(path)?;
    let mut hasher = blake3::Hasher::new();
    let mut buffer = [0u8; 65_536];
    loop {
        checkpoint().map_err(EvidenceError::from)?;
        let n = input.read(&mut buffer)?;
        if n == 0 {
            break;
        }
        hasher.update(&buffer[..n]);
    }
    Ok(hasher.finalize().to_hex().to_string())
}

fn is_cram(path: &Path) -> Result<bool> {
    let mut magic = [0u8; 4];
    let read = File::open(path)?.read(&mut magic)?;
    Ok(read == 4 && &magic == b"CRAM")
}

fn is_fasta(path: &Path) -> Result<bool> {
    let mut magic = [0u8; 1];
    Ok(File::open(path)?.read(&mut magic)? == 1 && magic[0] == b'>')
}

fn alignment_sidecar(path: &Path) -> Result<Option<PathBuf>> {
    let cram = is_cram(path)?;
    let extensions: &[&str] = if cram { &["crai"] } else { &["bai", "csi"] };
    Ok(extensions
        .iter()
        .flat_map(|extension| {
            [
                PathBuf::from(format!("{}.{}", path.display(), extension)),
                path.with_extension(extension),
            ]
        })
        .find(|candidate| candidate.is_file()))
}

fn fasta_sidecar(path: &Path) -> Option<PathBuf> {
    let candidate = PathBuf::from(format!("{}.fai", path.display()));
    candidate.is_file().then_some(candidate)
}

fn structured_digest(domain: &str, fields: &BTreeMap<String, String>) -> String {
    let mut hash = blake3::Hasher::new();
    hash.update(domain.as_bytes());
    hash.update(&[0]);
    for (key, value) in fields {
        hash.update(&(key.len() as u64).to_le_bytes());
        hash.update(key.as_bytes());
        hash.update(&(value.len() as u64).to_le_bytes());
        hash.update(value.as_bytes());
    }
    hash.finalize().to_hex().to_string()
}

fn resolved_destination(path: &Path) -> Result<PathBuf> {
    if path.exists() {
        if !path.is_file() {
            bail!("destination must be a file: {}", path.display());
        }
        return Ok(std::fs::canonicalize(path)?);
    }
    let parent = path
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let name = path
        .file_name()
        .ok_or_else(|| anyhow::anyhow!("destination has no filename: {}", path.display()))?;
    Ok(std::fs::canonicalize(parent)?.join(name))
}

fn same_file(left: &Path, right: &Path) -> bool {
    if left == right {
        return true;
    }
    #[cfg(unix)]
    {
        use std::os::unix::fs::MetadataExt;
        if let (Ok(left), Ok(right)) = (std::fs::metadata(left), std::fs::metadata(right)) {
            return left.dev() == right.dev() && left.ino() == right.ino();
        }
    }
    false
}

fn validate_destinations(command: &EvidenceCommand, receipt: Option<&Path>) -> Result<()> {
    let artifacts: Vec<_> = command
        .output
        .iter()
        .chain(command.options.position_output.iter())
        .collect();
    let partials: Vec<_> = artifacts
        .iter()
        .map(|p| PathBuf::from(format!("{}.partial", p.display())))
        .collect();
    let inputs: Vec<_> = [
        Some(&command.alignments),
        command.reference.as_ref(),
        command.options.alignment_index.as_ref(),
        command.options.reference_fai.as_ref(),
        command.options.cram_reference.as_ref(),
        command.options.cram_reference_fai.as_ref(),
        command.options.sites.as_ref(),
        command.selection.regions.as_ref(),
    ]
    .into_iter()
    .flatten()
    .map(std::fs::canonicalize)
    .collect::<std::io::Result<_>>()?;
    let mut destinations = Vec::<PathBuf>::new();
    for path in artifacts
        .into_iter()
        .map(PathBuf::as_path)
        .chain(receipt)
        .chain(partials.iter().map(PathBuf::as_path))
    {
        let resolved = resolved_destination(path)?;
        if inputs.iter().any(|input| same_file(&resolved, input)) {
            bail!(
                "output, receipt, and partial destinations must not overwrite an input: {}",
                path.display()
            );
        }
        if destinations.iter().any(|other| same_file(&resolved, other)) {
            bail!(
                "output, receipt, and partial destinations must use different paths: {}",
                path.display()
            );
        }
        ensure_destination(path, command.force)?;
        destinations.push(resolved);
    }
    Ok(())
}

fn run_inner(mut command: EvidenceCommand) -> Result<()> {
    let started = Instant::now();
    let initial_rss = peak_rss_bytes();
    let budget = command
        .memory_budget_mb
        .map(|mb| {
            mb.checked_mul(1 << 20)
                .ok_or_else(|| anyhow::anyhow!("memory budget is too large"))
        })
        .transpose()?;
    if command.enforce && budget.is_none() {
        bail!("--enforce requires --memory-budget-mb");
    }
    let os_limit = if command.require_os_limit {
        let budget =
            budget.ok_or_else(|| anyhow::anyhow!("--require-os-limit requires a budget"))?;
        match rosalind::contract::detected_os_memory_limit_bytes() {
            Some(limit) if limit <= budget => Some(limit),
            _ => return Err(EvidenceError::InvalidRequest("OS enforcement requires a Linux cgroup-v2 memory.max at or below the declared budget".into()).into()),
        }
    } else {
        None
    };
    let _governor = budget
        .map(|bytes| MemoryGovernor::start(bytes, Duration::from_millis(100), peak_rss_bytes))
        .transpose()?;
    checkpoint().map_err(EvidenceError::from)?;
    if command.selection.region.is_some() || command.selection.shard_count.is_some() {
        bail!("exact evidence requires --sites VCF or --regions BED; execution tiles do not select scientific loci");
    }
    if command.options.sites.is_some() == command.selection.regions.is_some() {
        bail!("supply exactly one of --sites VCF or --regions BED");
    }
    if command.panel && command.options.sites.is_some() {
        bail!("panel-qc requires --regions BED");
    }
    if !command.panel && command.reference.is_none() {
        bail!("evidence requires --reference, --reference-pack, or --index");
    }
    if !command.panel && command.options.position_output.is_some() {
        bail!("--position-output is a panel-qc option");
    }
    if command.panel && command.options.format != FeatureFormat::Tsv {
        bail!("panel summaries use --format tsv; optional per-position evidence uses Arrow");
    }
    if command.options.workers == 0 {
        bail!("--workers must be positive");
    }
    let receipt_path = command.manifest.clone().or_else(|| {
        command
            .output
            .as_ref()
            .map(|p| PathBuf::from(format!("{}.manifest.json", p.display())))
    });
    if !command.panel && command.options.min_callable_depth != 10 {
        bail!("--min-callable-depth is a panel-qc option");
    }
    let cram = is_cram(&command.alignments)?;
    if !cram
        && (command.options.cram_reference.is_some()
            || command.options.cram_reference_fai.is_some())
    {
        bail!("--cram-reference and --cram-reference-fai require CRAM input");
    }
    if command.options.cram_reference_fai.is_some() && command.options.cram_reference.is_none() {
        bail!("--cram-reference-fai requires --cram-reference; use --reference-fai for the analysis FASTA");
    }
    let fasta = command
        .reference
        .as_deref()
        .map(is_fasta)
        .transpose()?
        .unwrap_or(false);
    if command.options.reference_fai.is_some() && !fasta {
        bail!("--reference-fai requires an indexed FASTA analysis reference");
    }
    if command.options.alignment_index.is_none() {
        command.options.alignment_index = alignment_sidecar(&command.alignments)?;
    }
    if command.options.alignment_index.is_none() {
        bail!("indexed evidence requires a BAI/CSI beside BAM or CRAI beside CRAM; create it with samtools index or pass --alignment-index");
    }
    if fasta {
        let reference = command.reference.as_ref().unwrap();
        command.options.reference_fai = command
            .options
            .reference_fai
            .or_else(|| fasta_sidecar(reference));
    }
    if let Some(reference) = &command.options.cram_reference {
        command.options.cram_reference_fai = command
            .options
            .cram_reference_fai
            .or_else(|| fasta_sidecar(reference));
    }
    validate_destinations(&command, receipt_path.as_deref())?;
    let label = if command.panel {
        "analyze panel-qc"
    } else {
        "analyze evidence"
    };
    let mut capture = CommandCapture::new(label);
    let mut hashes = BTreeMap::<PathBuf, String>::new();
    let mut identities = BTreeMap::<String, String>::new();
    let inputs = [
        ("--alignments", Some(&command.alignments)),
        ("--reference", command.reference.as_ref()),
        (
            "--alignment-index",
            command.options.alignment_index.as_ref(),
        ),
        ("--reference-fai", command.options.reference_fai.as_ref()),
        ("--cram-reference", command.options.cram_reference.as_ref()),
        (
            "--cram-reference-fai",
            command.options.cram_reference_fai.as_ref(),
        ),
        ("--sites", command.options.sites.as_ref()),
        ("--regions", command.selection.regions.as_ref()),
    ];
    let input_snapshot = InputSnapshot::capture(
        inputs
            .iter()
            .filter_map(|(_, path)| path.map(|path| path.to_path_buf())),
    )?;
    for (flag, path) in inputs {
        if let Some(path) = path {
            let key = std::fs::canonicalize(path)?;
            let hash = match hashes.get(&key) {
                Some(hash) => hash.clone(),
                None => {
                    let hash = hash_file(path)?;
                    hashes.insert(key.clone(), hash.clone());
                    hash
                }
            };
            capture.input_hashed(flag, &key.display().to_string(), &hash);
            identities.insert(flag.trim_start_matches('-').to_string(), hash);
        }
    }
    input_snapshot.verify()?;
    let hashing_ms = started.elapsed().as_millis();
    let mut request = if let Some(reference) = &command.reference {
        EvidenceRequest::new(&command.alignments, reference)
    } else {
        EvidenceRequest::coverage(&command.alignments)
    };
    request.cram_reference = command.options.cram_reference.clone();
    request.alignment_index = command.options.alignment_index.clone();
    request.reference_fai = command.options.reference_fai.clone();
    request.cram_reference_fai = command.options.cram_reference_fai.clone();
    request.profile.min_mapq = command.mapq_threshold.unwrap_or(20);
    request.profile.min_base_quality = command.options.base_quality_threshold;
    request.execution.memory_budget_bytes = budget;
    request.execution.max_microtile_bases = command.options.tile_bases;
    request.execution.max_read_len = command.max_read_len as usize;
    request.execution.max_record_bytes = command.options.max_record_bytes;
    let mut engine = EvidenceEngine::open(request)?;
    let selection = if let Some(sites) = &command.options.sites {
        EvidenceSelection::from_vcf(sites, engine.contigs())?
    } else {
        EvidenceSelection::from_bed(
            command.selection.regions.as_ref().unwrap(),
            engine.contigs(),
        )?
    };
    engine.set_selection(selection)?;
    // The canonical selection digest is independent of BED ordering/overlap and
    // VCF row ordering. File identities remain separately byte-verifiable.
    let selection_digest = engine.selection_digest();
    let profile = &engine.request().profile;
    let mut science = BTreeMap::from([
        (
            "extractor_semantics".to_string(),
            EVIDENCE_SEMANTICS_VERSION.to_string(),
        ),
        (
            "package_version".to_string(),
            env!("CARGO_PKG_VERSION").to_string(),
        ),
        ("schema".to_string(), EVIDENCE_SCHEMA_VERSION.to_string()),
        (
            "fields_version".to_string(),
            EvidenceFields::VERSION.to_string(),
        ),
        (
            "fields".to_string(),
            engine.request().fields.bits().to_string(),
        ),
        ("profile".to_string(), EvidenceProfile::ID.to_string()),
        ("counting_unit".to_string(), "read".to_string()),
        ("selection".to_string(), selection_digest),
        ("mapq".to_string(), profile.min_mapq.to_string()),
        (
            "base_quality".to_string(),
            profile.min_base_quality.to_string(),
        ),
        (
            "exclude_secondary".to_string(),
            profile.exclude_secondary.to_string(),
        ),
        (
            "exclude_supplementary".to_string(),
            profile.exclude_supplementary.to_string(),
        ),
        (
            "exclude_qc_fail".to_string(),
            profile.exclude_qc_fail.to_string(),
        ),
        (
            "exclude_duplicates".to_string(),
            profile.exclude_duplicates.to_string(),
        ),
        ("unavailable_qualities".to_string(), "exclude".to_string()),
    ]);
    // Index identities are conservative compatibility constraints: a stale or
    // different index must never be trusted merely because data bytes match.
    // BED/VCF bytes are provenance only; normalized loci/ALT annotations above
    // define selection, independent of source ordering and redundant rows.
    for key in [
        "alignments",
        "reference",
        "cram-reference",
        "alignment-index",
        "reference-fai",
        "cram-reference-fai",
    ] {
        science.insert(
            format!("input.{key}"),
            identities
                .get(key)
                .cloned()
                .unwrap_or_else(|| "absent".into()),
        );
    }
    let science_digest = structured_digest("rosalind-evidence-science-v1", &science);
    execute(
        command,
        &mut engine,
        capture,
        receipt_path,
        started,
        initial_rss,
        hashing_ms,
        os_limit,
        science_digest,
        input_snapshot,
    )
}

// Kept separate from preflight so every failure below owns transactional output.
#[allow(clippy::too_many_arguments)]
fn execute(
    command: EvidenceCommand,
    engine: &mut EvidenceEngine,
    capture: CommandCapture,
    receipt_path: Option<PathBuf>,
    started: Instant,
    initial_rss: u64,
    hashing_ms: u128,
    os_limit: Option<u64>,
    science_digest: String,
    input_snapshot: InputSnapshot,
) -> Result<()> {
    // Implemented together with the first-party encoder interfaces below.
    execute_outputs(
        command,
        engine,
        capture,
        receipt_path,
        started,
        initial_rss,
        hashing_ms,
        os_limit,
        science_digest,
        input_snapshot,
    )
}
