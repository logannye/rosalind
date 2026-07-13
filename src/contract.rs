//! Public, analyzer-agnostic execution of Rosalind's bounded whole-genome
//! column stream. Downstream Rust binaries use this module to inherit planning,
//! refusal, the live RSS governor, deterministic output, and a sealed receipt
//! without copying CLI orchestration.

// The public error shape deliberately carries the complete breached outcome by
// value so callers can inspect its sealed receipt and telemetry without a second
// allocation or a breaking `Box` in the documented API.
#![allow(clippy::result_large_err)]

use std::collections::{BTreeMap, BTreeSet};
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::time::Duration;

use crate::call::{
    estimate_variants_working_set, predicted_peak_rss_bytes, run_bounded_selected_bam,
    run_bounded_whole_genome, ColumnAnalyzer,
};
use crate::core::governor::{GovernorError, MemoryGovernor};
use crate::core::{CoreError, MemoryBudget, WorkingSet, PILEUP_IO_RSS_OVERHEAD};
use crate::genomics::{AnalysisReference, ReferenceProvider};
use crate::io::bam::StreamingBamSource;
use crate::pileup::{PileupParams, SkipCounts};
use crate::provenance::{CommandCapture, RunManifest, MEASUREMENT_KEYS};
use crate::selection::AnalysisSelection;
use crate::util::atomic::{write_atomic, AtomicFile};
use crate::util::rss::peak_rss_bytes;

/// Identity of the binary that owns a contract run.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProducerIdentity {
    /// Product or fork name.
    pub name: String,
    /// Product version written to `tool_version` and `producer.version`.
    pub version: String,
    /// Source repository, when public.
    pub repository: Option<String>,
    /// Expected executable name for humans. Receipts never execute it implicitly.
    pub binary: String,
}

impl ProducerIdentity {
    /// Identity of the upstream Rosalind binary.
    pub fn rosalind() -> Self {
        Self {
            name: "rosalind".to_string(),
            version: env!("CARGO_PKG_VERSION").to_string(),
            repository: Some("https://github.com/logannye/rosalind".to_string()),
            binary: "rosalind".to_string(),
        }
    }
}

/// Identity and claim namespace of one column analyzer.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AnalyzerIdentity {
    /// Stable analyzer identifier.
    pub id: String,
    /// Analyzer version, independently evolvable from the host binary.
    pub version: String,
    /// Prefix applied to [`ColumnAnalyzer::params`]. Upstream `features` uses an
    /// empty prefix for receipt compatibility; third-party analyzers should use
    /// `analyzer.` or `x.<reverse-dns>.`.
    pub param_prefix: String,
}

impl AnalyzerIdentity {
    /// A namespaced analyzer identity suitable for downstream binaries.
    pub fn new(id: impl Into<String>, version: impl Into<String>) -> Self {
        Self {
            id: id.into(),
            version: version.into(),
            param_prefix: "analyzer.".to_string(),
        }
    }

    /// Override the claim prefix used for analyzer-provided parameters.
    pub fn with_param_prefix(mut self, prefix: impl Into<String>) -> Self {
        self.param_prefix = prefix.into();
        self
    }
}

/// Replayable, tokenized command prefix plus analyzer-specific options.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReplayInvocation {
    /// Tokens before standard contract operands, such as `["analyze", "coverage"]`.
    pub argv_prefix: Vec<String>,
    /// Analyzer-specific `--flag value` pairs.
    pub options: BTreeMap<String, String>,
    /// Analyzer-specific bare flags.
    pub flags: BTreeSet<String>,
}

impl ReplayInvocation {
    /// Start a replay recipe from an already-tokenized command prefix.
    pub fn new(prefix: impl IntoIterator<Item = impl Into<String>>) -> Self {
        Self {
            argv_prefix: prefix.into_iter().map(Into::into).collect(),
            options: BTreeMap::new(),
            flags: BTreeSet::new(),
        }
    }

    /// Add an analyzer-specific option.
    pub fn option(mut self, flag: impl Into<String>, value: impl ToString) -> Self {
        self.options.insert(flag.into(), value.to_string());
        self
    }

    /// Add an analyzer-specific bare flag.
    pub fn flag(mut self, flag: impl Into<String>) -> Self {
        self.flags.insert(flag.into());
        self
    }
}

/// Primary artifact destination.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OutputTarget {
    /// Write a hashable artifact to this path.
    File(PathBuf),
    /// Stream to stdout. A receipt can still be written explicitly, but it has
    /// no output artifact and therefore cannot be byte-reproduced.
    Stdout,
}

/// Maximum memory retained by an analyzer after the runner measures its baseline.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AnalyzerMemoryModel {
    /// No pre-run analyzer bound is available. Record-only execution remains valid.
    Unknown,
    /// A versioned upper bound on additional retained bytes after runner startup.
    Fixed {
        /// Stable name for the estimator used by the analyzer producer.
        model_id: String,
        /// Maximum additional retained bytes after the process baseline is measured.
        max_additional_bytes: u64,
    },
}

impl AnalyzerMemoryModel {
    fn additional_bytes(&self) -> Option<u64> {
        match self {
            Self::Unknown => None,
            Self::Fixed {
                max_additional_bytes,
                ..
            } => Some(*max_additional_bytes),
        }
    }
}

/// How strongly a run asks Rosalind to honor its memory declaration.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EnforcementMode {
    /// Measure and report only; an unknown analyzer model is allowed.
    RecordOnly,
    /// Predict up front and use Rosalind's cooperative RSS governor.
    Cooperative,
    /// Require cooperative enforcement plus an existing Linux cgroup-v2 hard limit.
    RequireOsLimit,
}

impl EnforcementMode {
    fn is_enforced(self) -> bool {
        self != Self::RecordOnly
    }
}

/// Evidence supporting the resource-contract result.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EnforcementAssurance {
    /// The run was measured without a declared analyzer bound.
    ObservedOnly,
    /// A declared bound and the cooperative process RSS governor were active.
    DeclaredBoundCooperative,
    /// Cooperative enforcement ran inside a matching Linux cgroup-v2 hard limit.
    CgroupV2,
}

impl EnforcementAssurance {
    fn as_str(self) -> &'static str {
        match self {
            Self::ObservedOnly => "observed-only",
            Self::DeclaredBoundCooperative => "declared-bound-cooperative",
            Self::CgroupV2 => "cgroup-v2",
        }
    }
}

/// Transaction policy for a file output and its receipt.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputPolicy {
    /// Refuse existing destinations and atomically create a new artifact.
    CreateNewAtomic,
    /// Atomically replace an existing destination (`--force` at the CLI).
    ReplaceAtomic,
}

impl OutputPolicy {
    fn replace(self) -> bool {
        self == Self::ReplaceAtomic
    }
}

impl OutputTarget {
    fn path(&self) -> Option<&Path> {
        match self {
            Self::File(path) => Some(path),
            Self::Stdout => None,
        }
    }
}

/// Complete configuration for one bounded column analysis.
#[derive(Debug, Clone)]
pub struct ContractRunSpec {
    /// Identity of the binary producing the output and receipt.
    pub producer: ProducerIdentity,
    /// Identity of the column analyzer implementation.
    pub analyzer: AnalyzerIdentity,
    /// Analyzer-owned memory retained after the process baseline is measured.
    pub analyzer_memory: AnalyzerMemoryModel,
    /// Tokenized replay command and analyzer-specific options.
    pub invocation: ReplayInvocation,
    /// Persisted Rosalind reference index.
    pub index: PathBuf,
    /// Coordinate-sorted BAM input.
    pub alignments: PathBuf,
    /// Primary artifact destination.
    pub output: OutputTarget,
    /// Whether file destinations must be new or may be replaced atomically.
    pub output_policy: OutputPolicy,
    /// Explicit receipt destination, or the output sidecar when omitted.
    pub manifest: Option<PathBuf>,
    /// Minimum accepted read mapping quality.
    pub mapq_threshold: u8,
    /// Maximum active reads retained at a locus; zero means uncapped.
    pub max_depth: u32,
    /// Maximum admitted read length when enforcement is enabled.
    pub max_read_len: u32,
    /// Optional declared process RSS budget in MiB.
    pub memory_budget_mb: Option<u64>,
    /// Requested resource-enforcement tier.
    pub enforcement: EnforcementMode,
}

/// Recorded resource-contract verdict.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ContractVerdict {
    /// No budget was declared.
    Unset,
    /// Realized RSS stayed within the declared budget.
    Within,
    /// Realized RSS exceeded the declared budget.
    Over,
}

impl ContractVerdict {
    fn as_str(self) -> &'static str {
        match self {
            Self::Unset => "unset",
            Self::Within => "within",
            Self::Over => "over",
        }
    }
}

/// Runtime-governor state recorded for the run.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GovernorState {
    /// No live enforcement; measurements are recorded after the run.
    RecordOnly,
    /// Live enforcement remained active and did not trip.
    Enforced,
    /// The live governor observed RSS above the budget.
    Tripped,
}

impl GovernorState {
    fn as_str(self) -> &'static str {
        match self {
            Self::RecordOnly => "record-only",
            Self::Enforced => "enforced",
            Self::Tripped => "tripped",
        }
    }
}

/// Why a contract run was refused before creating its output.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RefusalReport {
    /// Declared process RSS budget in MiB.
    pub budget_mb: u64,
    /// Predicted process RSS high-water mark.
    pub predicted_peak_rss_bytes: u64,
    /// Largest indexed contig used in the prediction.
    pub largest_contig_bytes: u64,
    /// Measured RSS before streaming begins.
    pub baseline_rss_bytes: u64,
    /// Configured active-read cap.
    pub max_depth: u32,
    /// Configured admitted read-length cap.
    pub max_read_len: u32,
}

/// Successful or breached run telemetry. A breached outcome still names the
/// partial output and sealed receipt when they were written.
#[derive(Debug, Clone)]
pub struct ContractRunOutcome {
    /// Self-hash of the sealed receipt, when a receipt was requested.
    pub claim_hash: Option<String>,
    /// Written receipt path, when present.
    pub manifest_path: Option<PathBuf>,
    /// Final successful primary output path, when file-backed.
    pub output_path: Option<PathBuf>,
    /// Preserved partial artifact after a governed breach.
    pub partial_output_path: Option<PathBuf>,
    /// Up-front process RSS prediction.
    pub predicted_peak_rss_bytes: u64,
    /// Realized or governor-observed process RSS high-water mark.
    pub peak_rss_bytes: u64,
    /// Modeled analyzer working-set high-water mark.
    pub max_working_set_bytes: u64,
    /// Analyzer contribution included in the up-front prediction, when known.
    pub analyzer_predicted_bytes: Option<u64>,
    /// Reads skipped for each bounded-ingest reason.
    pub skips: SkipCounts,
    /// Result of comparing realized RSS with the declared budget.
    pub verdict: ContractVerdict,
    /// Live-governor state for the run.
    pub governor: GovernorState,
    /// Strength of the evidence behind the resource-contract result.
    pub assurance: EnforcementAssurance,
}

/// Typed failure modes for library callers. The CLI alone maps refusal/breach
/// to process exit codes 3 and 4.
#[derive(Debug)]
pub enum ContractRunError {
    /// The requested contract cannot provide its stated guarantees.
    InvalidConfiguration(String),
    /// Enforced execution cannot proceed without an analyzer memory bound.
    UnknownAnalyzerBound,
    /// Safe output creation refused an existing destination.
    OutputExists(PathBuf),
    /// The requested OS-level enforcement tier is unavailable or insufficient.
    OsEnforcementUnavailable(String),
    /// Up-front prediction exceeded the enforced budget; no primary output exists.
    Refused(RefusalReport),
    /// Live RSS crossed the enforced budget after output creation.
    Breached(ContractRunOutcome),
    /// Index, alignment, or analyzer-stream failure.
    Input(CoreError),
    /// Filesystem or output-stream failure.
    Io(io::Error),
}

impl std::fmt::Display for ContractRunError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidConfiguration(message) => write!(f, "invalid contract run: {message}"),
            Self::UnknownAnalyzerBound => write!(
                f,
                "enforced contract run requires a declared analyzer memory model"
            ),
            Self::OutputExists(path) => write!(
                f,
                "output already exists: {} (choose another path or enable replacement)",
                path.display()
            ),
            Self::OsEnforcementUnavailable(message) => {
                write!(f, "OS memory enforcement unavailable: {message}")
            }
            Self::Refused(report) => write!(
                f,
                "contract refused: predicted peak {} MiB exceeds budget {} MiB",
                report.predicted_peak_rss_bytes / (1 << 20),
                report.budget_mb
            ),
            Self::Breached(outcome) => write!(
                f,
                "contract breached: realized peak {} MiB",
                outcome.peak_rss_bytes / (1 << 20)
            ),
            Self::Input(error) => write!(f, "analysis input failed: {error}"),
            Self::Io(error) => write!(f, "analysis I/O failed: {error}"),
        }
    }
}

impl std::error::Error for ContractRunError {}

impl From<CoreError> for ContractRunError {
    fn from(value: CoreError) -> Self {
        Self::Input(value)
    }
}

impl From<io::Error> for ContractRunError {
    fn from(value: io::Error) -> Self {
        Self::Io(value)
    }
}

/// Execute one analyzer over the bounded whole-genome pileup stream.
pub fn run_column_analysis(
    analyzer: &mut dyn ColumnAnalyzer,
    spec: ContractRunSpec,
) -> Result<ContractRunOutcome, ContractRunError> {
    run_column_analysis_selected(analyzer, spec, AnalysisSelection::WholeGenome)
}

/// Execute one analyzer over a whole-genome, interval, or deterministic shard selection.
pub fn run_column_analysis_selected(
    analyzer: &mut dyn ColumnAnalyzer,
    spec: ContractRunSpec,
    selection: AnalysisSelection,
) -> Result<ContractRunOutcome, ContractRunError> {
    validate_spec(&spec)?;
    validate_destinations(&spec)?;

    let analyzer_predicted_bytes = spec.analyzer_memory.additional_bytes();
    if spec.enforcement.is_enforced() && analyzer_predicted_bytes.is_none() {
        return Err(ContractRunError::UnknownAnalyzerBound);
    }
    let os_limit_bytes = if spec.enforcement == EnforcementMode::RequireOsLimit {
        let budget_bytes = MemoryBudget::from_mb(
            spec.memory_budget_mb
                .expect("validated OS-enforced run has a budget"),
        )
        .bytes;
        match effective_cgroup_v2_limit_bytes() {
            Some(limit) if limit <= budget_bytes => Some(limit),
            Some(limit) => {
                return Err(ContractRunError::OsEnforcementUnavailable(format!(
                    "active cgroup limit {} MiB exceeds declared budget {} MiB",
                    limit / (1 << 20),
                    budget_bytes / (1 << 20)
                )))
            }
            None => {
                return Err(ContractRunError::OsEnforcementUnavailable(
                    "run inside a cgroup-v2 container/systemd scope whose memory.max is at or below the declared budget"
                        .to_string(),
                ))
            }
        }
    } else {
        None
    };
    let assurance = match spec.enforcement {
        EnforcementMode::RecordOnly => EnforcementAssurance::ObservedOnly,
        EnforcementMode::Cooperative => EnforcementAssurance::DeclaredBoundCooperative,
        EnforcementMode::RequireOsLimit => EnforcementAssurance::CgroupV2,
    };

    let loaded = AnalysisReference::open(&spec.index).map_err(|error| {
        ContractRunError::InvalidConfiguration(format!(
            "failed to open analysis reference {}: {error}",
            spec.index.display()
        ))
    })?;
    let contigs = loaded.contigs();
    if !matches!(selection, AnalysisSelection::WholeGenome) {
        crate::io::bam::find_bai(&spec.alignments).ok_or_else(|| {
            ContractRunError::InvalidConfiguration(format!(
                "sparse analysis requires {}.bai or {}",
                spec.alignments.display(),
                spec.alignments.with_extension("bai").display()
            ))
        })?;
    }
    let largest = selection.largest_reference_span(contigs);
    let baseline = peak_rss_bytes();
    let predicted_ws = estimate_variants_working_set(largest, spec.max_depth, spec.max_read_len)
        .bytes
        .saturating_add(analyzer_predicted_bytes.unwrap_or(0));
    let predicted_peak =
        predicted_peak_rss_bytes(largest, spec.max_depth, spec.max_read_len, baseline)
            .saturating_add(analyzer_predicted_bytes.unwrap_or(0));

    if let Some(budget_mb) = spec
        .memory_budget_mb
        .filter(|_| spec.enforcement.is_enforced())
    {
        if !MemoryBudget::from_mb(budget_mb).admits(predicted_peak) {
            return Err(ContractRunError::Refused(RefusalReport {
                budget_mb,
                predicted_peak_rss_bytes: predicted_peak,
                largest_contig_bytes: largest,
                baseline_rss_bytes: baseline,
                max_depth: spec.max_depth,
                max_read_len: spec.max_read_len,
            }));
        }
    }

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
    let _governor = if spec.enforcement.is_enforced() {
        let budget = spec
            .memory_budget_mb
            .expect("validated enforced run has a budget");
        Some(
            MemoryGovernor::start(
                MemoryBudget::from_mb(budget).bytes,
                Duration::from_millis(poll_ms),
                live_rss,
            )
            .map_err(governor_error)?,
        )
    } else {
        None
    };

    let params = PileupParams {
        min_mapq: spec.mapq_threshold,
        max_depth: (spec.max_depth != 0).then_some(spec.max_depth),
        max_read_len: spec.enforcement.is_enforced().then_some(spec.max_read_len),
        ..PileupParams::default()
    };
    let mut atomic_output = match &spec.output {
        OutputTarget::File(path) => Some(AtomicFile::create(path)?),
        OutputTarget::Stdout => None,
    };
    let drive_result = match &mut atomic_output {
        Some(file) => {
            let mut writer = io::BufWriter::new(file.file_mut());
            let result = drive_selection(
                analyzer,
                &spec.alignments,
                &loaded,
                &selection,
                params,
                &mut writer,
            );
            if result.is_ok() {
                writer.flush()?;
            } else {
                let _ = writer.flush();
            }
            result
        }
        None => {
            let stdout = io::stdout();
            let mut writer = stdout.lock();
            let result = drive_selection(
                analyzer,
                &spec.alignments,
                &loaded,
                &selection,
                params,
                &mut writer,
            );
            if result.is_ok() {
                writer.flush()?;
            } else {
                let _ = writer.flush();
            }
            result
        }
    };
    let (max_ws, skips, breached, breach_peak) = match drive_result {
        Ok((working_set, skips)) => (working_set, skips, false, 0),
        Err(CoreError::BudgetExceeded { needed, .. }) => {
            (WorkingSet { bytes: 0 }, SkipCounts::default(), true, needed)
        }
        Err(error) => return Err(ContractRunError::Input(error)),
    };
    let peak = if breached {
        breach_peak
    } else {
        std::env::var("ROSALIND_FORCE_PEAK_RSS_BYTES")
            .ok()
            .and_then(|v| v.parse::<u64>().ok())
            .unwrap_or_else(peak_rss_bytes)
    };
    let verdict = match spec
        .memory_budget_mb
        .map(|mb| MemoryBudget::from_mb(mb).admits(peak))
    {
        None => ContractVerdict::Unset,
        Some(true) => ContractVerdict::Within,
        Some(false) => ContractVerdict::Over,
    };
    let governor = if breached {
        GovernorState::Tripped
    } else if spec.enforcement.is_enforced() {
        GovernorState::Enforced
    } else {
        GovernorState::RecordOnly
    };
    let is_breach =
        breached || (spec.enforcement.is_enforced() && verdict == ContractVerdict::Over);
    let (output_path, partial_output_path) = match (atomic_output, spec.output.path()) {
        (Some(file), Some(requested)) if is_breach => {
            let partial = partial_path(requested);
            let written = file.commit_as(&partial, spec.output_policy.replace())?;
            (None, Some(written))
        }
        (Some(file), Some(_)) => {
            let written = file.commit(spec.output_policy.replace())?;
            (Some(written), None)
        }
        _ => (None, None),
    };
    let manifest_path = receipt_destination(&spec);
    let mut receipt_spec = spec.clone();
    if let Some(path) = output_path.as_ref().or(partial_output_path.as_ref()) {
        receipt_spec.output = OutputTarget::File(path.clone());
    }
    let claim_hash = if let Some(path) = &manifest_path {
        Some(write_receipt(
            analyzer,
            &receipt_spec,
            path,
            predicted_ws,
            predicted_peak,
            baseline,
            peak,
            max_ws,
            skips,
            verdict,
            governor,
            assurance,
            os_limit_bytes,
            is_breach,
            &selection,
        )?)
    } else {
        None
    };
    let outcome = ContractRunOutcome {
        claim_hash,
        manifest_path,
        output_path,
        partial_output_path,
        predicted_peak_rss_bytes: predicted_peak,
        peak_rss_bytes: peak,
        max_working_set_bytes: max_ws.bytes,
        analyzer_predicted_bytes,
        skips,
        verdict,
        governor,
        assurance,
    };
    if is_breach {
        Err(ContractRunError::Breached(outcome))
    } else {
        Ok(outcome)
    }
}

/// Assertions intended for downstream analyzer integration tests.
#[cfg(feature = "contract-testkit")]
pub mod testkit {
    use std::path::Path;

    use crate::provenance::RunManifest;

    use super::{ContractRunError, ContractRunOutcome, ContractVerdict, GovernorState};

    /// Assert that two emitted artifacts are byte-identical.
    pub fn assert_byte_identical(a: &Path, b: &Path) {
        let left = std::fs::read(a).unwrap_or_else(|e| panic!("cannot read {}: {e}", a.display()));
        let right = std::fs::read(b).unwrap_or_else(|e| panic!("cannot read {}: {e}", b.display()));
        assert_eq!(
            left,
            right,
            "artifacts differ: {} vs {}",
            a.display(),
            b.display()
        );
    }

    /// Parse a receipt and assert both local integrity hashes are intact.
    pub fn assert_receipt_intact(path: &Path) -> RunManifest {
        let text = std::fs::read_to_string(path)
            .unwrap_or_else(|e| panic!("cannot read receipt {}: {e}", path.display()));
        let manifest = RunManifest::from_canonical_json(&text)
            .unwrap_or_else(|e| panic!("cannot parse receipt {}: {e}", path.display()));
        assert_eq!(
            manifest.self_hash_ok(),
            Some(true),
            "claim self-hash mismatch"
        );
        if manifest.claims_measurements() {
            assert_eq!(
                manifest.measurement_hash_ok(),
                Some(true),
                "measurement self-hash mismatch"
            );
        }
        manifest
    }

    /// Assert that relocating recorded file paths does not change a modern
    /// receipt's content-addressed claim.
    pub fn assert_paths_are_portable(manifest: &RunManifest) {
        let before = manifest.content_hash();
        let mut relocated = manifest.clone();
        for (index, file) in relocated.inputs.iter_mut().enumerate() {
            file.path = format!("/relocated/input-{index}");
        }
        for (index, file) in relocated.outputs.iter_mut().enumerate() {
            file.path = format!("/relocated/output-{index}");
        }
        assert_eq!(before, relocated.content_hash());
    }

    /// Extract a successful outcome or panic with the typed contract error.
    pub fn assert_completed(
        result: Result<ContractRunOutcome, ContractRunError>,
    ) -> ContractRunOutcome {
        result.unwrap_or_else(|e| panic!("contract run did not complete: {e}"))
    }

    /// Assert that a contract result was refused before execution.
    pub fn assert_refused(result: Result<ContractRunOutcome, ContractRunError>) {
        assert!(
            matches!(result, Err(ContractRunError::Refused(_))),
            "expected Refused, got {result:?}"
        );
    }

    /// Assert that an enforced run completed within its declared budget.
    pub fn assert_enforced_fit(
        result: Result<ContractRunOutcome, ContractRunError>,
    ) -> ContractRunOutcome {
        let outcome = assert_completed(result);
        assert_eq!(outcome.verdict, ContractVerdict::Within);
        assert_eq!(outcome.governor, GovernorState::Enforced);
        outcome
    }

    /// Assert that `reproduce --binary` performed a byte comparison with the
    /// explicitly selected executable and reported `REPRODUCED`.
    pub fn assert_reproduced_with_binary(
        report: &crate::reproduce::ReproReport,
        expected_binary: &Path,
    ) {
        assert_eq!(report.verdict_label, "REPRODUCED");
        assert_eq!(report.exit_code, 0);
        assert!(report.compared, "reproduction did not compare output bytes");
        assert_eq!(report.execution_binary, expected_binary);
    }
}

fn validate_spec(spec: &ContractRunSpec) -> Result<(), ContractRunError> {
    if spec.invocation.argv_prefix.is_empty() {
        return Err(ContractRunError::InvalidConfiguration(
            "replay argv prefix must not be empty".to_string(),
        ));
    }
    if spec.enforcement.is_enforced() && spec.memory_budget_mb.is_none() {
        return Err(ContractRunError::InvalidConfiguration(
            "--enforce requires a memory budget".to_string(),
        ));
    }
    if spec.enforcement.is_enforced() && spec.max_depth == 0 {
        return Err(ContractRunError::InvalidConfiguration(
            "--enforce requires max_depth > 0".to_string(),
        ));
    }
    const STANDARD_FLAGS: &[&str] = &[
        "--index",
        "--alignments",
        "--mapq-threshold",
        "--max-depth",
        "--max-read-len",
        "--memory-budget-mb",
        "--enforce",
        "--require-os-limit",
        "--force",
        "-o",
        "--output",
        "--manifest",
    ];
    for flag in spec
        .invocation
        .options
        .keys()
        .chain(spec.invocation.flags.iter())
    {
        if !flag.starts_with('-') {
            return Err(ContractRunError::InvalidConfiguration(format!(
                "replay option {flag:?} must begin with '-'"
            )));
        }
        if STANDARD_FLAGS.contains(&flag.as_str()) {
            return Err(ContractRunError::InvalidConfiguration(format!(
                "analyzer replay option {flag} collides with the contract runner"
            )));
        }
    }
    let is_bam = spec
        .alignments
        .extension()
        .and_then(|e| e.to_str())
        .is_some_and(|e| e.eq_ignore_ascii_case("bam"));
    if !is_bam {
        return Err(ContractRunError::InvalidConfiguration(
            "bounded index analysis requires a coordinate-sorted BAM".to_string(),
        ));
    }
    Ok(())
}

fn validate_destinations(spec: &ContractRunSpec) -> Result<(), ContractRunError> {
    if spec.output_policy == OutputPolicy::ReplaceAtomic {
        return Ok(());
    }
    if let Some(output) = spec.output.path() {
        for path in [output.to_path_buf(), partial_path(output)] {
            if path.exists() {
                return Err(ContractRunError::OutputExists(path));
            }
        }
    }
    if let Some(receipt) = receipt_destination(spec) {
        if receipt.exists() {
            return Err(ContractRunError::OutputExists(receipt));
        }
    }
    Ok(())
}

fn partial_path(path: &Path) -> PathBuf {
    let mut name = path.as_os_str().to_os_string();
    name.push(".partial");
    PathBuf::from(name)
}

fn effective_cgroup_v2_limit_bytes() -> Option<u64> {
    if let Ok(value) = std::env::var("ROSALIND_TEST_CGROUP_MEMORY_MAX") {
        return (value != "max").then(|| value.parse().ok()).flatten();
    }
    #[cfg(target_os = "linux")]
    {
        let cgroup = std::fs::read_to_string("/proc/self/cgroup").ok()?;
        let relative = cgroup.lines().find_map(|line| {
            let mut fields = line.splitn(3, ':');
            let hierarchy = fields.next()?;
            let controllers = fields.next()?;
            let path = fields.next()?;
            (hierarchy == "0" && controllers.is_empty()).then_some(path)
        })?;
        let path = Path::new("/sys/fs/cgroup")
            .join(relative.trim_start_matches('/'))
            .join("memory.max");
        let value = std::fs::read_to_string(path).ok()?;
        let value = value.trim();
        (value != "max").then(|| value.parse().ok()).flatten()
    }
    #[cfg(not(target_os = "linux"))]
    {
        None
    }
}

/// Detect the effective Linux cgroup-v2 `memory.max`, when finite and readable.
pub fn detected_os_memory_limit_bytes() -> Option<u64> {
    effective_cgroup_v2_limit_bytes()
}

fn governor_error(error: GovernorError) -> ContractRunError {
    ContractRunError::InvalidConfiguration(format!("failed to start memory governor: {error}"))
}

fn drive<S: crate::pileup::ReadSource>(
    analyzer: &mut dyn ColumnAnalyzer,
    source: S,
    ref_view: &dyn crate::genomics::ReferenceSequence,
    contigs: &crate::core::ContigSet,
    params: PileupParams,
    writer: &mut dyn Write,
) -> Result<(WorkingSet, SkipCounts), CoreError> {
    run_bounded_whole_genome(analyzer, source, ref_view, contigs, params, writer)
}

fn drive_selection(
    analyzer: &mut dyn ColumnAnalyzer,
    bam_path: &Path,
    reference: &AnalysisReference,
    selection: &AnalysisSelection,
    params: PileupParams,
    writer: &mut dyn Write,
) -> Result<(WorkingSet, SkipCounts), CoreError> {
    match selection {
        AnalysisSelection::WholeGenome => {
            let source = StreamingBamSource::new(bam_path, reference.contigs())?;
            drive(
                analyzer,
                source,
                reference,
                reference.contigs(),
                params,
                writer,
            )
        }
        AnalysisSelection::Intervals(_) | AnalysisSelection::Shard { .. } => {
            run_bounded_selected_bam(analyzer, bam_path, reference, selection, params, writer)
        }
    }
}

fn receipt_destination(spec: &ContractRunSpec) -> Option<PathBuf> {
    spec.manifest.clone().or_else(|| {
        spec.output.path().map(|path| {
            let mut name = path.as_os_str().to_os_string();
            name.push(".manifest.json");
            PathBuf::from(name)
        })
    })
}

#[allow(clippy::too_many_arguments)]
fn write_receipt(
    analyzer: &dyn ColumnAnalyzer,
    spec: &ContractRunSpec,
    destination: &Path,
    predicted_ws: u64,
    predicted_peak: u64,
    baseline: u64,
    peak: u64,
    max_ws: WorkingSet,
    skips: SkipCounts,
    verdict: ContractVerdict,
    governor: GovernorState,
    assurance: EnforcementAssurance,
    os_limit_bytes: Option<u64>,
    breached: bool,
    selection: &AnalysisSelection,
) -> Result<String, ContractRunError> {
    let mut manifest = RunManifest::new(spec.invocation.argv_prefix.join(" "));
    manifest.tool_version = spec.producer.version.clone();
    let mut command = CommandCapture::from_argv_prefix(spec.invocation.argv_prefix.clone());
    let reference_flag = if AnalysisReference::path_is_pack(&spec.index).unwrap_or(false) {
        "--reference-pack"
    } else {
        "--index"
    };
    command.input(reference_flag, &spec.index)?;
    command.input("--alignments", &spec.alignments)?;
    match selection {
        AnalysisSelection::WholeGenome => {}
        AnalysisSelection::Intervals(intervals) => {
            if let Some(region) = intervals.region_origin() {
                command.opt("--region", region);
            } else if let Some(path) = intervals.bed_origin() {
                command.input("--regions", path)?;
            }
        }
        AnalysisSelection::Shard { count, index, .. } => {
            command.opt("--shard-count", *count);
            command.opt("--shard-index", *index);
        }
    }
    command.opt("--mapq-threshold", spec.mapq_threshold);
    command.opt("--max-depth", spec.max_depth);
    command.opt("--max-read-len", spec.max_read_len);
    command.flag_if(spec.enforcement.is_enforced(), "--enforce");
    command.flag_if(
        spec.enforcement == EnforcementMode::RequireOsLimit,
        "--require-os-limit",
    );
    command.flag_if(spec.output_policy == OutputPolicy::ReplaceAtomic, "--force");
    let analyzer_params = analyzer.params();
    if spec.invocation.argv_prefix.first().map(String::as_str) == Some("features") {
        if let Some(format) = analyzer_params.get("artifact.format") {
            command.opt(
                "--format",
                if format == "arrow-ipc" {
                    "arrow-ipc"
                } else {
                    "tsv"
                },
            );
        }
    }
    if let Some(memory_mb) = spec.memory_budget_mb {
        command.opt("--memory-budget-mb", memory_mb);
    }
    for (flag, value) in &spec.invocation.options {
        command.opt(flag, value);
    }
    for flag in &spec.invocation.flags {
        command.flag(flag);
    }
    if let Some(output) = spec.output.path() {
        command.output("-o", output)?;
    }
    command.record_into(&mut manifest);
    manifest.params.insert(
        "artifact.input.0.role".to_string(),
        if AnalysisReference::path_is_pack(&spec.index).unwrap_or(false) {
            "analysis-reference-pack"
        } else {
            "reference-index"
        }
        .to_string(),
    );
    manifest.params.insert(
        "artifact.input.1.role".to_string(),
        "sorted-alignments".to_string(),
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
    manifest
        .params
        .insert("partition.kind".to_string(), selection.kind().to_string());
    if let Some(intervals) = selection.intervals() {
        manifest.params.insert(
            "partition.interval_count".to_string(),
            intervals.intervals().len().to_string(),
        );
        manifest.params.insert(
            "partition.total_bases".to_string(),
            intervals.total_bases().to_string(),
        );
        manifest
            .params
            .insert("partition.intervals_blake3".to_string(), intervals.blake3());
    }
    if let AnalysisSelection::Shard { count, index, .. } = selection {
        manifest.params.insert(
            "partition.algorithm".to_string(),
            "reference-span-v1".to_string(),
        );
        manifest
            .params
            .insert("partition.shard_count".to_string(), count.to_string());
        manifest
            .params
            .insert("partition.shard_index".to_string(), index.to_string());
    }
    if !manifest.outputs.is_empty() {
        manifest.params.insert(
            "artifact.output.0.role".to_string(),
            if breached {
                "partial-analyzer-output"
            } else {
                "analyzer-output"
            }
            .to_string(),
        );
        if let Some(format) = analyzer_params.get("artifact.format") {
            manifest
                .params
                .insert("artifact.output.0.format".to_string(), format.clone());
        }
    }
    manifest.params.insert(
        "replay.kind".to_string(),
        if spec.producer.name == "rosalind" {
            "rosalind"
        } else {
            "external-analyzer"
        }
        .to_string(),
    );

    manifest
        .params
        .insert("producer.name".to_string(), spec.producer.name.clone());
    manifest.params.insert(
        "producer.version".to_string(),
        spec.producer.version.clone(),
    );
    manifest
        .params
        .insert("producer.binary".to_string(), spec.producer.binary.clone());
    if let Some(repository) = &spec.producer.repository {
        manifest
            .params
            .insert("producer.repository".to_string(), repository.clone());
    }
    manifest
        .params
        .insert("analyzer.id".to_string(), spec.analyzer.id.clone());
    manifest.params.insert(
        "analyzer.version".to_string(),
        spec.analyzer.version.clone(),
    );
    for (key, value) in analyzer_params {
        let key = format!("{}{}", spec.analyzer.param_prefix, key);
        if MEASUREMENT_KEYS.contains(&key.as_str()) {
            return Err(ContractRunError::InvalidConfiguration(format!(
                "analyzer parameter {key} collides with a measurement key"
            )));
        }
        if manifest.params.contains_key(&key) {
            return Err(ContractRunError::InvalidConfiguration(format!(
                "analyzer parameter {key} collides with a reserved claim field"
            )));
        }
        manifest.params.insert(key, value);
    }
    manifest.params.insert(
        "predicted_working_set_bytes".to_string(),
        predicted_ws.to_string(),
    );
    match &spec.analyzer_memory {
        AnalyzerMemoryModel::Unknown => {
            manifest
                .params
                .insert("analyzer.memory_model".to_string(), "unknown".to_string());
        }
        AnalyzerMemoryModel::Fixed {
            model_id,
            max_additional_bytes,
        } => {
            manifest
                .params
                .insert("analyzer.memory_model".to_string(), model_id.clone());
            manifest.params.insert(
                "analyzer.max_additional_bytes".to_string(),
                max_additional_bytes.to_string(),
            );
        }
    }
    manifest.params.insert(
        "contract.assurance".to_string(),
        assurance.as_str().to_string(),
    );
    manifest.params.insert(
        "run_status".to_string(),
        if breached { "breached" } else { "completed" }.to_string(),
    );
    if let Some(limit) = os_limit_bytes {
        manifest
            .params
            .insert("os.memory_limit_bytes".to_string(), limit.to_string());
    }
    manifest
        .params
        .insert("peak_rss_bytes".to_string(), peak.to_string());
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
        .insert("governor".to_string(), governor.as_str().to_string());
    manifest
        .params
        .insert("baseline_rss_bytes".to_string(), baseline.to_string());
    manifest.params.insert(
        "rss_residual_bytes".to_string(),
        peak.saturating_sub(max_ws.bytes)
            .saturating_sub(baseline)
            .to_string(),
    );
    manifest.params.insert(
        "io_rss_overhead_assumed_bytes".to_string(),
        PILEUP_IO_RSS_OVERHEAD.to_string(),
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
        .insert("contract_verdict".to_string(), verdict.as_str().to_string());
    manifest.finalize();
    write_atomic(
        destination,
        manifest.to_canonical_json().as_bytes(),
        spec.output_policy.replace(),
    )?;
    Ok(manifest.content_hash())
}
