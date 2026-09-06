//! Managed, deterministic artifacts produced by independent evidence analyzers.
//!
//! The runner owns input verification, cooperative resource/cancellation checks,
//! canonical batch boundaries, staged files and receipt publication. Factories
//! own scientific reducers and encoders, and borrow the staged output sink.

use super::*;
use crate::contract::{
    AnalyzerIdentity, EnforcementMode, OutputPolicy, ProducerIdentity, ReplayInvocation,
};
use crate::core::cancellation::{CancellationScope, CancellationToken, SignalCancellationGuard};
use crate::core::governor::{checkpoint, MemoryGovernor};
use crate::core::CoreError;
use crate::dataset::{
    canonical_dataset_query, dataset_query_digest, DatasetError, DatasetQuery, DatasetReadLimits,
    DescriptorSelection, VerifiedEvidenceDataset, VerifiedInputSession,
};
use crate::provenance::{CommandCapture, FileHash, RunManifest};
use crate::util::atomic::{commit_group, AtomicFile};
use crate::variant_io::VariantLimits;
use serde::Deserialize;
use std::collections::{BTreeMap, BTreeSet};
use std::fs::{self, File};
use std::io::{self, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};

/// Maximum replayable inline canonical query; larger selections use BED/VCF.
pub const MAX_ARTIFACT_QUERY_BYTES: usize = 32 << 10;
/// Default serialized receipt envelope. Its actual reservation is size-derived.
pub const DEFAULT_ARTIFACT_RECEIPT_BYTES: usize = 32 << 20;

/// Scientific selection resolver, independent of execution tiling.
#[derive(Debug, Clone)]
pub enum ArtifactSelection {
    /// Use the supplied EvidenceRequest or DatasetQuery selection. Non-whole
    /// requests must fit the 32 KiB inline replay envelope.
    Request,
    /// Use exactly the stored dataset selection; valid only for dataset inputs.
    Stored,
    /// Read and verify a local BED, replacing the supplied selection.
    Bed(PathBuf),
    /// Read and verify local VCF/VCF.gz/BCF with explicit parsing envelopes.
    Variants(PathBuf, VariantLimits),
}

/// Raw indexed inputs or portable persisted evidence for one external analyzer.
#[derive(Debug, Clone)]
pub enum EvidenceArtifactSource {
    /// Hash verified inputs once before opening the indexed engine.
    Native {
        /// Scientific profile, fields, sample and resource request.
        request: EvidenceRequest,
        /// Typed selection source; file forms explicitly replace request.selection.
        selection: ArtifactSelection,
    },
    /// Read verified persisted rows without opening original alignment/reference files.
    Dataset {
        /// Portable evidence-dataset.manifest.json location.
        manifest: PathBuf,
        /// Requested physical projection and, for Request, selection.
        query: DatasetQuery,
        /// File forms replace query.selection; Stored uses the stored denominator.
        selection: ArtifactSelection,
        /// Consumer and process resource declarations.
        execution: EvidenceExecution,
        /// Untrusted portable metadata envelopes.
        limits: DatasetReadLimits,
        /// Optional replay proof operands; each must match an actual dependency.
        artifacts: Vec<PathBuf>,
    },
}

/// A factory supplies scientific state and a sink-borrowing batch consumer.
/// The complete peak state of reducer, encoder and flush transients belongs in
/// requirements.retained_bytes. Returned requirements must match exactly.
pub trait EvidenceArtifactFactory {
    /// Required fields, reference capability and peak additional retained bytes.
    fn requirements(&self) -> EvidenceRequirements;
    /// Deterministic scientific parameters, under the spec's analyzer namespace.
    fn params(&self) -> BTreeMap<String, String>;
    /// Construct after admission. The sink cannot escape the returned analyzer.
    fn create<'a>(
        &'a mut self,
        output: &'a mut dyn Write,
    ) -> Result<Box<dyn EvidenceAnalyzer + 'a>, EvidenceError>;
}

/// Complete managed artifact request. Existing ColumnAnalyzer APIs are unchanged.
#[derive(Debug, Clone)]
pub struct EvidenceArtifactSpec {
    /// Executable producer identity, independent of the Rosalind dependency version.
    pub producer: ProducerIdentity,
    /// Independent scientific analyzer identity and parameter namespace.
    pub analyzer: AnalyzerIdentity,
    /// Tokenized command prefix and analyzer-owned options; standard flags are reserved.
    pub invocation: ReplayInvocation,
    /// Verified extraction or persisted source.
    pub source: EvidenceArtifactSource,
    /// Primary persisted artifact path.
    pub output: PathBuf,
    /// Receipt path; defaults to OUTPUT.manifest.json.
    pub manifest: Option<PathBuf>,
    /// Transactional creation or explicit replacement.
    pub output_policy: OutputPolicy,
    /// Observation, cooperative enforcement, or required cgroup enforcement.
    pub enforcement: EnforcementMode,
    /// Cooperative cancellation request, also checked during input hashing.
    pub cancellation: Option<CancellationToken>,
    /// Opt in to temporary SIGINT/SIGTERM handling during this run.
    pub handle_signals: bool,
    /// Maximum emitted receipt bytes, checked before publication.
    pub max_receipt_bytes: usize,
}
impl EvidenceArtifactSpec {
    /// Construct an observed, create-new artifact request with default envelopes.
    pub fn new(
        source: EvidenceArtifactSource,
        output: impl Into<PathBuf>,
        producer: ProducerIdentity,
        analyzer: AnalyzerIdentity,
        invocation: ReplayInvocation,
    ) -> Self {
        Self {
            producer,
            analyzer,
            invocation,
            source,
            output: output.into(),
            manifest: None,
            output_policy: OutputPolicy::CreateNewAtomic,
            enforcement: EnforcementMode::RecordOnly,
            cancellation: None,
            handle_signals: false,
            max_receipt_bytes: DEFAULT_ARTIFACT_RECEIPT_BYTES,
        }
    }
}

/// A committed result or explicitly incomplete resource-failure artifact.
#[derive(Debug, Clone)]
pub struct EvidenceArtifactOutcome {
    /// Committed primary artifact, or OUTPUT.partial after resource failure.
    pub output: PathBuf,
    /// Committed receipt; failed resources use RECEIPT.partial.
    pub manifest: PathBuf,
    /// Sealed scientific claim hash.
    pub claim_hash: String,
    /// True only after complete execution, flush, source guards and publication.
    pub completed: bool,
    /// Exact emitted loci and observed extraction counters.
    pub stats: EvidenceRunStats,
    /// Preflight process reservation including managed artifact state.
    pub predicted_peak_rss_bytes: u64,
    /// Observed process high-water memory through receipt construction.
    pub peak_rss_bytes: u64,
}

/// Typed artifact failure. Cancellation and arbitrary failures never publish a
/// success or a partial artifact; resource failures may publish identified partials.
#[derive(Debug, thiserror::Error)]
pub enum EvidenceArtifactError {
    /// Inconsistent settings or replay recipe.
    #[error("invalid evidence artifact configuration: {0}")]
    InvalidConfiguration(String),
    /// Required input could not be read or interpreted.
    #[error("invalid evidence artifact input: {0}")]
    Input(String),
    /// A recorded persisted identity failed verification.
    #[error("evidence artifact integrity failed: {0}")]
    Integrity(String),
    /// Unknown analyzer state is not admissible under enforcement.
    #[error("enforced evidence artifact requires a declared analyzer memory bound")]
    UnknownAnalyzerBound,
    /// A destination exists and replacement was not requested.
    #[error("artifact destination already exists: {0}")]
    OutputExists(PathBuf),
    /// No adequate process-level OS boundary exists.
    #[error("OS memory enforcement unavailable: {0}")]
    OsEnforcementUnavailable(String),
    /// Prediction refused before creating output.
    #[error("evidence artifact refused: needs {needed} bytes, budget is {budget} bytes")]
    Refused {
        /// Minimum predicted process bytes.
        needed: u64,
        /// Declared process budget.
        budget: u64,
    },
    /// Actual resource/record limits stopped execution.
    #[error("evidence artifact resource failure: {message}")]
    Resource {
        /// Exact underlying resource failure.
        message: String,
        /// Explicit incomplete files when bytes had been staged.
        partial: Option<Box<EvidenceArtifactOutcome>>,
    },
    /// Caller or signal requested cooperative cancellation.
    #[error("evidence artifact cancelled")]
    Cancelled,
    /// Scientific consumer, encoder or factory failed.
    #[error("evidence artifact analyzer failed: {0}")]
    Analyzer(String),
    /// Local file or output-stream failure.
    #[error("evidence artifact I/O failed: {0}")]
    Io(#[from] io::Error),
}
impl EvidenceArtifactError {
    /// CLI exit convention; library callers should match the typed variant.
    pub fn exit_code(&self) -> i32 {
        match self {
            Self::Refused { .. }
            | Self::UnknownAnalyzerBound
            | Self::OsEnforcementUnavailable(_) => 3,
            Self::Resource { .. } => 4,
            Self::Integrity(_) => 5,
            Self::Cancelled => 130,
            _ => 2,
        }
    }
}
impl From<EvidenceError> for EvidenceArtifactError {
    fn from(error: EvidenceError) -> Self {
        match error {
            EvidenceError::Refused { needed, budget } => Self::Refused { needed, budget },
            EvidenceError::Core(CoreError::Cancelled) => Self::Cancelled,
            EvidenceError::Core(CoreError::BudgetExceeded { .. })
            | EvidenceError::RecordLimit(_) => Self::Resource {
                message: error.to_string(),
                partial: None,
            },
            EvidenceError::Io(error) => Self::Io(error),
            EvidenceError::Analyzer(error) => Self::Analyzer(error),
            EvidenceError::InvalidRequest(error) => Self::InvalidConfiguration(error),
            other => Self::Input(other.to_string()),
        }
    }
}
impl From<DatasetError> for EvidenceArtifactError {
    fn from(error: DatasetError) -> Self {
        match error {
            DatasetError::Evidence(error) => error.into(),
            DatasetError::Io(error) => Self::Io(error),
            DatasetError::Corrupt(error) => Self::Integrity(error),
            other => Self::Input(other.to_string()),
        }
    }
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct InlineQuery {
    version: u32,
    fields: u32,
    selection: DescriptorSelection,
}
/// Parse the bounded, versioned --query-json recipe operand. The selected input's
/// dictionary performs coordinate/reference validation when the runner opens it.
pub fn parse_artifact_query(text: &str) -> Result<DatasetQuery, EvidenceArtifactError> {
    if text.len() > MAX_ARTIFACT_QUERY_BYTES {
        return Err(invalid(
            "inline query exceeds32KiB; use a file-backed BED/VCF selection",
        ));
    }
    let query: InlineQuery =
        serde_json::from_str(text).map_err(|e| invalid(&format!("invalid inline query: {e}")))?;
    if query.version != 1 {
        return Err(invalid("unsupported inline query version"));
    }
    // The descriptor form carries both denominators and optional SNV labels.
    // Validate their agreement before choosing the engine representation.
    let intervals = &query.selection.intervals;
    if intervals.iter().any(|i| i.start >= i.end)
        || intervals.windows(2).any(|pair| {
            pair[0].contig > pair[1].contig
                || (pair[0].contig == pair[1].contig && pair[0].end >= pair[1].start)
        })
    {
        return Err(invalid("inline query intervals must be normalized"));
    }
    if let Some(sites) = &query.selection.sites {
        let rows: u64 = intervals.iter().map(|i| u64::from(i.end - i.start)).sum();
        if rows != sites.len() as u64
            || sites.windows(2).any(|pair| {
                (pair[0].contig, pair[0].position) >= (pair[1].contig, pair[1].position)
            })
            || sites.iter().any(|s| {
                !b"ACGT".contains(&s.reference)
                    || s.alternates.is_empty()
                    || s.alternates.len() > 3
                    || s.alternates.windows(2).any(|pair| pair[0] >= pair[1])
                    || s.alternates
                        .iter()
                        .any(|b| !b"ACGT".contains(b) || *b == s.reference)
            })
            || intervals
                .iter()
                .flat_map(|i| (i.start..i.end).map(move |p| (i.contig, p)))
                .zip(sites)
                .any(|(locus, site)| locus != (site.contig, site.position))
        {
            return Err(invalid(
                "inline query sites must exactly match canonical intervals and SNV alleles",
            ));
        }
    }
    Ok(DatasetQuery {
        fields: EvidenceFields::from_bits(query.fields)?,
        selection: query.selection.to_evidence_selection(),
    })
}

/// Run a factory under the managed evidence artifact lifecycle. The governor and
/// cancellation scope cover hashing, opening, extraction, encoding and receipt
/// preparation. Successful artifact/receipt files publish as one checked group.
pub fn run_evidence_artifact(
    factory: &mut dyn EvidenceArtifactFactory,
    mut spec: EvidenceArtifactSpec,
) -> Result<EvidenceArtifactOutcome, EvidenceArtifactError> {
    let started = Instant::now();
    let token = spec.cancellation.clone().unwrap_or_default();
    let _cancellation =
        CancellationScope::start(token.clone()).map_err(|e| invalid(&e.to_string()))?;
    let _signals = if spec.handle_signals {
        Some(SignalCancellationGuard::install(&_cancellation)?)
    } else {
        None
    };
    check_runtime(None)?;
    let requirements = factory.requirements();
    let declared_budget = execution(&spec.source).memory_budget_bytes;
    let enforced = spec.enforcement != EnforcementMode::RecordOnly;
    if enforced && requirements.retained_bytes.is_none() {
        return Err(EvidenceArtifactError::UnknownAnalyzerBound);
    }
    if enforced && declared_budget.is_none() {
        return Err(invalid("enforcement requires a memory budget"));
    }
    let budget = declared_budget.filter(|_| enforced);
    if let Some(budget) = budget {
        let needed = rss();
        if needed > budget {
            return Err(EvidenceArtifactError::Refused { needed, budget });
        }
    }
    let os_limit = if spec.enforcement == EnforcementMode::RequireOsLimit {
        let limit = crate::contract::detected_os_memory_limit_bytes().ok_or_else(|| {
            EvidenceArtifactError::OsEnforcementUnavailable(
                "require an existing Linux cgroup-v2 memory.max".into(),
            )
        })?;
        if limit > budget.unwrap() {
            return Err(EvidenceArtifactError::OsEnforcementUnavailable(
                "cgroup memory.max exceeds declared budget".into(),
            ));
        }
        Some(limit)
    } else {
        None
    };
    let _governor = budget
        .map(|bytes| MemoryGovernor::start(bytes, Duration::from_millis(100), rss))
        .transpose()
        .map_err(|e| invalid(&e.to_string()))?;
    let params = factory.params();
    check_runtime(budget)?;
    validate_spec(&spec, &params)?;
    execution_mut(&mut spec.source).memory_budget_bytes = budget;
    let receipt_path = spec
        .manifest
        .clone()
        .unwrap_or_else(|| suffix(&spec.output, ".manifest.json"));
    let mut capture = CommandCapture::from_argv_prefix(spec.invocation.argv_prefix.clone());
    capture_common(&mut capture, &spec, declared_budget);
    let mut opened = open_source(&spec.source, &mut capture)?;
    opened.verify()?;
    let inputs = opened.input_paths();
    validate_paths(&spec, &receipt_path, &inputs, opened.dataset_root())?;
    let fields = opened.fields();
    if !fields.contains(requirements.fields) {
        return Err(invalid("source projection omits factory-required fields"));
    }
    let batch_bytes = (fields.storage_bytes_per_locus() + 40)
        .checked_mul(EVIDENCE_ARROW_BATCH_ROWS as u64)
        .and_then(|n| n.checked_add(64 << 10))
        .ok_or_else(|| invalid("canonical buffer reservation overflow"))?;
    let retained = requirements
        .retained_bytes
        .unwrap_or(0)
        .max(execution(&spec.source).analyzer_bytes);
    let base_metadata = opened.receipt_memory_bytes().saturating_add(2 << 20);
    // Query normalization, compatibility serialization and receipt construction
    // can duplicate input metadata. Admit their explicit bound before those
    // allocations, even for a huge file-backed selection with a tiny output.
    let preliminary = retained
        .checked_add(base_metadata)
        .and_then(|n| n.checked_add(batch_bytes))
        .ok_or_else(|| invalid("artifact metadata reservation overflow"))?;
    opened.plan(&Declared {
        requirements: EvidenceRequirements {
            retained_bytes: Some(preliminary),
            ..requirements.clone()
        },
    })?;
    check_runtime(budget)?;
    let source_params = opened.params()?;
    let serialized_params =
        serde_json::to_vec(&(&params, &source_params)).map_err(|e| invalid(&e.to_string()))?;
    if serialized_params.len() > spec.max_receipt_bytes {
        return Err(invalid("scientific parameters exceed receipt envelope"));
    }
    let metadata = base_metadata
        .checked_add((serialized_params.len() as u64).saturating_mul(16))
        .ok_or_else(|| invalid("artifact metadata reservation overflow"))?;
    let reserve = retained
        .checked_add(metadata)
        .and_then(|n| n.checked_add(batch_bytes))
        .ok_or_else(|| invalid("factory/runner reservation overflow"))?;
    let planner = Declared {
        requirements: EvidenceRequirements {
            retained_bytes: Some(reserve),
            ..requirements.clone()
        },
    };
    let predicted = opened.plan(&planner)?;
    check_runtime(budget)?;
    let mut pending = AtomicFile::create(&spec.output)?;
    let mut emitted = 0;
    let mut stats = EvidenceRunStats::default();
    let result = (|| -> Result<(), EvidenceArtifactError> {
        let mut writer = BufWriter::new(pending.file_mut());
        let result = (|| -> Result<EvidenceRunStats, EvidenceArtifactError> {
            let mut analyzer = factory.create(&mut writer)?;
            if analyzer.requirements() != requirements {
                return Err(EvidenceArtifactError::Analyzer(
                    "factory-created analyzer requirements differ from the admitted declaration"
                        .into(),
                ));
            }
            check_runtime(budget)?;
            let mut canonical = CanonicalConsumer::new(analyzer.as_mut(), fields, reserve, budget);
            let result = opened.run(&mut canonical);
            emitted = canonical.emitted;
            result
        })();
        let flushed = writer.flush();
        if matches!(result, Err(EvidenceArtifactError::Resource { .. })) {
            flushed?;
            return result.map(|_| ());
        }
        result.map(|observed| stats = observed)?;
        flushed?;
        check_runtime(budget)?;
        opened.verify()?;
        Ok(())
    })();
    stats.emitted_loci = emitted;
    match result {
        Ok(()) => publish(
            pending,
            &receipt_path,
            &spec,
            &requirements,
            &token,
            &mut opened,
            capture,
            params,
            source_params,
            stats,
            predicted,
            started,
            os_limit,
            None,
            budget,
        ),
        Err(EvidenceArtifactError::Resource { message, .. }) => {
            if token.is_cancelled() {
                return Err(EvidenceArtifactError::Cancelled);
            }
            let partial = publish(
                pending,
                &suffix(&receipt_path, ".partial"),
                &spec,
                &requirements,
                &token,
                &mut opened,
                capture,
                params,
                source_params,
                stats,
                predicted,
                started,
                os_limit,
                Some(&message),
                None,
            )?;
            Err(EvidenceArtifactError::Resource {
                message,
                partial: Some(Box::new(partial)),
            })
        }
        Err(error) => Err(error),
    }
}

struct Declared {
    requirements: EvidenceRequirements,
}
impl EvidenceAnalyzer for Declared {
    fn requirements(&self) -> EvidenceRequirements {
        self.requirements.clone()
    }
    fn on_batch(&mut self, _: &EvidenceBatch) -> Result<(), EvidenceError> {
        Ok(())
    }
}
struct CanonicalConsumer<'a> {
    analyzer: &'a mut dyn EvidenceAnalyzer,
    batch: EvidenceBatch,
    reserve: u64,
    budget: Option<u64>,
    emitted: u64,
    last: Option<(u32, u32)>,
}
impl<'a> CanonicalConsumer<'a> {
    fn new(
        analyzer: &'a mut dyn EvidenceAnalyzer,
        fields: EvidenceFields,
        reserve: u64,
        budget: Option<u64>,
    ) -> Self {
        let mut batch = EvidenceBatch::new(0, "", 0, fields, Vec::new());
        batch.reserve_exact(EVIDENCE_ARROW_BATCH_ROWS);
        Self {
            analyzer,
            batch,
            reserve,
            budget,
            emitted: 0,
            last: None,
        }
    }
    fn flush(&mut self) -> Result<(), EvidenceError> {
        if !self.batch.is_empty() {
            checkpoint()?;
            self.analyzer.on_batch(&self.batch)?;
            self.emitted += self.batch.len() as u64;
            self.batch.clear();
            check_evidence_runtime(self.budget)?;
        }
        Ok(())
    }
}
impl EvidenceAnalyzer for CanonicalConsumer<'_> {
    fn requirements(&self) -> EvidenceRequirements {
        EvidenceRequirements {
            retained_bytes: Some(self.reserve),
            ..self.analyzer.requirements()
        }
    }
    fn on_batch(&mut self, batch: &EvidenceBatch) -> Result<(), EvidenceError> {
        check_evidence_runtime(self.budget)?;
        if !self.batch.is_empty()
            && (self.batch.contig_id, self.batch.canonical_tile_start)
                != (batch.contig_id, batch.canonical_tile_start)
        {
            self.flush()?;
        }
        self.batch.contig_id = batch.contig_id;
        self.batch.contig.clone_from(&batch.contig);
        self.batch.canonical_tile_start = batch.canonical_tile_start;
        for row in batch.rows() {
            let key = (batch.contig_id, row.position);
            if self.last.is_some_and(|last| last >= key) {
                return Err(EvidenceError::InvalidInput(
                    "evidence stream repeats or reorders loci".into(),
                ));
            }
            self.last = Some(key);
            self.batch.push_row(row)?;
            if self.batch.len() == EVIDENCE_ARROW_BATCH_ROWS {
                self.flush()?;
            }
        }
        Ok(())
    }
    fn finish(&mut self) -> Result<(), EvidenceError> {
        self.flush()?;
        checkpoint()?;
        self.analyzer.finish()?;
        check_evidence_runtime(self.budget)
    }
}

enum Opened {
    Native {
        engine: EvidenceEngine,
        session: VerifiedInputSession,
    },
    Dataset {
        dataset: VerifiedEvidenceDataset,
        query: DatasetQuery,
        execution: EvidenceExecution,
        selection_session: Option<VerifiedInputSession>,
        artifacts: Vec<PathBuf>,
    },
}
impl Opened {
    fn fields(&self) -> EvidenceFields {
        match self {
            Self::Native { engine, .. } => engine.request().fields,
            Self::Dataset { query, .. } => query.fields,
        }
    }
    fn verify(&self) -> Result<(), EvidenceArtifactError> {
        match self {
            Self::Native { session, .. } => session.verify()?,
            Self::Dataset {
                dataset,
                selection_session,
                ..
            } => {
                dataset.verify_unchanged()?;
                if let Some(s) = selection_session {
                    s.verify()?;
                }
            }
        }
        Ok(())
    }
    fn input_paths(&self) -> Vec<PathBuf> {
        match self {
            Self::Native { session, .. } => session
                .identities()
                .iter()
                .map(|s| PathBuf::from(&s.path))
                .collect(),
            Self::Dataset {
                dataset,
                selection_session,
                ..
            } => dataset
                .source_hashes()
                .iter()
                .map(|s| PathBuf::from(&s.path))
                .chain(
                    selection_session
                        .iter()
                        .flat_map(|s| s.identities().iter().map(|s| PathBuf::from(&s.path))),
                )
                .collect(),
        }
    }
    fn dataset_root(&self) -> Option<&Path> {
        match self {
            Self::Dataset { dataset, .. } => Path::new(&dataset.source_hashes()[0].path).parent(),
            _ => None,
        }
    }
    fn plan(&mut self, analyzer: &dyn EvidenceAnalyzer) -> Result<u64, EvidenceArtifactError> {
        Ok(match self {
            Self::Native { engine, .. } => {
                engine.plan_for_analyzer(analyzer)?.predicted_peak_rss_bytes
            }
            Self::Dataset {
                dataset,
                query,
                execution,
                ..
            } => {
                dataset
                    .plan(query, analyzer, execution)?
                    .predicted_peak_rss_bytes
            }
        })
    }
    fn run(
        &mut self,
        analyzer: &mut dyn EvidenceAnalyzer,
    ) -> Result<EvidenceRunStats, EvidenceArtifactError> {
        Ok(match self {
            Self::Native { engine, .. } => engine.run(analyzer)?,
            Self::Dataset {
                dataset,
                query,
                execution,
                ..
            } => dataset.visit_batches(query, analyzer, execution)?,
        })
    }
    fn receipt_memory_bytes(&self) -> u64 {
        match self {
            Self::Native { engine, session } => engine
                .plan()
                .sample_scope_bytes
                .saturating_mul(8)
                .saturating_add(engine.plan().selection_bytes.saturating_mul(8))
                .saturating_add(
                    engine
                        .contigs()
                        .iter()
                        .map(|c| c.name.len() as u64 + 64)
                        .sum::<u64>()
                        .saturating_mul(64),
                )
                .saturating_add(
                    session
                        .identities()
                        .iter()
                        .map(|s| s.path.len() as u64 + 256)
                        .sum::<u64>()
                        .saturating_mul(16),
                ),
            Self::Dataset { dataset, query, .. } => {
                let scope = &dataset.descriptor().sample_scope;
                let scope_bytes = scope
                    .declared_samples
                    .iter()
                    .map(|s| s.len() as u64 + 64)
                    .sum::<u64>()
                    .saturating_add(
                        scope
                            .read_groups
                            .iter()
                            .map(|g| {
                                g.id.len() as u64
                                    + g.sample.as_ref().map_or(0, |s| s.len()) as u64
                                    + 128
                            })
                            .sum::<u64>(),
                    )
                    .saturating_add(scope.selected_sample.as_ref().map_or(0, |s| s.len()) as u64);
                let query_bytes = match &query.selection {
                    EvidenceSelection::Sites(sites) => (sites.len() as u64).saturating_mul(4096),
                    EvidenceSelection::Intervals(intervals) => {
                        (intervals.len() as u64).saturating_mul(1024)
                    }
                    EvidenceSelection::WholeGenome => {
                        (dataset.contigs().len() as u64).saturating_mul(1024)
                    }
                };
                (dataset.descriptor().partitions.len() as u64)
                    .saturating_mul(24 << 10)
                    .saturating_add(scope_bytes.saturating_mul(64))
                    .saturating_add(query_bytes)
            }
        }
    }
    fn params(&self) -> Result<BTreeMap<String, String>, EvidenceArtifactError> {
        let (compatibility, query, sample, profile) = match self {
            Self::Native { engine, session } => (
                session.compatibility_key(engine)?,
                DatasetQuery {
                    selection: engine.request().selection.clone(),
                    fields: engine.request().fields,
                },
                engine.sample_scope().canonical_json(),
                engine.request().profile.clone(),
            ),
            Self::Dataset { dataset, query, .. } => (
                dataset.descriptor().compatibility_blake3.clone(),
                query.clone(),
                dataset.descriptor().sample_scope.canonical_json()?,
                dataset.descriptor().profile.to_profile(),
            ),
        };
        let contigs = match self {
            Self::Native { engine, .. } => engine.contigs(),
            Self::Dataset { dataset, .. } => dataset.contigs(),
        };
        Ok(BTreeMap::from([
            ("evidence.compatibility_blake3".into(), compatibility),
            (
                "evidence.query_blake3".into(),
                dataset_query_digest(&query, contigs)?,
            ),
            ("evidence.fields".into(), query.fields.bits().to_string()),
            (
                "evidence.fields_version".into(),
                query.fields.mask_version().to_string(),
            ),
            (
                "evidence.schema".into(),
                query.fields.schema_version().to_string(),
            ),
            (
                "evidence.semantics".into(),
                EVIDENCE_SEMANTICS_VERSION.into(),
            ),
            ("evidence.counting_unit".into(), "read".into()),
            ("evidence.sampling".into(), "none".into()),
            ("evidence.sample_scope".into(), sample),
            ("evidence.min_mapq".into(), profile.min_mapq.to_string()),
            (
                "evidence.min_base_quality".into(),
                profile.min_base_quality.to_string(),
            ),
            (
                "evidence.exclude_secondary".into(),
                profile.exclude_secondary.to_string(),
            ),
            (
                "evidence.exclude_supplementary".into(),
                profile.exclude_supplementary.to_string(),
            ),
            (
                "evidence.exclude_qc_fail".into(),
                profile.exclude_qc_fail.to_string(),
            ),
            (
                "evidence.exclude_duplicates".into(),
                profile.exclude_duplicates.to_string(),
            ),
        ]))
    }
    fn dependencies(&self, completed: bool) -> Result<Vec<FileHash>, EvidenceArtifactError> {
        Ok(match self {
            Self::Native { .. } => Vec::new(),
            Self::Dataset { dataset, .. } => dataset
                .source_hashes()
                .iter()
                .skip(1)
                .cloned()
                .chain(if completed {
                    dataset.verified_partition_hashes()?
                } else {
                    Vec::new()
                })
                .collect(),
        })
    }
    fn validate_artifacts(&self, dependencies: &[FileHash]) -> Result<(), EvidenceArtifactError> {
        if let Self::Dataset {
            dataset, artifacts, ..
        } = self
        {
            for path in artifacts {
                let hash = hash_file(path, true)?;
                if !dependencies
                    .iter()
                    .chain(dataset.source_hashes())
                    .any(|known| known.blake3 == hash)
                {
                    return Err(invalid(
                        "--dataset-artifact does not match a verified dependency of this query",
                    ));
                }
            }
        }
        Ok(())
    }
}

fn open_source(
    source: &EvidenceArtifactSource,
    capture: &mut CommandCapture,
) -> Result<Opened, EvidenceArtifactError> {
    match source {
        EvidenceArtifactSource::Native { request, selection } => {
            if matches!(selection, ArtifactSelection::Stored) {
                return Err(invalid("Stored selection requires a persisted dataset"));
            }
            let whole_genome = matches!(request.selection, EvidenceSelection::WholeGenome);
            let mut request = request.clone();
            if request.alignment_index.is_none() {
                request.alignment_index = ["bai", "csi", "crai"]
                    .iter()
                    .flat_map(|extension| {
                        [
                            PathBuf::from(format!(
                                "{}.{}",
                                request.alignments.display(),
                                extension
                            )),
                            request.alignments.with_extension(extension),
                        ]
                    })
                    .find(|path| path.is_file());
            }
            if request.alignment_index.is_none() {
                return Err(invalid("indexed alignments require BAI/CSI/CRAI"));
            }
            for (reference, fai) in [
                (&request.reference, &mut request.reference_fai),
                (&request.cram_reference, &mut request.cram_reference_fai),
            ] {
                if fai.is_none() {
                    *fai = reference
                        .as_ref()
                        .map(|p| suffix(p, ".fai"))
                        .filter(|p| p.is_file());
                }
            }
            let mut inputs = vec![("alignments".into(), request.alignments.clone())];
            for (role, path) in [
                ("reference", &request.reference),
                ("alignment-index", &request.alignment_index),
                ("reference-fai", &request.reference_fai),
                ("cram-reference", &request.cram_reference),
                ("cram-reference-fai", &request.cram_reference_fai),
            ] {
                if let Some(path) = path {
                    inputs.push((role.into(), path.clone()));
                }
            }
            if let Some((role, path)) = selection_file(selection) {
                inputs.push((role.into(), path.to_path_buf()));
            }
            let session = VerifiedInputSession::open(inputs)?;
            for source in session.identities() {
                capture.input_hashed(&format!("--{}", source.role), &source.path, &source.blake3);
            }
            capture_native(capture, &request);
            let mut engine = EvidenceEngine::open(request)?;
            let resolved = resolve_selection(
                selection,
                &engine.request().selection,
                engine.contigs(),
                None,
            )?;
            engine.set_selection(resolved)?;
            capture_selection(
                capture,
                selection,
                if whole_genome {
                    &EvidenceSelection::WholeGenome
                } else {
                    &engine.request().selection
                },
                engine.request().fields,
                engine.contigs(),
            )?;
            session.dataset_namespace(&engine)?;
            Ok(Opened::Native { engine, session })
        }
        EvidenceArtifactSource::Dataset {
            manifest,
            query,
            selection,
            execution,
            limits,
            artifacts,
        } => {
            let mut limits = *limits;
            limits.memory_budget_bytes = execution.memory_budget_bytes;
            let dataset = VerifiedEvidenceDataset::open(manifest, limits)?;
            let selection_session = selection_file(selection)
                .map(|(role, path)| VerifiedInputSession::open([(role.into(), path.to_path_buf())]))
                .transpose()?;
            let mut query = query.clone();
            query.selection = resolve_selection(
                selection,
                &query.selection,
                dataset.contigs(),
                Some(dataset.selection()),
            )?;
            let parent = &dataset.source_hashes()[0];
            capture.input_hashed("--dataset", &parent.path, &parent.blake3);
            if let Some(session) = &selection_session {
                for source in session.identities() {
                    capture.input_hashed(
                        &format!("--{}", source.role),
                        &source.path,
                        &source.blake3,
                    );
                }
            }
            capture.opt("--fields", query.fields.bits()).opt(
                "--max-dataset-metadata-bytes",
                limits.max_descriptor_bytes.min(limits.max_manifest_bytes),
            );
            capture_selection(
                capture,
                selection,
                &query.selection,
                query.fields,
                dataset.contigs(),
            )?;
            Ok(Opened::Dataset {
                dataset,
                query,
                execution: execution.clone(),
                selection_session,
                artifacts: artifacts.clone(),
            })
        }
    }
}
fn resolve_selection(
    selection: &ArtifactSelection,
    requested: &EvidenceSelection,
    contigs: &crate::core::ContigSet,
    stored: Option<EvidenceSelection>,
) -> Result<EvidenceSelection, EvidenceArtifactError> {
    Ok(match selection {
        ArtifactSelection::Request => requested.clone(),
        ArtifactSelection::Stored => {
            stored.ok_or_else(|| invalid("Stored selection requires dataset input"))?
        }
        ArtifactSelection::Bed(path) => EvidenceSelection::from_bed(path, contigs)?,
        ArtifactSelection::Variants(path, limits) => {
            EvidenceSelection::from_variants(path, contigs, *limits)?
        }
    })
}
fn selection_file(selection: &ArtifactSelection) -> Option<(&'static str, &Path)> {
    match selection {
        ArtifactSelection::Bed(path) => Some(("regions", path)),
        ArtifactSelection::Variants(path, _) => Some(("sites", path)),
        _ => None,
    }
}
fn capture_selection(
    capture: &mut CommandCapture,
    selection: &ArtifactSelection,
    requested: &EvidenceSelection,
    fields: EvidenceFields,
    contigs: &crate::core::ContigSet,
) -> Result<(), EvidenceArtifactError> {
    match selection {
        ArtifactSelection::Stored => {
            capture.flag("--whole-dataset");
        }
        ArtifactSelection::Request if matches!(requested, EvidenceSelection::WholeGenome) => {
            capture.flag("--whole-genome");
        }
        ArtifactSelection::Request => {
            let json = canonical_dataset_query(
                &DatasetQuery {
                    selection: requested.clone(),
                    fields,
                },
                contigs,
            )?;
            if json.len() > MAX_ARTIFACT_QUERY_BYTES {
                return Err(invalid(
                    "inline query exceeds32KiB; use a file-backed BED/VCF selection",
                ));
            }
            capture.opt("--query-json", json);
        }
        ArtifactSelection::Variants(_, limits) => {
            capture
                .opt("--max-variant-header-bytes", limits.max_header_bytes)
                .opt("--max-variant-record-bytes", limits.max_record_bytes);
        }
        ArtifactSelection::Bed(_) => {}
    }
    Ok(())
}
fn capture_native(capture: &mut CommandCapture, request: &EvidenceRequest) {
    capture
        .opt("--fields", request.fields.bits())
        .opt("--mapq-threshold", request.profile.min_mapq)
        .opt("--base-quality-threshold", request.profile.min_base_quality);
    for (flag, excluded) in [
        ("--include-secondary", request.profile.exclude_secondary),
        (
            "--include-supplementary",
            request.profile.exclude_supplementary,
        ),
        ("--include-qc-fail", request.profile.exclude_qc_fail),
        ("--include-duplicates", request.profile.exclude_duplicates),
    ] {
        capture.flag_if(!excluded, flag);
    }
    match &request.sample_selection {
        EvidenceSampleSelection::Auto => {}
        EvidenceSampleSelection::Named(sample) => {
            capture.opt("--sample", sample);
        }
        EvidenceSampleSelection::Pool => {
            capture.flag("--pool-samples");
        }
    }
}
fn capture_common(capture: &mut CommandCapture, spec: &EvidenceArtifactSpec, budget: Option<u64>) {
    let execution = execution(&spec.source);
    capture
        .opt("--tile-bases", execution.max_microtile_bases)
        .opt("--max-read-len", execution.max_read_len)
        .opt("--max-record-bytes", execution.max_record_bytes)
        .opt("--max-receipt-bytes", spec.max_receipt_bytes);
    if let Some(bytes) = budget {
        capture.opt("--memory-budget-bytes", bytes);
    }
    capture
        .flag_if(spec.enforcement != EnforcementMode::RecordOnly, "--enforce")
        .flag_if(
            spec.enforcement == EnforcementMode::RequireOsLimit,
            "--require-os-limit",
        )
        .flag_if(spec.output_policy == OutputPolicy::ReplaceAtomic, "--force");
    for (flag, value) in &spec.invocation.options {
        capture.opt(flag, value);
    }
    for flag in &spec.invocation.flags {
        capture.flag(flag);
    }
}

#[allow(clippy::too_many_arguments)]
fn publish(
    mut pending: AtomicFile,
    receipt_path: &Path,
    spec: &EvidenceArtifactSpec,
    requirements: &EvidenceRequirements,
    token: &CancellationToken,
    opened: &mut Opened,
    mut capture: CommandCapture,
    params: BTreeMap<String, String>,
    source_params: BTreeMap<String, String>,
    mut stats: EvidenceRunStats,
    predicted: u64,
    started: Instant,
    os_limit: Option<u64>,
    failure: Option<&str>,
    budget: Option<u64>,
) -> Result<EvidenceArtifactOutcome, EvidenceArtifactError> {
    let initial_completed = failure.is_none();
    let mut failure = failure.map(str::to_owned);
    pending.file_mut().flush()?;
    pending.file_mut().sync_all()?;
    let hash = match hash_file_with_token(pending.temporary_path(), failure.is_none(), Some(token))
    {
        Ok(hash) => hash,
        Err(EvidenceArtifactError::Resource { message, .. }) => {
            failure = Some(message);
            hash_file_with_token(pending.temporary_path(), false, Some(token))?
        }
        Err(error) => return Err(error),
    };
    let dependencies = match opened.dependencies(failure.is_none()) {
        Ok(dependencies) => dependencies,
        Err(EvidenceArtifactError::Resource { message, .. }) => {
            failure = Some(message);
            opened.dependencies(false)?
        }
        Err(error) => return Err(error),
    };
    if failure.is_none() {
        match opened.validate_artifacts(&dependencies) {
            Ok(()) => {}
            Err(EvidenceArtifactError::Resource { message, .. }) => failure = Some(message),
            Err(error) => return Err(error),
        }
    }
    let mut completed = failure.is_none();
    let mut output = if completed {
        spec.output.clone()
    } else {
        suffix(&spec.output, ".partial")
    };
    let mut final_receipt_path = if initial_completed && !completed {
        suffix(receipt_path, ".partial")
    } else {
        receipt_path.to_path_buf()
    };
    for dependency in dependencies {
        capture.input_hashed("--dataset-artifact", &dependency.path, &dependency.blake3);
    }
    capture.output_hashed("-o", &output.display().to_string(), &hash);
    let mut receipt = RunManifest::new(spec.invocation.argv_prefix.join(" "));
    capture.record_into(&mut receipt);
    receipt.tool_version = spec.producer.version.clone();
    let science = serde_json::to_vec(&(
        &spec.analyzer.id,
        &spec.analyzer.version,
        &params,
        &source_params,
    ))
    .map_err(|e| invalid(&e.to_string()))?;
    receipt.params.extend(source_params);
    for (key, value) in params {
        receipt
            .params
            .insert(format!("{}{key}", spec.analyzer.param_prefix), value);
    }
    for (key, value) in [
        (
            "run_status",
            if completed {
                "completed"
            } else {
                "resource-failure"
            }
            .to_string(),
        ),
        ("producer.name", spec.producer.name.clone()),
        ("producer.version", spec.producer.version.clone()),
        ("producer.binary", spec.producer.binary.clone()),
        ("analyzer.id", spec.analyzer.id.clone()),
        ("analyzer.version", spec.analyzer.version.clone()),
        ("replay.kind", "external-analyzer".into()),
        (
            "evidence.artifact_semantics",
            "managed-evidence-artifact-v1".into(),
        ),
        (
            "science.blake3",
            blake3::hash(&science).to_hex().to_string(),
        ),
        (
            "artifact.output.0.role",
            if completed {
                "analyzer-output"
            } else {
                "partial-analyzer-output"
            }
            .into(),
        ),
        (
            "contract.assurance",
            match spec.enforcement {
                EnforcementMode::RecordOnly => "observed-only",
                EnforcementMode::Cooperative => "declared-bound-cooperative",
                EnforcementMode::RequireOsLimit => "cgroup-v2",
            }
            .into(),
        ),
        (
            "contract.canonical_batch_rows",
            EVIDENCE_ARROW_BATCH_ROWS.to_string(),
        ),
        ("outcome.rows", stats.emitted_loci.to_string()),
    ] {
        receipt.params.insert(key.into(), value);
    }
    receipt.params.extend([
        (
            "analyzer.required_fields".into(),
            requirements.fields.bits().to_string(),
        ),
        (
            "analyzer.requires_reference".into(),
            requirements.requires_reference.to_string(),
        ),
        (
            "analyzer.context_bases".into(),
            requirements.context_bases.to_string(),
        ),
        (
            "analyzer.memory_model".into(),
            if requirements.retained_bytes.is_some() {
                "declared-evidence-factory-v1"
            } else {
                "unknown"
            }
            .into(),
        ),
    ]);
    if let Some(bytes) = requirements.retained_bytes {
        receipt
            .measurements
            .insert("analyzer.max_additional_bytes".into(), bytes.to_string());
    }
    if let Some(repository) = &spec.producer.repository {
        receipt
            .params
            .insert("producer.repository".into(), repository.clone());
    }
    if let Some(failure) = &failure {
        receipt
            .params
            .insert("resource.failure".into(), failure.into());
    }
    if let Some(limit) = os_limit {
        receipt
            .measurements
            .insert("resource.os_limit_bytes".into(), limit.to_string());
    }
    if let Opened::Native { engine, .. } = opened {
        receipt
            .measurements
            .extend(engine.plan().decoder_measurements());
    }
    receipt.measurements.extend([
        ("predicted_peak_rss_bytes".into(), predicted.to_string()),
        ("peak_rss_bytes".into(), rss().to_string()),
        (
            "execution.elapsed_ms".into(),
            started.elapsed().as_millis().to_string(),
        ),
        (
            "execution.alignment_record_visits".into(),
            stats.record_visits.to_string(),
        ),
        ("execution.microtiles".into(), stats.microtiles.to_string()),
        (
            "execution.max_returned_record_bytes".into(),
            stats.max_record_bytes.to_string(),
        ),
        (
            "execution.max_returned_read_length".into(),
            stats.max_read_length.to_string(),
        ),
    ]);
    let declared = receipt
        .params
        .get("memory_budget_bytes")
        .and_then(|value| value.parse::<u64>().ok());
    if let Some(bytes) = declared.filter(|bytes| bytes % (1 << 20) == 0) {
        receipt
            .params
            .insert("memory_budget_mb".into(), (bytes / (1 << 20)).to_string());
    }
    receipt.params.insert(
        "governor".into(),
        match (spec.enforcement, completed) {
            (EnforcementMode::RecordOnly, _) => "record-only",
            (_, true) => "enforced",
            (_, false) => "tripped",
        }
        .into(),
    );
    let mut receipt_file = AtomicFile::create(receipt_path)?;
    // A self-referential RSS fixed point is not a completion condition: input
    // revalidation and serialization may move the high-water mark every time.
    // Allocate canonical JSON once, observe after encoding/sync/validation, then
    // patch only the measurement values and digest in that existing buffer.
    let peak = loop {
        token.check().map_err(EvidenceError::from)?;
        receipt.params.insert(
            "resource.peak_sampling_phase".into(),
            "post-encoding-sync-and-source-validation-before-atomic-commit".into(),
        );
        receipt.params.insert(
            "contract_verdict".into(),
            if declared.is_some() {
                "within"
            } else {
                "unset"
            }
            .into(),
        );
        receipt
            .measurements
            .insert("peak_rss_bytes".into(), "00000000000000000000".into());
        receipt.finalize();
        let mut bytes = receipt.to_canonical_json();
        if bytes.len() > spec.max_receipt_bytes {
            return Err(invalid("receipt exceeds declared metadata envelope"));
        }
        write_staged_receipt(&mut receipt_file, &bytes)?;
        pending.file_mut().sync_all()?;
        let checked = if completed {
            opened.verify().and_then(|()| check_runtime(budget))
        } else {
            token
                .check()
                .map_err(EvidenceError::from)
                .map_err(Into::into)
        };
        let checked = checked.and_then(|()| {
            let peak = rss();
            reseal_observed_measurements(&mut bytes, peak, declared)?;
            write_staged_receipt(&mut receipt_file, &bytes)?;
            // No scientific work remains. A late breach still changes the
            // transaction to explicit partials; cancellation publishes nothing.
            if completed {
                check_runtime(budget)?;
            }
            token.check().map_err(EvidenceError::from)?;
            Ok(peak)
        });
        match checked {
            Ok(peak) => break peak,
            Err(EvidenceArtifactError::Resource { message, .. }) if completed => {
                completed = false;
                failure = Some(message);
                output = suffix(&spec.output, ".partial");
                final_receipt_path = suffix(receipt_path, ".partial");
                receipt.outputs[0].path = output.display().to_string();
                receipt
                    .params
                    .insert("run_status".into(), "resource-failure".into());
                receipt.params.insert(
                    "artifact.output.0.role".into(),
                    "partial-analyzer-output".into(),
                );
                receipt
                    .params
                    .insert("resource.failure".into(), failure.as_ref().unwrap().clone());
                receipt.params.insert("governor".into(), "tripped".into());
                // Drop the old canonical buffer before allocating the partial
                // receipt. This state transition can occur only once.
            }
            Err(error) => return Err(error),
        }
    };
    token.check().map_err(EvidenceError::from)?;
    commit_group(
        vec![
            (pending, output.clone()),
            (receipt_file, final_receipt_path.clone()),
        ],
        spec.output_policy == OutputPolicy::ReplaceAtomic,
    )?;
    stats.peak_rss_bytes = peak;
    let outcome = EvidenceArtifactOutcome {
        output,
        manifest: final_receipt_path,
        claim_hash: receipt.params["manifest_blake3"].clone(),
        completed,
        stats,
        predicted_peak_rss_bytes: predicted,
        peak_rss_bytes: peak,
    };
    if initial_completed && !completed {
        Err(EvidenceArtifactError::Resource {
            message: failure.unwrap(),
            partial: Some(Box::new(outcome)),
        })
    } else {
        Ok(outcome)
    }
}

fn write_staged_receipt(file: &mut AtomicFile, bytes: &str) -> Result<(), EvidenceArtifactError> {
    file.file_mut().set_len(0)?;
    std::io::Seek::rewind(file.file_mut())?;
    file.file_mut().write_all(bytes.as_bytes())?;
    file.file_mut().sync_all()?;
    Ok(())
}

// All edits retain or shorten the existing String; hashing borrows its slices.
// Measurement hashes omit their own entry, exactly as RunManifest does. Claim
// bytes are untouched because measurements are excluded from schema-3+ claims.
fn reseal_observed_measurements(
    json: &mut String,
    peak: u64,
    declared: Option<u64>,
) -> Result<(), EvidenceArtifactError> {
    let marker = "\"peak_rss_bytes\":\"";
    let start = json
        .find(marker)
        .ok_or_else(|| invalid("missing peak measurement"))?
        + marker.len();
    let mut digits = [b'0'; 20];
    let mut remaining = peak;
    for digit in digits.iter_mut().rev() {
        *digit += (remaining % 10) as u8;
        remaining /= 10;
    }
    // Replacing the 20-byte ASCII decimal placeholder keeps the existing
    // allocation and valid UTF-8 intact.
    json.replace_range(start..start + 20, std::str::from_utf8(&digits).unwrap());
    let verdict = match declared {
        None => "unset",
        Some(limit) if peak <= limit => "within",
        Some(_) => "over",
    };
    let verdict_marker = "\"contract_verdict\":\"";
    let verdict_start = json
        .find(verdict_marker)
        .ok_or_else(|| invalid("missing resource verdict"))?
        + verdict_marker.len();
    let verdict_end = json[verdict_start..]
        .find('"')
        .ok_or_else(|| invalid("invalid resource verdict"))?
        + verdict_start;
    json.replace_range(verdict_start..verdict_end, verdict);
    let map_marker = ",\"measurements\":";
    let map_start = json
        .find(map_marker)
        .ok_or_else(|| invalid("missing measurements"))?
        + map_marker.len();
    let map_end = json[map_start..]
        .find("},\"outputs\":")
        .ok_or_else(|| invalid("missing measurement terminator"))?
        + map_start
        + 1;
    let hash_marker = "\"measurement_blake3\":\"";
    let entry_start = json[map_start..map_end]
        .find(hash_marker)
        .ok_or_else(|| invalid("missing measurement digest"))?
        + map_start;
    let hash_start = entry_start + hash_marker.len();
    let entry_end = hash_start + 64 + 1;
    let (omit_start, omit_end) = if json.as_bytes().get(entry_end) == Some(&b',') {
        (entry_start, entry_end + 1)
    } else {
        (entry_start - 1, entry_end)
    };
    let mut hasher = blake3::Hasher::new();
    hasher.update(&json.as_bytes()[map_start..omit_start]);
    hasher.update(&json.as_bytes()[omit_end..map_end]);
    let digest = hasher.finalize().to_hex();
    json.replace_range(hash_start..hash_start + 64, digest.as_str());
    Ok(())
}

fn validate_spec(
    spec: &EvidenceArtifactSpec,
    params: &BTreeMap<String, String>,
) -> Result<(), EvidenceArtifactError> {
    if spec.max_receipt_bytes == 0 {
        return Err(invalid("receipt envelope must be positive"));
    }
    if spec.producer.name.is_empty()
        || spec.producer.name == "rosalind"
        || spec.producer.version.is_empty()
        || spec.producer.binary.is_empty()
        || spec.analyzer.id.is_empty()
        || spec.analyzer.version.is_empty()
    {
        return Err(invalid(
            "external producer/analyzer identity must be nonempty and cannot impersonate rosalind",
        ));
    }
    if spec.analyzer.param_prefix != "analyzer."
        && !(spec.analyzer.param_prefix.starts_with("x.")
            && spec.analyzer.param_prefix.ends_with('.'))
    {
        return Err(invalid(
            "external analyzer parameters require analyzer. or x.*. namespace",
        ));
    }
    if spec.invocation.argv_prefix.is_empty()
        || spec.invocation.argv_prefix.iter().any(|token| {
            token.is_empty()
                || token.bytes().any(|b| b.is_ascii_whitespace() || b == 0)
                || token.starts_with('-')
        })
    {
        return Err(invalid(
            "replay prefix must contain nonempty plain command tokens",
        ));
    }
    let reserved: BTreeSet<&str> = [
        "--alignments",
        "--reference",
        "--alignment-index",
        "--reference-fai",
        "--cram-reference",
        "--cram-reference-fai",
        "--dataset",
        "--dataset-artifact",
        "--regions",
        "--sites",
        "--whole-genome",
        "--whole-dataset",
        "--query-json",
        "--fields",
        "--mapq-threshold",
        "--base-quality-threshold",
        "--include-secondary",
        "--include-supplementary",
        "--include-qc-fail",
        "--include-duplicates",
        "--sample",
        "--pool-samples",
        "--tile-bases",
        "--max-read-len",
        "--max-record-bytes",
        "--max-variant-header-bytes",
        "--max-variant-record-bytes",
        "--max-receipt-bytes",
        "--max-dataset-metadata-bytes",
        "--memory-budget-mb",
        "--memory-budget-bytes",
        "--enforce",
        "--require-os-limit",
        "--force",
        "--output",
        "-o",
        "--manifest",
    ]
    .into_iter()
    .collect();
    for flag in spec.invocation.options.keys().chain(&spec.invocation.flags) {
        if !flag.starts_with("--")
            || flag.len() < 3
            || !flag
                .bytes()
                .all(|b| b.is_ascii_lowercase() || b.is_ascii_digit() || b == b'-')
            || reserved.contains(flag.as_str())
            || spec.invocation.options.contains_key(flag) && spec.invocation.flags.contains(flag)
        {
            return Err(invalid(
                "analyzer replay flag is malformed, duplicate or reserved",
            ));
        }
    }
    if spec.invocation.options.values().any(|value| {
        value == "true"
            || value.starts_with("@in:")
            || value.starts_with("@out:")
            || value.contains('\0')
    }) {
        return Err(invalid("analyzer option values cannot use reserved recipe markers or the bare-flag true representation"));
    }
    if params.keys().any(|key| {
        key.is_empty()
            || matches!(
                key.as_str(),
                "id" | "version"
                    | "required_fields"
                    | "requires_reference"
                    | "context_bases"
                    | "memory_model"
                    | "max_additional_bytes"
            )
            || key.bytes().any(|b| b.is_ascii_control())
    }) {
        return Err(invalid("analyzer parameters have empty/reserved keys"));
    }
    let text = serde_json::to_vec(&(
        params,
        &spec.invocation.argv_prefix,
        &spec.invocation.options,
        &spec.invocation.flags,
        &spec.producer.name,
        &spec.producer.version,
        &spec.producer.binary,
        &spec.producer.repository,
        &spec.analyzer.id,
        &spec.analyzer.version,
        &spec.analyzer.param_prefix,
    ))
    .map_err(|e| invalid(&e.to_string()))?;
    if text.len() > 64 << 10 || text.len() > spec.max_receipt_bytes / 4 {
        return Err(invalid(
            "analyzer identity, parameters or replay recipe exceed metadata envelope",
        ));
    }
    Ok(())
}
fn validate_paths(
    spec: &EvidenceArtifactSpec,
    receipt: &Path,
    inputs: &[PathBuf],
    dataset_root: Option<&Path>,
) -> Result<(), EvidenceArtifactError> {
    let paths = [
        spec.output.clone(),
        receipt.to_path_buf(),
        suffix(&spec.output, ".partial"),
        suffix(receipt, ".partial"),
    ];
    let mut destinations = Vec::<PathBuf>::new();
    for path in paths {
        if path.exists() && spec.output_policy == OutputPolicy::CreateNewAtomic {
            return Err(EvidenceArtifactError::OutputExists(path));
        }
        let resolved = resolve_destination(&path)?;
        if dataset_root.is_some_and(|root| resolved.starts_with(root)) {
            return Err(invalid(
                "write derived artifacts outside the immutable dataset directory",
            ));
        }
        if inputs.iter().any(|input| same_file(&resolved, input))
            || destinations.iter().any(|other| same_file(&resolved, other))
        {
            return Err(invalid(
                "output, receipt and partial paths must be distinct and cannot replace inputs",
            ));
        }
        destinations.push(resolved);
    }
    Ok(())
}
fn resolve_destination(path: &Path) -> Result<PathBuf, EvidenceArtifactError> {
    if path.exists() {
        if !path.is_file() {
            return Err(invalid("artifact destination is not a regular file"));
        }
        return Ok(fs::canonicalize(path)?);
    }
    let parent = path
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    Ok(fs::canonicalize(parent)?.join(
        path.file_name()
            .ok_or_else(|| invalid("artifact destination has no filename"))?,
    ))
}
fn same_file(left: &Path, right: &Path) -> bool {
    if left == right {
        return true;
    }
    #[cfg(unix)]
    {
        use std::os::unix::fs::MetadataExt;
        if let (Ok(a), Ok(b)) = (fs::metadata(left), fs::metadata(right)) {
            return a.dev() == b.dev() && a.ino() == b.ino();
        }
    }
    false
}
fn execution(source: &EvidenceArtifactSource) -> &EvidenceExecution {
    match source {
        EvidenceArtifactSource::Native { request, .. } => &request.execution,
        EvidenceArtifactSource::Dataset { execution, .. } => execution,
    }
}
fn execution_mut(source: &mut EvidenceArtifactSource) -> &mut EvidenceExecution {
    match source {
        EvidenceArtifactSource::Native { request, .. } => &mut request.execution,
        EvidenceArtifactSource::Dataset { execution, .. } => execution,
    }
}
fn suffix(path: &Path, ending: &str) -> PathBuf {
    let mut path = path.as_os_str().to_os_string();
    path.push(ending);
    PathBuf::from(path)
}
fn rss() -> u64 {
    crate::util::rss::peak_rss_bytes()
}
fn check_evidence_runtime(budget: Option<u64>) -> Result<(), EvidenceError> {
    checkpoint()?;
    if let Some(budget) = budget {
        let needed = rss();
        if needed > budget {
            return Err(CoreError::BudgetExceeded { needed, budget }.into());
        }
    }
    Ok(())
}
fn check_runtime(budget: Option<u64>) -> Result<(), EvidenceArtifactError> {
    Ok(check_evidence_runtime(budget)?)
}
fn hash_file(path: &Path, checked: bool) -> Result<String, EvidenceArtifactError> {
    hash_file_with_token(path, checked, None)
}
fn hash_file_with_token(
    path: &Path,
    checked: bool,
    token: Option<&CancellationToken>,
) -> Result<String, EvidenceArtifactError> {
    let mut file = File::open(path)?;
    let mut hash = blake3::Hasher::new();
    let mut bytes = [0u8; 65536];
    loop {
        if let Some(token) = token {
            token.check().map_err(EvidenceError::from)?;
        }
        if checked {
            check_runtime(None)?;
        }
        let n = file.read(&mut bytes)?;
        if n == 0 {
            break;
        }
        hash.update(&bytes[..n]);
    }
    Ok(hash.finalize().to_hex().to_string())
}
fn invalid(message: &str) -> EvidenceArtifactError {
    EvidenceArtifactError::InvalidConfiguration(message.into())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn in_place_measurement_reseal_preserves_claim_hash_and_allocation() {
        for (peak, declared, verdict) in [
            (0, None, "unset"),
            (42, Some(42), "within"),
            (43, Some(42), "over"),
            (u64::MAX, Some(u64::MAX), "within"),
            (u64::MAX, Some(42), "over"),
        ] {
            let mut receipt = RunManifest::new("example run");
            receipt
                .params
                .insert("science.example".into(), "stable".into());
            receipt.params.insert(
                "contract_verdict".into(),
                if declared.is_some() {
                    "within"
                } else {
                    "unset"
                }
                .into(),
            );
            receipt
                .measurements
                .insert("peak_rss_bytes".into(), "00000000000000000000".into());
            receipt.finalize();
            let claim = receipt.params["manifest_blake3"].clone();
            let mut json = receipt.to_canonical_json();
            let capacity = json.capacity();
            let pointer = json.as_ptr();
            reseal_observed_measurements(&mut json, peak, declared).unwrap();
            assert_eq!(json.capacity(), capacity);
            assert_eq!(json.as_ptr(), pointer);
            let decoded = RunManifest::from_canonical_json(&json).unwrap();
            assert_eq!(decoded.params["manifest_blake3"], claim);
            assert_eq!(decoded.self_hash_ok(), Some(true));
            assert_eq!(decoded.measurement_hash_ok(), Some(true));
            assert_eq!(
                decoded
                    .get_recorded("peak_rss_bytes")
                    .unwrap()
                    .parse::<u64>()
                    .unwrap(),
                peak
            );
            assert_eq!(decoded.get_recorded("contract_verdict").unwrap(), verdict);
            assert_eq!(decoded.to_canonical_json(), json);
        }
    }

    #[test]
    fn resource_failure_during_publication_seals_partials() {
        // The governor is process-wide: exercise its publication-time transition
        // in a child process so concurrent ordinary library tests stay unarmed.
        const CHILD: &str = "ROSALIND_ARTIFACT_PUBLICATION_CHILD";
        if std::env::var_os(CHILD).is_none() {
            let output = std::process::Command::new(std::env::current_exe().unwrap())
                .args([
                    "--exact",
                    "evidence::artifact::tests::resource_failure_during_publication_seals_partials",
                    "--nocapture",
                ])
                .env(CHILD, "1")
                .output()
                .unwrap();
            assert!(
                output.status.success(),
                "{}{}",
                String::from_utf8_lossy(&output.stdout),
                String::from_utf8_lossy(&output.stderr)
            );
            return;
        }
        let root =
            std::env::temp_dir().join(format!("rosalind-publish-resource-{}", std::process::id()));
        fs::create_dir_all(&root).unwrap();
        let reference = root.join("reference.fa");
        fs::write(&reference, ">chr1\nAA\n").unwrap();
        fs::write(root.join("reference.fa.fai"), "chr1\t2\t6\t2\t3\n").unwrap();
        let bam_path = root.join("reads.bam");
        let mut header = rust_htslib::bam::Header::new();
        header.push_record(
            rust_htslib::bam::header::HeaderRecord::new(b"SQ")
                .push_tag(b"SN", "chr1")
                .push_tag(b"LN", 2),
        );
        drop(
            rust_htslib::bam::Writer::from_path(&bam_path, &header, rust_htslib::bam::Format::Bam)
                .unwrap(),
        );
        rust_htslib::bam::index::build(
            &bam_path,
            None::<&PathBuf>,
            rust_htslib::bam::index::Type::Bai,
            1,
        )
        .unwrap();
        let mut request = EvidenceRequest::new(&bam_path, &reference);
        request.fields = EvidenceFields::DEPTHS;
        let spec = EvidenceArtifactSpec::new(
            EvidenceArtifactSource::Native {
                request,
                selection: ArtifactSelection::Request,
            },
            root.join("result.tsv"),
            ProducerIdentity {
                name: "publication-test".into(),
                version: "1".into(),
                binary: "publication-test".into(),
                repository: None,
            },
            AnalyzerIdentity::new("publication-test", "1"),
            ReplayInvocation::new(["run"]),
        );
        let mut capture = CommandCapture::new("run");
        let mut opened = open_source(&spec.source, &mut capture).unwrap();
        let params = opened.params().unwrap();
        let mut pending = AtomicFile::create(&spec.output).unwrap();
        pending
            .file_mut()
            .write_all(b"complete staged custom output\n")
            .unwrap();
        let governor = MemoryGovernor::start(1, Duration::from_millis(100), || 2).unwrap();
        let requirements = EvidenceRequirements {
            fields: EvidenceFields::DEPTHS,
            requires_reference: false,
            context_bases: 0,
            retained_bytes: Some(0),
        };
        let error = publish(
            pending,
            &root.join("result.tsv.manifest.json"),
            &spec,
            &requirements,
            &CancellationToken::new(),
            &mut opened,
            capture,
            BTreeMap::new(),
            params,
            EvidenceRunStats::default(),
            0,
            Instant::now(),
            None,
            None,
            Some(1),
        )
        .unwrap_err();
        drop(governor);
        let EvidenceArtifactError::Resource {
            partial: Some(partial),
            ..
        } = error
        else {
            panic!("late resource failure did not preserve identified partials");
        };
        assert!(!partial.completed);
        assert!(!spec.output.exists());
        assert!(!root.join("result.tsv.manifest.json").exists());
        let receipt =
            RunManifest::from_canonical_json(&fs::read_to_string(&partial.manifest).unwrap())
                .unwrap();
        assert_eq!(receipt.params["run_status"], "resource-failure");
        assert_eq!(
            receipt.params["artifact.output.0.role"],
            "partial-analyzer-output"
        );
        assert_eq!(
            receipt.outputs[0].blake3,
            hash_file(&partial.output, true).unwrap()
        );
        assert_eq!(receipt.self_hash_ok(), Some(true));
        assert_eq!(receipt.measurement_hash_ok(), Some(true));
        fs::remove_dir_all(root).unwrap();
    }
}
