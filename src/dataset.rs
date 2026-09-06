//! Canonical evidence partition caching and bounded first-party parallel work.
//!
//! Cache identity depends on verified scientific inputs and fixed coordinate
//! ownership, never on worker count or execution microtiles. Workers only produce
//! first-party Arrow evidence; the caller's analyzer consumes it serially.

mod lookup;
pub use lookup::VerifiedEvidenceLookup;
mod persisted;
pub use persisted::{canonical_dataset_query, dataset_query_digest};
pub use persisted::{
    DatasetCoverage, DatasetQuery, DatasetReadLimits, DatasetReadPlan, VerifiedEvidenceDataset,
};
mod parquet;
pub use parquet::plan_parquet_export_with_inputs;
pub use parquet::{
    export_parquet_dataset, export_parquet_dataset_with_inputs, plan_parquet_export,
    ParquetExportOptions, ParquetExportOutcome, PARQUET_EXPORT_MANIFEST_NAME,
    PARQUET_EXPORT_SEMANTICS,
};
mod reuse;
pub use reuse::*;
mod descriptor;
pub use descriptor::*;

use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::{BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{mpsc, Arc, Mutex};

use thiserror::Error;

use crate::evidence::{
    evidence_reader_memory_bytes, read_evidence_batches_expected_fields, EvidenceAnalyzer,
    EvidenceArrowWriter, EvidenceEngine, EvidenceError, EvidenceExecution, EvidenceFields,
    EvidenceRunStats, EvidenceSelection, CANONICAL_TILE_BASES, EVIDENCE_SEMANTICS_VERSION,
};
use crate::provenance::{blake3_file, FileHash, RunManifest};
use crate::selection::GenomicInterval;
use crate::util::atomic::write_atomic;

const MAX_WORKERS: usize = 64;
const MANIFEST_NAME: &str = "dataset.manifest.json";
static TEMP_SEQUENCE: AtomicU64 = AtomicU64::new(0);

/// Cheap input-change detector for an immutable, content-verified input session.
/// This detects ordinary replacement, size and timestamp changes; it is not a
/// cryptographic identity or protection against deliberate metadata restoration.
#[derive(Debug, Clone)]
pub struct InputSnapshot {
    files: Arc<Vec<(PathBuf, InputStamp)>>,
}
#[derive(Debug, PartialEq, Eq)]
struct InputStamp {
    length: u64,
    modified: Option<std::time::SystemTime>,
    #[cfg(unix)]
    unix_identity: (u64, u64, i64, i64, i64, i64),
}
impl InputStamp {
    fn read(path: &Path) -> Result<Self, EvidenceError> {
        let metadata = fs::metadata(path)?;
        if !metadata.is_file() {
            return Err(EvidenceError::InvalidInput(format!(
                "input must remain a regular file: {}",
                path.display()
            )));
        }
        Ok(Self {
            length: metadata.len(),
            modified: metadata.modified().ok(),
            #[cfg(unix)]
            unix_identity: {
                use std::os::unix::fs::MetadataExt;
                (
                    metadata.dev(),
                    metadata.ino(),
                    metadata.mtime(),
                    metadata.mtime_nsec(),
                    metadata.ctime(),
                    metadata.ctime_nsec(),
                )
            },
        })
    }
}
impl InputSnapshot {
    /// Capture before hashing/opening, then verify after hashing and before
    /// publication. Callers must keep every source and index immutable throughout.
    pub fn capture(paths: impl IntoIterator<Item = PathBuf>) -> Result<Self, EvidenceError> {
        let files = paths
            .into_iter()
            .map(|path| {
                let stamp = InputStamp::read(&path)?;
                Ok((path, stamp))
            })
            .collect::<Result<Vec<_>, EvidenceError>>()?;
        Ok(Self {
            files: Arc::new(files),
        })
    }
    /// Refuse an ordinary mutation of any captured input.
    pub fn verify(&self) -> Result<(), EvidenceError> {
        for (path, before) in self.files.iter() {
            let after = InputStamp::read(path).map_err(|error| {
                EvidenceError::InvalidInput(format!(
                    "input changed during analysis: {}: {error}",
                    path.display()
                ))
            })?;
            if before != &after {
                return Err(EvidenceError::InvalidInput(format!(
                    "input changed during analysis: {}; all inputs must remain immutable",
                    path.display()
                )));
            }
        }
        Ok(())
    }
    fn for_engine(engine: &EvidenceEngine) -> Result<Self, EvidenceError> {
        let request = engine.request();
        let mut paths = vec![request.alignments.clone()];
        for path in [
            &request.alignment_index,
            &request.reference,
            &request.reference_fai,
            &request.cram_reference,
            &request.cram_reference_fai,
        ]
        .into_iter()
        .flatten()
        {
            paths.push(path.clone());
        }
        if request.alignment_index.is_none() {
            for extension in ["bai", "csi", "crai"] {
                for path in [
                    PathBuf::from(format!("{}.{}", request.alignments.display(), extension)),
                    request.alignments.with_extension(extension),
                ] {
                    if path.is_file() {
                        paths.push(path);
                    }
                }
            }
        }
        for path in [&request.reference, &request.cram_reference]
            .into_iter()
            .flatten()
        {
            let fai = PathBuf::from(format!("{}.fai", path.display()));
            if fai.is_file() {
                paths.push(fai);
            }
        }
        paths.sort();
        paths.dedup();
        Self::capture(paths)
    }
}

fn request_digest(engine: &EvidenceEngine) -> String {
    let request = engine.request();
    let mut hash = blake3::Hasher::new();
    let mut field = |value: &[u8]| {
        hash.update(&(value.len() as u64).to_le_bytes());
        hash.update(value);
    };
    field(b"rosalind-dataset-request-v1");
    field(EVIDENCE_SEMANTICS_VERSION.as_bytes());
    field(env!("CARGO_PKG_VERSION").as_bytes());
    field(&request.fields.schema_version().to_le_bytes());
    field(&request.fields.mask_version().to_le_bytes());
    field(&request.fields.bits().to_le_bytes());
    field(crate::evidence::EvidenceProfile::ID.as_bytes());
    field(engine.sample_scope().canonical_json().as_bytes());
    let profile = &request.profile;
    field(&[
        profile.min_mapq,
        profile.min_base_quality,
        u8::from(profile.exclude_secondary),
        u8::from(profile.exclude_supplementary),
        u8::from(profile.exclude_qc_fail),
        u8::from(profile.exclude_duplicates),
    ]);
    field(engine.selection_digest().as_bytes());
    field(&[u8::from(
        request.reference.is_some() || request.cram_reference.is_some(),
    )]);
    for contig in engine.contigs().iter() {
        field(&contig.id.to_le_bytes());
        field(contig.name.as_bytes());
        field(&contig.length.to_le_bytes());
    }
    hash.finalize().to_hex().to_string()
}

/// Cache location and bounded worker controls.
#[derive(Debug, Clone)]
pub struct DatasetOptions {
    /// Parent directory; scientific identities occupy separate subdirectories.
    pub cache_dir: PathBuf,
    /// Permit reuse of complete partitions from an interrupted dataset.
    pub resume: bool,
    /// Requested first-party workers, from 1 through 64.
    pub workers: usize,
}

/// Conservative aggregate admission model for parallel extraction and reduction.
#[derive(Debug, Clone)]
pub struct DatasetPlan {
    /// Fixed coordinate partitions, independent of execution settings.
    pub partition_count: usize,
    /// Active worker count, limited by the number of partitions.
    pub worker_count: usize,
    /// Per-worker computation width selected by the aggregate budget.
    pub microtile_bases: u32,
    /// One captured whole-process baseline, counted once.
    pub baseline_rss_bytes: u64,
    /// Each worker's decoder, tile and Arrow encoder envelope.
    pub worker_bytes: u64,
    /// Serial Arrow decoding and final-analyzer envelope.
    pub reducer_bytes: u64,
    /// Partition descriptors, receipt bookkeeping and bounded queue envelope.
    pub metadata_queue_bytes: u64,
    /// Aggregate predicted peak, including all worker and reducer reservations.
    pub predicted_peak_rss_bytes: u64,
}

/// Results of a complete cached/parallel dataset execution.
#[derive(Debug, Clone)]
pub struct DatasetOutcome {
    /// Admitted aggregate execution plan.
    pub plan: DatasetPlan,
    /// Existing complete partitions checked and consumed.
    pub reused_partitions: usize,
    /// Newly computed and transactionally published partitions.
    pub computed_partitions: usize,
    /// Current-run computation statistics; reused rows add no record visits.
    pub stats: EvidenceRunStats,
    /// Content identity supplied by the verified input session.
    pub science_digest: String,
    /// Complete cache dataset receipt.
    pub dataset_manifest: PathBuf,
}

/// Cache integrity, execution, or compatibility failure.
#[derive(Debug, Error)]
pub enum DatasetError {
    /// Underlying evidence admission/input/analyzer failure.
    #[error(transparent)]
    Evidence(#[from] EvidenceError),
    /// Local cache I/O failure.
    #[error("dataset cache I/O failed: {0}")]
    Io(#[from] std::io::Error),
    /// Existing cached bytes or their claim failed verification.
    #[error("corrupt dataset cache: {0}")]
    Corrupt(String),
    /// Request cannot safely use the selected cache.
    #[error("incompatible dataset request: {0}")]
    Incompatible(String),
    /// A first-party worker panicked before completing its partition.
    #[error("dataset worker terminated unexpectedly")]
    WorkerPanic,
}

impl DatasetError {
    /// CLI convention: invalid requests/admission 3, integrity 5, I/O 2, breach 4.
    pub fn exit_code(&self) -> i32 {
        match self {
            Self::Corrupt(_) => 5,
            Self::Incompatible(_) => 3,
            Self::Evidence(EvidenceError::Refused { .. })
            | Self::Evidence(EvidenceError::InvalidRequest(_)) => 3,
            Self::Evidence(EvidenceError::Core(crate::core::CoreError::BudgetExceeded {
                ..
            }))
            | Self::Evidence(EvidenceError::RecordLimit(_)) => 4,
            _ => 2,
        }
    }
}

#[derive(Debug, Clone)]
struct Partition {
    contig: u32,
    start: u32,
    intervals: Vec<GenomicInterval>,
    selection: EvidenceSelection,
}

impl Partition {
    fn name(&self) -> String {
        format!("c{:08}-p{:010}", self.contig, self.start)
    }
    fn row_count(&self) -> u64 {
        self.intervals
            .iter()
            .map(|interval| u64::from(interval.end - interval.start))
            .sum()
    }
    fn selection_hash(&self) -> String {
        let mut hash = blake3::Hasher::new();
        for interval in &self.intervals {
            hash.update(&interval.contig.to_le_bytes());
            hash.update(&interval.start.to_le_bytes());
            hash.update(&interval.end.to_le_bytes());
        }
        if let EvidenceSelection::Sites(sites) = &self.selection {
            for site in sites {
                hash.update(&site.position.to_le_bytes());
                hash.update(&[site.reference]);
                hash.update(&(site.alternates.len() as u64).to_le_bytes());
                hash.update(&site.alternates);
            }
        }
        hash.finalize().to_hex().to_string()
    }
}

fn partitions(engine: &EvidenceEngine) -> Vec<Partition> {
    let mut grouped: BTreeMap<(u32, u32), Vec<GenomicInterval>> = BTreeMap::new();
    for interval in engine.intervals() {
        let mut start = interval.start;
        while start < interval.end {
            let canonical = start / CANONICAL_TILE_BASES * CANONICAL_TILE_BASES;
            let end = interval
                .end
                .min(canonical.saturating_add(CANONICAL_TILE_BASES));
            grouped
                .entry((interval.contig, canonical))
                .or_default()
                .push(GenomicInterval {
                    contig: interval.contig,
                    start,
                    end,
                });
            start = end;
        }
    }
    let mut site_groups = BTreeMap::new();
    if let EvidenceSelection::Sites(sites) = &engine.request().selection {
        for site in sites {
            site_groups
                .entry((
                    site.contig,
                    site.position / CANONICAL_TILE_BASES * CANONICAL_TILE_BASES,
                ))
                .or_insert_with(Vec::new)
                .push(site.clone());
        }
    }
    grouped
        .into_iter()
        .map(|((contig, start), intervals)| {
            let selection = site_groups
                .remove(&(contig, start))
                .map(EvidenceSelection::Sites)
                .unwrap_or_else(|| EvidenceSelection::Intervals(intervals.clone()));
            Partition {
                contig,
                start,
                intervals,
                selection,
            }
        })
        .collect()
}

/// Plan all first-party workers, queues, encoders and the serial consumer under
/// one process baseline. No cache files or output records are created.
pub fn plan_dataset(
    engine: &mut EvidenceEngine,
    options: &DatasetOptions,
    analyzer: &dyn EvidenceAnalyzer,
) -> Result<DatasetPlan, DatasetError> {
    if !(1..=MAX_WORKERS).contains(&options.workers) {
        return Err(DatasetError::Incompatible(format!(
            "workers must be in 1..={MAX_WORKERS}"
        )));
    }
    // Count canonical ownership without allocating descriptors before admission.
    let mut partition_count = 0usize;
    let mut last = None;
    for interval in engine.intervals() {
        if interval.start == interval.end {
            continue;
        }
        let first = interval.start / CANONICAL_TILE_BASES;
        let final_tile = (interval.end - 1) / CANONICAL_TILE_BASES;
        let shared = usize::from(last == Some((interval.contig, first)));
        partition_count = partition_count
            .checked_add((final_tile - first + 1) as usize - shared)
            .ok_or_else(|| DatasetError::Incompatible("too many partitions".into()))?;
        last = Some((interval.contig, final_tile));
    }
    let workers = options.workers.min(partition_count.max(1));
    engine.plan_for_analyzer(analyzer)?;
    let execution = &engine.request().execution;
    let final_analyzer = match analyzer.requirements().retained_bytes {
        Some(bytes) => bytes.max(execution.analyzer_bytes),
        None if execution.memory_budget_bytes.is_some() => {
            return Err(DatasetError::Incompatible(
                "budgeted reduction requires a declared analyzer bound".into(),
            ))
        }
        None => execution.analyzer_bytes,
    };
    let arrow_bytes = EvidenceArrowWriter::with_fields(std::io::sink(), engine.request().fields)
        .additional_memory_bytes()
        .expect("first-party encoder bound");
    let bytes_per_locus = engine.plan().bytes_per_locus;
    let decoder_bytes = engine
        .plan()
        .fixed_bytes
        .saturating_sub(engine.plan().analyzer_bytes)
        .saturating_sub(engine.plan().selection_bytes);
    let worker_fixed = decoder_bytes.saturating_add(arrow_bytes);
    let reducer_bytes =
        evidence_reader_memory_bytes(engine.request().fields).saturating_add(final_analyzer);
    let site_count = match &engine.request().selection {
        EvidenceSelection::Sites(sites) => sites.len(),
        _ => 0,
    };
    // Descriptors and the eventual top-level receipt are O(partitions), while
    // each queued selection owns at most one canonical tile's site annotations.
    let metadata_queue_bytes = (partition_count as u64)
        .saturating_mul(4096)
        .saturating_add(engine.plan().selection_bytes)
        .saturating_add((site_count as u64).saturating_mul(256))
        .saturating_add((engine.intervals().len() as u64).saturating_mul(256))
        .saturating_add((workers as u64).saturating_mul(u64::from(CANONICAL_TILE_BASES) * 512));
    let baseline = engine.plan().baseline_rss_bytes;
    let fixed = baseline
        .saturating_add(worker_fixed.saturating_mul(workers as u64))
        .saturating_add(reducer_bytes)
        .saturating_add(metadata_queue_bytes);
    let maximum = execution.max_microtile_bases.clamp(1, CANONICAL_TILE_BASES);
    let microtile = if let Some(budget) = execution.memory_budget_bytes {
        let minimum = fixed.saturating_add(bytes_per_locus.saturating_mul(workers as u64));
        if minimum > budget {
            return Err(EvidenceError::Refused {
                needed: minimum,
                budget,
            }
            .into());
        }
        u64::from(maximum).min((budget - fixed) / (bytes_per_locus * workers as u64)) as u32
    } else {
        maximum
    };
    let worker_bytes =
        worker_fixed.saturating_add(bytes_per_locus.saturating_mul(u64::from(microtile)));
    Ok(DatasetPlan {
        partition_count,
        worker_count: workers,
        microtile_bases: microtile,
        baseline_rss_bytes: baseline,
        worker_bytes,
        reducer_bytes,
        metadata_queue_bytes,
        predicted_peak_rss_bytes: fixed
            .saturating_add(bytes_per_locus * u64::from(microtile) * workers as u64),
    })
}

/// Compute or resume exact first-party Arrow partitions, then feed the final
/// analyzer in canonical order. `science_digest` is the caller's verified input
/// content/producer namespace. Request semantics are additionally bound and
/// validated internally, so reusing that namespace with different filters or
/// selection is refused. Inputs must be immutable from engine-open to completion.
/// Completed datasets can be reused without `resume`; incomplete ones require it.
pub fn run_dataset(
    engine: &mut EvidenceEngine,
    science_digest: &str,
    options: &DatasetOptions,
    analyzer: &mut dyn EvidenceAnalyzer,
) -> Result<DatasetOutcome, DatasetError> {
    let inputs = InputSnapshot::for_engine(engine)?;
    run_dataset_with_snapshot(engine, science_digest, options, analyzer, &inputs)
}

/// Execute using the same snapshot captured before the caller hashed and opened
/// its inputs. Every worker checks it before publishing a partition; reduction
/// checks again before publishing the complete dataset receipt.
pub fn run_dataset_with_snapshot(
    engine: &mut EvidenceEngine,
    science_digest: &str,
    options: &DatasetOptions,
    analyzer: &mut dyn EvidenceAnalyzer,
    inputs: &InputSnapshot,
) -> Result<DatasetOutcome, DatasetError> {
    inputs.verify()?;
    let request_digest = request_digest(engine);
    let fields = engine.request().fields;
    if science_digest.len() != 64 || !science_digest.bytes().all(|byte| byte.is_ascii_hexdigit()) {
        return Err(DatasetError::Incompatible(
            "science digest must be a 64-character BLAKE3 hex string".into(),
        ));
    }
    let plan = plan_dataset(engine, options, analyzer)?;
    let parts = partitions(engine);
    let root = options.cache_dir.join(science_digest);
    let dataset_manifest = root.join(MANIFEST_NAME);
    if root.exists()
        && !dataset_manifest.is_file()
        && !options.resume
        && fs::read_dir(&root)?.next().is_some()
    {
        return Err(DatasetError::Incompatible(
            "incomplete cached dataset; pass --resume to reuse its complete partitions".into(),
        ));
    }
    let completed = if dataset_manifest.exists() {
        let manifest = read_receipt(&dataset_manifest)?;
        if manifest.params.get("dataset.request_blake3") != Some(&request_digest) {
            return Err(DatasetError::Incompatible("cache namespace was used with different request semantics; choose a new verified input/science key".into()));
        }
        if manifest
            .params
            .get("evidence.science_blake3")
            .map(String::as_str)
            != Some(science_digest)
            || manifest
                .params
                .get("dataset.partition_count")
                .and_then(|value| value.parse::<usize>().ok())
                != Some(parts.len())
            || manifest.inputs.len() != parts.len()
            || manifest.outputs.len() != parts.len()
            || manifest.params.get("run_status").map(String::as_str) != Some("completed")
            || manifest
                .params
                .get("evidence.schema")
                .and_then(|value| value.parse::<u32>().ok())
                != Some(fields.schema_version())
            || receipt_fields(&manifest)? != fields
            || manifest.params.get("pileup.semantics").map(String::as_str)
                != Some("exact-or-fail-v1")
        {
            return Err(DatasetError::Corrupt(
                "complete dataset identity or partition count differs".into(),
            ));
        }
        Some(manifest)
    } else {
        None
    };
    fs::create_dir_all(&root)?;
    let mut missing = Vec::new();
    for (index, part) in parts.iter().enumerate() {
        let path = root.join(part.name());
        if !path.exists() {
            if completed.is_some() {
                return Err(DatasetError::Corrupt(format!(
                    "complete dataset is missing partition {}",
                    part.name()
                )));
            }
            missing.push(index);
        } else {
            let receipt = verify_partition(&path, part, science_digest, &request_digest, fields)?;
            if let Some(complete) = &completed {
                if blake3_file(&path.join("manifest.json"))? != complete.inputs[index].blake3
                    || receipt.outputs[0].blake3 != complete.outputs[index].blake3
                {
                    return Err(DatasetError::Corrupt(format!(
                        "partition differs from completed dataset: {}",
                        part.name()
                    )));
                }
            }
        }
    }
    let reused_partitions = parts.len() - missing.len();
    let computed_partitions = missing.len();
    let factory = engine.worker_factory();
    let mut worker_execution: EvidenceExecution = engine.request().execution.clone();
    worker_execution.memory_budget_bytes = None;
    worker_execution.max_microtile_bases = plan.microtile_bases;
    worker_execution.analyzer_bytes =
        EvidenceArrowWriter::with_fields(std::io::sink(), engine.request().fields)
            .additional_memory_bytes()
            .unwrap();
    let mut stats = EvidenceRunStats::default();
    if !missing.is_empty() {
        let cancelled = Arc::new(AtomicBool::new(false));
        let (sender, receiver) = mpsc::sync_channel::<usize>(plan.worker_count);
        let receiver = Arc::new(Mutex::new(receiver));
        let results = std::thread::scope(|scope| {
            let mut handles = Vec::new();
            for _ in 0..plan.worker_count {
                let receiver = Arc::clone(&receiver);
                let cancelled = Arc::clone(&cancelled);
                let factory = factory.clone();
                let execution = worker_execution.clone();
                let root = &root;
                let parts = &parts;
                let request_digest = &request_digest;
                handles.push(scope.spawn(move || {
                    let result = (|| {
                        let mut total = EvidenceRunStats::default();
                        let mut worker = None;
                        loop {
                            if cancelled.load(Ordering::Acquire) {
                                break;
                            }
                            let job = receiver
                                .lock()
                                .map_err(|_| DatasetError::WorkerPanic)?
                                .recv();
                            let Ok(index) = job else {
                                break;
                            };
                            let part = &parts[index];
                            crate::core::governor::checkpoint().map_err(EvidenceError::from)?;
                            if worker.is_none() {
                                worker = Some(factory.open_with_execution(
                                    part.selection.clone(),
                                    execution.clone(),
                                )?);
                            }
                            let engine = worker.as_mut().unwrap();
                            engine.set_selection(part.selection.clone())?;
                            let computed = write_partition(
                                engine,
                                root,
                                part,
                                science_digest,
                                request_digest,
                                inputs,
                            )?;
                            add_stats(&mut total, &computed)?;
                        }
                        Ok::<_, DatasetError>(total)
                    })();
                    if result.is_err() {
                        cancelled.store(true, Ordering::Release);
                    }
                    result
                }));
            }
            // The submitting thread must not retain a receiver: if every worker
            // fails while send() is blocked on a full queue, disconnect wakes it.
            drop(receiver);
            for index in missing {
                if cancelled.load(Ordering::Acquire) || sender.send(index).is_err() {
                    break;
                }
            }
            drop(sender);
            handles
                .into_iter()
                .map(|handle| handle.join().unwrap_or(Err(DatasetError::WorkerPanic)))
                .collect::<Vec<_>>()
        });
        for result in results {
            add_stats(&mut stats, &result?)?;
        }
    }
    let mut dataset = RunManifest::new("dataset evidence");
    dataset.params.extend([
        ("evidence.science_blake3".into(), science_digest.into()),
        ("dataset.request_blake3".into(), request_digest.clone()),
        (
            "evidence.schema".into(),
            engine.request().fields.schema_version().to_string(),
        ),
        (
            "evidence.fields".into(),
            engine.request().fields.bits().to_string(),
        ),
        (
            "evidence.fields_version".into(),
            engine.request().fields.mask_version().to_string(),
        ),
        ("pileup.semantics".into(), "exact-or-fail-v1".into()),
        ("dataset.partition_count".into(), parts.len().to_string()),
        (
            "dataset.partition_bases".into(),
            CANONICAL_TILE_BASES.to_string(),
        ),
        ("run_status".into(), "completed".into()),
    ]);
    let mut emitted = 0u64;
    for part in &parts {
        crate::core::governor::checkpoint().map_err(EvidenceError::from)?;
        let path = root.join(part.name());
        let receipt = verify_partition(&path, part, science_digest, &request_digest, fields)?;
        let mut interval_index = 0usize;
        let mut next_position = part
            .intervals
            .first()
            .map(|interval| interval.start)
            .unwrap_or(0);
        let mut consumer_failed = false;
        let read_result = read_evidence_batches_expected_fields(
            BufReader::new(File::open(path.join("evidence.arrow"))?),
            engine.contigs(),
            fields,
            |batch| {
                if batch.fields() != fields {
                    return Err(EvidenceError::InvalidInput(
                        "cached batch field mask differs from request".into(),
                    ));
                }
                if batch.contig_id != part.contig {
                    return Err(EvidenceError::InvalidInput(
                        "cached batch contig differs from ownership".into(),
                    ));
                }
                for row in batch.rows() {
                    let interval = part.intervals.get(interval_index).ok_or_else(|| {
                        EvidenceError::InvalidInput("cached partition contains extra loci".into())
                    })?;
                    if row.position != next_position {
                        return Err(EvidenceError::InvalidInput(
                            "cached partition omits or duplicates a selected locus".into(),
                        ));
                    }
                    next_position += 1;
                    if next_position == interval.end {
                        interval_index += 1;
                        if let Some(next) = part.intervals.get(interval_index) {
                            next_position = next.start;
                        }
                    }
                }
                emitted = emitted
                    .checked_add(batch.len() as u64)
                    .ok_or(EvidenceError::CounterOverflow)?;
                let result = analyzer.on_batch(batch);
                consumer_failed = result.is_err();
                result
            },
        );
        let metadata = read_result.map_err(|error| {
            if consumer_failed {
                error.into()
            } else {
                DatasetError::Corrupt(format!("invalid cached Arrow evidence: {error}"))
            }
        })?;
        if metadata.fields != fields || metadata.schema_version != fields.schema_version() {
            return Err(DatasetError::Corrupt(
                "cached Arrow metadata differs from requested fields/schema".into(),
            ));
        }
        if interval_index != part.intervals.len() {
            return Err(DatasetError::Corrupt(
                "cached partition is missing selected loci".into(),
            ));
        }
        dataset.inputs.push(FileHash {
            path: path.join("manifest.json").display().to_string(),
            blake3: blake3_file(&path.join("manifest.json"))?,
        });
        dataset.outputs.push(FileHash {
            path: path.join("evidence.arrow").display().to_string(),
            blake3: receipt.outputs[0].blake3.clone(),
        });
    }
    analyzer.finish()?;
    inputs.verify()?;
    stats.emitted_loci = emitted;
    stats.peak_rss_bytes = crate::util::rss::peak_rss_bytes();
    dataset
        .params
        .insert("outcome.rows".into(), emitted.to_string());
    dataset.finalize();
    inputs.verify()?;
    write_atomic(
        &dataset_manifest,
        dataset.to_canonical_json().as_bytes(),
        true,
    )?;
    Ok(DatasetOutcome {
        plan,
        reused_partitions,
        computed_partitions,
        stats,
        science_digest: science_digest.into(),
        dataset_manifest,
    })
}

fn add_stats(total: &mut EvidenceRunStats, part: &EvidenceRunStats) -> Result<(), DatasetError> {
    for (sum, value) in [
        (&mut total.emitted_loci, part.emitted_loci),
        (&mut total.record_visits, part.record_visits),
        (
            &mut total.filtered_record_visits,
            part.filtered_record_visits,
        ),
        (
            &mut total.sample_filtered_record_visits,
            part.sample_filtered_record_visits,
        ),
        (&mut total.microtiles, part.microtiles),
    ] {
        *sum = sum
            .checked_add(value)
            .ok_or(EvidenceError::CounterOverflow)?;
    }
    total.max_record_bytes = total.max_record_bytes.max(part.max_record_bytes);
    total.max_read_length = total.max_read_length.max(part.max_read_length);
    total.peak_rss_bytes = total.peak_rss_bytes.max(part.peak_rss_bytes);
    Ok(())
}

struct StagingDirectory(PathBuf);
impl Drop for StagingDirectory {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn write_partition(
    engine: &mut EvidenceEngine,
    root: &Path,
    part: &Partition,
    science: &str,
    request_digest: &str,
    inputs: &InputSnapshot,
) -> Result<EvidenceRunStats, DatasetError> {
    let destination = root.join(part.name());
    let temporary = StagingDirectory(root.join(format!(
        ".{}-{}-{}",
        part.name(),
        std::process::id(),
        TEMP_SEQUENCE.fetch_add(1, Ordering::Relaxed)
    )));
    fs::create_dir(&temporary.0)?;
    let output = temporary.0.join("evidence.arrow");
    let mut writer = EvidenceArrowWriter::with_fields(
        BufWriter::new(File::create(&output)?),
        engine.request().fields,
    );
    let stats = engine.run(&mut writer)?;
    let mut file = writer.into_inner()?;
    file.flush()?;
    file.get_ref().sync_all()?;
    drop(file);
    if stats.emitted_loci != part.row_count() {
        return Err(DatasetError::Corrupt(
            "worker emitted incomplete ownership span".into(),
        ));
    }
    let mut receipt = RunManifest::new("dataset partition");
    receipt.params.extend([
        ("evidence.science_blake3".into(), science.into()),
        ("dataset.request_blake3".into(), request_digest.into()),
        (
            "evidence.schema".into(),
            engine.request().fields.schema_version().to_string(),
        ),
        (
            "evidence.fields".into(),
            engine.request().fields.bits().to_string(),
        ),
        (
            "evidence.fields_version".into(),
            engine.request().fields.mask_version().to_string(),
        ),
        ("pileup.semantics".into(), "exact-or-fail-v1".into()),
        ("dataset.contig".into(), part.contig.to_string()),
        ("dataset.partition_start".into(), part.start.to_string()),
        (
            "dataset.partition_bases".into(),
            CANONICAL_TILE_BASES.to_string(),
        ),
        ("dataset.selection_blake3".into(), part.selection_hash()),
        ("outcome.rows".into(), stats.emitted_loci.to_string()),
        ("run_status".into(), "completed".into()),
    ]);
    receipt.outputs.push(FileHash {
        path: destination.join("evidence.arrow").display().to_string(),
        blake3: blake3_file(&output)?,
    });
    receipt.finalize();
    write_atomic(
        &temporary.0.join("manifest.json"),
        receipt.to_canonical_json().as_bytes(),
        false,
    )?;
    if destination.exists() {
        return Err(DatasetError::Incompatible(
            "a concurrent writer published this partition; retry with --resume".into(),
        ));
    }
    inputs.verify()?;
    fs::rename(&temporary.0, &destination)?;
    Ok(stats)
}

fn read_receipt(path: &Path) -> Result<RunManifest, DatasetError> {
    let receipt = RunManifest::from_canonical_json(&fs::read_to_string(path)?)
        .map_err(|error| DatasetError::Corrupt(error.to_string()))?;
    if receipt.self_hash_ok() != Some(true) || receipt.measurement_hash_ok() == Some(false) {
        return Err(DatasetError::Corrupt(format!(
            "receipt integrity failed: {}",
            path.display()
        )));
    }
    Ok(receipt)
}

fn receipt_fields(receipt: &RunManifest) -> Result<EvidenceFields, DatasetError> {
    let schema = receipt.params.get("evidence.schema").map(String::as_str);
    let fields = match receipt.params.get("evidence.fields") {
        None if schema == Some("1") => EvidenceFields::ALL,
        Some(bits) => bits
            .parse::<u32>()
            .ok()
            .and_then(|bits| EvidenceFields::from_bits(bits).ok())
            .ok_or_else(|| DatasetError::Corrupt("invalid evidence field mask".into()))?,
        None => {
            return Err(DatasetError::Corrupt(
                "missing projected evidence field mask".into(),
            ))
        }
    };
    match receipt.params.get("evidence.fields_version") {
        None if schema == Some("1") => (),
        Some(version) if version == &fields.mask_version().to_string() => (),
        _ => {
            return Err(DatasetError::Corrupt(
                "unsupported evidence field mask version".into(),
            ))
        }
    }
    if schema != Some(fields.schema_version().to_string().as_str()) {
        return Err(DatasetError::Corrupt(
            "inconsistent evidence schema and field mask".into(),
        ));
    }
    Ok(fields)
}

fn verify_partition(
    path: &Path,
    part: &Partition,
    science: &str,
    request_digest: &str,
    fields: EvidenceFields,
) -> Result<RunManifest, DatasetError> {
    if !path.join("manifest.json").is_file() || !path.join("evidence.arrow").is_file() {
        return Err(DatasetError::Corrupt(format!(
            "partition is incomplete: {}",
            path.display()
        )));
    }
    let receipt = read_receipt(&path.join("manifest.json"))?;
    if receipt
        .params
        .get("dataset.request_blake3")
        .map(String::as_str)
        != Some(request_digest)
    {
        return Err(DatasetError::Incompatible(
            "cache partition was produced with different request semantics".into(),
        ));
    }
    if receipt_fields(&receipt)? != fields {
        return Err(DatasetError::Corrupt(
            "partition field mask differs from request".into(),
        ));
    }
    let expected = [
        ("evidence.science_blake3", science.to_string()),
        ("evidence.schema", fields.schema_version().to_string()),
        ("pileup.semantics", "exact-or-fail-v1".into()),
        ("dataset.contig", part.contig.to_string()),
        ("dataset.partition_start", part.start.to_string()),
        ("dataset.partition_bases", CANONICAL_TILE_BASES.to_string()),
        ("dataset.selection_blake3", part.selection_hash()),
        ("outcome.rows", part.row_count().to_string()),
        ("run_status", "completed".into()),
    ];
    if expected
        .iter()
        .any(|(key, value)| receipt.params.get(*key) != Some(value))
        || receipt.outputs.len() != 1
    {
        return Err(DatasetError::Corrupt(format!(
            "partition claim differs: {}",
            path.display()
        )));
    }
    if blake3_file(&path.join("evidence.arrow"))? != receipt.outputs[0].blake3 {
        return Err(DatasetError::Corrupt(format!(
            "Arrow artifact hash differs: {}",
            path.display()
        )));
    }
    Ok(receipt)
}
