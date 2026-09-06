//! Offline, bounded queries over a portable, integrity-checked evidence dataset.

use super::*;
use crate::core::ContigSet;
use crate::evidence::{
    EvidenceBatch, EvidenceProfile, EvidenceSampleScope, SnvSite, EVIDENCE_ARROW_BATCH_ROWS,
};
use std::io::Read;

/// Envelopes for untrusted portable metadata. Original BAM/reference files are
/// never opened: their recorded identities remain provenance claims.
#[derive(Debug, Clone, Copy)]
pub struct DatasetReadLimits {
    /// Largest accepted parent receipt, before parsing.
    pub max_manifest_bytes: usize,
    /// Largest accepted descriptor, before parsing.
    pub max_descriptor_bytes: usize,
    /// Largest accepted partition receipt, before parsing.
    pub max_partition_receipt_bytes: usize,
    /// Optional whole-process admission budget during metadata opening.
    pub memory_budget_bytes: Option<u64>,
}
impl Default for DatasetReadLimits {
    fn default() -> Self {
        Self {
            max_manifest_bytes: 32 << 20,
            max_descriptor_bytes: 32 << 20,
            max_partition_receipt_bytes: 64 << 10,
            memory_budget_bytes: None,
        }
    }
}

/// Scientific selection and physical output groups for one offline query.
#[derive(Debug, Clone)]
pub struct DatasetQuery {
    /// WholeGenome means the full dictionary; use dataset.selection() for all stored loci.
    pub selection: EvidenceSelection,
    /// Must be a subset of the groups physically present in the dataset.
    pub fields: EvidenceFields,
}

/// Canonical normalized scientific query, including exact selected coordinates,
/// SNV REF/ALT annotations and physical fields. Parent compatibility must also
/// participate in the identity of a derived scientific artifact.
pub fn canonical_dataset_query(
    query: &DatasetQuery,
    contigs: &ContigSet,
) -> Result<String, DatasetError> {
    let (intervals, sites) = query.selection.normalize(contigs)?;
    let selection = DescriptorSelection {
        intervals: intervals
            .iter()
            .map(|i| DescriptorInterval {
                contig: i.contig,
                start: i.start,
                end: i.end,
            })
            .collect(),
        sites: matches!(query.selection, EvidenceSelection::Sites(_)).then(|| {
            sites
                .values()
                .map(|site| DescriptorSite {
                    contig: site.contig,
                    position: site.position,
                    reference: site.reference,
                    alternates: site.alternates.clone(),
                })
                .collect()
        }),
    };
    serde_json::to_string(
        &serde_json::json!({"version":1,"fields":query.fields.bits(),"selection":selection}),
    )
    .map_err(|error| {
        DatasetError::Incompatible(format!("cannot serialize normalized query: {error}"))
    })
}

/// Hash the canonical query alone. Combine with source compatibility and analyzer
/// semantics when identifying a derived scientific result.
pub fn dataset_query_digest(
    query: &DatasetQuery,
    contigs: &ContigSet,
) -> Result<String, DatasetError> {
    let mut hash = blake3::Hasher::new();
    hash.update(b"rosalind-persisted-query-v1\0");
    hash.update(canonical_dataset_query(query, contigs)?.as_bytes());
    Ok(hash.finalize().to_hex().to_string())
}

/// Coverage membership, independent of query ALT annotations. Missing loci must
/// be computed from alignments or rejected; they never become zero-depth rows.
#[derive(Debug, Clone)]
pub struct DatasetCoverage {
    /// Requested loci represented by stored rows, retaining query REF/ALT.
    pub covered: EvidenceSelection,
    /// Requested loci absent from the dataset, retaining query REF/ALT.
    pub missing: EvidenceSelection,
    /// Number of represented unique requested loci.
    pub covered_loci: u64,
    /// Number of absent unique requested loci.
    pub missing_loci: u64,
}

/// Conservative plan for serial persisted reads. Physical source IPC buffers
/// retain the source field set even when output requests a smaller projection.
#[derive(Debug, Clone)]
pub struct DatasetReadPlan {
    /// Stable working-set model identifier.
    pub model_id: &'static str,
    /// Process RSS captured before opening portable metadata, counted once.
    pub baseline_rss_bytes: u64,
    /// Retained metadata, query normalization and receipt parsing reservation.
    pub metadata_bytes: u64,
    /// Bounded decoder workspace for the source field mask.
    pub source_decoder_bytes: u64,
    /// One projected canonical output batch of at most 1,024 rows.
    pub projection_bytes: u64,
    /// Consumer-retained state plus the explicitly requested additional reserve.
    pub analyzer_bytes: u64,
    /// Predicted process high-water memory for a successful read.
    pub predicted_peak_rss_bytes: u64,
    /// Unique selected positions; all must be represented.
    pub selected_loci: u64,
    /// Physical groups decoded from source partitions.
    pub source_fields: EvidenceFields,
    /// Physical groups passed to the consumer.
    pub output_fields: EvidenceFields,
    /// Fixed output row envelope, independent of source batch boundaries.
    pub output_batch_rows: usize,
}

/// A relocated dataset reader. The descriptor and parent are checked at open;
/// each consumed partition is hash-checked before decoding and its exact stored
/// coordinate/REF/ALT contract is checked while streaming. Failed callbacks may
/// have consumed earlier batches, so artifact consumers should stage outputs.
/// Files must remain immutable for this session. Metadata guards detect ordinary
/// replacement or mutation, not deliberate restoration of file metadata.
#[derive(Debug)]
pub struct VerifiedEvidenceDataset {
    root: PathBuf,
    descriptor: DatasetDescriptor,
    contigs: ContigSet,
    fields: EvidenceFields,
    limits: DatasetReadLimits,
    baseline: u64,
    metadata_bytes: u64,
    identity: Vec<FileHash>,
    snapshot: InputSnapshot,
    observed: BTreeMap<usize, InputSnapshot>,
}

impl VerifiedEvidenceDataset {
    /// Open only the portable parent and descriptor. Payloads are verified on
    /// demand; opening does not imply every partition body has been decoded.
    pub fn open(
        manifest: impl AsRef<Path>,
        limits: DatasetReadLimits,
    ) -> Result<Self, DatasetError> {
        if limits.max_manifest_bytes == 0
            || limits.max_descriptor_bytes == 0
            || limits.max_partition_receipt_bytes == 0
        {
            return Err(DatasetError::Incompatible(
                "metadata envelopes must be positive".into(),
            ));
        }
        let baseline = crate::util::rss::peak_rss_bytes();
        let manifest = fs::canonicalize(manifest)?;
        let root = manifest
            .parent()
            .ok_or_else(|| DatasetError::Incompatible("portable manifest has no directory".into()))?
            .to_path_buf();
        let descriptor_path = local_path(&root, DATASET_DESCRIPTOR_NAME)?;
        let snapshot = InputSnapshot::capture([manifest.clone(), descriptor_path.clone()])?;
        let manifest_len = limited_size(&manifest, limits.max_manifest_bytes)?;
        let descriptor_len = limited_size(&descriptor_path, limits.max_descriptor_bytes)?;
        // Bound decoding, validation clones, maps, and canonical reserialization.
        let metadata_bytes = (manifest_len as u64)
            .saturating_add(descriptor_len as u64)
            .saturating_mul(16)
            .saturating_add(64 << 10);
        admit(
            baseline.saturating_add(metadata_bytes),
            limits.memory_budget_bytes,
        )?;
        let manifest_bytes = read_limited(&manifest, manifest_len)?;
        let parent = parse_receipt(&manifest_bytes)?;
        let descriptor_bytes = read_limited(&descriptor_path, descriptor_len)?;
        let descriptor_hash = blake3::hash(&descriptor_bytes).to_hex().to_string();
        let descriptor = DatasetDescriptor::from_bytes(
            &descriptor_bytes,
            DescriptorLimits {
                max_bytes: limits.max_descriptor_bytes,
            },
        )?;
        let metadata_bytes = descriptor
            .partitions
            .iter()
            .fold(metadata_bytes, |bytes, part| {
                bytes
                    .saturating_add((root.as_os_str().len() as u64).saturating_mul(2))
                    .saturating_add(
                        part.arrow.path.len() as u64 + part.receipt.path.len() as u64 + 1024,
                    )
            });
        admit(
            baseline.saturating_add(metadata_bytes),
            limits.memory_budget_bytes,
        )?;
        validate_parent(&parent, &descriptor, &descriptor_hash)?;
        let contigs = descriptor.contig_set();
        let fields = EvidenceFields::from_bits(descriptor.fields)?;
        let observed = BTreeMap::new();
        let identity = vec![
            FileHash {
                path: manifest.display().to_string(),
                blake3: blake3::hash(&manifest_bytes).to_hex().to_string(),
            },
            FileHash {
                path: descriptor_path.display().to_string(),
                blake3: descriptor_hash,
            },
        ];
        snapshot.verify()?;
        crate::core::governor::checkpoint().map_err(EvidenceError::from)?;
        Ok(Self {
            root,
            descriptor,
            contigs,
            fields,
            limits,
            baseline,
            metadata_bytes,
            identity,
            snapshot,
            observed,
        })
    }

    /// Verified descriptor, including original source identities as claims only.
    pub fn descriptor(&self) -> &DatasetDescriptor {
        &self.descriptor
    }
    /// Ordered coordinate dictionary; no alignment or reference is reopened.
    pub fn contigs(&self) -> &ContigSet {
        &self.contigs
    }
    /// Actual local parent/descriptor hashes for derived artifact lineage.
    pub fn source_hashes(&self) -> &[FileHash] {
        &self.identity
    }
    /// Locally verified consumed partition receipts and Arrow artifacts for a
    /// derived receipt. Unvisited partitions are not falsely marked verified.
    pub fn verified_partition_hashes(&self) -> Result<Vec<FileHash>, DatasetError> {
        self.verify_unchanged()?;
        let mut files = Vec::new();
        for index in self.observed.keys() {
            let part = &self.descriptor.partitions[*index];
            for artifact in [&part.receipt, &part.arrow] {
                files.push(FileHash {
                    path: self.root.join(&artifact.path).display().to_string(),
                    blake3: artifact.blake3.clone(),
                });
            }
        }
        Ok(files)
    }
    /// All stored loci, including their original selection annotations.
    pub fn selection(&self) -> EvidenceSelection {
        self.descriptor.selection.to_evidence_selection()
    }
    /// Physically available groups.
    pub fn fields(&self) -> EvidenceFields {
        self.fields
    }
    pub(super) fn execution_budget(&self, execution: &EvidenceExecution) -> Option<u64> {
        execution
            .memory_budget_bytes
            .or(self.limits.memory_budget_bytes)
    }
    pub(crate) fn memory_budget_bytes(&self) -> Option<u64> {
        self.limits.memory_budget_bytes
    }

    /// Profile changes require raw observations and cannot be reconstructed from
    /// these summaries. Reject even an apparently stricter threshold.
    pub fn require_profile(&self, profile: &EvidenceProfile) -> Result<(), DatasetError> {
        if self.descriptor.profile.to_profile() != *profile {
            return Err(DatasetError::Incompatible(
                "stored evidence uses a different filter profile; recompute from alignments".into(),
            ));
        }
        Ok(())
    }
    /// Sample scopes are exact scientific identities, never a post-hoc row filter.
    pub fn require_sample_scope(&self, scope: &EvidenceSampleScope) -> Result<(), DatasetError> {
        if self.descriptor.sample_scope.canonical_json()? != scope.canonical_json() {
            return Err(DatasetError::Incompatible(
                "stored evidence uses a different sample scope; recompute from alignments".into(),
            ));
        }
        Ok(())
    }

    /// Split a normalized query by stored coverage without reading Arrow. This
    /// checks membership only; requested SNV REF is validated against stored rows
    /// when visited, and its ALT set need not match the source selection's ALT.
    pub fn coverage_split(
        &self,
        selection: &EvidenceSelection,
    ) -> Result<DatasetCoverage, DatasetError> {
        let (intervals, sites) = selection.normalize(&self.contigs)?;
        let (covered, missing) = split_intervals(&intervals, &self.descriptor.selection.intervals);
        let covered_loci = interval_count(&covered);
        let missing_loci = interval_count(&missing);
        let select = |intervals: Vec<GenomicInterval>| {
            if matches!(selection, EvidenceSelection::Sites(_)) {
                EvidenceSelection::Sites(
                    sites
                        .values()
                        .filter(|site| contains(&intervals, site.contig, site.position))
                        .cloned()
                        .collect(),
                )
            } else {
                EvidenceSelection::Intervals(intervals)
            }
        };
        Ok(DatasetCoverage {
            covered: select(covered),
            missing: select(missing),
            covered_loci,
            missing_loci,
        })
    }

    /// Admit the source decoder, projection, metadata and declared consumer state.
    /// Only execution.memory_budget_bytes/analyzer_bytes affect offline reads;
    /// BAM record/read envelopes and alignment microtiles do not apply to Arrow.
    pub fn plan(
        &self,
        query: &DatasetQuery,
        analyzer: &dyn EvidenceAnalyzer,
        execution: &EvidenceExecution,
    ) -> Result<DatasetReadPlan, DatasetError> {
        if !self.fields.contains(query.fields) {
            return Err(DatasetError::Incompatible(
                "query requires fields absent from the stored dataset".into(),
            ));
        }
        let requirements = analyzer.requirements();
        if !query.fields.contains(requirements.fields) {
            return Err(DatasetError::Incompatible(
                "query omits fields required by its analyzer".into(),
            ));
        }
        if requirements.context_bases != 0
            || (requirements.requires_reference && !self.descriptor.has_reference)
        {
            return Err(DatasetError::Incompatible(
                "stored evidence does not provide the requested reference context".into(),
            ));
        }
        let budget = self.execution_budget(execution);
        let analyzer_bytes = match requirements.retained_bytes {
            Some(bytes) => bytes.max(execution.analyzer_bytes),
            None if budget.is_some() => {
                return Err(DatasetError::Incompatible(
                    "budgeted persisted analyzer requires a declared memory bound".into(),
                ))
            }
            None => execution.analyzer_bytes,
        };
        let query_bytes = query_memory_bytes(&query.selection, self.contigs.len());
        let metadata_bytes = self
            .metadata_bytes
            .saturating_add(query_bytes)
            .saturating_add((self.limits.max_partition_receipt_bytes as u64).saturating_mul(16))
            .saturating_add((self.descriptor.partitions.len() as u64).saturating_mul(8));
        let source_decoder_bytes = evidence_reader_memory_bytes(self.fields);
        let max_name = self
            .contigs
            .iter()
            .map(|c| c.name.len() as u64)
            .max()
            .unwrap_or(0);
        let projection_bytes = (query.fields.storage_bytes_per_locus() + 40)
            .saturating_mul(EVIDENCE_ARROW_BATCH_ROWS as u64)
            .saturating_add(max_name);
        let predicted_peak_rss_bytes = self
            .baseline
            .saturating_add(metadata_bytes)
            .saturating_add(source_decoder_bytes)
            .saturating_add(projection_bytes)
            .saturating_add(analyzer_bytes);
        admit(predicted_peak_rss_bytes, budget)?;
        let split = self.coverage_split(&query.selection)?;
        if split.missing_loci != 0 {
            return Err(DatasetError::Incompatible(format!("query includes {} loci absent from the dataset; use coverage_split and compute missing loci from alignments", split.missing_loci)));
        }
        Ok(DatasetReadPlan {
            model_id: "persisted-evidence-reader-v1",
            baseline_rss_bytes: self.baseline,
            metadata_bytes,
            source_decoder_bytes,
            projection_bytes,
            analyzer_bytes,
            predicted_peak_rss_bytes,
            selected_loci: split.covered_loci,
            source_fields: self.fields,
            output_fields: query.fields,
            output_batch_rows: EVIDENCE_ARROW_BATCH_ROWS,
        })
    }

    /// Stream complete selected evidence in dictionary order. Each touched source
    /// partition is fully validated, even when only one of its rows is requested.
    /// Output batches never cross canonical ownership and contain at most 1,024
    /// rows. Fresh query ALT annotations replace the source annotations.
    pub fn visit_batches(
        &mut self,
        query: &DatasetQuery,
        analyzer: &mut dyn EvidenceAnalyzer,
        execution: &EvidenceExecution,
    ) -> Result<EvidenceRunStats, DatasetError> {
        let plan = self.plan(query, analyzer, execution)?;
        runtime_check(self.execution_budget(execution))?;
        self.snapshot.verify()?;
        let (intervals, sites) = query.selection.normalize(&self.contigs)?;
        let mut stats = EvidenceRunStats::default();
        let mut selected = Vec::new();
        for interval in &intervals {
            let first = self.descriptor.partitions.partition_point(|part| {
                (part.contig, part.start.saturating_add(CANONICAL_TILE_BASES))
                    <= (interval.contig, interval.start)
            });
            let end = self.descriptor.partitions.partition_point(|part| {
                (part.contig, part.start) < (interval.contig, interval.end)
            });
            for index in first..end {
                if selected.last() != Some(&index) {
                    selected.push(index);
                }
            }
        }
        let budget = self.execution_budget(execution);
        for &index in &selected {
            let part = &self.descriptor.partitions[index];
            if !overlaps(
                &intervals,
                part.contig,
                part.start,
                part.start.saturating_add(CANONICAL_TILE_BASES),
            ) {
                continue;
            }
            self.snapshot.verify()?;
            if let Some(observed) = self.observed.get(&index) {
                observed.verify()?;
            }
            let receipt_path = local_path(&self.root, &part.receipt.path)?;
            let arrow_path = local_path(&self.root, &part.arrow.path)?;
            let snapshot = InputSnapshot::capture([receipt_path.clone(), arrow_path.clone()])?;
            let receipt_bytes = read_limited(
                &receipt_path,
                limited_size(&receipt_path, self.limits.max_partition_receipt_bytes)?,
            )?;
            if blake3::hash(&receipt_bytes).to_hex().as_str() != part.receipt.blake3 {
                return Err(DatasetError::Corrupt(
                    "partition receipt hash differs from descriptor".into(),
                ));
            }
            let receipt = parse_receipt(&receipt_bytes)?;
            validate_partition_receipt(&receipt, part, &self.descriptor, self.fields)?;
            if hash_budgeted(&arrow_path, budget)? != part.arrow.blake3 {
                return Err(DatasetError::Corrupt(
                    "partition Arrow hash differs from descriptor".into(),
                ));
            }
            snapshot.verify()?;
            let mut output = EvidenceBatch::new(
                part.contig,
                self.contigs
                    .by_id(part.contig)
                    .expect("validated dictionary")
                    .name
                    .to_string(),
                part.start,
                query.fields,
                Vec::new(),
            );
            output.reserve_exact(EVIDENCE_ARROW_BATCH_ROWS);
            let mut interval_index = 0usize;
            let mut position = part.intervals[0].start;
            let mut row_count = 0usize;
            let mut consumer_failed = false;
            let mut query_failed = false;
            let result = read_evidence_batches_expected_fields(
                BufReader::new(File::open(&arrow_path)?),
                &self.contigs,
                self.fields,
                |batch| {
                    runtime_check(budget)?;
                    if batch.contig_id != part.contig || batch.canonical_tile_start != part.start {
                        return Err(EvidenceError::InvalidInput(
                            "source batch has incorrect canonical ownership".into(),
                        ));
                    }
                    for mut row in batch.rows() {
                        let owned = part.intervals.get(interval_index).ok_or_else(|| {
                            EvidenceError::InvalidInput("source partition has extra rows".into())
                        })?;
                        if row.position != position {
                            return Err(EvidenceError::InvalidInput(
                                "source partition omits or duplicates selected loci".into(),
                            ));
                        }
                        if let Some(source_sites) = &part.sites {
                            let source = source_sites.get(row_count).ok_or_else(|| {
                                EvidenceError::InvalidInput(
                                    "source partition has an unselected site".into(),
                                )
                            })?;
                            if source.position != row.position
                                || source.reference != row.reference
                                || source.alternates != row.requested_alts
                            {
                                return Err(EvidenceError::InvalidInput(
                                    "stored row REF/ALT differs from its descriptor".into(),
                                ));
                            }
                        } else if !row.requested_alts.is_empty() {
                            return Err(EvidenceError::InvalidInput(
                                "interval source has unexpected ALT annotations".into(),
                            ));
                        }
                        row_count += 1;
                        position += 1;
                        if position == owned.end {
                            interval_index += 1;
                            if let Some(next) = part.intervals.get(interval_index) {
                                position = next.start;
                            }
                        }
                        if !contains(&intervals, part.contig, row.position) {
                            continue;
                        }
                        if let Some(site) = sites.get(&(part.contig, row.position)) {
                            if site.reference != row.reference {
                                query_failed = true;
                                return Err(EvidenceError::InvalidInput(
                                    "query SNV REF differs from stored reference base".into(),
                                ));
                            }
                            row.requested_alts = &site.alternates;
                        } else {
                            row.requested_alts = &[];
                        }
                        output.push_row(row)?;
                        if output.len() == EVIDENCE_ARROW_BATCH_ROWS {
                            let result = analyzer.on_batch(&output);
                            consumer_failed = result.is_err();
                            result?;
                            runtime_check(budget)?;
                            stats.emitted_loci += output.len() as u64;
                            output.clear();
                        }
                    }
                    Ok(())
                },
            );
            let metadata = result.map_err(|error| {
                if consumer_failed || query_failed || matches!(error, EvidenceError::Core(_)) {
                    DatasetError::Evidence(error)
                } else {
                    DatasetError::Corrupt(format!("invalid persisted Arrow evidence: {error}"))
                }
            })?;
            if interval_index != part.intervals.len()
                || row_count as u64 != part.expected_rows
                || metadata.fields != self.fields
                || metadata.schema_version != self.fields.schema_version()
            {
                return Err(DatasetError::Corrupt(
                    "source partition is incomplete or has inconsistent fields".into(),
                ));
            }
            snapshot.verify()?;
            self.snapshot.verify()?;
            if !output.is_empty() {
                analyzer.on_batch(&output)?;
                runtime_check(budget)?;
                stats.emitted_loci += output.len() as u64;
            }
            self.observed.insert(index, snapshot);
            stats.microtiles += 1;
        }
        if stats.emitted_loci != plan.selected_loci {
            return Err(DatasetError::Corrupt(
                "query did not emit every selected locus".into(),
            ));
        }
        analyzer.finish()?;
        runtime_check(budget)?;
        self.snapshot.verify()?;
        for index in selected {
            self.observed[&index].verify()?;
        }
        stats.peak_rss_bytes = crate::util::rss::peak_rss_bytes();
        if let Some(budget) = execution
            .memory_budget_bytes
            .or(self.limits.memory_budget_bytes)
        {
            if stats.peak_rss_bytes > budget {
                return Err(EvidenceError::Core(crate::core::CoreError::BudgetExceeded {
                    needed: stats.peak_rss_bytes,
                    budget,
                })
                .into());
            }
        }
        Ok(stats)
    }

    /// Verify ordinary mutation of the parent, descriptor and every previously
    /// consumed partition. Call again immediately before artifact publication.
    pub fn verify_unchanged(&self) -> Result<(), DatasetError> {
        crate::core::governor::checkpoint().map_err(EvidenceError::from)?;
        self.snapshot.verify()?;
        for snapshot in self.observed.values() {
            snapshot.verify()?;
        }
        Ok(())
    }
}

fn admit(needed: u64, budget: Option<u64>) -> Result<(), DatasetError> {
    if let Some(budget) = budget {
        if needed > budget {
            return Err(EvidenceError::Refused { needed, budget }.into());
        }
    }
    Ok(())
}
pub(super) fn runtime_check(budget: Option<u64>) -> Result<(), EvidenceError> {
    crate::core::governor::checkpoint()?;
    if let Some(budget) = budget {
        let needed = crate::util::rss::peak_rss_bytes();
        if needed > budget {
            return Err(crate::core::CoreError::BudgetExceeded { needed, budget }.into());
        }
    }
    Ok(())
}
pub(super) fn hash_budgeted(path: &Path, budget: Option<u64>) -> Result<String, EvidenceError> {
    let mut reader = File::open(path)?;
    let mut buffer = [0u8; 64 << 10];
    let mut hash = blake3::Hasher::new();
    loop {
        runtime_check(budget)?;
        let count = reader.read(&mut buffer)?;
        if count == 0 {
            break;
        }
        hash.update(&buffer[..count]);
    }
    Ok(hash.finalize().to_hex().to_string())
}
fn limited_size(path: &Path, maximum: usize) -> Result<usize, DatasetError> {
    let metadata = fs::metadata(path)?;
    if !metadata.is_file() || metadata.len() > maximum as u64 {
        return Err(DatasetError::Corrupt(format!(
            "metadata exceeds its {maximum}-byte regular-file envelope: {}",
            path.display()
        )));
    }
    Ok(metadata.len() as usize)
}
fn read_limited(path: &Path, size: usize) -> Result<Vec<u8>, DatasetError> {
    let mut bytes = Vec::with_capacity(size);
    File::open(path)?
        .take(size as u64 + 1)
        .read_to_end(&mut bytes)?;
    if bytes.len() != size {
        return Err(DatasetError::Corrupt(
            "metadata length changed while reading".into(),
        ));
    }
    Ok(bytes)
}
fn parse_receipt(bytes: &[u8]) -> Result<RunManifest, DatasetError> {
    let text = std::str::from_utf8(bytes)
        .map_err(|_| DatasetError::Corrupt("receipt is not UTF-8".into()))?;
    let receipt = RunManifest::from_canonical_json(text)
        .map_err(|error| DatasetError::Corrupt(error.to_string()))?;
    if receipt.self_hash_ok() != Some(true) || receipt.measurement_hash_ok() == Some(false) {
        return Err(DatasetError::Corrupt(
            "portable receipt integrity failed".into(),
        ));
    }
    Ok(receipt)
}
fn local_path(root: &Path, relative: &str) -> Result<PathBuf, DatasetError> {
    let path = Path::new(relative);
    if path
        .components()
        .any(|c| !matches!(c, std::path::Component::Normal(_)))
    {
        return Err(DatasetError::Corrupt(
            "dataset artifact path is not a relative descendant".into(),
        ));
    }
    let path = root.join(path);
    if path.as_os_str().len() > 4096 || !fs::canonicalize(&path)?.starts_with(root) {
        return Err(DatasetError::Corrupt(
            "dataset artifact escapes its relocated directory or path envelope".into(),
        ));
    }
    Ok(path)
}
fn validate_parent(
    receipt: &RunManifest,
    descriptor: &DatasetDescriptor,
    hash: &str,
) -> Result<(), DatasetError> {
    let expected = [
        ("dataset.descriptor_version", descriptor.version.to_string()),
        ("dataset.namespace", descriptor.dataset_namespace.clone()),
        (
            "dataset.compatibility_blake3",
            descriptor.compatibility_blake3.clone(),
        ),
        ("dataset.request_blake3", descriptor.request_blake3.clone()),
        ("evidence.schema", descriptor.schema_version.to_string()),
        ("evidence.fields", descriptor.fields.to_string()),
        (
            "evidence.fields_version",
            descriptor.fields_version.to_string(),
        ),
        ("evidence.semantics", descriptor.semantics.clone()),
        (
            "dataset.partition_count",
            descriptor.partitions.len().to_string(),
        ),
        ("run_status", "completed".into()),
        (
            "outcome.rows",
            descriptor
                .partitions
                .iter()
                .map(|part| part.expected_rows)
                .sum::<u64>()
                .to_string(),
        ),
    ];
    if receipt.subcommand != "dataset evidence portable"
        || expected
            .iter()
            .any(|(key, value)| receipt.params.get(*key) != Some(value))
    {
        return Err(DatasetError::Corrupt(
            "portable receipt disagrees with its descriptor".into(),
        ));
    }
    let mut outputs = vec![(DATASET_DESCRIPTOR_NAME, hash)];
    outputs.extend(
        descriptor
            .partitions
            .iter()
            .map(|part| (part.arrow.path.as_str(), part.arrow.blake3.as_str())),
    );
    outputs.sort_unstable();
    let mut actual: Vec<_> = receipt
        .outputs
        .iter()
        .map(|file| (file.path.as_str(), file.blake3.as_str()))
        .collect();
    actual.sort_unstable();
    if actual != outputs {
        return Err(DatasetError::Corrupt(
            "portable receipt output ownership differs from descriptor".into(),
        ));
    }
    let mut expected_inputs: Vec<_> = descriptor
        .sources
        .iter()
        .map(|source| (source.path.as_str(), source.blake3.as_str()))
        .collect();
    expected_inputs.extend(
        descriptor
            .partitions
            .iter()
            .map(|part| (part.receipt.path.as_str(), part.receipt.blake3.as_str())),
    );
    expected_inputs.sort_unstable();
    // A source may have multiple provenance roles but one receipt input identity.
    expected_inputs.dedup();
    let mut actual: Vec<_> = receipt
        .inputs
        .iter()
        .map(|file| (file.path.as_str(), file.blake3.as_str()))
        .collect();
    actual.sort_unstable();
    actual.dedup();
    if actual != expected_inputs {
        return Err(DatasetError::Corrupt(
            "portable receipt input identities differ from descriptor".into(),
        ));
    }
    Ok(())
}
fn validate_partition_receipt(
    receipt: &RunManifest,
    part: &DescriptorPartition,
    descriptor: &DatasetDescriptor,
    fields: EvidenceFields,
) -> Result<(), DatasetError> {
    let ownership = Partition {
        contig: part.contig,
        start: part.start,
        intervals: part
            .intervals
            .iter()
            .map(|i| GenomicInterval {
                contig: i.contig,
                start: i.start,
                end: i.end,
            })
            .collect(),
        selection: match &part.sites {
            Some(sites) => EvidenceSelection::Sites(
                sites
                    .iter()
                    .map(|s| SnvSite {
                        contig: s.contig,
                        position: s.position,
                        reference: s.reference,
                        alternates: s.alternates.clone(),
                    })
                    .collect(),
            ),
            None => EvidenceSelection::Intervals(Vec::new()),
        },
    };
    let expected = [
        ("dataset.request_blake3", descriptor.request_blake3.clone()),
        (
            "evidence.science_blake3",
            descriptor.dataset_namespace.clone(),
        ),
        ("pileup.semantics", "exact-or-fail-v1".into()),
        ("dataset.contig", part.contig.to_string()),
        ("dataset.partition_start", part.start.to_string()),
        ("dataset.partition_bases", CANONICAL_TILE_BASES.to_string()),
        ("dataset.selection_blake3", ownership.selection_hash()),
        ("outcome.rows", part.expected_rows.to_string()),
        ("run_status", "completed".into()),
    ];
    if expected
        .iter()
        .any(|(key, value)| receipt.params.get(*key) != Some(value))
        || receipt_fields(receipt)? != fields
        || receipt.outputs.len() != 1
        || receipt.outputs[0].blake3 != part.arrow.blake3
    {
        return Err(DatasetError::Corrupt(
            "partition receipt disagrees with portable ownership/identity".into(),
        ));
    }
    Ok(())
}
fn interval_count(intervals: &[GenomicInterval]) -> u64 {
    intervals.iter().map(|i| u64::from(i.end - i.start)).sum()
}
fn contains(intervals: &[GenomicInterval], contig: u32, position: u32) -> bool {
    let index = intervals.partition_point(|i| (i.contig, i.end) <= (contig, position));
    intervals
        .get(index)
        .is_some_and(|i| i.contig == contig && i.start <= position)
}
fn overlaps(intervals: &[GenomicInterval], contig: u32, start: u32, end: u32) -> bool {
    let index = intervals.partition_point(|i| (i.contig, i.end) <= (contig, start));
    intervals
        .get(index)
        .is_some_and(|i| i.contig == contig && i.start < end)
}
fn split_intervals(
    query: &[GenomicInterval],
    source: &[DescriptorInterval],
) -> (Vec<GenomicInterval>, Vec<GenomicInterval>) {
    let mut covered = Vec::new();
    let mut missing = Vec::new();
    for query in query {
        let mut position = query.start;
        let cursor =
            source.partition_point(|part| (part.contig, part.end) <= (query.contig, position));
        let mut current = cursor;
        while position < query.end {
            let part = source
                .get(current)
                .filter(|s| s.contig == query.contig && s.start < query.end);
            let Some(part) = part else {
                missing.push(GenomicInterval {
                    contig: query.contig,
                    start: position,
                    end: query.end,
                });
                break;
            };
            if part.start > position {
                let end = part.start.min(query.end);
                missing.push(GenomicInterval {
                    contig: query.contig,
                    start: position,
                    end,
                });
                position = end;
            }
            if position < part.end {
                let end = part.end.min(query.end);
                covered.push(GenomicInterval {
                    contig: query.contig,
                    start: position,
                    end,
                });
                position = end;
            }
            current += 1;
        }
    }
    (covered, missing)
}
fn query_memory_bytes(selection: &EvidenceSelection, contigs: usize) -> u64 {
    match selection {
        EvidenceSelection::WholeGenome => (contigs as u64).saturating_mul(256),
        EvidenceSelection::Intervals(intervals) => (intervals.len() as u64).saturating_mul(256),
        EvidenceSelection::Sites(sites) => sites.iter().fold(0u64, |sum, site| {
            sum.saturating_add(1024 + site.alternates.len() as u64 * 4)
        }),
    }
}

#[cfg(test)]
pub(super) mod tests {
    use super::*;
    use crate::evidence::{
        EvidenceAlleleQuality, EvidenceCallback, EvidenceDepths, EvidenceRequest,
    };
    use rust_htslib::bam::{self, header::HeaderRecord, record::Cigar, record::CigarString};

    type Row = (
        u32,
        u32,
        u8,
        Vec<u8>,
        EvidenceDepths,
        [u64; 4],
        EvidenceAlleleQuality,
    );
    fn fields() -> EvidenceFields {
        EvidenceFields::DEPTHS
            .union(EvidenceFields::ALLELES)
            .union(EvidenceFields::ALLELE_QUALITY)
    }
    fn rows(batch: &EvidenceBatch, output: &mut Vec<Row>) {
        output.extend(batch.rows().map(|row| {
            (
                batch.contig_id,
                row.position,
                row.reference,
                row.requested_alts.to_vec(),
                *row.depths.unwrap(),
                row.alleles.unwrap().allele_counts,
                *row.allele_quality.unwrap(),
            )
        }));
    }
    pub(in crate::dataset) struct Fixture {
        pub(in crate::dataset) root: PathBuf,
        engine: EvidenceEngine,
        pub(in crate::dataset) manifest: PathBuf,
    }
    impl Fixture {
        pub(in crate::dataset) fn new(sites: bool) -> Self {
            let root = std::env::temp_dir().join(format!(
                "rosalind-persisted-{}-{}",
                std::process::id(),
                TEMP_SEQUENCE.fetch_add(1, Ordering::Relaxed)
            ));
            fs::create_dir(&root).unwrap();
            let length = CANONICAL_TILE_BASES + 6;
            let reference = root.join("reference.fa");
            fs::write(
                &reference,
                format!(">chr1\n{}\n", "A".repeat(length as usize)),
            )
            .unwrap();
            let fai = root.join("reference.fa.fai");
            fs::write(
                &fai,
                format!("chr1\t{length}\t6\t{length}\t{}\n", length + 1),
            )
            .unwrap();
            let alignment = root.join("reads.bam");
            let mut header = bam::Header::new();
            header.push_record(HeaderRecord::new(b"HD").push_tag(b"SO", "coordinate"));
            header.push_record(
                HeaderRecord::new(b"SQ")
                    .push_tag(b"SN", "chr1")
                    .push_tag(b"LN", length),
            );
            let mut writer = bam::Writer::from_path(&alignment, &header, bam::Format::Bam).unwrap();
            for (position, base, quality, mapq) in [
                (1, b'C', 30, 60),
                (1025, b'G', 35, 42),
                (CANONICAL_TILE_BASES + 1, b'T', 40, 55),
            ] {
                let mut record = bam::Record::new();
                record.set(
                    format!("r{position}").as_bytes(),
                    Some(&CigarString(vec![Cigar::Match(1)])),
                    &[base],
                    &[quality],
                );
                record.set_tid(0);
                record.set_pos(i64::from(position));
                record.set_flags(0);
                record.set_mapq(mapq);
                writer.write(&record).unwrap();
            }
            drop(writer);
            bam::index::build(&alignment, None::<&PathBuf>, bam::index::Type::Bai, 1).unwrap();
            let index = root.join("reads.bam.bai");
            let mut request = EvidenceRequest::new(&alignment, &reference);
            request.alignment_index = Some(index.clone());
            request.reference_fai = Some(fai.clone());
            request.fields = EvidenceFields::ALL_SUPPORTED;
            request.selection = if sites {
                EvidenceSelection::Sites(
                    [1, 2, 1025, CANONICAL_TILE_BASES + 1]
                        .into_iter()
                        .map(|position| SnvSite {
                            contig: 0,
                            position,
                            reference: b'A',
                            alternates: vec![b'C'],
                        })
                        .collect(),
                )
            } else {
                EvidenceSelection::Intervals(vec![
                    GenomicInterval {
                        contig: 0,
                        start: 0,
                        end: 2300,
                    },
                    GenomicInterval {
                        contig: 0,
                        start: CANONICAL_TILE_BASES,
                        end: CANONICAL_TILE_BASES + 5,
                    },
                ])
            };
            let mut engine = EvidenceEngine::open(request).unwrap();
            let session = VerifiedInputSession::open([
                ("alignments".into(), alignment),
                ("alignment-index".into(), index),
                ("reference".into(), reference),
                ("reference-fai".into(), fai),
            ])
            .unwrap();
            let namespace = session.dataset_namespace(&engine).unwrap();
            let mut analyzer = EvidenceCallback::with_fields(
                |_: &EvidenceBatch| Ok(()),
                0,
                EvidenceFields::ALL_SUPPORTED,
            );
            let outcome = run_dataset_with_snapshot(
                &mut engine,
                &namespace,
                &DatasetOptions {
                    cache_dir: root.join("cache"),
                    resume: false,
                    workers: 1,
                },
                &mut analyzer,
                session.snapshot(),
            )
            .unwrap();
            let manifest =
                publish_evidence_dataset(&engine, &outcome, &session, DescriptorLimits::default())
                    .unwrap();
            Self {
                root,
                engine,
                manifest,
            }
        }
        pub(in crate::dataset) fn open(&self) -> VerifiedEvidenceDataset {
            VerifiedEvidenceDataset::open(&self.manifest, DatasetReadLimits::default()).unwrap()
        }
        fn fresh(&self, query: &DatasetQuery) -> Vec<Row> {
            let mut request = self.engine.request().clone();
            request.selection = query.selection.clone();
            request.fields = query.fields;
            let mut engine = EvidenceEngine::open(request).unwrap();
            let mut output = Vec::new();
            let mut analyzer = EvidenceCallback::with_fields(
                |batch: &EvidenceBatch| {
                    rows(batch, &mut output);
                    Ok(())
                },
                1 << 20,
                query.fields,
            );
            engine.run(&mut analyzer).unwrap();
            output
        }
        fn reseal_descriptor(&self, update: impl FnOnce(&mut DatasetDescriptor)) {
            let path = self
                .manifest
                .parent()
                .unwrap()
                .join(DATASET_DESCRIPTOR_NAME);
            let mut descriptor = DatasetDescriptor::from_bytes(
                &fs::read(&path).unwrap(),
                DescriptorLimits::default(),
            )
            .unwrap();
            update(&mut descriptor);
            // Intentionally bypass validate when constructing corruption cases.
            fs::write(&path, serde_json::to_vec(&descriptor).unwrap()).unwrap();
            let mut parent = read_receipt(&self.manifest).unwrap();
            parent
                .outputs
                .iter_mut()
                .find(|f| f.path == DATASET_DESCRIPTOR_NAME)
                .unwrap()
                .blake3 = blake3_file(&path).unwrap();
            parent.finalize();
            fs::write(&self.manifest, parent.to_canonical_json()).unwrap();
        }
    }
    impl Drop for Fixture {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.root);
        }
    }
    fn collect(
        dataset: &mut VerifiedEvidenceDataset,
        query: &DatasetQuery,
    ) -> Result<Vec<Row>, DatasetError> {
        let mut output = Vec::new();
        let mut analyzer = EvidenceCallback::with_fields(
            |batch: &EvidenceBatch| {
                assert!(batch.len() <= 1024);
                assert_eq!(batch.fields(), query.fields);
                rows(batch, &mut output);
                Ok(())
            },
            1 << 20,
            query.fields,
        );
        dataset.visit_batches(query, &mut analyzer, &EvidenceExecution::default())?;
        Ok(output)
    }

    #[test]
    fn projected_subset_and_new_alts_equal_fresh_after_relocation_without_sources() {
        let fixture = Fixture::new(false);
        let query = DatasetQuery {
            selection: EvidenceSelection::Sites(vec![
                SnvSite {
                    contig: 0,
                    position: 1,
                    reference: b'A',
                    alternates: vec![b'T', b'G'],
                },
                SnvSite {
                    contig: 0,
                    position: 2,
                    reference: b'A',
                    alternates: vec![b'C'],
                },
                SnvSite {
                    contig: 0,
                    position: CANONICAL_TILE_BASES + 1,
                    reference: b'A',
                    alternates: vec![b'C', b'T'],
                },
            ]),
            fields: fields(),
        };
        let expected = fixture.fresh(&query);
        let original = fixture.manifest.parent().unwrap();
        let relocated = fixture.root.join("relocated");
        fs::rename(original, &relocated).unwrap();
        for name in [
            "reads.bam",
            "reads.bam.bai",
            "reference.fa",
            "reference.fa.fai",
        ] {
            fs::remove_file(fixture.root.join(name)).unwrap();
        }
        let mut dataset = VerifiedEvidenceDataset::open(
            relocated.join(EVIDENCE_DATASET_MANIFEST_NAME),
            DatasetReadLimits::default(),
        )
        .unwrap();
        assert_eq!(collect(&mut dataset, &query).unwrap(), expected);
        assert_eq!(dataset.source_hashes().len(), 2);
        assert_eq!(dataset.verified_partition_hashes().unwrap().len(), 4);
        assert_eq!(expected[1].4.callable_depth, 0);
        assert_eq!(expected[0].3, b"GT");
    }

    #[test]
    fn source_alt_annotations_do_not_limit_new_queries_but_reference_must_match() {
        let fixture = Fixture::new(true);
        let mut dataset = fixture.open();
        let mut query = DatasetQuery {
            selection: EvidenceSelection::Sites(vec![SnvSite {
                contig: 0,
                position: 1,
                reference: b'A',
                alternates: vec![b'G', b'T'],
            }]),
            fields: fields(),
        };
        assert_eq!(
            collect(&mut dataset, &query).unwrap(),
            fixture.fresh(&query)
        );
        if let EvidenceSelection::Sites(sites) = &mut query.selection {
            sites[0].reference = b'C';
        }
        assert!(matches!(
            collect(&mut dataset, &query),
            Err(DatasetError::Evidence(EvidenceError::InvalidInput(_)))
        ));
    }

    #[test]
    fn coverage_split_preserves_query_annotations_and_never_synthesizes_holes() {
        let fixture = Fixture::new(true);
        let dataset = fixture.open();
        let query = DatasetQuery {
            selection: EvidenceSelection::Sites(vec![
                SnvSite {
                    contig: 0,
                    position: 1,
                    reference: b'A',
                    alternates: vec![b'G'],
                },
                SnvSite {
                    contig: 0,
                    position: 3,
                    reference: b'A',
                    alternates: vec![b'T'],
                },
            ]),
            fields: fields(),
        };
        let split = dataset.coverage_split(&query.selection).unwrap();
        assert_eq!((split.covered_loci, split.missing_loci), (1, 1));
        let EvidenceSelection::Sites(missing) = split.missing else {
            panic!("site annotations must survive")
        };
        assert_eq!(missing[0].alternates, b"T");
        let analyzer = EvidenceCallback::with_fields(|_: &EvidenceBatch| Ok(()), 0, fields());
        assert!(dataset
            .plan(&query, &analyzer, &EvidenceExecution::default())
            .is_err());
        let coverage = dataset
            .coverage_split(&EvidenceSelection::Intervals(vec![GenomicInterval {
                contig: 0,
                start: 0,
                end: 5,
            }]))
            .unwrap();
        assert_eq!((coverage.covered_loci, coverage.missing_loci), (2, 3));
    }

    #[test]
    fn output_batches_are_canonical_and_planner_reserves_source_projection_and_analyzer() {
        let fixture = Fixture::new(false);
        let mut dataset = fixture.open();
        let query = DatasetQuery {
            selection: dataset.selection(),
            fields: EvidenceFields::DEPTHS,
        };
        let mut batches = Vec::new();
        let mut analyzer = EvidenceCallback::with_fields(
            |batch: &EvidenceBatch| {
                assert!(batch.rows().all(|row| row.quality_histograms.is_none()));
                batches.push((batch.canonical_tile_start, batch.len()));
                Ok(())
            },
            777,
            query.fields,
        );
        let plan = dataset
            .plan(&query, &analyzer, &EvidenceExecution::default())
            .unwrap();
        assert_eq!(
            plan.source_decoder_bytes,
            evidence_reader_memory_bytes(EvidenceFields::ALL_SUPPORTED)
        );
        assert!(plan.source_decoder_bytes > evidence_reader_memory_bytes(query.fields));
        assert_eq!(plan.analyzer_bytes, 777);
        let mut execution = EvidenceExecution {
            memory_budget_bytes: Some(plan.predicted_peak_rss_bytes - 1),
            ..Default::default()
        };
        assert!(matches!(
            dataset.plan(&query, &analyzer, &execution),
            Err(DatasetError::Evidence(EvidenceError::Refused { .. }))
        ));
        execution.memory_budget_bytes = None;
        dataset
            .visit_batches(&query, &mut analyzer, &execution)
            .unwrap();
        assert_eq!(
            batches,
            [(0, 1024), (0, 1024), (0, 252), (CANONICAL_TILE_BASES, 5)]
        );
        let different = EvidenceProfile {
            min_base_quality: 21,
            ..Default::default()
        };
        assert!(dataset.require_profile(&different).is_err());
        dataset
            .require_profile(&EvidenceProfile::default())
            .unwrap();
    }

    #[test]
    fn corrupt_hash_unknown_versions_and_mutation_are_rejected() {
        let fixture = Fixture::new(true);
        let mut dataset = fixture.open();
        let query = DatasetQuery {
            selection: dataset.selection(),
            fields: fields(),
        };
        collect(&mut dataset, &query).unwrap();
        let path = dataset
            .root
            .join(&dataset.descriptor.partitions[0].arrow.path);
        let mut bytes = fs::read(&path).unwrap();
        bytes[0] ^= 1;
        fs::write(&path, bytes).unwrap();
        assert!(dataset.verify_unchanged().is_err());
        let mut reopened = fixture.open();
        assert!(collect(&mut reopened, &query).is_err());
        fixture.reseal_descriptor(|descriptor| descriptor.version = 999);
        assert!(
            VerifiedEvidenceDataset::open(&fixture.manifest, DatasetReadLimits::default()).is_err()
        );
    }

    #[test]
    fn resealed_omitted_or_duplicate_source_rows_still_fail_exact_ownership() {
        for missing in [true, false] {
            let fixture = Fixture::new(true);
            let dataset = fixture.open();
            let part = &dataset.descriptor.partitions[0];
            let path = dataset.root.join(&part.arrow.path);
            let mut batches = Vec::new();
            read_evidence_batches_expected_fields(
                BufReader::new(File::open(&path).unwrap()),
                dataset.contigs(),
                dataset.fields(),
                |batch| {
                    batches.push(batch.clone());
                    Ok(())
                },
            )
            .unwrap();
            let mut changed = EvidenceBatch::new(0, "chr1", 0, dataset.fields(), Vec::new());
            for (index, row) in batches.iter().flat_map(EvidenceBatch::rows).enumerate() {
                if missing && index == 0 {
                    continue;
                }
                changed.push_row(row).unwrap();
                if !missing && index == 0 {
                    changed.push_row(row).unwrap();
                }
            }
            let mut writer =
                EvidenceArrowWriter::with_fields(File::create(&path).unwrap(), dataset.fields());
            writer.on_batch(&changed).unwrap();
            writer.finish().unwrap();
            let arrow_hash = blake3_file(&path).unwrap();
            let receipt_path = dataset.root.join(&part.receipt.path);
            let mut receipt = read_receipt(&receipt_path).unwrap();
            receipt.outputs[0].blake3 = arrow_hash.clone();
            receipt.finalize();
            fs::write(&receipt_path, receipt.to_canonical_json()).unwrap();
            let receipt_hash = blake3_file(&receipt_path).unwrap();
            fixture.reseal_descriptor(|descriptor| {
                descriptor.partitions[0].arrow.blake3 = arrow_hash.clone();
                descriptor.partitions[0].receipt.blake3 = receipt_hash.clone();
            });
            let mut parent = read_receipt(&fixture.manifest).unwrap();
            parent
                .outputs
                .iter_mut()
                .find(|f| f.path == part.arrow.path)
                .unwrap()
                .blake3 = arrow_hash;
            parent
                .inputs
                .iter_mut()
                .find(|f| f.path == part.receipt.path)
                .unwrap()
                .blake3 = receipt_hash;
            parent.finalize();
            fs::write(&fixture.manifest, parent.to_canonical_json()).unwrap();
            let mut dataset = fixture.open();
            let query = DatasetQuery {
                selection: dataset.selection(),
                fields: fields(),
            };
            assert!(collect(&mut dataset, &query).is_err());
        }
    }

    #[test]
    fn metadata_open_refuses_before_parsing_when_budget_or_envelope_is_too_small() {
        let fixture = Fixture::new(true);
        assert!(matches!(
            VerifiedEvidenceDataset::open(
                &fixture.manifest,
                DatasetReadLimits {
                    memory_budget_bytes: Some(1),
                    ..Default::default()
                }
            ),
            Err(DatasetError::Evidence(EvidenceError::Refused { .. }))
        ));
        assert!(VerifiedEvidenceDataset::open(
            &fixture.manifest,
            DatasetReadLimits {
                max_descriptor_bytes: 1,
                ..Default::default()
            }
        )
        .is_err());
    }

    #[test]
    fn local_budget_stops_after_first_breaching_callback_without_a_global_governor() {
        const CHILD: &str = "ROSALIND_PERSISTED_BUDGET_CHILD";
        if std::env::var_os(CHILD).is_none() {
            let result=std::process::Command::new(std::env::current_exe().unwrap())
                .args(["--exact","dataset::persisted::tests::local_budget_stops_after_first_breaching_callback_without_a_global_governor","--nocapture"])
                .env(CHILD,"1").output().unwrap();
            assert!(
                result.status.success(),
                "{}{}",
                String::from_utf8_lossy(&result.stdout),
                String::from_utf8_lossy(&result.stderr)
            );
            return;
        }
        let fixture = Fixture::new(false);
        let mut dataset = fixture.open();
        let query = DatasetQuery {
            selection: dataset.selection(),
            fields: EvidenceFields::DEPTHS,
        };
        let mut calls = 0;
        let mut retained = Vec::new();
        let mut analyzer = EvidenceCallback::with_fields(
            |_: &EvidenceBatch| {
                calls += 1;
                retained.resize(128 << 20, 1u8);
                std::hint::black_box(&retained);
                Ok(())
            },
            0,
            query.fields,
        );
        let plan = dataset
            .plan(&query, &analyzer, &EvidenceExecution::default())
            .unwrap();
        let execution = EvidenceExecution {
            memory_budget_bytes: Some(plan.predicted_peak_rss_bytes + (16 << 20)),
            ..Default::default()
        };
        let result = dataset.visit_batches(&query, &mut analyzer, &execution);
        assert!(matches!(
            result,
            Err(DatasetError::Evidence(EvidenceError::Core(
                crate::core::CoreError::BudgetExceeded { .. }
            )))
        ));
        assert_eq!(
            calls, 1,
            "local checks must stop before delivering the next batch"
        );
    }
}
