//! Source-compatible reuse with one projected canonical partition of merge state.
//!
//! The source is read once per touched partition. All missing loci are passed to
//! one indexed extraction run; output ownership and order do not depend on which
//! source supplied a row. Inputs and persisted files remain immutable throughout.

use super::{DatasetError, DatasetQuery, VerifiedEvidenceDataset, VerifiedInputSession};
use crate::evidence::{
    EvidenceAnalyzer, EvidenceBatch, EvidenceCallback, EvidenceEngine, EvidenceError,
    EvidenceExecution, EvidenceFields, EvidenceRequirements, EvidenceRunStats, EvidenceSelection,
    SnvSite, CANONICAL_TILE_BASES, EVIDENCE_ARROW_BATCH_ROWS,
};
use crate::selection::GenomicInterval;

/// Execution controls for explicit persisted-source reuse.
#[derive(Debug, Clone, Copy)]
pub struct ReuseOptions {
    /// Currently exactly one; unsupported parallel execution is refused.
    pub workers: usize,
}
impl Default for ReuseOptions {
    fn default() -> Self {
        Self { workers: 1 }
    }
}

/// Joint source decoder, native extractor and projected merge reservation.
#[derive(Debug, Clone)]
pub struct ReusePlan {
    /// Versioned reservation identity.
    pub model_id: &'static str,
    /// Unique requested positions supplied by verified persisted evidence.
    pub reused_loci: u64,
    /// Unique requested positions requiring indexed extraction.
    pub computed_loci: u64,
    /// Maximum requested positions retained in one ownership partition.
    pub retained_loci: usize,
    /// Source metadata and additional normalized selection state.
    pub metadata_bytes: u64,
    /// Additional original, split and per-partition selection state.
    pub selection_bytes: u64,
    /// Effective whole-process limit, inherited from the source when needed.
    pub memory_budget_bytes: Option<u64>,
    /// Decoder reservation for physical source fields, even for smaller output.
    pub source_decoder_bytes: u64,
    /// Projected merge state, ordering indices and canonical output buffer.
    pub merge_bytes: u64,
    /// Declared downstream consumer reservation.
    pub analyzer_bytes: u64,
    /// Native extraction workspace; zero when every locus is reusable.
    pub native_bytes: u64,
    /// Whole-process admission estimate including the measured baseline once.
    pub predicted_peak_rss_bytes: u64,
    /// Fields physically read from the persisted source.
    pub source_fields: EvidenceFields,
    /// Exact projected groups emitted to the consumer.
    pub output_fields: EvidenceFields,
}

/// Successful reuse measurements. Fetch visits count only newly computed loci;
/// emitted_loci counts both sources, with each requested coordinate owned once.
#[derive(Debug, Clone)]
pub struct ReuseOutcome {
    /// Aggregate extraction measurements plus total emitted rows and final RSS.
    pub stats: EvidenceRunStats,
    /// Positions obtained from the verified source dataset.
    pub reused_loci: u64,
    /// Positions newly computed from alignments.
    pub computed_loci: u64,
}

/// Admit source-compatible reuse before starting an artifact consumer. Source
/// identities, profile, sample scope and dictionary must match exactly; source
/// selection and physical fields may be supersets of the requested capabilities.
pub fn plan_reuse(
    engine: &mut EvidenceEngine,
    session: &VerifiedInputSession,
    source: &VerifiedEvidenceDataset,
    analyzer: &dyn EvidenceAnalyzer,
    options: &ReuseOptions,
) -> Result<ReusePlan, DatasetError> {
    if options.workers != 1 {
        return Err(incompatible(
            "persisted-source reuse currently requires exactly one worker",
        ));
    }
    if session.compatibility_key(engine)? != source.descriptor().compatibility_blake3 {
        return Err(incompatible(
            "persisted source identities, filters, sample scope or reference dictionary differ",
        ));
    }
    source.verify_unchanged()?;
    let fields = engine.request().fields;
    if !source.fields().contains(fields) {
        return Err(incompatible(
            "persisted source is missing requested physical evidence groups",
        ));
    }
    let requirements = analyzer.requirements();
    if !fields.contains(requirements.fields)
        || requirements.context_bases != 0
        || (requirements.requires_reference && !source.descriptor().has_reference)
    {
        return Err(incompatible(
            "reuse request does not satisfy analyzer scientific requirements",
        ));
    }
    let execution = engine.request().execution.clone();
    let budget = execution
        .memory_budget_bytes
        .or(source.memory_budget_bytes());
    let consumer = match requirements.retained_bytes {
        Some(bytes) => bytes.max(execution.analyzer_bytes),
        None if budget.is_some() => {
            return Err(incompatible(
                "budgeted reuse requires a declared consumer memory bound",
            ));
        }
        None => execution.analyzer_bytes,
    };
    let split = source.coverage_split(&engine.request().selection)?;
    let retained_loci = PartitionSelections::new(engine.intervals(), &engine.request().selection)
        .map(|part| part.count())
        .max()
        .unwrap_or(0);
    let merge_bytes = checked_mul(fields.storage_bytes_per_locus() + 64, retained_loci as u64)?
        .checked_add(checked_mul(
            fields.storage_bytes_per_locus() + 40,
            EVIDENCE_ARROW_BATCH_ROWS as u64,
        )?)
        .and_then(|n| n.checked_add(64 << 10))
        .ok_or_else(|| incompatible("reuse merge reservation overflow"))?;
    let query = DatasetQuery {
        selection: split.covered,
        fields,
    };
    let mut read_execution = execution.clone();
    read_execution.memory_budget_bytes = None;
    read_execution.analyzer_bytes = 0;
    let read_plan = source.plan(
        &query,
        &EvidenceCallback::with_fields(|_: &EvidenceBatch| Ok(()), 0, fields),
        &read_execution,
    )?;
    // Original request, split selections and the temporarily replaced engine
    // index coexist. Preserve the full original selection reservation even when
    // the missing subset is smaller. Per-partition metadata is also bounded.
    let selection_bytes = selection_bytes(&engine.request().selection, engine.intervals().len())
        .saturating_add(selection_bytes(&query.selection, 0))
        .saturating_add(selection_bytes(&split.missing, 0))
        .saturating_mul(4)
        .saturating_add((retained_loci as u64).saturating_mul(512));
    let metadata_bytes = read_plan.metadata_bytes.saturating_add(selection_bytes);
    let baseline = engine
        .plan()
        .baseline_rss_bytes
        .max(read_plan.baseline_rss_bytes);
    let baseline_adjustment = baseline.saturating_sub(engine.plan().baseline_rss_bytes);
    let extra = consumer
        .saturating_add(metadata_bytes)
        .saturating_add(read_plan.source_decoder_bytes)
        .saturating_add(read_plan.projection_bytes)
        .saturating_add(merge_bytes)
        .saturating_add(baseline_adjustment);
    let (native_bytes, predicted_peak_rss_bytes) = if split.missing_loci == 0 {
        (
            0,
            baseline.saturating_add(extra.saturating_sub(baseline_adjustment)),
        )
    } else {
        let plan = engine.plan_for_analyzer(&EvidenceCallback::with_fields(
            |_: &EvidenceBatch| Ok(()),
            extra,
            fields,
        ))?;
        (
            plan.predicted_peak_rss_bytes
                .saturating_sub(plan.baseline_rss_bytes)
                .saturating_sub(extra),
            plan.predicted_peak_rss_bytes,
        )
    };
    if let Some(budget) = budget {
        if predicted_peak_rss_bytes > budget {
            return Err(EvidenceError::Refused {
                needed: predicted_peak_rss_bytes,
                budget,
            }
            .into());
        }
    }
    Ok(ReusePlan {
        model_id: "source-compatible-reuse-v1",
        reused_loci: split.covered_loci,
        computed_loci: split.missing_loci,
        retained_loci,
        metadata_bytes,
        selection_bytes,
        memory_budget_bytes: budget,
        source_decoder_bytes: read_plan.source_decoder_bytes,
        merge_bytes,
        analyzer_bytes: consumer,
        native_bytes,
        predicted_peak_rss_bytes,
        source_fields: source.fields(),
        output_fields: fields,
    })
}

/// Stream the same requested evidence as fresh extraction while computing only
/// uncovered loci. The original engine selection is restored on success or
/// failure. Consumer.finish is called once, after all source guards pass.
/// Artifact consumers must stage output until this function succeeds.
pub fn run_reusing_dataset(
    engine: &mut EvidenceEngine,
    session: &VerifiedInputSession,
    source: &mut VerifiedEvidenceDataset,
    analyzer: &mut dyn EvidenceAnalyzer,
    options: &ReuseOptions,
) -> Result<ReuseOutcome, DatasetError> {
    let plan = plan_reuse(engine, session, source, analyzer, options)?;
    let original = engine.request().selection.clone();
    let original_intervals = engine.intervals().to_vec();
    let mut execution = engine.request().execution.clone();
    execution.memory_budget_bytes = plan.memory_budget_bytes;
    let split = source.coverage_split(&original)?;
    // Keep the same aggregate reservation when the engine changes selection.
    // Any extra original-selection reservation is retained conservatively.
    let extra = plan
        .predicted_peak_rss_bytes
        .saturating_sub(engine.plan().baseline_rss_bytes)
        .saturating_sub(plan.native_bytes);
    let mut merge = Merger::new(
        source,
        analyzer,
        session,
        &original_intervals,
        &original,
        &execution,
        &plan,
        extra,
    );
    let result = (|| {
        let mut stats = if plan.computed_loci == 0 {
            EvidenceRunStats::default()
        } else {
            engine.set_selection(split.missing)?;
            let result = engine.run(&mut merge);
            if let Some(error) = merge.error.take() {
                return Err(error);
            }
            result?
        };
        merge.drain()?;
        session.verify()?;
        merge.source.verify_unchanged()?;
        if merge.emitted != plan.reused_loci.saturating_add(plan.computed_loci) {
            return Err(incompatible(
                "reuse output denominator differs from requested selection",
            ));
        }
        runtime_check(plan.memory_budget_bytes)?;
        merge.analyzer.finish()?;
        runtime_check(plan.memory_budget_bytes)?;
        session.verify()?;
        merge.source.verify_unchanged()?;
        stats.emitted_loci = merge.emitted;
        stats.peak_rss_bytes = crate::util::rss::peak_rss_bytes();
        Ok(ReuseOutcome {
            stats,
            reused_loci: plan.reused_loci,
            computed_loci: plan.computed_loci,
        })
    })();
    drop(merge);
    let restored = engine.set_selection(original).map_err(DatasetError::from);
    match result {
        Ok(outcome) => restored.map(|()| outcome),
        Err(error) => Err(error),
    }
}

struct Merger<'a> {
    source: &'a mut VerifiedEvidenceDataset,
    analyzer: &'a mut dyn EvidenceAnalyzer,
    session: &'a VerifiedInputSession,
    parts: PartitionSelections<'a>,
    current: Option<SelectionPart>,
    batch: EvidenceBatch,
    output: EvidenceBatch,
    order: Vec<usize>,
    execution: EvidenceExecution,
    collector_bytes: u64,
    extra: u64,
    emitted: u64,
    error: Option<DatasetError>,
}
impl<'a> Merger<'a> {
    #[allow(clippy::too_many_arguments)]
    fn new(
        source: &'a mut VerifiedEvidenceDataset,
        analyzer: &'a mut dyn EvidenceAnalyzer,
        session: &'a VerifiedInputSession,
        intervals: &'a [GenomicInterval],
        selection: &'a EvidenceSelection,
        execution: &EvidenceExecution,
        plan: &ReusePlan,
        extra: u64,
    ) -> Self {
        let mut batch = EvidenceBatch::new(0, "", 0, plan.output_fields, Vec::new());
        batch.reserve_exact(plan.retained_loci);
        let mut output = EvidenceBatch::new(0, "", 0, plan.output_fields, Vec::new());
        output.reserve_exact(EVIDENCE_ARROW_BATCH_ROWS);
        // Source reads execute inside the native callback, so both native and
        // merge buffers must remain reserved during source decoding as well.
        let collector_bytes = plan
            .merge_bytes
            .saturating_add(plan.native_bytes)
            .saturating_add(plan.analyzer_bytes)
            .saturating_add(plan.selection_bytes);
        Self {
            source,
            analyzer,
            session,
            parts: PartitionSelections::new(intervals, selection),
            current: None,
            batch,
            output,
            order: Vec::with_capacity(plan.retained_loci),
            execution: execution.clone(),
            collector_bytes,
            extra,
            emitted: 0,
            error: None,
        }
    }
    fn load_next(&mut self) -> Result<bool, DatasetError> {
        let Some(part) = self.parts.next() else {
            return Ok(false);
        };
        self.batch.clear();
        self.batch.contig_id = part.contig;
        self.batch.contig = self
            .source
            .contigs()
            .by_id(part.contig)
            .expect("verified dictionary")
            .name
            .to_string();
        self.batch.canonical_tile_start = part.start;
        let split = self.source.coverage_split(&part.selection)?;
        if split.covered_loci != 0 {
            let query = DatasetQuery {
                selection: split.covered,
                fields: self.batch.fields(),
            };
            let max_rows = part.count();
            let collector = &mut self.batch;
            self.source.visit_batches(
                &query,
                &mut EvidenceCallback::with_fields(
                    |batch: &EvidenceBatch| {
                        if batch.contig_id != collector.contig_id
                            || batch.canonical_tile_start != collector.canonical_tile_start
                            || collector.len().saturating_add(batch.len()) > max_rows
                        {
                            return Err(EvidenceError::InvalidInput(
                                "reused partition exceeds canonical ownership or row envelope"
                                    .into(),
                            ));
                        }
                        for row in batch.rows() {
                            collector.push_row(row)?;
                        }
                        Ok(())
                    },
                    self.collector_bytes,
                    query.fields,
                ),
                &self.execution,
            )?;
        }
        self.current = Some(part);
        Ok(true)
    }
    fn flush(&mut self) -> Result<(), DatasetError> {
        let Some(part) = self.current.take() else {
            return Ok(());
        };
        if self.batch.len() != part.count() {
            return Err(incompatible(
                "reused and computed rows do not cover the canonical partition",
            ));
        }
        self.order.clear();
        self.order.extend(0..self.batch.len());
        self.order
            .sort_unstable_by_key(|&index| self.batch.loci()[index].position);
        self.output.contig_id = part.contig;
        self.output.contig.clone_from(&self.batch.contig);
        self.output.canonical_tile_start = part.start;
        let mut expected = part
            .intervals
            .iter()
            .flat_map(|interval| interval.start..interval.end);
        for &index in &self.order {
            let row = self.batch.row(index).expect("owned index");
            if Some(row.position) != expected.next() {
                return Err(incompatible(
                    "reused and computed loci overlap or omit requested coordinates",
                ));
            }
            self.output.push_row(row)?;
            if self.output.len() == EVIDENCE_ARROW_BATCH_ROWS {
                self.analyzer.on_batch(&self.output)?;
                runtime_check(self.execution.memory_budget_bytes)?;
                self.emitted += self.output.len() as u64;
                self.output.clear();
            }
        }
        if !self.output.is_empty() {
            self.analyzer.on_batch(&self.output)?;
            runtime_check(self.execution.memory_budget_bytes)?;
            self.emitted += self.output.len() as u64;
            self.output.clear();
        }
        self.session.verify()?;
        crate::core::governor::checkpoint().map_err(EvidenceError::from)?;
        Ok(())
    }
    fn add_missing(&mut self, batch: &EvidenceBatch) -> Result<(), DatasetError> {
        runtime_check(self.execution.memory_budget_bytes)?;
        let key = (batch.contig_id, batch.canonical_tile_start);
        loop {
            if self.current.is_none() && !self.load_next()? {
                return Err(incompatible(
                    "native extraction emitted an unrequested partition",
                ));
            }
            let part = self.current.as_ref().expect("loaded partition");
            let current = (part.contig, part.start);
            if current == key {
                break;
            }
            if current > key {
                return Err(incompatible("native extraction partition order regressed"));
            }
            self.flush()?;
        }
        if self.batch.len().saturating_add(batch.len())
            > self.current.as_ref().expect("loaded partition").count()
        {
            return Err(incompatible(
                "native extraction exceeds requested partition denominator",
            ));
        }
        for row in batch.rows() {
            self.batch.push_row(row)?;
        }
        self.session.verify()?;
        Ok(())
    }
    fn drain(&mut self) -> Result<(), DatasetError> {
        self.flush()?;
        while self.load_next()? {
            self.flush()?;
        }
        Ok(())
    }
}
impl EvidenceAnalyzer for Merger<'_> {
    fn requirements(&self) -> EvidenceRequirements {
        EvidenceRequirements {
            fields: self.batch.fields(),
            requires_reference: false,
            context_bases: 0,
            retained_bytes: Some(self.extra),
        }
    }
    fn on_batch(&mut self, batch: &EvidenceBatch) -> Result<(), EvidenceError> {
        if let Err(error) = self.add_missing(batch) {
            self.error = Some(error);
            return Err(EvidenceError::InvalidInput(
                "verified dataset reuse failed".into(),
            ));
        }
        Ok(())
    }
}

struct SelectionPart {
    contig: u32,
    start: u32,
    intervals: Vec<GenomicInterval>,
    selection: EvidenceSelection,
}
impl SelectionPart {
    fn count(&self) -> usize {
        self.intervals
            .iter()
            .map(|i| (i.end - i.start) as usize)
            .sum()
    }
}
struct PartitionSelections<'a> {
    intervals: &'a [GenomicInterval],
    sites: Option<&'a [SnvSite]>,
    index: usize,
    position: u32,
}
impl<'a> PartitionSelections<'a> {
    fn new(intervals: &'a [GenomicInterval], selection: &'a EvidenceSelection) -> Self {
        Self {
            intervals,
            sites: match selection {
                EvidenceSelection::Sites(sites) => Some(sites),
                _ => None,
            },
            index: 0,
            position: intervals.first().map_or(0, |i| i.start),
        }
    }
}
impl Iterator for PartitionSelections<'_> {
    type Item = SelectionPart;
    fn next(&mut self) -> Option<Self::Item> {
        let interval = self.intervals.get(self.index)?;
        let contig = interval.contig;
        let start = self.position / CANONICAL_TILE_BASES * CANONICAL_TILE_BASES;
        let end = start.saturating_add(CANONICAL_TILE_BASES);
        let mut intervals = Vec::new();
        while let Some(interval) = self.intervals.get(self.index) {
            if interval.contig != contig || self.position >= end {
                break;
            }
            let stop = interval.end.min(end);
            intervals.push(GenomicInterval {
                contig,
                start: self.position,
                end: stop,
            });
            self.position = stop;
            if self.position == interval.end {
                self.index += 1;
                self.position = self.intervals.get(self.index).map_or(0, |i| i.start);
            } else {
                break;
            }
        }
        let selection = match self.sites {
            Some(sites) => {
                let from =
                    sites.partition_point(|site| (site.contig, site.position) < (contig, start));
                let to = sites.partition_point(|site| (site.contig, site.position) < (contig, end));
                EvidenceSelection::Sites(sites[from..to].to_vec())
            }
            None => EvidenceSelection::Intervals(intervals.clone()),
        };
        Some(SelectionPart {
            contig,
            start,
            intervals,
            selection,
        })
    }
}
fn runtime_check(budget: Option<u64>) -> Result<(), DatasetError> {
    crate::core::governor::checkpoint().map_err(EvidenceError::from)?;
    if let Some(budget) = budget {
        let needed = crate::util::rss::peak_rss_bytes();
        if needed > budget {
            return Err(EvidenceError::Core(crate::core::CoreError::BudgetExceeded {
                needed,
                budget,
            })
            .into());
        }
    }
    Ok(())
}

fn selection_bytes(selection: &EvidenceSelection, interval_count: usize) -> u64 {
    match selection {
        EvidenceSelection::Sites(sites) => (sites.len() as u64).saturating_mul(256),
        EvidenceSelection::Intervals(intervals) => (intervals.len() as u64).saturating_mul(64),
        EvidenceSelection::WholeGenome => (interval_count as u64).saturating_mul(64),
    }
}
fn checked_mul(left: u64, right: u64) -> Result<u64, DatasetError> {
    left.checked_mul(right)
        .ok_or_else(|| incompatible("reuse memory envelope overflow"))
}
fn incompatible(message: &str) -> DatasetError {
    DatasetError::Incompatible(message.into())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dataset::{
        publish_evidence_dataset, run_dataset_with_snapshot, DatasetOptions, DatasetReadLimits,
        DescriptorLimits,
    };
    use crate::evidence::{EvidenceArrowWriter, EvidenceRequest};
    use rust_htslib::bam::{
        self,
        record::{Cigar, CigarString},
    };
    use std::{
        fs,
        path::PathBuf,
        sync::atomic::{AtomicU64, Ordering},
    };
    static NEXT: AtomicU64 = AtomicU64::new(0);
    struct Fixture {
        root: PathBuf,
        request: EvidenceRequest,
        session: VerifiedInputSession,
    }
    impl Fixture {
        fn new() -> Self {
            let root = std::env::temp_dir().join(format!(
                "rosalind-reuse-{}-{}",
                std::process::id(),
                NEXT.fetch_add(1, Ordering::Relaxed)
            ));
            fs::create_dir_all(&root).unwrap();
            let reference = root.join("reference.fa");
            fs::write(&reference, format!(">chr1\n{}\n", "A".repeat(40000))).unwrap();
            fs::write(
                root.join("reference.fa.fai"),
                "chr1\t40000\t6\t40000\t40001\n",
            )
            .unwrap();
            let bam_path = root.join("reads.bam");
            let mut header = bam::Header::new();
            header.push_record(
                bam::header::HeaderRecord::new(b"SQ")
                    .push_tag(b"SN", "chr1")
                    .push_tag(b"LN", 40000),
            );
            let mut writer = bam::Writer::from_path(&bam_path, &header, bam::Format::Bam).unwrap();
            for (index, position) in [0, 16380, 18000, 32760].into_iter().enumerate() {
                let mut record = bam::Record::new();
                record.set(
                    format!("read{index}").as_bytes(),
                    Some(&CigarString(vec![Cigar::Match(12)])),
                    b"AACAGAAAACTA",
                    &[30; 12],
                );
                record.set_tid(0);
                record.set_pos(position);
                record.set_mapq(60);
                record.set_flags(if index % 2 == 0 { 0 } else { 16 });
                writer.write(&record).unwrap();
            }
            drop(writer);
            bam::index::build(&bam_path, None::<&PathBuf>, bam::index::Type::Bai, 1).unwrap();
            let session = VerifiedInputSession::open(vec![
                ("alignments".into(), bam_path.clone()),
                ("alignment-index".into(), root.join("reads.bam.bai")),
                ("reference".into(), reference.clone()),
                ("reference-fai".into(), root.join("reference.fa.fai")),
            ])
            .unwrap();
            let mut request = EvidenceRequest::new(bam_path, reference);
            request.alignment_index = Some(root.join("reads.bam.bai"));
            request.reference_fai = Some(root.join("reference.fa.fai"));
            request.fields = EvidenceFields::ALL_SUPPORTED;
            Self {
                root,
                request,
                session,
            }
        }
        fn engine(
            &self,
            selection: EvidenceSelection,
            fields: EvidenceFields,
            width: u32,
        ) -> EvidenceEngine {
            let mut request = self.request.clone();
            request.selection = selection;
            request.fields = fields;
            request.execution.max_microtile_bases = width;
            EvidenceEngine::open(request).unwrap()
        }
        fn source(
            &self,
            selection: EvidenceSelection,
            fields: EvidenceFields,
        ) -> VerifiedEvidenceDataset {
            let mut engine = self.engine(selection, fields, 256);
            let namespace = self.session.dataset_namespace(&engine).unwrap();
            let outcome = run_dataset_with_snapshot(
                &mut engine,
                &namespace,
                &DatasetOptions {
                    cache_dir: self.root.join("cache"),
                    resume: false,
                    workers: 1,
                },
                &mut EvidenceCallback::with_fields(|_: &EvidenceBatch| Ok(()), 0, fields),
                self.session.snapshot(),
            )
            .unwrap();
            let manifest = publish_evidence_dataset(
                &engine,
                &outcome,
                &self.session,
                DescriptorLimits::default(),
            )
            .unwrap();
            VerifiedEvidenceDataset::open(manifest, DatasetReadLimits::default()).unwrap()
        }
    }
    impl Drop for Fixture {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.root);
        }
    }
    fn intervals(ranges: &[(u32, u32)]) -> EvidenceSelection {
        EvidenceSelection::Intervals(
            ranges
                .iter()
                .map(|&(start, end)| GenomicInterval {
                    contig: 0,
                    start,
                    end,
                })
                .collect(),
        )
    }
    fn sites(positions: &[u32], alternate: u8) -> EvidenceSelection {
        EvidenceSelection::Sites(
            positions
                .iter()
                .map(|&position| SnvSite {
                    contig: 0,
                    position,
                    reference: b'A',
                    alternates: vec![alternate],
                })
                .collect(),
        )
    }
    fn fresh(engine: &mut EvidenceEngine) -> Vec<u8> {
        let mut writer = EvidenceArrowWriter::with_fields(Vec::new(), engine.request().fields);
        engine.run(&mut writer).unwrap();
        writer.into_inner().unwrap()
    }

    #[test]
    fn partial_overlap_and_projection_equal_fresh_across_microtiles() {
        let f = Fixture::new();
        let mut source = f.source(
            intervals(&[
                (0, 4),
                (6, 9),
                (1024, 1030),
                (16380, 16386),
                (18000, 18004),
                (32760, 32766),
            ]),
            EvidenceFields::ALL_SUPPORTED,
        );
        for fields in [
            EvidenceFields::DEPTHS.union(EvidenceFields::ALLELE_QUALITY),
            EvidenceFields::ALL_SUPPORTED,
        ] {
            for width in [1, 256, 16384] {
                let query = intervals(&[(0, 2051), (16378, 16400), (18000, 18012), (32759, 32775)]);
                let mut engine = f.engine(query, fields, width);
                let expected = fresh(&mut engine);
                let original_digest = engine.selection_digest();
                let mut writer = EvidenceArrowWriter::with_fields(Vec::new(), fields);
                let outcome = run_reusing_dataset(
                    &mut engine,
                    &f.session,
                    &mut source,
                    &mut writer,
                    &ReuseOptions::default(),
                )
                .unwrap();
                assert_eq!(writer.into_inner().unwrap(), expected);
                assert_eq!(engine.selection_digest(), original_digest);
                assert_eq!(outcome.reused_loci, 29);
                assert_eq!(outcome.computed_loci, 2072);
                assert_eq!(outcome.stats.emitted_loci, 2101);
                assert!(outcome.stats.microtiles > 0);
                assert!(outcome.stats.record_visits > 0);
            }
        }
    }

    #[test]
    fn subset_changes_query_alt_and_reuses_zero_rows_without_alignment_visits() {
        let f = Fixture::new();
        let mut source = f.source(
            sites(&[1, 2, 5, 1024, 16384, 18000, 33000], b'C'),
            EvidenceFields::ALL_SUPPORTED,
        );
        let mut engine = f.engine(sites(&[2, 5, 16384, 33000], b'T'), EvidenceFields::ALL, 1);
        let expected = fresh(&mut engine);
        let mut writer = EvidenceArrowWriter::with_fields(Vec::new(), EvidenceFields::ALL);
        let plan = plan_reuse(
            &mut engine,
            &f.session,
            &source,
            &writer,
            &ReuseOptions::default(),
        )
        .unwrap();
        assert_eq!(plan.native_bytes, 0);
        let outcome = run_reusing_dataset(
            &mut engine,
            &f.session,
            &mut source,
            &mut writer,
            &ReuseOptions::default(),
        )
        .unwrap();
        assert_eq!(writer.into_inner().unwrap(), expected);
        assert_eq!(
            (
                outcome.reused_loci,
                outcome.computed_loci,
                outcome.stats.record_visits
            ),
            (4, 0, 0)
        );
        assert_eq!(outcome.stats.microtiles, 0);
    }

    #[test]
    fn missing_only_and_empty_selection_finish_once_and_match_fresh() {
        let f = Fixture::new();
        let mut source = f.source(intervals(&[(0, 12)]), EvidenceFields::DEPTHS);
        for query in [intervals(&[(18000, 18012)]), intervals(&[])] {
            let mut engine = f.engine(query, EvidenceFields::DEPTHS, 4);
            let expected = fresh(&mut engine);
            struct Consumer {
                writer: EvidenceArrowWriter<Vec<u8>>,
                finishes: usize,
            }
            impl EvidenceAnalyzer for Consumer {
                fn requirements(&self) -> EvidenceRequirements {
                    self.writer.requirements()
                }
                fn on_batch(&mut self, batch: &EvidenceBatch) -> Result<(), EvidenceError> {
                    self.writer.on_batch(batch)
                }
                fn finish(&mut self) -> Result<(), EvidenceError> {
                    self.finishes += 1;
                    self.writer.finish()
                }
            }
            let mut consumer = Consumer {
                writer: EvidenceArrowWriter::with_fields(Vec::new(), EvidenceFields::DEPTHS),
                finishes: 0,
            };
            let outcome = run_reusing_dataset(
                &mut engine,
                &f.session,
                &mut source,
                &mut consumer,
                &ReuseOptions::default(),
            )
            .unwrap();
            assert_eq!(consumer.finishes, 1);
            assert_eq!(consumer.writer.into_inner().unwrap(), expected);
            assert_eq!(outcome.reused_loci, 0);
        }
    }

    #[test]
    fn scientific_mismatch_missing_fields_and_parallel_options_are_refused() {
        let f = Fixture::new();
        let source = f.source(intervals(&[(0, 12)]), EvidenceFields::DEPTHS);
        let consumer =
            EvidenceCallback::with_fields(|_: &EvidenceBatch| Ok(()), 0, EvidenceFields::DEPTHS);
        let mut engine = f.engine(intervals(&[(0, 12)]), EvidenceFields::DEPTHS, 4);
        assert!(plan_reuse(
            &mut engine,
            &f.session,
            &source,
            &consumer,
            &ReuseOptions { workers: 2 }
        )
        .is_err());
        let mut request = engine.request().clone();
        request.profile.min_mapq = 30;
        let mut filtered = EvidenceEngine::open(request).unwrap();
        assert!(plan_reuse(
            &mut filtered,
            &f.session,
            &source,
            &consumer,
            &ReuseOptions::default()
        )
        .unwrap_err()
        .to_string()
        .contains("identities, filters"));
        let mut absent = f.engine(
            intervals(&[(0, 12)]),
            EvidenceFields::DEPTHS.union(EvidenceFields::ALLELE_QUALITY),
            4,
        );
        assert!(plan_reuse(
            &mut absent,
            &f.session,
            &source,
            &consumer,
            &ReuseOptions::default()
        )
        .unwrap_err()
        .to_string()
        .contains("missing requested"));
        let other = Fixture::new();
        let unrelated =
            VerifiedInputSession::open(vec![("alignments".into(), other.root.join("reads.bam"))])
                .unwrap();
        assert!(plan_reuse(
            &mut engine,
            &unrelated,
            &source,
            &consumer,
            &ReuseOptions::default()
        )
        .is_err());
        // Source content, not merely file names or public digest strings, must
        // match. Reopen a complete session after changing the second reference.
        fs::write(
            other.root.join("reference.fa"),
            format!(">chr1\n{}\n", "C".repeat(40000)),
        )
        .unwrap();
        let changed_session = VerifiedInputSession::open(vec![
            ("alignments".into(), other.root.join("reads.bam")),
            ("alignment-index".into(), other.root.join("reads.bam.bai")),
            ("reference".into(), other.root.join("reference.fa")),
            ("reference-fai".into(), other.root.join("reference.fa.fai")),
        ])
        .unwrap();
        let mut changed_engine = other.engine(intervals(&[(0, 12)]), EvidenceFields::DEPTHS, 4);
        assert!(plan_reuse(
            &mut changed_engine,
            &changed_session,
            &source,
            &consumer,
            &ReuseOptions::default()
        )
        .unwrap_err()
        .to_string()
        .contains("identities, filters"));
        let mut request = engine.request().clone();
        request.execution.memory_budget_bytes =
            Some(crate::util::rss::peak_rss_bytes() + (128 << 20));
        let oversized = EvidenceCallback::with_fields(
            |_: &EvidenceBatch| Ok(()),
            u64::MAX,
            EvidenceFields::DEPTHS,
        );
        let mut bounded = EvidenceEngine::open(request).unwrap();
        assert!(plan_reuse(
            &mut bounded,
            &f.session,
            &source,
            &oversized,
            &ReuseOptions::default()
        )
        .is_err());
    }

    #[test]
    fn consumer_failure_restores_engine_selection_and_does_not_finish() {
        let f = Fixture::new();
        let mut source = f.source(intervals(&[(0, 3)]), EvidenceFields::DEPTHS);
        let mut engine = f.engine(
            intervals(&[(0, 12), (18000, 18012)]),
            EvidenceFields::DEPTHS,
            4,
        );
        let before = engine.selection_digest();
        struct Failing {
            finishes: usize,
        }
        impl EvidenceAnalyzer for Failing {
            fn requirements(&self) -> EvidenceRequirements {
                EvidenceRequirements {
                    fields: EvidenceFields::DEPTHS,
                    requires_reference: false,
                    context_bases: 0,
                    retained_bytes: Some(0),
                }
            }
            fn on_batch(&mut self, _: &EvidenceBatch) -> Result<(), EvidenceError> {
                Err(EvidenceError::Analyzer("intentional test failure".into()))
            }
            fn finish(&mut self) -> Result<(), EvidenceError> {
                self.finishes += 1;
                Ok(())
            }
        }
        let mut consumer = Failing { finishes: 0 };
        assert!(run_reusing_dataset(
            &mut engine,
            &f.session,
            &mut source,
            &mut consumer,
            &ReuseOptions::default()
        )
        .unwrap_err()
        .to_string()
        .contains("intentional test failure"));
        assert_eq!(consumer.finishes, 0);
        assert_eq!(engine.selection_digest(), before);
    }
    #[test]
    fn joint_plan_counts_source_decoder_and_inherits_source_budget() {
        let f = Fixture::new();
        let source = f.source(
            intervals(&[(0, 3), (16380, 16384)]),
            EvidenceFields::ALL_SUPPORTED,
        );
        let manifest = PathBuf::from(&source.source_hashes()[0].path);
        drop(source);
        let budget = crate::util::rss::peak_rss_bytes() + (128 << 20);
        let mut source = VerifiedEvidenceDataset::open(
            manifest,
            DatasetReadLimits {
                memory_budget_bytes: Some(budget),
                ..DatasetReadLimits::default()
            },
        )
        .unwrap();
        let query = intervals(&[(0, 2051), (16380, 16400)]);
        let mut engine = f.engine(query, EvidenceFields::DEPTHS, 256);
        struct Unknown;
        impl EvidenceAnalyzer for Unknown {
            fn requirements(&self) -> EvidenceRequirements {
                EvidenceRequirements {
                    fields: EvidenceFields::DEPTHS,
                    requires_reference: false,
                    context_bases: 0,
                    retained_bytes: None,
                }
            }
            fn on_batch(&mut self, _: &EvidenceBatch) -> Result<(), EvidenceError> {
                Ok(())
            }
        }
        assert!(plan_reuse(
            &mut engine,
            &f.session,
            &source,
            &Unknown,
            &ReuseOptions::default()
        )
        .unwrap_err()
        .to_string()
        .contains("declared consumer memory"));
        let mut writer = EvidenceArrowWriter::with_fields(Vec::new(), EvidenceFields::DEPTHS);
        let plan = plan_reuse(
            &mut engine,
            &f.session,
            &source,
            &writer,
            &ReuseOptions::default(),
        )
        .unwrap();
        assert_eq!(
            plan.source_decoder_bytes,
            crate::evidence::evidence_reader_memory_bytes(EvidenceFields::ALL_SUPPORTED)
        );
        assert_eq!(plan.memory_budget_bytes, Some(budget));
        assert_eq!(plan.retained_loci, 2055);
        assert!(plan.native_bytes > 0 && plan.merge_bytes > 0);
        assert!(plan.predicted_peak_rss_bytes <= budget);
        let outcome = run_reusing_dataset(
            &mut engine,
            &f.session,
            &mut source,
            &mut writer,
            &ReuseOptions::default(),
        )
        .unwrap();
        assert_eq!((outcome.reused_loci, outcome.computed_loci), (7, 2064));
    }
}
