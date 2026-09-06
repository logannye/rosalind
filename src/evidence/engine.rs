use super::*;
use crate::core::ContigSet;
use crate::selection::GenomicInterval;
use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::{self, Read};
use std::collections::BTreeMap;
use std::path::PathBuf;
use std::sync::Arc;

const DECODER_SLACK_BYTES: u64 = 8 << 20;

/// Conservative model of retained summaries plus declared decoder/analyzer space.
/// Decoder allocations are governed cooperatively; CRAM internal workspaces and
/// htslib record allocation are not a hard preallocation guarantee.
#[derive(Debug, Clone)]
pub struct EvidencePlan {
    /// Versioned resource-estimator identity.
    pub model_id: &'static str,
    /// Measured process high-water RSS before planned allocations.
    pub baseline_rss_bytes: u64,
    /// Reserved decoder, I/O, allocator and analyzer bytes.
    pub fixed_bytes: u64,
    /// Conservative bytes reserved per live locus summary.
    pub bytes_per_locus: u64,
    /// Chosen execution window width in reference bases.
    pub microtile_bases: u32,
    /// Width of canonical ownership tiles, independent of execution budget.
    pub canonical_tile_bases: u32,
    /// Baseline plus the admitted decoder, summary and analyzer reservation.
    pub predicted_peak_rss_bytes: u64,
    /// Number of unique requested reference positions.
    pub selected_loci: u64,
    /// Reserved additional memory declared by the consumer.
    pub analyzer_bytes: u64,
    /// Conservative retained canonical selection and site annotation bytes.
    pub selection_bytes: u64,
    /// Conservative retained sample scope, read-group header, and worker setup bytes.
    pub sample_scope_bytes: u64,
}
#[derive(Debug, Default, Clone, PartialEq, Eq)]
/// Aggregate execution measurements; fetch-level visits may include the same read in multiple windows.
pub struct EvidenceRunStats {
    /// Number of exact rows emitted, including zero-depth positions.
    pub emitted_loci: u64,
    /// Fetch-level records, deliberately not a count of unique reads.
    pub record_visits: u64,
    /// Fetch-level records excluded by read-level filters; not unique reads.
    pub filtered_record_visits: u64,
    /// Fetch-level records belonging to another declared sample, before locus filters.
    pub sample_filtered_record_visits: u64,
    /// Number of indexed execution windows processed.
    pub microtiles: u64,
    /// Largest decoded BAM-record payload observed in bytes.
    pub max_record_bytes: u64,
    /// Largest decoded sequence length observed.
    pub max_read_length: u64,
    /// Measured whole-process high-water resident memory at completion.
    pub peak_rss_bytes: u64,
}

#[derive(Debug)]
/// Reusable exact evidence executor with one indexed decoder and shared reference windows.
pub struct EvidenceEngine {
    request: EvidenceRequest,
    reference: EvidenceReference,
    reader: bam::IndexedReader,
    contigs: ContigSet,
    tids: Vec<u32>,
    intervals: Vec<GenomicInterval>,
    sites: BTreeMap<(u32, u32), SnvSite>,
    sample_scope: Arc<EvidenceSampleScope>,
    plan: EvidencePlan,
}
impl EvidenceEngine {
    /// Open and validate local inputs, retaining only indexed reference access and metadata.
    pub fn open(mut request: EvidenceRequest) -> Result<Self, EvidenceError> {
        if request.fields != EvidenceFields::ALL {
            return Err(EvidenceError::InvalidRequest(
                "evidence schema v1 emits all fields; physical field projection is not supported"
                    .into(),
            ));
        }
        if request.execution.max_microtile_bases == 0
            || request.execution.max_read_len == 0
            || request.execution.max_record_bytes == 0
        {
            return Err(EvidenceError::InvalidRequest(
                "execution widths and record/read envelopes must be positive".into(),
            ));
        }
        if request.profile.min_base_quality > 93 || request.profile.min_mapq == 255 {
            return Err(EvidenceError::InvalidRequest(
                "quality thresholds must be BQ<=93 and MAPQ<=254".into(),
            ));
        }
        let opened = match &request.alignment_index {
            Some(index) => bam::IndexedReader::from_path_and_index(&request.alignments, index),
            None => bam::IndexedReader::from_path(&request.alignments),
        };
        let mut reader = opened.map_err(|error| {
            EvidenceError::InvalidInput(format!(
                "open indexed BAM/CRAM {} (BAI/CSI/CRAI required): {error}",
                request.alignments.display()
            ))
        })?;
        let sample_scope = Arc::new(EvidenceSampleScope::resolve(
            reader.header(),
            &request.sample_selection,
        )?);
        let sample_scope_bytes = sample_scope
            .memory_bytes()
            .saturating_add((reader.header().as_bytes().len() as u64).saturating_mul(2));
        let alignment_contigs = header_contigs(reader.header())?;
        let mut reference = match &request.reference {
            Some(path) => EvidenceReference::open_with_fai(path, request.reference_fai.as_deref())?,
            None => EvidenceReference::unavailable(alignment_contigs.clone()),
        };
        let is_cram = alignment_is_cram(&request.alignments)?;
        if is_cram {
            let cram_fasta: PathBuf = request.cram_reference.clone()
                .or_else(|| reference.fasta_path().map(PathBuf::from))
                .ok_or_else(|| EvidenceError::InvalidRequest("CRAM requires an explicit local FASTA with .fai; no remote reference lookup is permitted".into()))?;
            if let Some(index) = request
                .cram_reference_fai
                .as_ref()
                .or(request.reference_fai.as_ref())
            {
                let adjacent = PathBuf::from(format!("{}.fai", cram_fasta.display()));
                if std::fs::canonicalize(index)?
                    != std::fs::canonicalize(&adjacent).map_err(|_| {
                        EvidenceError::InvalidRequest(
                            "CRAM replay requires staging its verified FASTA and FAI adjacently"
                                .into(),
                        )
                    })?
                {
                    return Err(EvidenceError::InvalidRequest("CRAM FAI must be adjacent to its FASTA; stage verified reference inputs before replay".into()));
                }
            }
            let cram_reference = EvidenceReference::open(&cram_fasta)?;
            if cram_reference.fasta_path().is_none() {
                return Err(EvidenceError::InvalidRequest(
                    "CRAM decoder reference must be a local FASTA, not a pack".into(),
                ));
            }
            validate_dictionary(cram_reference.contigs(), &alignment_contigs)?;
            reader.set_reference(&cram_fasta).map_err(|error| {
                EvidenceError::InvalidInput(format!("set local CRAM reference: {error}"))
            })?;
            // An explicit decode reference is also the analysis sequence when no
            // separate reference was requested.
            if request.reference.is_none() {
                reference = cram_reference;
            }
        }
        let contigs = reference.contigs().clone();
        validate_dictionary(&contigs, &alignment_contigs)?;
        let tids = contigs
            .iter()
            .map(|contig| {
                reader.header().tid(contig.name.as_bytes()).ok_or_else(|| {
                    EvidenceError::InvalidInput(format!(
                        "reference contig {} absent from alignment header",
                        contig.name
                    ))
                })
            })
            .collect::<Result<Vec<_>, _>>()?;
        let (intervals, sites) = request.selection.normalize(&contigs)?;
        if !sites.is_empty() && !reference.has_sequence() {
            return Err(EvidenceError::InvalidRequest(
                "SNV selection requires reference sequence for REF validation".into(),
            ));
        }
        for site in sites.values() {
            let base = reference.read_window(site.contig, site.position, site.position + 1)?[0];
            if base != site.reference {
                return Err(EvidenceError::InvalidInput(format!(
                    "VCF REF mismatch at {}:{}: requested {}, reference {}",
                    contigs.by_id(site.contig).unwrap().name,
                    site.position + 1,
                    char::from(site.reference),
                    char::from(base)
                )));
            }
        }
        let selected_loci = intervals.iter().try_fold(0u64, |sum, interval| {
            sum.checked_add(interval.len() as u64)
                .ok_or(EvidenceError::CounterOverflow)
        })?;
        request.selection = canonical_selection(&request.selection, &intervals, &sites);
        let plan = make_plan(
            &request.execution,
            selected_loci,
            request.execution.analyzer_bytes,
            crate::util::rss::peak_rss_bytes(),
            selection_memory_bytes(&intervals, &sites),
            sample_scope_bytes,
        )?;
        Ok(Self {
            request,
            reference,
            reader,
            contigs,
            tids,
            intervals,
            sites,
            sample_scope,
            plan,
        })
    }
    /// Read the current resource plan; include the consumer before final admission.
    pub fn plan(&self) -> &EvidencePlan {
        &self.plan
    }
    /// Read the canonical reference contig dictionary.
    pub fn contigs(&self) -> &ContigSet {
        &self.contigs
    }
    /// Read the request with its normalized scientific selection.
    pub fn request(&self) -> &EvidenceRequest {
        &self.request
    }
    /// Resolved sample identity and assignment policy for every emitted row.
    pub fn sample_scope(&self) -> &EvidenceSampleScope {
        &self.sample_scope
    }
    /// Make a worker factory that shares already-open reference storage. The
    /// caller must preserve input immutability for the lifetime of the run.
    pub fn worker_factory(&self) -> EvidenceWorkerFactory {
        EvidenceWorkerFactory {
            // A factory retains no original selection: workers receive only
            // their canonical partition, so cloning cannot multiply site lists.
            request: EvidenceRequest {
                alignments: self.request.alignments.clone(),
                alignment_index: self.request.alignment_index.clone(),
                reference: self.request.reference.clone(),
                reference_fai: self.request.reference_fai.clone(),
                cram_reference: self.request.cram_reference.clone(),
                cram_reference_fai: self.request.cram_reference_fai.clone(),
                selection: EvidenceSelection::WholeGenome,
                sample_selection: self.request.sample_selection.clone(),
                fields: self.request.fields,
                profile: self.request.profile.clone(),
                execution: self.request.execution.clone(),
            },
            reference: self.reference.clone(),
            sample_scope: self.sample_scope.clone(),
        }
    }
    /// Stable, explicit coordinate/annotation digest for scientific cache keys.
    pub fn selection_digest(&self) -> String {
        let mut hash = blake3::Hasher::new();
        hash.update(b"rosalind-evidence-selection-v1\0");
        for interval in &self.intervals {
            hash.update(&interval.contig.to_le_bytes());
            hash.update(&interval.start.to_le_bytes());
            hash.update(&interval.end.to_le_bytes());
        }
        hash.update(b"\0sites\0");
        for site in self.sites.values() {
            hash.update(&site.contig.to_le_bytes());
            hash.update(&site.position.to_le_bytes());
            hash.update(&[site.reference, site.alternates.len() as u8]);
            hash.update(&site.alternates);
        }
        hash.finalize().to_hex().to_string()
    }
    /// Read the ordered disjoint intervals used for extraction.
    pub fn intervals(&self) -> &[GenomicInterval] {
        &self.intervals
    }
    /// Resolve a VCF/BED selection after inspecting this engine's dictionary,
    /// without opening the alignment or revalidating the reference a second time.
    pub fn set_selection(&mut self, selection: EvidenceSelection) -> Result<(), EvidenceError> {
        let (intervals, sites) = selection.normalize(&self.contigs)?;
        if !sites.is_empty() && !self.reference.has_sequence() {
            return Err(EvidenceError::InvalidRequest(
                "SNV selection requires reference sequence".into(),
            ));
        }
        for site in sites.values() {
            let base = self
                .reference
                .read_window(site.contig, site.position, site.position + 1)?[0];
            if base != site.reference {
                return Err(EvidenceError::InvalidInput(format!(
                    "VCF REF mismatch at {}:{}",
                    self.contigs.by_id(site.contig).unwrap().name,
                    site.position + 1
                )));
            }
        }
        let selected = intervals.iter().try_fold(0u64, |sum, interval| {
            sum.checked_add(interval.len() as u64)
                .ok_or(EvidenceError::CounterOverflow)
        })?;
        let plan = make_plan(
            &self.request.execution,
            selected,
            self.plan.analyzer_bytes,
            self.plan.baseline_rss_bytes,
            selection_memory_bytes(&intervals, &sites),
            self.plan.sample_scope_bytes,
        )?;
        self.request.selection = canonical_selection(&selection, &intervals, &sites);
        self.intervals = intervals;
        self.sites = sites;
        self.plan = plan;
        Ok(())
    }
    /// Include the consumer's declared peak retained bytes before opening any
    /// output. Call this before displaying or admitting the final plan.
    pub fn plan_for_analyzer(
        &mut self,
        analyzer: &dyn EvidenceAnalyzer,
    ) -> Result<&EvidencePlan, EvidenceError> {
        let requirements = analyzer.requirements();
        if !self.request.fields.contains(requirements.fields) {
            return Err(EvidenceError::InvalidRequest(
                "analyzer required fields are absent from the evidence request".into(),
            ));
        }
        if requirements.context_bases > 0 {
            return Err(EvidenceError::InvalidRequest(
                "evidence schema v1 does not provide flanking reference context".into(),
            ));
        }
        if requirements.requires_reference && !self.reference.has_sequence() {
            return Err(EvidenceError::InvalidRequest(
                "analyzer requires actual reference bases".into(),
            ));
        }
        let bytes = match requirements.retained_bytes {
            Some(bytes) => bytes.max(self.request.execution.analyzer_bytes),
            None if self.request.execution.memory_budget_bytes.is_some() => {
                return Err(EvidenceError::InvalidRequest(
                    "budgeted evidence analyzer requires a declared memory bound".into(),
                ))
            }
            None => self.request.execution.analyzer_bytes,
        };
        self.plan = make_plan(
            &self.request.execution,
            self.plan.selected_loci,
            bytes,
            self.plan.baseline_rss_bytes,
            self.plan.selection_bytes,
            self.plan.sample_scope_bytes,
        )?;
        Ok(&self.plan)
    }
    /// Execute in canonical coordinate order. Safe to execute once per engine;
    /// subsequent calls rerun the same request from indexed input.
    pub fn run(
        &mut self,
        analyzer: &mut dyn EvidenceAnalyzer,
    ) -> Result<EvidenceRunStats, EvidenceError> {
        self.plan_for_analyzer(analyzer)?;
        let mut stats = EvidenceRunStats::default();
        let mut record = bam::Record::new();
        for interval in &self.intervals {
            let mut start = interval.start;
            while start < interval.end {
                crate::core::governor::checkpoint()?;
                let canonical_start = start / CANONICAL_TILE_BASES * CANONICAL_TILE_BASES;
                let canonical_end = canonical_start.saturating_add(CANONICAL_TILE_BASES);
                let end = interval
                    .end
                    .min(canonical_end)
                    .min(start.saturating_add(self.plan.microtile_bases));
                let sequence = self.reference.read_window(interval.contig, start, end)?;
                let mut rows = Vec::with_capacity((end - start) as usize);
                for (offset, reference) in sequence.into_iter().enumerate() {
                    let position = start + offset as u32;
                    let requested_alts = self
                        .sites
                        .get(&(interval.contig, position))
                        .map_or_else(Vec::new, |site| site.alternates.clone());
                    rows.push(EvidenceRow {
                        position,
                        reference,
                        requested_alts,
                        ..EvidenceRow::default()
                    });
                }
                self.reader
                    .fetch((
                        self.tids[interval.contig as usize],
                        start as i64,
                        end as i64,
                    ))
                    .map_err(|error| {
                        EvidenceError::InvalidInput(format!("indexed interval fetch: {error}"))
                    })?;
                let mut previous = None;
                while let Some(result) = self.reader.read(&mut record) {
                    result.map_err(|error| {
                        EvidenceError::InvalidInput(format!("decode alignment: {error}"))
                    })?;
                    crate::core::governor::checkpoint()?;
                    checked_add(&mut stats.record_visits, 1)?;
                    let record_bytes = record.inner().l_data.max(0) as usize;
                    stats.max_record_bytes = stats.max_record_bytes.max(record_bytes as u64);
                    stats.max_read_length = stats.max_read_length.max(record.seq_len() as u64);
                    if record_bytes > self.request.execution.max_record_bytes {
                        return Err(EvidenceError::RecordLimit(format!(
                            "record has {record_bytes} bytes, maximum is {}",
                            self.request.execution.max_record_bytes
                        )));
                    }
                    if record.seq_len() > self.request.execution.max_read_len {
                        return Err(EvidenceError::RecordLimit(format!(
                            "read has {} bases, maximum is {}",
                            record.seq_len(),
                            self.request.execution.max_read_len
                        )));
                    }
                    if record.is_unmapped() || record.tid() < 0 {
                        continue;
                    }
                    let coordinate = (record.tid(), record.pos());
                    if previous.is_some_and(|previous| coordinate < previous) {
                        return Err(EvidenceError::InvalidInput(
                            "indexed alignment records are not coordinate sorted".into(),
                        ));
                    }
                    previous = Some(coordinate);
                    if record.tid() != self.tids[interval.contig as usize] as i32
                        || record.pos() < 0
                    {
                        return Err(EvidenceError::InvalidInput(
                            "indexed fetch returned an invalid contig/coordinate".into(),
                        ));
                    }
                    if !self.sample_scope.includes_record(&record)? {
                        checked_add(&mut stats.sample_filtered_record_visits, 1)?;
                        continue;
                    }
                    let rejected = read_filter(&record, &self.request.profile);
                    if rejected.is_some() {
                        checked_add(&mut stats.filtered_record_visits, 1)?;
                    }
                    accumulate_record(
                        &record,
                        &mut rows,
                        start,
                        end,
                        self.contigs.by_id(interval.contig).unwrap().length,
                        &self.request.profile,
                        rejected,
                        self.reference.has_sequence(),
                    )?;
                }
                checked_add(&mut stats.microtiles, 1)?;
                checked_add(&mut stats.emitted_loci, rows.len() as u64)?;
                let batch = EvidenceBatch {
                    contig_id: interval.contig,
                    contig: self
                        .contigs
                        .by_id(interval.contig)
                        .unwrap()
                        .name
                        .to_string(),
                    canonical_tile_start: canonical_start,
                    rows,
                };
                analyzer.on_batch(&batch)?;
                self.check_budget()?;
                start = end;
            }
        }
        analyzer.finish()?;
        self.check_budget()?;
        stats.peak_rss_bytes = crate::util::rss::peak_rss_bytes();
        Ok(stats)
    }
    fn check_budget(&self) -> Result<(), EvidenceError> {
        crate::core::governor::checkpoint()?;
        if let Some(budget) = self.request.execution.memory_budget_bytes {
            let needed = crate::util::rss::peak_rss_bytes();
            if needed > budget {
                return Err(crate::core::CoreError::BudgetExceeded { needed, budget }.into());
            }
        }
        Ok(())
    }
}

fn make_plan(
    execution: &EvidenceExecution,
    selected_loci: u64,
    analyzer_bytes: u64,
    baseline: u64,
    selection_bytes: u64,
    sample_scope_bytes: u64,
) -> Result<EvidencePlan, EvidenceError> {
    // Capacity, site annotation and one-byte reference window included. One
    // reusable decoded BAM record and one CIGAR view each fit this envelope.
    let bytes_per_locus = std::mem::size_of::<EvidenceRow>() as u64 + 16;
    let fixed = DECODER_SLACK_BYTES
        .checked_add(
            (execution.max_record_bytes as u64)
                .checked_mul(3)
                .ok_or(EvidenceError::CounterOverflow)?,
        )
        .and_then(|bytes| bytes.checked_add(analyzer_bytes))
        .and_then(|bytes| bytes.checked_add(selection_bytes))
        .and_then(|bytes| bytes.checked_add(sample_scope_bytes))
        .ok_or(EvidenceError::CounterOverflow)?;
    let maximum = execution.max_microtile_bases.clamp(1, CANONICAL_TILE_BASES) as u64;
    let microtile = if let Some(budget) = execution.memory_budget_bytes {
        let minimum = baseline
            .saturating_add(fixed)
            .saturating_add(bytes_per_locus);
        if budget < minimum {
            return Err(EvidenceError::Refused {
                needed: minimum,
                budget,
            });
        }
        maximum.min((budget - baseline - fixed) / bytes_per_locus)
    } else {
        maximum
    };
    Ok(EvidencePlan {
        model_id: "exact-summary-tiles-v2",
        baseline_rss_bytes: baseline,
        fixed_bytes: fixed,
        bytes_per_locus,
        microtile_bases: microtile as u32,
        canonical_tile_bases: CANONICAL_TILE_BASES,
        predicted_peak_rss_bytes: baseline
            .saturating_add(fixed)
            .saturating_add(microtile * bytes_per_locus),
        selected_loci,
        analyzer_bytes,
        selection_bytes,
        sample_scope_bytes,
    })
}

fn selection_memory_bytes(
    intervals: &[GenomicInterval],
    sites: &BTreeMap<(u32, u32), SnvSite>,
) -> u64 {
    // Both the normalized request and execution indexes are retained. Reserve
    // B-tree/node/ALT allocation overhead conservatively for each SNV. This is
    // also charged after set_selection without resetting the original baseline.
    (intervals.len() as u64)
        .saturating_mul(2 * std::mem::size_of::<GenomicInterval>() as u64)
        .saturating_add((sites.len() as u64).saturating_mul(256))
}

fn validate_dictionary(reference: &ContigSet, alignment: &ContigSet) -> Result<(), EvidenceError> {
    if reference.len() != alignment.len() {
        return Err(EvidenceError::InvalidInput(
            "reference and alignment contig dictionaries differ".into(),
        ));
    }
    for contig in reference.iter() {
        if alignment
            .by_name(&contig.name)
            .is_none_or(|other| other.length != contig.length)
        {
            return Err(EvidenceError::InvalidInput(format!(
                "reference and alignment dictionary disagree for {}",
                contig.name
            )));
        }
    }
    Ok(())
}
#[derive(Clone, Copy)]
enum ReadFilter {
    Secondary,
    Supplementary,
    QcFail,
    Duplicate,
    UnavailableMapq,
    LowMapq,
}
fn read_filter(record: &bam::Record, profile: &EvidenceProfile) -> Option<ReadFilter> {
    if profile.exclude_secondary && record.is_secondary() {
        Some(ReadFilter::Secondary)
    } else if profile.exclude_supplementary && record.is_supplementary() {
        Some(ReadFilter::Supplementary)
    } else if profile.exclude_qc_fail && record.is_quality_check_failed() {
        Some(ReadFilter::QcFail)
    } else if profile.exclude_duplicates && record.is_duplicate() {
        Some(ReadFilter::Duplicate)
    } else if record.mapq() == 255 {
        Some(ReadFilter::UnavailableMapq)
    } else if record.mapq() < profile.min_mapq {
        Some(ReadFilter::LowMapq)
    } else {
        None
    }
}
fn increment_filter(row: &mut EvidenceRow, reason: ReadFilter) -> Result<(), EvidenceError> {
    let counter = match reason {
        ReadFilter::Secondary => &mut row.filters.secondary,
        ReadFilter::Supplementary => &mut row.filters.supplementary,
        ReadFilter::QcFail => &mut row.filters.qc_fail,
        ReadFilter::Duplicate => &mut row.filters.duplicate,
        ReadFilter::UnavailableMapq => &mut row.filters.unavailable_mapq,
        ReadFilter::LowMapq => &mut row.filters.low_mapq,
    };
    checked_add(counter, 1)
}
#[allow(clippy::too_many_arguments)]
fn accumulate_record(
    record: &bam::Record,
    rows: &mut [EvidenceRow],
    start: u32,
    end: u32,
    contig_len: u32,
    profile: &EvidenceProfile,
    rejected: Option<ReadFilter>,
    has_reference: bool,
) -> Result<(), EvidenceError> {
    let mut reference_position = record.pos() as u64;
    let mut query_position = 0usize;
    let sequence = record.seq();
    let qualities = record.qual();
    for operation in record.cigar().iter() {
        match *operation {
            Cigar::Match(length) | Cigar::Equal(length) | Cigar::Diff(length) => {
                let reference_end = reference_position
                    .checked_add(length as u64)
                    .ok_or(EvidenceError::CounterOverflow)?;
                let query_end = query_position
                    .checked_add(length as usize)
                    .ok_or(EvidenceError::CounterOverflow)?;
                if query_end > record.seq_len() || reference_end > contig_len as u64 {
                    return Err(EvidenceError::InvalidInput(
                        "CIGAR exceeds read or reference bounds".into(),
                    ));
                }
                let overlap_start = reference_position.max(start as u64);
                let overlap_end = reference_end.min(end as u64);
                for position in overlap_start..overlap_end {
                    let offset = query_position + (position - reference_position) as usize;
                    let row = &mut rows[(position - start as u64) as usize];
                    checked_add(&mut row.prefilter_depth, 1)?;
                    if let Some(reason) = rejected {
                        increment_filter(row, reason)?;
                        continue;
                    }
                    checked_add(&mut row.aligned_depth, 1)?;
                    let quality = qualities[offset];
                    if quality == 255 {
                        checked_add(&mut row.filters.unavailable_base_quality, 1)?;
                        continue;
                    }
                    if quality > 93 {
                        return Err(EvidenceError::InvalidInput(
                            "base quality outside SAM range 0..93 (or missing 255)".into(),
                        ));
                    }
                    if quality < profile.min_base_quality {
                        checked_add(&mut row.filters.low_base_quality, 1)?;
                        continue;
                    }
                    // SAMv1 §1.4 SEQ: '=' means identical to the reference,
                    // including on reverse-strand records already stored in
                    // forward-reference orientation. It is not an ambiguous base.
                    let base = if sequence[offset] == b'=' {
                        if !has_reference {
                            return Err(EvidenceError::InvalidInput("reference-equality SEQ '=' requires a local analysis reference for callable evidence; supply --reference".into()));
                        }
                        row.reference
                    } else {
                        sequence[offset]
                    };
                    let allele = match base {
                        b'A' => 0,
                        b'C' => 1,
                        b'G' => 2,
                        b'T' => 3,
                        _ => {
                            checked_add(&mut row.filters.ambiguous_base, 1)?;
                            continue;
                        }
                    };
                    checked_add(&mut row.callable_depth, 1)?;
                    checked_add(&mut row.allele_counts[allele], 1)?;
                    checked_add(
                        &mut row.strand_counts[allele][usize::from(record.is_reverse())],
                        1,
                    )?;
                    checked_add(&mut row.base_quality_sum, quality as u64)?;
                    checked_add(&mut row.mapping_quality_sum, record.mapq() as u64)?;
                    let cycle = if record.is_reverse() {
                        record.seq_len() - 1 - offset
                    } else {
                        offset
                    };
                    checked_add(&mut row.read_position_sum, cycle as u64)?;
                    checked_add(&mut row.read_length_sum, record.seq_len() as u64)?;
                    checked_add(&mut row.base_quality_histogram[quality as usize], 1)?;
                    checked_add(
                        &mut row.mapping_quality_histogram[record.mapq() as usize],
                        1,
                    )?;
                }
                reference_position = reference_end;
                query_position = query_end;
            }
            Cigar::Ins(length) | Cigar::SoftClip(length) => {
                query_position = query_position
                    .checked_add(length as usize)
                    .ok_or(EvidenceError::CounterOverflow)?;
                if query_position > record.seq_len() {
                    return Err(EvidenceError::InvalidInput(
                        "CIGAR consumes more bases than read".into(),
                    ));
                }
            }
            Cigar::Del(length) | Cigar::RefSkip(length) => {
                reference_position = reference_position
                    .checked_add(length as u64)
                    .ok_or(EvidenceError::CounterOverflow)?;
                if reference_position > contig_len as u64 {
                    return Err(EvidenceError::InvalidInput(
                        "CIGAR exceeds reference bounds".into(),
                    ));
                }
            }
            Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }
    if query_position != record.seq_len() {
        return Err(EvidenceError::InvalidInput(
            "CIGAR read consumption does not match sequence length".into(),
        ));
    }
    Ok(())
}

fn canonical_selection(
    original: &EvidenceSelection,
    intervals: &[GenomicInterval],
    sites: &BTreeMap<(u32, u32), SnvSite>,
) -> EvidenceSelection {
    match original {
        EvidenceSelection::WholeGenome => EvidenceSelection::WholeGenome,
        EvidenceSelection::Intervals(_) => EvidenceSelection::Intervals(intervals.to_vec()),
        EvidenceSelection::Sites(_) => EvidenceSelection::Sites(sites.values().cloned().collect()),
    }
}

/// Cloneable, thread-safe access to verified reference storage. Each worker owns
/// its own indexed decoder; reference bytes and FASTA file handles are shared.
#[derive(Debug, Clone)]
pub struct EvidenceWorkerFactory {
    request: EvidenceRequest,
    reference: EvidenceReference,
    sample_scope: Arc<EvidenceSampleScope>,
}
impl EvidenceWorkerFactory {
    /// Open a worker for one canonical selection without hashing the reference
    /// again. Dataset immutability is a caller-owned session invariant.
    pub fn open(&self, selection: EvidenceSelection) -> Result<EvidenceEngine, EvidenceError> {
        self.open_with_execution(selection, self.request.execution.clone())
    }
    /// Open with a worker share of an already admitted aggregate execution plan.
    /// A parent governor may govern total RSS while the worker budget is None.
    pub fn open_with_execution(
        &self,
        selection: EvidenceSelection,
        execution: EvidenceExecution,
    ) -> Result<EvidenceEngine, EvidenceError> {
        let mut request = self.request.clone();
        request.selection = selection;
        request.execution = execution;
        let mut reader = match &request.alignment_index {
            Some(index) => bam::IndexedReader::from_path_and_index(&request.alignments, index),
            None => bam::IndexedReader::from_path(&request.alignments),
        }
        .map_err(|error| EvidenceError::InvalidInput(format!("open worker alignment: {error}")))?;
        let sample_scope =
            EvidenceSampleScope::resolve(reader.header(), &request.sample_selection)?;
        if sample_scope != *self.sample_scope {
            return Err(EvidenceError::InvalidInput(
                "worker alignment sample metadata differs from the opened session".into(),
            ));
        }
        let sample_scope_bytes = sample_scope
            .memory_bytes()
            .saturating_add((reader.header().as_bytes().len() as u64).saturating_mul(2));
        drop(sample_scope);
        if alignment_is_cram(&request.alignments)? {
            let fasta = request
                .cram_reference
                .as_deref()
                .or_else(|| self.reference.fasta_path())
                .ok_or_else(|| {
                    EvidenceError::InvalidRequest(
                        "worker CRAM requires verified local FASTA".into(),
                    )
                })?;
            reader.set_reference(fasta).map_err(|error| {
                EvidenceError::InvalidInput(format!("set worker CRAM reference: {error}"))
            })?;
        }
        let contigs = self.reference.contigs().clone();
        validate_dictionary(&contigs, &header_contigs(reader.header())?)?;
        let tids = contigs
            .iter()
            .map(|contig| {
                reader.header().tid(contig.name.as_bytes()).ok_or_else(|| {
                    EvidenceError::InvalidInput("worker contig dictionary mismatch".into())
                })
            })
            .collect::<Result<Vec<_>, _>>()?;
        let plan = make_plan(
            &request.execution,
            0,
            request.execution.analyzer_bytes,
            crate::util::rss::peak_rss_bytes(),
            0,
            sample_scope_bytes,
        )?;
        let mut engine = EvidenceEngine {
            request,
            reference: self.reference.clone(),
            reader,
            contigs,
            tids,
            intervals: Vec::new(),
            sites: BTreeMap::new(),
            sample_scope: self.sample_scope.clone(),
            plan,
        };
        engine.set_selection(engine.request.selection.clone())?;
        Ok(engine)
    }
}

fn alignment_is_cram(path: &std::path::Path) -> Result<bool, EvidenceError> {
    let mut file = std::fs::File::open(path)?;
    let mut magic = [0u8; 4];
    let count = std::io::Read::read(&mut file, &mut magic)?;
    Ok(count == 4 && &magic == b"CRAM")
}
