use super::cram::CramPreflight;
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
    /// Codec-specific decoder reservation, included in fixed_bytes.
    pub decoder_bytes: u64,
    /// Codec-specific admission contract.
    pub decoder_model: &'static str,
    /// Checked CRAM container maxima; absent for BAM.
    pub cram_envelope: Option<CramEnvelope>,
}
impl EvidencePlan {
    /// Additive execution diagnostics; these do not change scientific identity.
    pub fn decoder_measurements(&self) -> BTreeMap<String, String> {
        let mut values = BTreeMap::from([
            ("execution.decoder_model".into(), self.decoder_model.into()),
            (
                "execution.decoder_bytes".into(),
                self.decoder_bytes.to_string(),
            ),
        ]);
        if let Some(cram) = &self.cram_envelope {
            for (key, value) in [
                ("containers", cram.containers),
                ("max_container_records", cram.max_container_records),
                ("max_container_bases", cram.max_container_bases),
                (
                    "base_reservation_bases",
                    cram.max_container_records
                        .saturating_mul(cram.validated_max_read_length)
                        .max(cram.max_container_bases),
                ),
                ("max_container_slices", cram.max_container_slices),
                ("max_container_blocks", cram.max_container_blocks),
                ("max_compressed_bytes", cram.max_compressed_bytes),
                ("max_uncompressed_bytes", cram.max_uncompressed_bytes),
                (
                    "max_compression_header_bytes",
                    cram.max_compression_header_bytes,
                ),
                ("header_bytes", cram.header_bytes),
                ("reference_bytes", cram.reference_bytes),
                ("index_bytes", cram.index_bytes),
                ("validated_records", cram.validated_records),
                (
                    "validated_max_record_bytes",
                    cram.validated_max_record_bytes,
                ),
                ("validated_max_read_length", cram.validated_max_read_length),
                ("validation_wall_micros", cram.validation_wall_micros),
                ("validation_peak_rss_bytes", cram.validation_peak_rss_bytes),
                ("declared_records", cram.declared_records),
                ("validated_bases", cram.validated_bases),
                ("declared_bases", cram.declared_bases),
                ("metadata_cap_bytes", 8 << 20),
                ("block_cap_bytes", 64 << 20),
                ("container_cap_bytes", 256 << 20),
                ("records_cap", 1_000_000),
                ("blocks_cap", 4096),
                ("slices_cap", 1024),
                ("codec_depth_cap", 16),
            ] {
                values.insert(format!("execution.cram.{key}"), value.to_string());
            }
            values.insert(
                "execution.cram.profile".into(),
                "3.0;coordinate;single-reference;RAW,GZIP,rANS4;checked-metadata-v1".into(),
            );
        }
        values
    }
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
    cram: Option<Arc<CramPreflight>>,
}
impl EvidenceEngine {
    /// Open and validate local inputs, retaining only indexed reference access and metadata.
    pub fn open(mut request: EvidenceRequest) -> Result<Self, EvidenceError> {
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
        let cram = if alignment_is_cram(&request.alignments)? {
            let checked = CramPreflight::inspect(&request)?;
            request.alignment_index = Some(checked.index.clone());
            Some(Arc::new(checked))
        } else {
            None
        };
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
        let is_cram = cram.is_some();
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
            request.fields,
            selected_loci,
            request.execution.analyzer_bytes,
            crate::util::rss::peak_rss_bytes(),
            selection_memory_bytes(&intervals, &sites),
            sample_scope_bytes,
            cram.as_deref().map(|p| &p.envelope),
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
            cram,
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
            cram: self.cram.clone(),
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
            self.request.fields,
            selected,
            self.plan.analyzer_bytes,
            self.plan.baseline_rss_bytes,
            selection_memory_bytes(&intervals, &sites),
            self.plan.sample_scope_bytes,
            self.cram.as_deref().map(|p| &p.envelope),
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
            self.request.fields,
            self.plan.selected_loci,
            bytes,
            self.plan.baseline_rss_bytes,
            self.plan.selection_bytes,
            self.plan.sample_scope_bytes,
            self.cram.as_deref().map(|p| &p.envelope),
        )?;
        Ok(&self.plan)
    }
    /// Execute in canonical coordinate order. Safe to execute once per engine;
    /// subsequent calls rerun the same request from indexed input.
    pub fn run(
        &mut self,
        analyzer: &mut dyn EvidenceAnalyzer,
    ) -> Result<EvidenceRunStats, EvidenceError> {
        if let Some(cram) = &self.cram {
            cram.verify()?;
        }
        self.plan_for_analyzer(analyzer)?;
        let mut stats = EvidenceRunStats::default();
        let mut record = bam::Record::new();
        // Determine local selected pressure from normalized intervals alone.
        // Allocate each summary group once, avoiding allocator-retained full
        // tiles followed by differently sized boundary tiles. Sparse requests
        // reserve their largest actual selected batch, not a dense chromosome.
        let mut interval_index = 0;
        let mut start = self.intervals.first().map_or(0, |interval| interval.start);
        let mut max_batch_loci = 0;
        while interval_index < self.intervals.len() {
            crate::core::governor::checkpoint()?;
            let window = coalesced_window(
                &self.intervals,
                interval_index,
                start,
                self.plan.microtile_bases,
            );
            max_batch_loci = max_batch_loci.max(window.selected_loci);
            interval_index = window.next_interval;
            start = window.next_start;
        }
        let mut batch = EvidenceBatch::new(0, "", 0, self.request.fields, Vec::new());
        batch.reserve_exact(max_batch_loci);
        let mut loci = Vec::with_capacity(max_batch_loci);
        interval_index = 0;
        start = self.intervals.first().map_or(0, |interval| interval.start);
        while interval_index < self.intervals.len() {
            if let Some(cram) = &self.cram {
                cram.verify()?;
            }
            let interval = &self.intervals[interval_index];
            crate::core::governor::checkpoint()?;
            let canonical_start = start / CANONICAL_TILE_BASES * CANONICAL_TILE_BASES;
            let window = coalesced_window(
                &self.intervals,
                interval_index,
                start,
                self.plan.microtile_bases,
            );
            let end = window.end;
            let sequence = self.reference.read_window(interval.contig, start, end)?;
            loci.clear();
            for selected in &self.intervals[interval_index..] {
                if selected.contig != interval.contig || selected.start >= end {
                    break;
                }
                for position in selected.start.max(start)..selected.end.min(end) {
                    let requested_alts = self
                        .sites
                        .get(&(interval.contig, position))
                        .map_or_else(Vec::new, |site| site.alternates.clone());
                    loci.push(EvidenceLocus {
                        position,
                        reference: sequence[(position - start) as usize],
                        requested_alts,
                    });
                }
            }
            batch.contig_id = interval.contig;
            batch.contig.clear();
            batch
                .contig
                .push_str(&self.contigs.by_id(interval.contig).unwrap().name);
            batch.canonical_tile_start = canonical_start;
            loci = batch.reset_loci(loci);
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
                if record.tid() != self.tids[interval.contig as usize] as i32 || record.pos() < 0 {
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
                    &mut batch,
                    start,
                    end,
                    self.contigs.by_id(interval.contig).unwrap().length,
                    &self.request.profile,
                    rejected,
                    self.reference.has_sequence(),
                )?;
            }
            checked_add(&mut stats.microtiles, 1)?;
            checked_add(&mut stats.emitted_loci, batch.len() as u64)?;
            analyzer.on_batch(&batch)?;
            self.check_budget()?;
            interval_index = window.next_interval;
            start = window.next_start;
        }
        analyzer.finish()?;
        if let Some(cram) = &self.cram {
            cram.verify()?;
        }
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

struct CoalescedWindow {
    end: u32,
    selected_loci: usize,
    next_interval: usize,
    next_start: u32,
}

// Fetch nearby selected intervals together, but allocate and count only their
// requested loci. Neither a canonical boundary nor the admitted coordinate span
// is crossed. This keeps a sparse VCF from fetching once for every single SNV.
fn coalesced_window(
    intervals: &[GenomicInterval],
    interval_index: usize,
    start: u32,
    width: u32,
) -> CoalescedWindow {
    let contig = intervals[interval_index].contig;
    let canonical_end =
        (start / CANONICAL_TILE_BASES * CANONICAL_TILE_BASES).saturating_add(CANONICAL_TILE_BASES);
    let limit = start.saturating_add(width).min(canonical_end);
    let mut next_interval = interval_index;
    let mut next_start = start;
    let mut selected_loci = 0;
    let mut end = start;
    while let Some(interval) = intervals.get(next_interval) {
        if interval.contig != contig || next_start >= limit {
            break;
        }
        end = interval.end.min(limit);
        selected_loci += (end - next_start) as usize;
        if end < interval.end {
            next_start = end;
            break;
        }
        next_interval += 1;
        next_start = intervals
            .get(next_interval)
            .map_or(0, |interval| interval.start);
    }
    CoalescedWindow {
        end,
        selected_loci,
        next_interval,
        next_start,
    }
}

// Keep independently auditable retained-memory terms explicit at the few callers.
#[allow(clippy::too_many_arguments)]
fn make_plan(
    execution: &EvidenceExecution,
    fields: EvidenceFields,
    selected_loci: u64,
    analyzer_bytes: u64,
    baseline: u64,
    selection_bytes: u64,
    sample_scope_bytes: u64,
    cram: Option<&CramEnvelope>,
) -> Result<EvidencePlan, EvidenceError> {
    // Two reusable identity vectors coexist while selecting the next window.
    // Reserve two small ALT allocations, reference-window bytes, and alignment
    // slack per coordinate in addition to both identity/vector element layouts.
    let bytes_per_locus =
        fields.storage_bytes_per_locus() + std::mem::size_of::<EvidenceLocus>() as u64 + 48;
    let cram_bytes = cram
        .map(|value| value.additional_bytes(execution))
        .transpose()?
        .unwrap_or(0);
    let decoder_bytes = DECODER_SLACK_BYTES
        .checked_add(
            (execution.max_record_bytes as u64)
                .checked_mul(3)
                .ok_or(EvidenceError::CounterOverflow)?,
        )
        .and_then(|bytes| bytes.checked_add(cram_bytes))
        .ok_or(EvidenceError::CounterOverflow)?;
    let fixed = decoder_bytes
        .checked_add(analyzer_bytes)
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
        model_id: "exact-summary-tiles-v4",
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
        decoder_bytes,
        decoder_model: cram.map_or("bam-record-envelope-v1", |value| value.model),
        cram_envelope: cram.cloned(),
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
fn increment_filter(row: &mut EvidenceDepths, reason: ReadFilter) -> Result<(), EvidenceError> {
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
    batch: &mut EvidenceBatch,
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
                let first = batch
                    .loci()
                    .partition_point(|locus| u64::from(locus.position) < overlap_start);
                for index in first..batch.len() {
                    let mut row = batch.row_mut(index).expect("selected locus index");
                    let position = u64::from(row.position);
                    if position >= overlap_end {
                        break;
                    }
                    let offset = query_position + (position - reference_position) as usize;
                    if let Some(depths) = row.depths.as_deref_mut() {
                        checked_add(&mut depths.prefilter_depth, 1)?;
                    }
                    if let Some(reason) = rejected {
                        if let Some(depths) = row.depths.as_deref_mut() {
                            increment_filter(depths, reason)?;
                        }
                        continue;
                    }
                    if let Some(depths) = row.depths.as_deref_mut() {
                        checked_add(&mut depths.aligned_depth, 1)?;
                    }
                    let quality = qualities[offset];
                    if quality == 255 {
                        if let Some(depths) = row.depths.as_deref_mut() {
                            checked_add(&mut depths.filters.unavailable_base_quality, 1)?;
                        }
                        continue;
                    }
                    if quality > 93 {
                        return Err(EvidenceError::InvalidInput(
                            "base quality outside SAM range 0..93 (or missing 255)".into(),
                        ));
                    }
                    if quality < profile.min_base_quality {
                        if let Some(depths) = row.depths.as_deref_mut() {
                            checked_add(&mut depths.filters.low_base_quality, 1)?;
                        }
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
                            if let Some(depths) = row.depths.as_deref_mut() {
                                checked_add(&mut depths.filters.ambiguous_base, 1)?;
                            }
                            continue;
                        }
                    };
                    if let Some(depths) = row.depths {
                        checked_add(&mut depths.callable_depth, 1)?;
                    }
                    if let Some(alleles) = row.alleles {
                        checked_add(&mut alleles.allele_counts[allele], 1)?;
                    }
                    if let Some(strands) = row.strands {
                        checked_add(
                            &mut strands.strand_counts[allele][usize::from(record.is_reverse())],
                            1,
                        )?;
                    }
                    if let Some(sums) = row.quality_sums {
                        checked_add(&mut sums.base_quality_sum, quality as u64)?;
                        checked_add(&mut sums.mapping_quality_sum, record.mapq() as u64)?;
                    }
                    let cycle = if record.is_reverse() {
                        record.seq_len() - 1 - offset
                    } else {
                        offset
                    };
                    if let Some(position) = row.read_position {
                        checked_add(&mut position.read_position_sum, cycle as u64)?;
                        checked_add(&mut position.read_length_sum, record.seq_len() as u64)?;
                    }
                    if let Some(quality_sums) = row.allele_quality {
                        checked_add(&mut quality_sums.base_quality_sum[allele], quality as u64)?;
                        checked_add(
                            &mut quality_sums.mapping_quality_sum[allele],
                            record.mapq() as u64,
                        )?;
                        checked_add(&mut quality_sums.read_position_sum[allele], cycle as u64)?;
                        checked_add(
                            &mut quality_sums.read_length_sum[allele],
                            record.seq_len() as u64,
                        )?;
                    }
                    if let Some(histograms) = row.quality_histograms {
                        checked_add(&mut histograms.base_quality_histogram[quality as usize], 1)?;
                        checked_add(
                            &mut histograms.mapping_quality_histogram[record.mapq() as usize],
                            1,
                        )?;
                    }
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
    cram: Option<Arc<CramPreflight>>,
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
        if let Some(cram) = &self.cram {
            cram.verify()?;
            cram.admit_decoder(&request.execution)?;
        }
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
            request.fields,
            0,
            request.execution.analyzer_bytes,
            crate::util::rss::peak_rss_bytes(),
            0,
            sample_scope_bytes,
            self.cram.as_deref().map(|p| &p.envelope),
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
            cram: self.cram.clone(),
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

#[cfg(test)]
mod projection_tests {
    use super::*;

    fn selected_loci() -> Vec<EvidenceLocus> {
        [0, 2, 3, 5, 9, 20, 24, 30]
            .into_iter()
            .map(|position| EvidenceLocus {
                position,
                reference: b'A',
                requested_alts: Vec::new(),
            })
            .collect()
    }

    #[test]
    fn every_physical_projection_matches_full_selected_observations() {
        let mut records = Vec::new();
        for (flags, mapq) in [(0, 60), (16, 60), (1024, 60), (0, 255)] {
            let mut record = bam::Record::new();
            let mut sequence = b"ACGTACGTACGTACGTACGTACGT".to_vec();
            sequence[5] = b'N';
            let mut qualities = vec![30; sequence.len()];
            qualities[2] = 255;
            qualities[3] = 10;
            record.set(
                b"read",
                Some(&bam::record::CigarString(vec![Cigar::Match(
                    sequence.len() as u32,
                )])),
                &sequence,
                &qualities,
            );
            record.set_tid(0);
            record.set_pos(0);
            record.set_flags(flags);
            record.set_mapq(mapq);
            records.push(record);
        }
        let extract = |fields| {
            let mut batch = EvidenceBatch::new(0, "chr1", 0, fields, selected_loci());
            for record in &records {
                let profile = EvidenceProfile::default();
                accumulate_record(
                    record,
                    &mut batch,
                    0,
                    31,
                    100,
                    &profile,
                    read_filter(record, &profile),
                    true,
                )
                .unwrap();
            }
            batch
        };
        let full = extract(EvidenceFields::ALL_SUPPORTED);
        assert_eq!(full.row(0).unwrap().depths.unwrap().prefilter_depth, 4);
        assert_eq!(full.row(0).unwrap().depths.unwrap().callable_depth, 2);
        assert_eq!(
            full.row(1)
                .unwrap()
                .depths
                .unwrap()
                .filters
                .unavailable_base_quality,
            2
        );
        assert_eq!(
            full.row(2)
                .unwrap()
                .depths
                .unwrap()
                .filters
                .low_base_quality,
            2
        );
        assert_eq!(
            full.row(3).unwrap().depths.unwrap().filters.ambiguous_base,
            2
        );
        assert_eq!(full.row(6).unwrap().depths.unwrap().prefilter_depth, 0);
        for bits in 0..=EvidenceFields::ALL_SUPPORTED.bits() {
            let fields = EvidenceFields::from_bits(bits).unwrap();
            let projected = extract(fields);
            assert_eq!(projected.len(), selected_loci().len());
            for (row, expected) in projected.rows().zip(full.rows()) {
                assert_eq!(
                    (row.position, row.reference),
                    (expected.position, expected.reference)
                );
                assert_eq!(
                    row.depths,
                    expected
                        .depths
                        .filter(|_| fields.contains(EvidenceFields::DEPTHS))
                );
                assert_eq!(
                    row.alleles,
                    expected
                        .alleles
                        .filter(|_| fields.contains(EvidenceFields::ALLELES))
                );
                assert_eq!(
                    row.strands,
                    expected
                        .strands
                        .filter(|_| fields.contains(EvidenceFields::STRANDS))
                );
                assert_eq!(
                    row.quality_sums,
                    expected
                        .quality_sums
                        .filter(|_| fields.contains(EvidenceFields::QUALITY_SUMS))
                );
                assert_eq!(
                    row.quality_histograms,
                    expected
                        .quality_histograms
                        .filter(|_| fields.contains(EvidenceFields::QUALITY_HISTOGRAMS))
                );
                assert_eq!(
                    row.allele_quality,
                    expected
                        .allele_quality
                        .filter(|_| fields.contains(EvidenceFields::ALLELE_QUALITY))
                );
                assert_eq!(
                    row.read_position,
                    expected
                        .read_position
                        .filter(|_| fields.contains(EvidenceFields::READ_POSITION))
                );
            }
        }
    }

    #[test]
    fn sparse_windows_obey_span_ownership_and_exact_selection() {
        let intervals = vec![
            GenomicInterval {
                contig: 0,
                start: 1,
                end: 2,
            },
            GenomicInterval {
                contig: 0,
                start: 10,
                end: 12,
            },
            GenomicInterval {
                contig: 0,
                start: 16_380,
                end: 16_390,
            },
            GenomicInterval {
                contig: 1,
                start: 3,
                end: 5,
            },
        ];
        for width in [1, 2, 10, 128, 16_384] {
            let mut index = 0;
            let mut start = intervals[0].start;
            let mut observed = Vec::new();
            while index < intervals.len() {
                let contig = intervals[index].contig;
                let window = coalesced_window(&intervals, index, start, width);
                assert!(window.end > start);
                assert!(window.end - start <= width);
                assert_eq!(
                    start / CANONICAL_TILE_BASES,
                    (window.end - 1) / CANONICAL_TILE_BASES
                );
                let before = observed.len();
                for interval in &intervals[index..] {
                    if interval.contig != contig || interval.start >= window.end {
                        break;
                    }
                    observed.extend(
                        (interval.start.max(start)..interval.end.min(window.end))
                            .map(|position| (contig, position)),
                    );
                }
                assert_eq!(observed.len() - before, window.selected_loci);
                index = window.next_interval;
                start = window.next_start;
            }
            let expected: Vec<_> = intervals
                .iter()
                .flat_map(|interval| {
                    (interval.start..interval.end).map(|position| (interval.contig, position))
                })
                .collect();
            assert_eq!(observed, expected);
        }
        let combined = coalesced_window(&intervals, 0, 1, 128);
        assert_eq!(
            (combined.end, combined.selected_loci, combined.next_interval),
            (12, 3, 2)
        );
    }

    #[test]
    fn projected_memory_model_admits_larger_windows_without_changing_fields() {
        let mut execution = EvidenceExecution::default();
        let full = make_plan(&execution, EvidenceFields::ALL, 16_384, 0, 0, 0, 0, None).unwrap();
        execution.memory_budget_bytes = Some(full.fixed_bytes + full.bytes_per_locus * 128);
        let full = make_plan(&execution, EvidenceFields::ALL, 16_384, 0, 0, 0, 0, None).unwrap();
        let fields = EvidenceFields::DEPTHS.union(EvidenceFields::QUALITY_SUMS);
        let projected = make_plan(&execution, fields, 16_384, 0, 0, 0, 0, None).unwrap();
        assert_eq!(full.microtile_bases, 128);
        assert!(projected.microtile_bases > full.microtile_bases * 10);
        assert!(projected.predicted_peak_rss_bytes <= execution.memory_budget_bytes.unwrap());
        assert_eq!(projected.model_id, "exact-summary-tiles-v4");
    }

    #[test]
    fn indexed_sparse_fetches_reduce_visits_and_preserve_full_output_bytes() {
        struct Fixture(PathBuf);
        impl Drop for Fixture {
            fn drop(&mut self) {
                let _ = std::fs::remove_dir_all(&self.0);
            }
        }
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let fixture = Fixture(std::env::temp_dir().join(format!(
            "rosalind-projection-{}-{nonce}",
            std::process::id()
        )));
        std::fs::create_dir(&fixture.0).unwrap();
        let fasta = fixture.0.join("reference.fa");
        std::fs::write(&fasta, format!(">chr1\n{}\n", "A".repeat(20_000))).unwrap();
        std::fs::write(
            fixture.0.join("reference.fa.fai"),
            "chr1\t20000\t6\t20000\t20001\n",
        )
        .unwrap();
        let alignment = fixture.0.join("reads.bam");
        let mut header = bam::Header::new();
        let mut hd = bam::header::HeaderRecord::new(b"HD");
        hd.push_tag(b"VN", "1.6").push_tag(b"SO", "coordinate");
        header.push_record(&hd);
        let mut sq = bam::header::HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr1").push_tag(b"LN", 20_000);
        header.push_record(&sq);
        let mut writer = bam::Writer::from_path(&alignment, &header, bam::Format::Bam).unwrap();
        for (index, position) in [0, 9000, 16_375].into_iter().enumerate() {
            let mut record = bam::Record::new();
            record.set(
                format!("read-{index}").as_bytes(),
                Some(&bam::record::CigarString(vec![Cigar::Match(100)])),
                &[b'T'; 100],
                &[30; 100],
            );
            record.set_tid(0);
            record.set_pos(position);
            record.set_mapq(60);
            record.set_flags(0);
            writer.write(&record).unwrap();
        }
        drop(writer);
        bam::index::build(&alignment, None::<&PathBuf>, bam::index::Type::Bai, 1).unwrap();
        let positions = [1, 10, 20, 9000, 16_383, 16_384, 19_999];
        let request = |width, fields| {
            let mut request = EvidenceRequest::new(&alignment, &fasta);
            request.fields = fields;
            request.execution.max_microtile_bases = width;
            request.selection = EvidenceSelection::Sites(
                positions
                    .into_iter()
                    .map(|position| SnvSite {
                        contig: 0,
                        position,
                        reference: b'A',
                        alternates: vec![b'T'],
                    })
                    .collect(),
            );
            request
        };
        let run = |width| {
            let mut engine = EvidenceEngine::open(request(width, EvidenceFields::ALL)).unwrap();
            let mut writer = EvidenceArrowWriter::new(Vec::new());
            let stats = engine.run(&mut writer).unwrap();
            let arrow = writer.into_inner().unwrap();
            let mut writer = EvidenceTsvWriter::new(Vec::new());
            engine.run(&mut writer).unwrap();
            (stats, arrow, writer.into_inner())
        };
        let narrow = run(1);
        let coalesced = run(CANONICAL_TILE_BASES);
        assert_eq!(narrow.1, coalesced.1);
        assert_eq!(narrow.2, coalesced.2);
        assert_eq!(narrow.0.emitted_loci, positions.len() as u64);
        assert_eq!(narrow.0.record_visits, 6);
        assert_eq!(coalesced.0.record_visits, 4);
        assert!(coalesced.0.microtiles < narrow.0.microtiles);
        let fields = EvidenceFields::DEPTHS.union(EvidenceFields::QUALITY_SUMS);
        let mut engine = EvidenceEngine::open(request(CANONICAL_TILE_BASES, fields)).unwrap();
        let mut observed = Vec::new();
        let mut callback = EvidenceCallback::with_fields(
            |batch: &EvidenceBatch| {
                assert_eq!(batch.fields(), fields);
                for row in batch.rows() {
                    assert!(row.quality_histograms.is_none());
                    assert!(row.alleles.is_none());
                    assert!(row.try_to_full_row().is_err());
                    observed.push((
                        row.position,
                        row.depths.unwrap().callable_depth,
                        row.quality_sums.unwrap().base_quality_sum,
                    ));
                }
                Ok(())
            },
            4096,
            fields,
        );
        engine.run(&mut callback).unwrap();
        assert_eq!(
            observed,
            positions
                .into_iter()
                .map(|position| (
                    position,
                    u64::from(position != 19_999),
                    if position == 19_999 { 0 } else { 30 }
                ))
                .collect::<Vec<_>>()
        );
    }
}
