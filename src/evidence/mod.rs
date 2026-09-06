//! Exact, budget-planned short-read evidence over indexed coordinate intervals.
//!
//! Every successful row includes every eligible observation. Memory selects an
//! execution microtile; it never selects biological observations. Read counts
//! intentionally count overlapping mates separately and do not collapse UMIs.

mod encoding;
mod engine;
mod panel;
mod reference;
mod selection;

pub use encoding::{
    read_evidence_batches, EvidenceArrowWriter, EvidenceTsvWriter, EVIDENCE_ARROW_BATCH_ROWS,
};
pub use engine::{EvidenceEngine, EvidencePlan, EvidenceRunStats, EvidenceWorkerFactory};
pub use panel::{FusedAnalyzers, PanelQcAnalyzer, PanelSummary, PanelTarget};
pub use reference::EvidenceReference;
pub use selection::{EvidenceSelection, SnvSite};

use crate::core::ContigSet;
use std::path::PathBuf;
use thiserror::Error;

/// Stable scientific schema; changes to fields or interpretation require a new version.
pub const EVIDENCE_SCHEMA_VERSION: u32 = 1;
/// Exact counting/filter semantics; increment when successful scientific results change.
pub const EVIDENCE_SEMANTICS_VERSION: &str = "exact-read-summary-v1";
/// Stable ownership boundaries, independent of the chosen microtile width.
pub const CANONICAL_TILE_BASES: u32 = 16_384;
/// Supported base quality values 0 through 93; 255 denotes missing quality.
pub const BASE_QUALITY_BINS: usize = 94;
/// Reported mapping qualities 0 through 254; 255 denotes unavailable MAPQ.
pub const MAPPING_QUALITY_BINS: usize = 255;

/// Failures distinguish admission, input envelope, arithmetic and I/O failures.
#[derive(Debug, Error)]
pub enum EvidenceError {
    #[error("invalid evidence request: {0}")]
    /// A request, selection or declared parameter is inconsistent.
    InvalidRequest(String),
    #[error("invalid evidence input: {0}")]
    /// Input data, coordinates or reference identity failed validation.
    InvalidInput(String),
    #[error("evidence plan refused: needs at least {needed} bytes, budget is {budget} bytes")]
    /// The declared process budget cannot admit the minimum execution window.
    Refused {
        /// Minimum predicted process bytes needed for one locus.
        needed: u64,
        /// Declared process memory budget in bytes.
        budget: u64,
    },
    #[error("evidence record exceeds declared envelope: {0}")]
    /// A decoded read or alignment record exceeds the declared envelope.
    RecordLimit(String),
    #[error("evidence integer counter overflow")]
    /// Checked integer accumulation exceeded its representable range.
    CounterOverflow,
    #[error(transparent)]
    /// An underlying local I/O operation failed.
    Io(#[from] std::io::Error),
    #[error(transparent)]
    /// The shared process governor or core input validation failed.
    Core(#[from] crate::core::CoreError),
    #[error("evidence analyzer failed: {0}")]
    /// A downstream consumer or encoder failed.
    Analyzer(String),
}

/// Explicit, versioned short-read DNA filter profile. Qualities marked unavailable
/// never satisfy a quality threshold, including a threshold of zero.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EvidenceProfile {
    /// Minimum reported mapping quality; unavailable MAPQ255 is excluded.
    pub min_mapq: u8,
    /// Minimum reported base quality; missing quality255 is excluded.
    pub min_base_quality: u8,
    /// Exclude SAM secondary alignments (0x100).
    pub exclude_secondary: bool,
    /// Exclude SAM supplementary alignments (0x800).
    pub exclude_supplementary: bool,
    /// Exclude SAM QC-failed reads (0x200).
    pub exclude_qc_fail: bool,
    /// Exclude duplicate-flagged reads (0x400); no UMI grouping is performed.
    pub exclude_duplicates: bool,
}
impl Default for EvidenceProfile {
    fn default() -> Self {
        Self {
            min_mapq: 20,
            min_base_quality: 20,
            exclude_secondary: true,
            exclude_supplementary: true,
            exclude_qc_fail: true,
            exclude_duplicates: true,
        }
    }
}
impl EvidenceProfile {
    /// Stable profile identifier; parameter overrides remain explicit recipe fields.
    pub const ID: &'static str = "shortread-dna-readcount-v1";
}

/// Execution controls do not alter the scientific result of a successful run.
#[derive(Debug, Clone)]
pub struct EvidenceExecution {
    /// Whole-process admission budget. None keeps planning advisory.
    pub memory_budget_bytes: Option<u64>,
    /// Maximum execution window width; the planner may choose a smaller width.
    pub max_microtile_bases: u32,
    /// Maximum decoded sequence length accepted by this execution.
    pub max_read_len: usize,
    /// Record length is checked immediately after htslib decode. This is a
    /// cooperative envelope, not a guarantee against the decoder's allocation.
    pub max_record_bytes: usize,
    /// Additional retained analyzer/encoder bytes reserved before admission.
    pub analyzer_bytes: u64,
}
impl Default for EvidenceExecution {
    fn default() -> Self {
        Self {
            memory_budget_bytes: None,
            max_microtile_bases: CANONICAL_TILE_BASES,
            max_read_len: 250,
            max_record_bytes: 1 << 20,
            analyzer_bytes: 0,
        }
    }
}

/// Versioned evidence capability set. Version 1 always emits ALL physically;
/// smaller consumer requirements are recorded without silently zeroing fields.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EvidenceFields(u32);
impl EvidenceFields {
    /// Depth and exclusive filter-reason counters.
    pub const DEPTHS: Self = Self(1);
    /// A/C/G/T callable counts.
    pub const ALLELES: Self = Self(2);
    /// Allele counts separated by read strand.
    pub const STRANDS: Self = Self(4);
    /// Exact callable base-quality and mapping-quality sums.
    pub const QUALITY_SUMS: Self = Self(8);
    /// Full callable BQ94 and MAPQ255 distributions.
    pub const QUALITY_HISTOGRAMS: Self = Self(16);
    /// Sequencing-cycle and read-length sums.
    pub const READ_POSITION: Self = Self(32);
    /// Complete evidence schema version1.
    pub const ALL: Self = Self(63);
    /// Stable field-mask version for receipts and cache keys.
    pub const VERSION: u32 = 1;
    /// Construct a validated field set from its portable bit representation.
    pub fn from_bits(bits: u32) -> Result<Self, EvidenceError> {
        if bits & !Self::ALL.0 != 0 {
            Err(EvidenceError::InvalidRequest(
                "unknown evidence field bits".into(),
            ))
        } else {
            Ok(Self(bits))
        }
    }
    /// Portable bit representation for scientific recipe identity.
    pub fn bits(self) -> u32 {
        self.0
    }
    /// Whether every required field is available in this set.
    pub fn contains(self, required: Self) -> bool {
        self.0 & required.0 == required.0
    }
    /// Union of two capability sets.
    pub fn union(self, other: Self) -> Self {
        Self(self.0 | other.0)
    }
    /// Canonical names in fixed schema order for human-readable receipts.
    pub fn names(self) -> Vec<&'static str> {
        [
            (Self::DEPTHS, "depths"),
            (Self::ALLELES, "alleles"),
            (Self::STRANDS, "strands"),
            (Self::QUALITY_SUMS, "quality-sums"),
            (Self::QUALITY_HISTOGRAMS, "quality-histograms"),
            (Self::READ_POSITION, "read-position"),
        ]
        .into_iter()
        .filter_map(|(field, name)| self.contains(field).then_some(name))
        .collect()
    }
}
impl Default for EvidenceFields {
    fn default() -> Self {
        Self::ALL
    }
}
/// Consumer capabilities and memory requirements validated before execution.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EvidenceRequirements {
    /// Scientific summaries consumed by the analyzer.
    pub fields: EvidenceFields,
    /// Whether N placeholders from reference-free coverage are unacceptable.
    pub requires_reference: bool,
    /// Flanking reference radius; only zero is supported by schema version1.
    pub context_bases: u32,
    /// Conservative additional retained bytes; None means unknown.
    pub retained_bytes: Option<u64>,
}

/// Cloneable request so a canonical partition can be executed in another process.
#[derive(Debug, Clone)]
pub struct EvidenceRequest {
    /// Local coordinate-sorted indexed BAM or CRAM path.
    pub alignments: PathBuf,
    /// Explicit BAI/CSI/CRAI, supporting content-addressed relocated replay.
    pub alignment_index: Option<PathBuf>,
    /// Explicit FAI for a relocated uncompressed analysis FASTA.
    pub reference_fai: Option<PathBuf>,
    /// Explicit CRAM FAI; must be adjacent to the FASTA for htslib.
    pub cram_reference_fai: Option<PathBuf>,
    /// None is supported for coverage/evidence counts, with reference bases N.
    pub reference: Option<PathBuf>,
    /// Required local FASTA for CRAM when `reference` is a pack or absent.
    pub cram_reference: Option<PathBuf>,
    /// Requested loci, normalized before execution.
    pub selection: EvidenceSelection,
    /// Physical output fields; schema v1 requires the complete field set.
    pub fields: EvidenceFields,
    /// Scientific filtering rules, independent of the memory budget.
    pub profile: EvidenceProfile,
    /// Resource limits and execution window preferences.
    pub execution: EvidenceExecution,
}
impl EvidenceRequest {
    /// Construct the configured value without starting evidence extraction.
    pub fn new(alignments: impl Into<PathBuf>, reference: impl Into<PathBuf>) -> Self {
        Self {
            reference: Some(reference.into()),
            ..Self::coverage(alignments)
        }
    }
    /// Create a coverage request using the alignment dictionary and unavailable reference bases.
    pub fn coverage(alignments: impl Into<PathBuf>) -> Self {
        Self {
            alignments: alignments.into(),
            alignment_index: None,
            reference_fai: None,
            cram_reference_fai: None,
            reference: None,
            cram_reference: None,
            selection: EvidenceSelection::WholeGenome,
            fields: EvidenceFields::ALL,
            profile: EvidenceProfile::default(),
            execution: EvidenceExecution::default(),
        }
    }
}

/// Exclusive first-failure counts at each matched reference base. The order is
/// secondary, supplementary, QC-fail, duplicate, unavailable/low MAPQ, then
/// unavailable/low base quality and ambiguous base.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct EvidenceFilterCounts {
    /// Matched observations excluded by the secondary-alignment filter.
    pub secondary: u64,
    /// Matched observations excluded by the supplementary-alignment filter.
    pub supplementary: u64,
    /// Matched observations excluded by the QC-failed flag.
    pub qc_fail: u64,
    /// Matched observations excluded by the duplicate flag.
    pub duplicate: u64,
    /// Matched observations whose mapping quality is unavailable (255).
    pub unavailable_mapq: u64,
    /// Matched observations below the mapping-quality threshold.
    pub low_mapq: u64,
    /// Eligible observations whose base quality is missing (255).
    pub unavailable_base_quality: u64,
    /// Eligible observations below the base-quality threshold.
    pub low_base_quality: u64,
    /// Quality-qualified observations whose base is not A/C/G/T.
    pub ambiguous_base: u64,
}

/// One exact row. Positions are zero-based internally; egress uses one-based POS.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EvidenceRow {
    /// Zero-based reference coordinate within the contig.
    pub position: u32,
    /// Normalized reference base at this locus, or N when unavailable.
    pub reference: u8,
    /// M/=/X observations before the profile's flag and quality filters.
    pub prefilter_depth: u64,
    /// Observations after read-level filtering, before base-level filtering.
    pub aligned_depth: u64,
    /// Number of A/C/G/T observations passing all declared filters.
    pub callable_depth: u64,
    /// Callable observation counts in A/C/G/T order.
    pub allele_counts: [u64; 4],
    /// Callable allele counts indexed by A/C/G/T then forward/reverse.
    pub strand_counts: [[u64; 2]; 4],
    /// Exact sum of callable-observation base qualities.
    pub base_quality_sum: u64,
    /// Exact sum of callable-observation mapping qualities.
    pub mapping_quality_sum: u64,
    /// Zero-based stored-SEQ offset sum, reversed back on reverse alignments.
    /// Includes soft-clipped offsets; hard-clipped bases are not in stored SEQ.
    pub read_position_sum: u64,
    /// Sum of full sequence lengths for callable observations, including soft clips.
    pub read_length_sum: u64,
    /// Exact callable-observation counts for base qualities0 through93.
    pub base_quality_histogram: [u64; BASE_QUALITY_BINS],
    /// Exact callable-observation counts for mapping qualities0 through254.
    pub mapping_quality_histogram: [u64; MAPPING_QUALITY_BINS],
    /// Exclusive first-failure counts under the declared profile.
    pub filters: EvidenceFilterCounts,
    /// Requested ALT alleles for VCF selection, in A/C/G/T order; otherwise empty.
    pub requested_alts: Vec<u8>,
}
impl Default for EvidenceRow {
    fn default() -> Self {
        Self {
            position: 0,
            reference: b'N',
            prefilter_depth: 0,
            aligned_depth: 0,
            callable_depth: 0,
            allele_counts: [0; 4],
            strand_counts: [[0; 2]; 4],
            base_quality_sum: 0,
            mapping_quality_sum: 0,
            read_position_sum: 0,
            read_length_sum: 0,
            base_quality_histogram: [0; BASE_QUALITY_BINS],
            mapping_quality_histogram: [0; MAPPING_QUALITY_BINS],
            filters: EvidenceFilterCounts::default(),
            requested_alts: Vec::new(),
        }
    }
}

/// Borrowed by analyzers during a callback. Retaining rows requires an explicit
/// copy, which the analyzer must include in its declared memory model.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EvidenceBatch {
    /// Stable contig identifier from the reference dictionary.
    pub contig_id: u32,
    /// Reference contig identifier or canonical contig name for a batch.
    pub contig: String,
    /// Zero-based start of the immutable16384-base ownership tile.
    pub canonical_tile_start: u32,
    /// Exact rows in strictly increasing coordinate order.
    pub rows: Vec<EvidenceRow>,
}

/// Batch consumer. The engine owns exact extraction; consumers own bounded
/// downstream state. Unknown memory models are refused for budgeted runs.
pub trait EvidenceAnalyzer {
    /// Declare needed scientific capabilities and a retained-memory bound.
    fn requirements(&self) -> EvidenceRequirements {
        EvidenceRequirements {
            fields: EvidenceFields::ALL,
            requires_reference: false,
            context_bases: 0,
            retained_bytes: self.additional_memory_bytes(),
        }
    }
    /// The fn value.
    fn on_batch(&mut self, batch: &EvidenceBatch) -> Result<(), EvidenceError>;
    /// The fn value.
    fn finish(&mut self) -> Result<(), EvidenceError> {
        Ok(())
    }
    /// The fn value.
    fn additional_memory_bytes(&self) -> Option<u64> {
        None
    }
}

/// Callback adapter for stateless/bounded SDK consumers.
#[derive(Debug)]
pub struct EvidenceCallback<F> {
    callback: F,
    bound: u64,
}
impl<F> EvidenceCallback<F> {
    /// Construct the configured value without starting evidence extraction.
    pub fn new(callback: F, additional_memory_bytes: u64) -> Self {
        Self {
            callback,
            bound: additional_memory_bytes,
        }
    }
}
impl<F: FnMut(&EvidenceBatch) -> Result<(), EvidenceError>> EvidenceAnalyzer
    for EvidenceCallback<F>
{
    fn on_batch(&mut self, batch: &EvidenceBatch) -> Result<(), EvidenceError> {
        (self.callback)(batch)
    }
    fn additional_memory_bytes(&self) -> Option<u64> {
        Some(self.bound)
    }
}

pub(crate) fn checked_add(value: &mut u64, increment: u64) -> Result<(), EvidenceError> {
    *value = value
        .checked_add(increment)
        .ok_or(EvidenceError::CounterOverflow)?;
    Ok(())
}
pub(crate) fn header_contigs(
    header: &rust_htslib::bam::HeaderView,
) -> Result<ContigSet, EvidenceError> {
    let mut contigs = ContigSet::new();
    for tid in 0..header.target_count() {
        let name = std::str::from_utf8(header.tid2name(tid))
            .map_err(|_| EvidenceError::InvalidInput("non-UTF8 alignment contig name".into()))?;
        let length = header
            .target_len(tid)
            .ok_or_else(|| EvidenceError::InvalidInput("missing contig length".into()))?;
        let length = u32::try_from(length)
            .map_err(|_| EvidenceError::InvalidInput("contig exceeds u32 coordinates".into()))?;
        if length == 0 || contigs.by_name(name).is_some() {
            return Err(EvidenceError::InvalidInput(
                "empty or duplicate alignment contig".into(),
            ));
        }
        contigs.push(name, length);
    }
    Ok(contigs)
}
