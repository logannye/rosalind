//! Genomics-specific utilities and data structures built on top of the O(√t)
//! engine.
//!
//! This module exposes foundational components used by higher-level genomic
//! algorithms (alignment, variant calling, etc.).

mod alignment;
mod block_alignment;
mod bwt_aligner;
mod compressed_dna;
mod eval;
mod fm_backing;
mod fm_index;
mod genome_index;
mod index;
mod io;
mod rank_select;
mod reference_pack;
mod sampled_sa;
mod sort;
mod suffix_array;
mod types;

pub use block_alignment::{align_within_block, AlignmentError, BlockAlignmentSummary, FMInterval};
pub use bwt_aligner::{AlignerError, AlignmentResult, BWTAligner};
pub use compressed_dna::{AmbiguityMask, CompressedDNA, CompressedDNAError};
pub use eval::{
    compare_callsets, normalize_variant, read_vcf_variants, BedIndex, BedParseError,
    ComparisonReport, NormalizeError, NormalizedVariant, VariantType, VcfParseError, VcfVariant,
    Zygosity,
};
pub use fm_index::{
    BWTBlock, BlockBoundary, BlockedFMIndex, CompressedBoundaries, FMIndexError, FmSymbol,
};
pub use genome_index::{GenomeIndex, GenomeIndexError, MAX_GENOME_LEN};
pub use index::{
    render_plan_line, BuildMemoryModel, FmIndexView, GenomeIndexView, IndexBuildReport,
    IndexHeader, IndexReader, IndexWriter, ReferenceIndex, ReferenceView,
};
pub use io::create_bam_writer;
pub(crate) use rank_select::popcount_range;
pub use rank_select::{
    BaseCode, RankSelectCheckpoint, RankSelectIndex, ALPHABET_SIZE, CHECKPOINT_STRIDE,
};
pub use reference_pack::{
    AnalysisReference, ReferencePackBuilder, ReferencePackError, ReferencePackMetadata,
    ReferencePackReader, ReferenceProvider, ReferenceSequence,
};
pub use sampled_sa::SampledSuffixArray;
pub(crate) use sampled_sa::RANK_STRIDE;
pub use sort::sort_bam_deterministic;
pub use types::{AlignedRead, CigarOp, CigarOpKind};
