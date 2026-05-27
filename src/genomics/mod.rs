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
mod pileup;
mod pileup_stream;
mod rank_select;
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
};
pub use fm_index::{
    BWTBlock, BlockBoundary, BlockedFMIndex, CompressedBoundaries, FMIndexError, FmSymbol,
};
pub use genome_index::{GenomeIndex, GenomeIndexError, MAX_GENOME_LEN};
pub use index::{IndexHeader, IndexReader, IndexWriter, ReferenceIndex};
pub use io::create_bam_writer;
pub use pileup::{PileupNode, PileupProcessor, PileupSummary, PileupWorkload};
pub use pileup_stream::BamPileupStream;
pub use rank_select::{
    BaseCode, RankSelectCheckpoint, RankSelectIndex, ALPHABET_SIZE, CHECKPOINT_STRIDE,
};
pub use sampled_sa::SampledSuffixArray;
pub use sort::sort_bam_deterministic;
pub use types::{AlignedRead, CigarOp, CigarOpKind};
