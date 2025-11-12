//! Genomics-specific utilities and data structures built on top of the O(√t)
//! engine.
//!
//! This module exposes foundational components used by higher-level genomic
//! algorithms (alignment, variant calling, etc.).

mod block_alignment;
mod bwt_aligner;
mod compressed_dna;
mod fm_index;
mod io;
mod pileup;
mod rank_select;
mod statistics;
mod types;
mod variant_caller;
mod vcf;

pub use block_alignment::{align_within_block, AlignmentError, BlockAlignmentSummary, FMInterval};
pub use bwt_aligner::{AlignerError, AlignmentResult, BWTAligner};
pub use compressed_dna::{AmbiguityMask, CompressedDNA, CompressedDNAError};
pub use fm_index::{
    BWTBlock, BlockBoundary, BlockedFMIndex, CompressedBoundaries, FMIndexError, FmSymbol,
};
pub use io::create_bam_writer;
pub use pileup::{PileupNode, PileupProcessor, PileupSummary, PileupWorkload};
pub use rank_select::{
    BaseCode, RankSelectCheckpoint, RankSelectIndex, ALPHABET_SIZE, CHECKPOINT_STRIDE,
};
pub use statistics::{bayesian_variant_caller, VariantCall};
pub use types::{AlignedRead, CigarOp, CigarOpKind};
pub use variant_caller::{StreamingVariantCaller, Variant, VariantCallerError};
pub use vcf::{render_vcf, write_vcf};
