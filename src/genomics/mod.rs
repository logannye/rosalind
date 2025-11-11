//! Genomics-specific utilities and data structures built on top of the O(√t)
//! engine.
//!
//! This module exposes foundational components used by higher-level genomic
//! algorithms (alignment, variant calling, etc.).

mod compressed_dna;
mod fm_index;
mod block_alignment;
mod rank_select;
mod bwt_aligner;
mod types;
mod pileup;
mod statistics;
mod variant_caller;

pub use compressed_dna::{AmbiguityMask, CompressedDNA, CompressedDNAError};
pub use fm_index::{
    BlockBoundary, BlockedFMIndex, CompressedBoundaries, FMIndexError, FmSymbol, BWTBlock,
};
pub use block_alignment::{
    align_within_block, AlignmentError, BlockAlignmentSummary, FMInterval,
};
pub use rank_select::{
    BaseCode, RankSelectCheckpoint, RankSelectIndex, ALPHABET_SIZE, CHECKPOINT_STRIDE,
};
pub use bwt_aligner::{BWTAligner, AlignerError, AlignmentResult};
pub use types::{AlignedRead, CigarOp, CigarOpKind};
pub use pileup::{PileupNode, PileupProcessor, PileupSummary, PileupWorkload};
pub use statistics::{VariantCall, bayesian_variant_caller};
pub use variant_caller::{StreamingVariantCaller, Variant, VariantCallerError};

