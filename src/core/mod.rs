//! Core types — the lingua franca shared by every Rosalind layer.
//!
//! Note: this module is named `core`; inside the crate always reference it as
//! `crate::core::…`. Reach the std `core` crate (rarely needed) as `::core::…`.

pub mod budget;
pub mod error;
pub mod locus;
pub mod record;
pub mod sequence;
pub use sequence::{allele_index, BaseCode};

pub use budget::{
    MemoryBudget, WorkingSet, PILEUP_ENGINE_OVERHEAD, PILEUP_IO_RSS_OVERHEAD,
    PILEUP_MAP_BYTES_PER_BASE, PILEUP_PER_READ_OVERHEAD, PILEUP_SEQQUAL_BYTES_PER_BASE,
};
pub use error::CoreError;
pub use locus::{Contig, ContigSet, Locus, Position};
pub use record::{AlignedRead, CigarOp, CigarOpKind, RefBase, SamFlags};
