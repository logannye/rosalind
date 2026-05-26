//! Core types — the lingua franca shared by every Rosalind layer.
//!
//! Note: this module is named `core`; inside the crate always reference it as
//! `crate::core::…`. Reach the std `core` crate (rarely needed) as `::core::…`.

pub mod budget;
pub mod error;
pub mod locus;

pub use budget::{MemoryBudget, WorkingSet};
pub use error::CoreError;
pub use locus::{Contig, ContigSet, Locus, Position};
