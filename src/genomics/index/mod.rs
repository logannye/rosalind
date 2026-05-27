//! Reference index IO (on-disk, deterministic, mmap-friendly).
//!
//! This module defines Rosalind’s long-lived artifact boundary for WGS: a
//! versioned binary index format that can be generated once and reused across
//! runs, with stable serialization and the ability to memory-map on consumer
//! hardware.

mod format;
mod io;
mod view;

pub use format::IndexHeader;
pub use io::{IndexReader, IndexWriter, ReferenceIndex};
pub use view::{FmIndexView, GenomeIndexView};
