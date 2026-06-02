//! The streaming pileup kernel.
//!
//! `PileupEngine` consumes coordinate-sorted reads and yields one `PileupColumn`
//! per covered reference position. It is the single bounded-memory substrate the
//! germline/somatic callers build on; build your own bounded per-locus analytics
//! over the same stream (see `examples/custom_pileup_analytics.rs`). The legacy
//! `GenomicPlugin`/`framework` lineage is separate and NOT memory-bounded.
//! Reference as `crate::pileup::…`.

pub mod column;
pub mod engine;
pub mod source;

pub use column::{Obs, PileupColumn};
pub use engine::{PileupEngine, PileupParams, SkipCounts};
pub use source::{ReadSource, SliceSource};
