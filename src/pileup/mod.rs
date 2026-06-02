//! The streaming pileup kernel.
//!
//! `PileupEngine` consumes coordinate-sorted reads and yields one `PileupColumn`
//! per covered reference position. It is the single bounded-memory substrate the
//! germline/somatic callers build on; build your own bounded per-locus analytics
//! over the same stream via the `ColumnAnalyzer` trait (see
//! `examples/columnkit_coverage.rs`) or by iterating it directly
//! (`examples/custom_pileup_analytics.rs`). Reference as `crate::pileup::…`.

pub mod column;
pub mod engine;
pub mod source;

pub use column::{Obs, PileupColumn};
pub use engine::{PileupEngine, PileupParams, SkipCounts};
pub use source::{ReadSource, SliceSource};
