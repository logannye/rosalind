//! The streaming pileup kernel.
//!
//! `PileupEngine` consumes coordinate-sorted reads and yields one `PileupColumn`
//! per covered reference position. It is the single substrate that variant
//! callers and plugins build on. Reference as `crate::pileup::…`.

pub mod column;
pub mod engine;
pub mod source;

pub use column::{Obs, PileupColumn};
pub use engine::{PileupEngine, PileupParams, SkipCounts};
pub use source::{ReadSource, SliceSource};
