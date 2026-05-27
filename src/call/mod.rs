//! The calling layer: turn the `PileupColumn` stream into calibrated,
//! abstention-aware variant calls. Built on `crate::core` + `crate::pileup`
//! only; no VCF writing or CLI wiring (those are later phases).

pub mod germline;
pub mod pipeline;
pub mod somatic;
pub mod types;

pub use germline::call_germline;
pub use pipeline::call_germline_region;
pub use somatic::call_somatic;
pub use types::{Filter, Genotype, GermlineCall, GermlineParams, SomaticCall, SomaticParams};
