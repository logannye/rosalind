//! The calling layer: turn the `PileupColumn` stream into calibrated,
//! abstention-aware variant calls. Built on `crate::core` + `crate::pileup`
//! only; no VCF writing or CLI wiring (those are later phases).

pub mod germline;
pub mod somatic;
pub mod types;

// pub use germline::call_germline;  // restored in Task 3
// pub use somatic::call_somatic;    // restored in Task 4
pub use types::{
    Filter, GermlineCall, GermlineParams, Genotype, SomaticCall, SomaticParams,
};
