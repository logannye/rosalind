//! The calling layer: turn the `PileupColumn` stream into calibrated,
//! abstention-aware variant calls. Built on `crate::core` + `crate::pileup`
//! only; no VCF writing or CLI wiring (those are later phases).

pub mod features;
pub mod germline;
pub mod pipeline;
pub mod plan;
pub mod somatic;
pub mod types;
pub mod whole_genome;

pub use features::{
    stream_features_region, stream_features_whole_genome, write_feature_header, write_feature_row,
};
pub use germline::call_germline;
pub use pipeline::{
    call_germline_region, call_germline_region_streaming, call_germline_region_tracked,
    call_somatic_region,
};
pub use plan::{estimate_variants_working_set, predicted_peak_rss_bytes, render_variants_plan};
pub use somatic::call_somatic;
pub use types::{Filter, Genotype, GermlineCall, GermlineParams, SomaticCall, SomaticParams};
pub use whole_genome::call_germline_whole_genome;
