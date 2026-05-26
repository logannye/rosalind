//! Somatic variant calling (tumor/normal) for SNVs.
//!
//! This module provides a deterministic, streaming-friendly somatic SNV caller.

mod model;
mod vcf;

pub use model::{SomaticCaller, SomaticCallerConfig, SomaticVariant};
pub use model::SomaticIndel;
pub use vcf::{render_somatic_vcf, write_somatic_vcf};


