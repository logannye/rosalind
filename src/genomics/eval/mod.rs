//! Clinical evaluation utilities (truth comparison and metrics).
//!
//! This module is used to evaluate end-to-end somatic calling quality on
//! clinically relevant benchmarks (e.g. tumor/normal WGS slices).

mod bed;
mod compare;
mod normalize;
mod vcf;

pub use bed::{BedIndex, BedParseError};
pub use compare::{compare_callsets, ComparisonReport, VariantType};
pub use normalize::{normalize_variant, NormalizeError, NormalizedVariant};
pub use vcf::{read_vcf_variants, VcfParseError, VcfVariant, Zygosity};
