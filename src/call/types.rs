//! Result and parameter types for the calling layer — pure data, no logic.

/// Diploid genotype over the biallelic {reference, alt} pair.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Genotype {
    /// Homozygous reference (0/0).
    HomRef,
    /// Heterozygous (0/1).
    Het,
    /// Homozygous alternate (1/1).
    HomAlt,
}

/// VCF FILTER status for an emitted call.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Filter {
    /// Passed all thresholds.
    Pass,
    /// QUAL below the configured minimum.
    LowQual,
    /// Depth below the configured minimum.
    LowDepth,
}

/// ACGT index (0..=3) → uppercase ASCII base. Shared by the calling functions.
pub(crate) const ACGT: [u8; 4] = [b'A', b'C', b'G', b'T'];

/// Tunable thresholds for germline genotyping and honest filtering.
#[derive(Debug, Clone)]
pub struct GermlineParams {
    /// Prior heterozygosity θ (default 1e-3).
    pub heterozygosity: f64,
    /// Calls with DP below this are flagged `LowDepth` (default 8).
    pub min_depth: u32,
    /// Calls with QUAL below this are flagged `LowQual` (default 30.0).
    pub min_qual: f64,
}

impl Default for GermlineParams {
    fn default() -> Self {
        Self {
            heterozygosity: 1e-3,
            min_depth: 8,
            min_qual: 30.0,
        }
    }
}

/// A germline variant call. Emitted only when GT ≠ HomRef. Locus-free — the
/// consumer pairs it with the originating column's `Locus`.
#[derive(Debug, Clone, PartialEq)]
pub struct GermlineCall {
    /// Called genotype.
    pub genotype: Genotype,
    /// Alternate base (uppercase ASCII).
    pub alt_base: u8,
    /// Phred site-is-variant quality, −10·log10 P(0/0 | data).
    pub qual: f64,
    /// Genotype quality (second-smallest PL, capped 99).
    pub gq: u8,
    /// Phred genotype likelihoods, ordered [0/0, 0/1, 1/1], min 0, cap 255.
    pub pl: [u32; 3],
    /// Allelic depths [ref, alt].
    pub ad: [u32; 2],
    /// Callable read depth.
    pub dp: u32,
    /// FILTER status.
    pub filter: Filter,
}

/// Tunable thresholds for tumor/normal somatic SNV calling.
#[derive(Debug, Clone)]
pub struct SomaticParams {
    /// Minimum tumor depth to consider a site (default 8).
    pub min_tumor_depth: u32,
    /// Minimum normal depth to consider a site (default 8).
    pub min_normal_depth: u32,
    /// Minimum tumor alt allele fraction (default 0.05).
    pub min_tumor_af: f32,
    /// Maximum tolerated normal alt allele fraction (default 0.01).
    pub max_normal_af: f32,
    /// Minimum LLR-derived quality to emit (default 10.0).
    pub min_quality: f64,
    /// Sequencing error rate for the noise model (default 1e-3).
    pub seq_error_rate: f64,
}

impl Default for SomaticParams {
    fn default() -> Self {
        Self {
            min_tumor_depth: 8,
            min_normal_depth: 8,
            min_tumor_af: 0.05,
            max_normal_af: 0.01,
            min_quality: 10.0,
            seq_error_rate: 1e-3,
        }
    }
}

/// A tumor/normal somatic SNV call. Locus-free.
#[derive(Debug, Clone, PartialEq)]
pub struct SomaticCall {
    /// Reference base (uppercase ASCII).
    pub ref_base: u8,
    /// Alternate base (uppercase ASCII).
    pub alt_base: u8,
    /// Tumor alt count / depth.
    pub tumor_alt: u32,
    /// Tumor depth.
    pub tumor_depth: u32,
    /// Normal alt count.
    pub normal_alt: u32,
    /// Normal depth.
    pub normal_depth: u32,
    /// Tumor alt allele fraction.
    pub tumor_af: f32,
    /// Normal alt allele fraction.
    pub normal_af: f32,
    /// LLR-derived Phred-scaled quality (clamped 0..200).
    pub quality: f64,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn defaults_match_spec() {
        let g = GermlineParams::default();
        assert_eq!(g.heterozygosity, 1e-3);
        assert_eq!(g.min_depth, 8);
        assert_eq!(g.min_qual, 30.0);

        let s = SomaticParams::default();
        assert_eq!(s.min_tumor_depth, 8);
        assert_eq!(s.min_tumor_af, 0.05);
        assert_eq!(s.max_normal_af, 0.01);
        assert_eq!(s.seq_error_rate, 1e-3);

        assert_eq!(ACGT[0], b'A');
        assert_eq!(ACGT[3], b'T');
        assert_ne!(Genotype::Het, Genotype::HomAlt);
        assert_ne!(Filter::Pass, Filter::LowDepth);
    }
}
