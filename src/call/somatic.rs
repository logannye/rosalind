//! The tumor/normal somatic SNV caller. The likelihood math
//! (`ln_factorial`/`ln_choose`/`ln_binom_pmf`/`somatic_llr`) is re-homed
//! verbatim from the legacy `genomics::somatic::model`; here it is fed corrected
//! `PileupColumn`s, so the legacy reverse-complement bug never enters.

use crate::call::types::{SomaticCall, SomaticParams, ACGT};
use crate::core::allele_index;
use crate::pileup::PileupColumn;

const LN10: f64 = std::f64::consts::LN_10;

/// Call a somatic SNV from a tumor column and a position-matched normal column.
/// Filter order matches the legacy caller: depth → alt-exists → AF gates →
/// LLR/quality. Returns `None` when any gate fails. Locus-free.
pub fn call_somatic(
    tumor: &PileupColumn,
    normal: &PileupColumn,
    params: &SomaticParams,
) -> Option<SomaticCall> {
    let t_depth = tumor.depth();
    let n_depth = normal.depth();
    if t_depth < params.min_tumor_depth || n_depth < params.min_normal_depth {
        return None;
    }

    let ref_idx = allele_index(tumor.ref_base)?;
    let t_counts = tumor.allele_counts();
    let n_counts = normal.allele_counts();

    // Alt = tumor's most-supported non-reference allele; ties → lowest index.
    let mut alt_idx: Option<usize> = None;
    let mut best = 0u32;
    for (i, &cnt) in t_counts.iter().enumerate() {
        if i == ref_idx {
            continue;
        }
        if cnt > best {
            best = cnt;
            alt_idx = Some(i);
        }
    }
    let alt_idx = alt_idx?;
    let t_alt = t_counts[alt_idx];
    if t_alt == 0 {
        return None;
    }
    let n_alt = n_counts[alt_idx];

    let t_af = t_alt as f32 / t_depth as f32;
    let n_af = n_alt as f32 / n_depth as f32;
    if t_af < params.min_tumor_af || n_af > params.max_normal_af {
        return None;
    }

    let llr = somatic_llr(t_alt, t_depth, n_alt, n_depth, params.seq_error_rate);
    let quality = (llr / LN10 * 10.0).clamp(0.0, 200.0);
    if quality < params.min_quality {
        return None;
    }

    Some(SomaticCall {
        ref_base: tumor.ref_base.to_ascii_uppercase(),
        alt_base: ACGT[alt_idx],
        tumor_alt: t_alt,
        tumor_depth: t_depth,
        normal_alt: n_alt,
        normal_depth: n_depth,
        tumor_af: t_af,
        normal_af: n_af,
        quality,
    })
}

// --- Likelihood math, re-homed verbatim from genomics::somatic::model ---

fn ln_factorial(n: u32) -> f64 {
    // For sequencing depths, n is small; compute exactly and deterministically.
    let mut acc = 0.0f64;
    for k in 2..=n {
        acc += (k as f64).ln();
    }
    acc
}

fn ln_choose(n: u32, k: u32) -> f64 {
    if k > n {
        return f64::NEG_INFINITY;
    }
    ln_factorial(n) - ln_factorial(k) - ln_factorial(n - k)
}

fn ln_binom_pmf(k: u32, n: u32, p: f64) -> f64 {
    let p = p.clamp(1e-12, 1.0 - 1e-12);
    ln_choose(n, k) + (k as f64) * p.ln() + ((n - k) as f64) * (1.0 - p).ln()
}

fn somatic_llr(t_alt: u32, t_depth: u32, n_alt: u32, n_depth: u32, err: f64) -> f64 {
    // Somatic vs noise-only: under somatic the tumor alt fraction is the MLE and
    // the normal alt fraction is ~err; under noise both are ~err.
    let t_p_hat = (t_alt as f64 / t_depth.max(1) as f64).clamp(err, 1.0 - 1e-6);
    let somatic = ln_binom_pmf(t_alt, t_depth, t_p_hat) + ln_binom_pmf(n_alt, n_depth, err);
    let noise = ln_binom_pmf(t_alt, t_depth, err) + ln_binom_pmf(n_alt, n_depth, err);
    somatic - noise
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::{Locus, Position};
    use crate::pileup::Obs;

    fn col(ref_base: u8, obs_spec: &[(u8, u8)]) -> PileupColumn {
        let obs: Vec<Obs> = obs_spec
            .iter()
            .map(|&(allele, base_qual)| Obs {
                allele,
                base_qual,
                mapq: 60,
                reverse: false,
            })
            .collect();
        PileupColumn {
            locus: Locus {
                contig: 0,
                pos: Position(100),
            },
            ref_base,
            raw_depth: obs.len() as u32,
            obs,
        }
    }

    /// `n` observations of a single allele index, all bq30.
    fn reads(allele: u8, n: usize) -> Vec<(u8, u8)> {
        (0..n).map(|_| (allele, 30u8)).collect()
    }

    #[test]
    fn filters_normal_contamination() {
        // Tumor 15A+5C (alt C, AF 0.25); normal 18A+2C (AF 0.10 > max 0.01) → None.
        let mut t = reads(0, 15);
        t.extend(reads(1, 5));
        let mut n = reads(0, 18);
        n.extend(reads(1, 2));
        let call = call_somatic(&col(b'A', &t), &col(b'A', &n), &SomaticParams::default());
        assert!(call.is_none());
    }

    #[test]
    fn calls_clean_somatic_snv() {
        // Tumor 10A+10C (AF 0.5); normal 20A (AF 0) → somatic C call.
        let mut t = reads(0, 10);
        t.extend(reads(1, 10));
        let n = reads(0, 20);
        let call = call_somatic(&col(b'A', &t), &col(b'A', &n), &SomaticParams::default()).unwrap();
        assert_eq!(call.alt_base, b'C');
        assert_eq!(call.tumor_alt, 10);
        assert_eq!(call.normal_alt, 0);
        assert!((call.tumor_af - 0.5).abs() < 1e-6);
        assert!(call.quality >= SomaticParams::default().min_quality);
    }

    #[test]
    fn respects_min_depth() {
        // Below min_tumor_depth → None even with a clean signal.
        let t = reads(1, 4); // 4 alt reads, depth 4 < 8
        let n = reads(0, 20);
        assert!(call_somatic(&col(b'A', &t), &col(b'A', &n), &SomaticParams::default()).is_none());
    }

    #[test]
    fn llr_math_matches_legacy_reference() {
        // Pin the re-homed math against hand-computed expectations.
        assert_eq!(ln_choose(5, 0), 0.0);
        assert!((ln_choose(5, 2) - 10f64.ln()).abs() < 1e-9); // C(5,2)=10
        assert_eq!(ln_choose(2, 5), f64::NEG_INFINITY);
        // Strong tumor signal, clean normal → positive LLR (somatic favored).
        assert!(somatic_llr(10, 20, 0, 20, 1e-3) > 0.0);
        // No tumor signal → non-positive LLR (noise favored).
        assert!(somatic_llr(0, 20, 0, 20, 1e-3) <= 0.0);
    }
}
