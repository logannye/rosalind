//! The germline diploid genotype-likelihood model: per-observation base-quality
//! error integrated into [0/0, 0/1, 1/1] log-likelihoods, then a prior, then a
//! calibrated, abstention-aware call.

use crate::call::types::{Filter, GermlineCall, GermlineParams, Genotype, ACGT};
use crate::core::allele_index;
use crate::pileup::PileupColumn;

const LN10: f64 = std::f64::consts::LN_10;

/// Per-genotype log-likelihoods plus the biallelic context for a column.
struct SiteLikelihoods {
    /// Natural-log likelihoods, ordered [0/0, 0/1, 1/1].
    log_l: [f64; 3],
    /// Index (0..=3) of the alt allele.
    alt_idx: usize,
    /// Allelic depths [ref, alt].
    ad: [u32; 2],
    /// Callable depth.
    dp: u32,
}

/// Accumulate diploid genotype log-likelihoods from per-observation base
/// qualities. Returns `None` when the reference base is non-callable or no
/// non-reference allele is observed (nothing to call).
fn site_likelihoods(column: &PileupColumn) -> Option<SiteLikelihoods> {
    let ref_idx = allele_index(column.ref_base)?;
    let counts = column.allele_counts();

    // Alt = most-supported non-reference allele; ties → lowest index.
    let mut alt_idx: Option<usize> = None;
    let mut best = 0u32;
    for i in 0..4 {
        if i == ref_idx {
            continue;
        }
        if counts[i] > best {
            best = counts[i];
            alt_idx = Some(i);
        }
    }
    let alt_idx = alt_idx?;
    if best == 0 {
        return None;
    }

    let mut log_l = [0.0f64; 3];
    for o in &column.obs {
        // ε from Phred base quality, capped at 0.75 so 1−ε stays positive (and
        // the log finite) even at q=0 — no base is worse than a random draw.
        let eps = (10f64.powf(-(o.base_qual as f64) / 10.0)).min(0.75);
        let a = o.allele as usize;
        let p_ref = if a == ref_idx { 1.0 - eps } else { eps / 3.0 };
        let p_alt = if a == alt_idx { 1.0 - eps } else { eps / 3.0 };
        log_l[0] += p_ref.ln();
        log_l[1] += (0.5 * p_ref + 0.5 * p_alt).ln();
        log_l[2] += p_alt.ln();
    }

    Some(SiteLikelihoods {
        log_l,
        alt_idx,
        ad: [counts[ref_idx], counts[alt_idx]],
        dp: column.depth(),
    })
}

/// Call a germline genotype from a pileup column. Returns `None` (abstains) when
/// the most-probable genotype is homozygous reference — honest by default: thin
/// or absent evidence is a no-call, not a confident reference assertion turned
/// variant.
pub fn call_germline(column: &PileupColumn, params: &GermlineParams) -> Option<GermlineCall> {
    let sl = site_likelihoods(column)?;

    let theta = params.heterozygosity;
    let log_prior = [
        (1.0 - 1.5 * theta).ln(),
        theta.ln(),
        (theta / 2.0).ln(),
    ];
    let log_post = [
        sl.log_l[0] + log_prior[0],
        sl.log_l[1] + log_prior[1],
        sl.log_l[2] + log_prior[2],
    ];

    // GT = argmax posterior; hom-ref → abstain.
    let gt_idx = argmax3(&log_post);
    if gt_idx == 0 {
        return None;
    }

    // PL: −10·log10 L(g), re-normalized so the max-likelihood genotype is 0,
    // each capped at 255. Order [0/0, 0/1, 1/1].
    let max_log_l = sl.log_l.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let mut pl = [0u32; 3];
    for g in 0..3 {
        let phred = -10.0 * (sl.log_l[g] - max_log_l) / LN10;
        pl[g] = phred.round().min(255.0) as u32;
    }

    // GQ = second-smallest PL, capped 99.
    let mut sorted = pl;
    sorted.sort_unstable();
    let gq = sorted[1].min(99) as u8;

    // QUAL = −10·log10 P(0/0 | data), via a numerically stable log-sum-exp.
    let max_lp = log_post.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let log_z = max_lp
        + log_post
            .iter()
            .map(|lp| (lp - max_lp).exp())
            .sum::<f64>()
            .ln();
    let log_p_homref = log_post[0] - log_z; // ≤ 0
    let qual = (-10.0 * log_p_homref / LN10).max(0.0);

    let genotype = if gt_idx == 1 {
        Genotype::Het
    } else {
        Genotype::HomAlt
    };
    let filter = if sl.dp < params.min_depth {
        Filter::LowDepth
    } else if qual < params.min_qual {
        Filter::LowQual
    } else {
        Filter::Pass
    };

    Some(GermlineCall {
        genotype,
        alt_base: ACGT[sl.alt_idx],
        qual,
        gq,
        pl,
        ad: sl.ad,
        dp: sl.dp,
        filter,
    })
}

/// Index of the maximum of three values; ties resolve to the lowest index.
fn argmax3(v: &[f64; 3]) -> usize {
    let mut best = 0;
    for i in 1..3 {
        if v[i] > v[best] {
            best = i;
        }
    }
    best
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::{Locus, Position};
    use crate::pileup::Obs;

    /// Build a column at chr0:100 with `ref_base` and observations given as
    /// `(allele_index, base_qual)` pairs (all forward strand, mapq 60).
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

    #[test]
    fn likelihoods_rank_genotypes_correctly() {
        // 15 ref (A) + 15 alt (C), all bq30 → het is the most-likely genotype.
        let het: Vec<(u8, u8)> = (0..15)
            .map(|_| (0u8, 30u8))
            .chain((0..15).map(|_| (1u8, 30u8)))
            .collect();
        let sl = site_likelihoods(&col(b'A', &het)).unwrap();
        assert_eq!(sl.alt_idx, 1); // C
        assert_eq!(sl.ad, [15, 15]);
        assert_eq!(sl.dp, 30);
        // het (index 1) has the largest log-likelihood
        assert!(sl.log_l[1] > sl.log_l[0]);
        assert!(sl.log_l[1] > sl.log_l[2]);

        // 20 alt only → hom-alt dominates.
        let homalt: Vec<(u8, u8)> = (0..20).map(|_| (1u8, 30u8)).collect();
        let sl = site_likelihoods(&col(b'A', &homalt)).unwrap();
        assert!(sl.log_l[2] > sl.log_l[1]);
        assert!(sl.log_l[2] > sl.log_l[0]);
    }

    #[test]
    fn no_alt_or_non_callable_ref_yields_none() {
        // All reference → nothing to call.
        let allref: Vec<(u8, u8)> = (0..20).map(|_| (0u8, 30u8)).collect();
        assert!(site_likelihoods(&col(b'A', &allref)).is_none());
        // Non-callable reference base.
        assert!(site_likelihoods(&col(b'N', &[(1, 30), (1, 30)])).is_none());
    }

    #[test]
    fn hom_ref_site_abstains() {
        // No alt observed → no variant record.
        let allref: Vec<(u8, u8)> = (0..20).map(|_| (0u8, 30u8)).collect();
        assert!(call_germline(&col(b'A', &allref), &GermlineParams::default()).is_none());
    }

    #[test]
    fn thin_alt_evidence_abstains() {
        // 3 ref + 1 alt at bq30: a lone alt read cannot overcome the θ=1e-3
        // prior, so the calibrated model abstains (no-call) rather than overcall.
        let thin = [(0u8, 30u8), (0, 30), (0, 30), (1, 30)];
        assert!(call_germline(&col(b'A', &thin), &GermlineParams::default()).is_none());
    }

    #[test]
    fn clear_het_passes() {
        let het: Vec<(u8, u8)> = (0..15)
            .map(|_| (0u8, 30u8))
            .chain((0..15).map(|_| (1u8, 30u8)))
            .collect();
        let call = call_germline(&col(b'A', &het), &GermlineParams::default()).unwrap();
        assert_eq!(call.genotype, Genotype::Het);
        assert_eq!(call.alt_base, b'C');
        assert_eq!(call.ad, [15, 15]);
        assert_eq!(call.dp, 30);
        assert_eq!(call.pl[1], 0); // het is the most-likely-by-likelihood genotype
        assert!(call.pl[0] > 0 && call.pl[2] > 0);
        assert_eq!(call.filter, Filter::Pass);
        assert!(call.gq > 0);
    }

    #[test]
    fn clear_hom_alt_passes() {
        let homalt: Vec<(u8, u8)> = (0..20).map(|_| (1u8, 30u8)).collect();
        let call = call_germline(&col(b'A', &homalt), &GermlineParams::default()).unwrap();
        assert_eq!(call.genotype, Genotype::HomAlt);
        assert_eq!(call.pl[2], 0);
        assert_eq!(call.filter, Filter::Pass);
    }

    #[test]
    fn shallow_variant_is_flagged_low_depth_not_dropped() {
        // 3 clean alt reads, 0 ref: unambiguously hom-alt, but DP=3 < 8 → emitted
        // with the LowDepth flag (a real variant, too shallow to fully trust).
        let shallow: Vec<(u8, u8)> = (0..3).map(|_| (1u8, 30u8)).collect();
        let call = call_germline(&col(b'A', &shallow), &GermlineParams::default()).unwrap();
        assert_eq!(call.genotype, Genotype::HomAlt);
        assert_eq!(call.dp, 3);
        assert_eq!(call.filter, Filter::LowDepth);
    }

    #[test]
    fn qual_increases_with_evidence() {
        let mk = |n: usize| -> f64 {
            let obs: Vec<(u8, u8)> = (0..n)
                .map(|_| (0u8, 30u8))
                .chain((0..n).map(|_| (1u8, 30u8)))
                .collect();
            call_germline(&col(b'A', &obs), &GermlineParams::default())
                .unwrap()
                .qual
        };
        // Same 0.5 alt fraction, deeper site → higher site-is-variant QUAL.
        assert!(mk(20) > mk(8));
    }
}
