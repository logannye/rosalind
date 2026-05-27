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
}
