use std::collections::BTreeMap;

use super::{normalize_variant, BedIndex, NormalizedVariant, VcfVariant, Zygosity};

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
/// Variant category used for stratified metrics.
pub enum VariantType {
    /// Single-nucleotide variant.
    Snv,
    /// Insertion or deletion (including 1bp).
    Indel,
}

/// Summary of truth comparison between a called set and a truth set.
#[derive(Debug, Clone)]
pub struct ComparisonReport {
    /// Total truth variants considered (after masking and normalization).
    pub total_truth: usize,
    /// Total called variants considered (after masking and normalization).
    pub total_calls: usize,
    /// True positives (calls present in truth).
    pub true_positive: usize,
    /// False positives (calls absent in truth).
    pub false_positive: usize,
    /// False negatives (truth variants absent in calls).
    pub false_negative: usize,
    /// Per-type counts (tp, fp, fn).
    pub by_type: BTreeMap<VariantType, (usize, usize, usize)>, // (tp, fp, fn)
    /// True positives whose genotype matched truth (both genotypes known).
    pub genotype_concordant: usize,
    /// True positives whose genotype disagreed with truth (both known) — a genotype error.
    pub genotype_discordant: usize,
    /// True positives where either genotype was absent (unscored for concordance).
    pub genotype_unknown: usize,
}

impl ComparisonReport {
    /// Genotype concordance = concordant / (concordant + discordant), over the scorable
    /// true positives. `0.0` when none are scorable.
    pub fn genotype_concordance(&self) -> f64 {
        let scorable = self.genotype_concordant + self.genotype_discordant;
        if scorable == 0 {
            0.0
        } else {
            self.genotype_concordant as f64 / scorable as f64
        }
    }

    /// Precision/PPV = TP / (TP + FP).
    pub fn precision(&self) -> f64 {
        if self.true_positive + self.false_positive == 0 {
            0.0
        } else {
            self.true_positive as f64 / (self.true_positive + self.false_positive) as f64
        }
    }

    /// Recall/sensitivity = TP / (TP + FN).
    pub fn recall(&self) -> f64 {
        if self.true_positive + self.false_negative == 0 {
            0.0
        } else {
            self.true_positive as f64 / (self.true_positive + self.false_negative) as f64
        }
    }
}

/// Compare two VCF callsets against a per-contig reference map, optionally masked
/// by BED regions. Each variant is normalized against the sequence of ITS OWN
/// contig (`references[v.chrom]`) — loading only the first FASTA record and
/// applying it to every contig silently miscompares (or crashes) on any
/// multi-contig benchmark, e.g. a real GIAB run.
pub fn compare_callsets(
    references: &BTreeMap<String, Vec<u8>>,
    calls: &[VcfVariant],
    truth: &[VcfVariant],
    bed: Option<&BedIndex>,
) -> Result<ComparisonReport, anyhow::Error> {
    // Key = normalized locus+allele (the detection identity); value = decomposed zygosity.
    let mut calls_map: BTreeMap<NormalizedVariant, Option<Zygosity>> = BTreeMap::new();
    let mut truth_map: BTreeMap<NormalizedVariant, Option<Zygosity>> = BTreeMap::new();

    for v in calls {
        if let Some(bed) = bed {
            if !bed.contains(&v.chrom, v.pos0) {
                continue;
            }
        }
        calls_map.insert(
            normalize_variant(reference_for(references, v)?, v)?,
            v.genotype,
        );
    }
    for v in truth {
        if let Some(bed) = bed {
            if !bed.contains(&v.chrom, v.pos0) {
                continue;
            }
        }
        truth_map.insert(
            normalize_variant(reference_for(references, v)?, v)?,
            v.genotype,
        );
    }

    let mut tp = 0usize;
    let mut fp = 0usize;
    let mut fn_ = 0usize;
    let mut gt_concordant = 0usize;
    let mut gt_discordant = 0usize;
    let mut gt_unknown = 0usize;

    let mut by_type: BTreeMap<VariantType, (usize, usize, usize)> = BTreeMap::new();

    for (v, call_gt) in calls_map.iter() {
        let vt = variant_type(v);
        if let Some(truth_gt) = truth_map.get(v) {
            tp += 1;
            by_type.entry(vt).or_insert((0, 0, 0)).0 += 1;
            match (call_gt, truth_gt) {
                (Some(c), Some(t)) if c == t => gt_concordant += 1,
                (Some(_), Some(_)) => gt_discordant += 1,
                _ => gt_unknown += 1,
            }
        } else {
            fp += 1;
            by_type.entry(vt).or_insert((0, 0, 0)).1 += 1;
        }
    }

    for v in truth_map.keys() {
        if !calls_map.contains_key(v) {
            fn_ += 1;
            by_type.entry(variant_type(v)).or_insert((0, 0, 0)).2 += 1;
        }
    }

    Ok(ComparisonReport {
        total_truth: truth_map.len(),
        total_calls: calls_map.len(),
        true_positive: tp,
        false_positive: fp,
        false_negative: fn_,
        by_type,
        genotype_concordant: gt_concordant,
        genotype_discordant: gt_discordant,
        genotype_unknown: gt_unknown,
    })
}

/// The reference sequence for a variant's contig, or a clear error naming the
/// missing contig and the available ones (the common contig-naming mismatch).
fn reference_for<'a>(
    references: &'a BTreeMap<String, Vec<u8>>,
    v: &VcfVariant,
) -> Result<&'a [u8], anyhow::Error> {
    references.get(&v.chrom).map(Vec::as_slice).ok_or_else(|| {
        let available: Vec<&str> = references.keys().map(String::as_str).collect();
        anyhow::anyhow!(
            "variant on contig '{}' (pos {}) has no matching sequence in the reference FASTA — \
             check the contig naming scheme (e.g. UCSC 'chr1' vs Ensembl '1'). \
             Reference contigs: [{}]",
            v.chrom,
            v.pos0 + 1,
            available.join(", ")
        )
    })
}

fn variant_type(v: &NormalizedVariant) -> VariantType {
    if v.reference.len() == 1 && v.alternate.len() == 1 {
        VariantType::Snv
    } else {
        VariantType::Indel
    }
}

#[cfg(test)]
mod gt_concordance_tests {
    use super::*;
    use crate::genomics::eval::Zygosity;

    fn vv(pos1: u32, refb: &str, alt: &str, gt: Option<Zygosity>) -> VcfVariant {
        VcfVariant {
            chrom: "chr1".to_string(),
            pos0: pos1 - 1,
            reference: refb.as_bytes().to_vec(),
            alternate: alt.as_bytes().to_vec(),
            qual: None,
            filter: ".".to_string(),
            info: BTreeMap::new(),
            genotype: gt,
        }
    }

    fn refs() -> BTreeMap<String, Vec<u8>> {
        BTreeMap::from([("chr1".to_string(), b"ACGTACGTAC".to_vec())])
    }

    #[test]
    fn het_called_as_hom_is_a_detection_tp_but_a_genotype_error() {
        let truth = vec![vv(5, "A", "C", Some(Zygosity::Het))];
        let calls = vec![vv(5, "A", "C", Some(Zygosity::Hom))]; // right site, wrong genotype
        let r = compare_callsets(&refs(), &calls, &truth, None).unwrap();
        assert_eq!(r.true_positive, 1, "detection still matches");
        assert_eq!(r.false_positive, 0);
        assert_eq!(r.false_negative, 0);
        assert_eq!(r.genotype_concordant, 0);
        assert_eq!(r.genotype_discordant, 1);
        assert_eq!(r.genotype_concordance(), 0.0);
    }

    #[test]
    fn matching_genotype_is_concordant() {
        let truth = vec![vv(5, "A", "C", Some(Zygosity::Het))];
        let calls = vec![vv(5, "A", "C", Some(Zygosity::Het))];
        let r = compare_callsets(&refs(), &calls, &truth, None).unwrap();
        assert_eq!((r.genotype_concordant, r.genotype_discordant), (1, 0));
        assert_eq!(r.genotype_concordance(), 1.0);
    }

    #[test]
    fn unknown_genotype_is_counted_separately_not_as_error() {
        let truth = vec![vv(5, "A", "C", None)]; // truth without GT
        let calls = vec![vv(5, "A", "C", Some(Zygosity::Het))];
        let r = compare_callsets(&refs(), &calls, &truth, None).unwrap();
        assert_eq!(r.true_positive, 1);
        assert_eq!(
            (
                r.genotype_concordant,
                r.genotype_discordant,
                r.genotype_unknown
            ),
            (0, 0, 1)
        );
    }
}
