use std::collections::{BTreeMap, BTreeSet};

use super::{normalize_variant, BedIndex, NormalizedVariant, VcfVariant};

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
}

impl ComparisonReport {
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
    let mut calls_set: BTreeSet<NormalizedVariant> = BTreeSet::new();
    let mut truth_set: BTreeSet<NormalizedVariant> = BTreeSet::new();

    for v in calls {
        if let Some(bed) = bed {
            if !bed.contains(&v.chrom, v.pos0) {
                continue;
            }
        }
        calls_set.insert(normalize_variant(reference_for(references, v)?, v)?);
    }
    for v in truth {
        if let Some(bed) = bed {
            if !bed.contains(&v.chrom, v.pos0) {
                continue;
            }
        }
        truth_set.insert(normalize_variant(reference_for(references, v)?, v)?);
    }

    let mut tp = 0usize;
    let mut fp = 0usize;
    let mut fn_ = 0usize;

    let mut by_type: BTreeMap<VariantType, (usize, usize, usize)> = BTreeMap::new();

    for v in calls_set.iter() {
        let vt = variant_type(v);
        if truth_set.contains(v) {
            tp += 1;
            let e = by_type.entry(vt).or_insert((0, 0, 0));
            e.0 += 1;
        } else {
            fp += 1;
            let e = by_type.entry(vt).or_insert((0, 0, 0));
            e.1 += 1;
        }
    }

    for v in truth_set.iter() {
        let vt = variant_type(v);
        if !calls_set.contains(v) {
            fn_ += 1;
            let e = by_type.entry(vt).or_insert((0, 0, 0));
            e.2 += 1;
        }
    }

    Ok(ComparisonReport {
        total_truth: truth_set.len(),
        total_calls: calls_set.len(),
        true_positive: tp,
        false_positive: fp,
        false_negative: fn_,
        by_type,
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
