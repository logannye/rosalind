use std::path::Path;
use std::sync::Arc;

use anyhow::Result;
use rust_htslib::bam;
use rust_htslib::bam::Read as BamRead;

use crate::genomics::{BamPileupStream, PileupNode};

/// Somatic variant call (SNV only for now).
#[derive(Debug, Clone, PartialEq)]
pub struct SomaticVariant {
    /// Contig name.
    pub chrom: Arc<str>,
    /// 0-based genomic coordinate.
    pub position: u32, // 0-based
    /// Reference base (uppercase ASCII).
    pub reference: u8,
    /// Alternate base (uppercase ASCII).
    pub alternate: u8,
    /// Tumor depth at the locus.
    pub tumor_depth: u32,
    /// Normal depth at the locus.
    pub normal_depth: u32,
    /// Tumor alternate allele count.
    pub tumor_alt_count: u32,
    /// Normal alternate allele count.
    pub normal_alt_count: u32,
    /// Tumor alt allele fraction.
    pub tumor_af: f32,
    /// Normal alt allele fraction.
    pub normal_af: f32,
    /// Phred-like quality score (heuristic).
    pub quality: f32,
    /// FILTER column value (currently `PASS` only).
    pub filter: &'static str,
}

/// Somatic indel call (simple insertion/deletion only).
#[derive(Debug, Clone, PartialEq)]
pub struct SomaticIndel {
    /// Contig name.
    pub chrom: Arc<str>,
    /// 0-based genomic coordinate (VCF POS-1) anchored at the left base.
    pub position: u32,
    /// Reference allele (ASCII A/C/G/T/N), including anchor base.
    pub reference: Vec<u8>,
    /// Alternate allele (ASCII A/C/G/T/N), including anchor base.
    pub alternate: Vec<u8>,
    /// Tumor supporting read count.
    pub tumor_support: u32,
    /// Normal supporting read count.
    pub normal_support: u32,
    /// Heuristic quality.
    pub quality: f32,
    /// FILTER column value.
    pub filter: &'static str,
}

/// Configuration for tumor/normal SNV calling.
#[derive(Debug, Clone)]
pub struct SomaticCallerConfig {
    /// Minimum tumor depth to consider a site.
    pub min_tumor_depth: u32,
    /// Minimum normal depth to consider a site.
    pub min_normal_depth: u32,
    /// Minimum tumor alt allele fraction.
    pub min_tumor_af: f32,
    /// Maximum allowed normal alt allele fraction.
    pub max_normal_af: f32,
    /// Minimum quality to emit a call.
    pub min_quality: f32,
    /// Sequencing error rate used for the noise model.
    pub sequencing_error_rate: f64,
    /// Minimum supporting reads for indel calls.
    pub min_indel_support: u32,
}

impl Default for SomaticCallerConfig {
    fn default() -> Self {
        Self {
            min_tumor_depth: 8,
            min_normal_depth: 8,
            min_tumor_af: 0.05,
            max_normal_af: 0.01,
            min_quality: 10.0,
            sequencing_error_rate: 1e-3,
            min_indel_support: 3,
        }
    }
}

/// A deterministic tumor/normal somatic SNV caller.
#[derive(Debug, Clone)]
pub struct SomaticCaller {
    cfg: SomaticCallerConfig,
}

impl SomaticCaller {
    /// Construct a new caller with the supplied configuration.
    pub fn new(cfg: SomaticCallerConfig) -> Self {
        Self { cfg }
    }

    /// Stream somatic SNVs from coordinate-sorted BAMs over a region.
    ///
    /// Current limitations:
    /// - SNVs only
    /// - pileup streamer currently supports simple match-only alignments
    pub fn call_snvs_from_sorted_bams(
        &self,
        chrom: Arc<str>,
        reference: Arc<[u8]>,
        region_start: u32,
        tumor_bam: impl AsRef<Path>,
        normal_bam: impl AsRef<Path>,
    ) -> Result<Vec<SomaticVariant>> {
        let region = region_start..(region_start + reference.len() as u32);
        let mut tumor = BamPileupStream::new(tumor_bam, Arc::clone(&chrom), region.clone())?;
        let mut normal = BamPileupStream::new(normal_bam, Arc::clone(&chrom), region.clone())?;

        let mut t_next = tumor.next().transpose()?;
        let mut n_next = normal.next().transpose()?;

        let mut out = Vec::new();
        while t_next.is_some() || n_next.is_some() {
            let pos = match (&t_next, &n_next) {
                (Some(t), Some(n)) => t.position.min(n.position),
                (Some(t), None) => t.position,
                (None, Some(n)) => n.position,
                (None, None) => break,
            };

            let t_node = if t_next.as_ref().map(|n| n.position) == Some(pos) {
                let node = t_next.take();
                t_next = tumor.next().transpose()?;
                node
            } else {
                None
            };
            let n_node = if n_next.as_ref().map(|n| n.position) == Some(pos) {
                let node = n_next.take();
                n_next = normal.next().transpose()?;
                node
            } else {
                None
            };

            let offset = (pos - region_start) as usize;
            if offset >= reference.len() {
                continue;
            }
            let ref_base = reference[offset];

            if let Some(call) =
                self.call_position(&chrom, pos, ref_base, t_node.as_ref(), n_node.as_ref())
            {
                out.push(call);
            }
        }

        // Deterministic ordering.
        out.sort_by(|a, b| {
            a.chrom
                .as_ref()
                .cmp(b.chrom.as_ref())
                .then_with(|| a.position.cmp(&b.position))
                .then_with(|| a.reference.cmp(&b.reference))
                .then_with(|| a.alternate.cmp(&b.alternate))
        });

        Ok(out)
    }

    /// Call simple indels (insertions/deletions) from tumor/normal BAMs.
    ///
    /// This is a conservative first implementation:
    /// - counts indel events implied by CIGAR operations
    /// - anchors alleles in VCF-style left-anchored representation
    /// - requires strong tumor support and near-zero normal support
    pub fn call_indels_from_sorted_bams(
        &self,
        chrom: Arc<str>,
        reference: Arc<[u8]>,
        region_start: u32,
        tumor_bam: impl AsRef<Path>,
        normal_bam: impl AsRef<Path>,
    ) -> Result<Vec<SomaticIndel>> {
        let tumor_counts =
            collect_indel_events(tumor_bam.as_ref(), &chrom, &reference, region_start)?;
        let normal_counts =
            collect_indel_events(normal_bam.as_ref(), &chrom, &reference, region_start)?;

        let mut out = Vec::new();
        for (key, t_support) in tumor_counts.iter() {
            let n_support = *normal_counts.get(key).unwrap_or(&0);
            if *t_support < self.cfg.min_indel_support {
                continue;
            }
            if n_support > 0 {
                continue;
            }
            out.push(SomaticIndel {
                chrom: Arc::clone(&chrom),
                position: key.pos,
                reference: key.ref_allele.clone(),
                alternate: key.alt_allele.clone(),
                tumor_support: *t_support,
                normal_support: n_support,
                quality: 60.0,
                filter: "PASS",
            });
        }

        out.sort_by(|a, b| {
            a.chrom
                .as_ref()
                .cmp(b.chrom.as_ref())
                .then_with(|| a.position.cmp(&b.position))
                .then_with(|| a.reference.cmp(&b.reference))
                .then_with(|| a.alternate.cmp(&b.alternate))
        });

        Ok(out)
    }

    fn call_position(
        &self,
        chrom: &Arc<str>,
        position: u32,
        reference_base: u8,
        tumor: Option<&PileupNode>,
        normal: Option<&PileupNode>,
    ) -> Option<SomaticVariant> {
        let (t_depth, t_counts) = pileup_counts(tumor);
        let (n_depth, n_counts) = pileup_counts(normal);

        if t_depth < self.cfg.min_tumor_depth || n_depth < self.cfg.min_normal_depth {
            return None;
        }

        let ref_idx = base_index(reference_base)?;
        let (alt_idx, t_alt) = best_alt(ref_idx, &t_counts)?;
        if t_alt == 0 {
            return None;
        }
        let n_alt = n_counts[alt_idx];

        let t_af = t_alt as f32 / t_depth as f32;
        let n_af = if n_depth > 0 {
            n_alt as f32 / n_depth as f32
        } else {
            0.0
        };

        if t_af < self.cfg.min_tumor_af {
            return None;
        }
        if n_af > self.cfg.max_normal_af {
            return None;
        }

        let llr = somatic_llr(
            t_alt,
            t_depth,
            n_alt,
            n_depth,
            self.cfg.sequencing_error_rate,
        );
        let qual = (llr / std::f64::consts::LN_10 * 10.0).max(0.0).min(200.0) as f32;
        if qual < self.cfg.min_quality {
            return None;
        }

        Some(SomaticVariant {
            chrom: Arc::clone(chrom),
            position,
            reference: reference_base.to_ascii_uppercase(),
            alternate: idx_to_base(alt_idx),
            tumor_depth: t_depth,
            normal_depth: n_depth,
            tumor_alt_count: t_alt,
            normal_alt_count: n_alt,
            tumor_af: t_af,
            normal_af: n_af,
            quality: qual,
            filter: "PASS",
        })
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct IndelKey {
    pos: u32,
    ref_allele: Vec<u8>,
    alt_allele: Vec<u8>,
}

fn collect_indel_events(
    bam_path: &Path,
    chrom: &Arc<str>,
    reference: &[u8],
    region_start: u32,
) -> Result<std::collections::BTreeMap<IndelKey, u32>> {
    let mut reader = bam::Reader::from_path(bam_path)?;
    let header = reader.header().to_owned();

    let mut counts: std::collections::BTreeMap<IndelKey, u32> = std::collections::BTreeMap::new();

    for rec in reader.records() {
        let rec = rec?;
        if rec.is_unmapped() {
            continue;
        }
        let tid = rec.tid();
        if tid < 0 {
            continue;
        }
        let name = header.tid2name(tid as u32);
        let name = std::str::from_utf8(name).unwrap_or("");
        if name != chrom.as_ref() {
            continue;
        }

        let start0 = rec.pos();
        if start0 < 0 {
            continue;
        }
        let mut ref_pos = start0 as u32;

        let mut seq: Vec<u8> = rec
            .seq()
            .as_bytes()
            .iter()
            .map(|b| b.to_ascii_uppercase())
            .collect();
        if rec.is_reverse() {
            seq = seq.into_iter().rev().map(complement).collect();
        }

        let mut read_i: usize = 0;
        for op in rec.cigar().iter() {
            match *op {
                bam::record::Cigar::Match(l)
                | bam::record::Cigar::Equal(l)
                | bam::record::Cigar::Diff(l) => {
                    ref_pos = ref_pos.saturating_add(l);
                    read_i = read_i.saturating_add(l as usize);
                }
                bam::record::Cigar::Ins(l) => {
                    let l_usize = l as usize;
                    if ref_pos == 0 {
                        read_i = read_i.saturating_add(l_usize);
                        continue;
                    }
                    let anchor_pos = ref_pos - 1;
                    let offset = (anchor_pos - region_start) as usize;
                    if offset >= reference.len() {
                        read_i = read_i.saturating_add(l_usize);
                        continue;
                    }
                    let anchor = reference[offset].to_ascii_uppercase();
                    let inserted = seq.get(read_i..read_i + l_usize).unwrap_or(&[]);
                    let ref_allele = vec![anchor];
                    let mut alt_allele = vec![anchor];
                    alt_allele.extend_from_slice(inserted);

                    *counts
                        .entry(IndelKey {
                            pos: anchor_pos,
                            ref_allele,
                            alt_allele,
                        })
                        .or_insert(0) += 1;

                    read_i = read_i.saturating_add(l_usize);
                }
                bam::record::Cigar::Del(l) => {
                    if ref_pos == 0 {
                        ref_pos = ref_pos.saturating_add(l);
                        continue;
                    }
                    let anchor_pos = ref_pos - 1;
                    let offset = (anchor_pos - region_start) as usize;
                    let del_len = l as usize;
                    if offset + 1 + del_len > reference.len() {
                        ref_pos = ref_pos.saturating_add(l);
                        continue;
                    }
                    let ref_allele = reference[offset..offset + 1 + del_len]
                        .iter()
                        .map(|b| b.to_ascii_uppercase())
                        .collect::<Vec<u8>>();
                    let alt_allele = vec![reference[offset].to_ascii_uppercase()];

                    *counts
                        .entry(IndelKey {
                            pos: anchor_pos,
                            ref_allele,
                            alt_allele,
                        })
                        .or_insert(0) += 1;

                    ref_pos = ref_pos.saturating_add(l);
                }
                bam::record::Cigar::SoftClip(l) => {
                    read_i = read_i.saturating_add(l as usize);
                }
                bam::record::Cigar::HardClip(_) => {}
                bam::record::Cigar::RefSkip(l) => {
                    ref_pos = ref_pos.saturating_add(l);
                }
                bam::record::Cigar::Pad(_) => {}
            }
        }
    }

    Ok(counts)
}

fn pileup_counts(node: Option<&PileupNode>) -> (u32, [u32; 4]) {
    match node {
        None => (0, [0; 4]),
        Some(n) => (n.depth, n.base_counts),
    }
}

fn best_alt(ref_idx: usize, counts: &[u32; 4]) -> Option<(usize, u32)> {
    let mut best_idx = None;
    let mut best = 0u32;
    for (i, &c) in counts.iter().enumerate() {
        if i == ref_idx {
            continue;
        }
        if c > best {
            best = c;
            best_idx = Some(i);
        }
    }
    best_idx.map(|i| (i, best))
}

fn base_index(base: u8) -> Option<usize> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

fn idx_to_base(idx: usize) -> u8 {
    match idx {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        3 => b'T',
        _ => b'N',
    }
}

fn complement(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => b'N',
    }
}

fn ln_factorial(n: u32) -> f64 {
    // For WGS depths, n is small; compute exactly and deterministically.
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
    // Compare somatic vs noise-only:
    // - under somatic: tumor alt fraction approximated by MLE, normal alt fraction ~ err
    // - under noise: both tumor and normal alt fraction ~ err
    let t_p_hat = (t_alt as f64 / t_depth.max(1) as f64).clamp(err, 1.0 - 1e-6);
    let somatic = ln_binom_pmf(t_alt, t_depth, t_p_hat) + ln_binom_pmf(n_alt, n_depth, err);
    let noise = ln_binom_pmf(t_alt, t_depth, err) + ln_binom_pmf(n_alt, n_depth, err);
    somatic - noise
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genomics::PileupNode;

    #[test]
    fn somatic_caller_filters_normal_contamination() {
        let cfg = SomaticCallerConfig {
            min_tumor_depth: 10,
            min_normal_depth: 10,
            min_tumor_af: 0.1,
            max_normal_af: 0.01,
            min_quality: 0.0,
            sequencing_error_rate: 1e-3,
            min_indel_support: 3,
        };
        let caller = SomaticCaller::new(cfg);
        let chrom: Arc<str> = Arc::from("chr1");

        let mut tumor = PileupNode::new(100);
        for _ in 0..15 {
            tumor.observe(0, 30); // A
        }
        for _ in 0..5 {
            tumor.observe(1, 30); // C
        }

        let mut normal = PileupNode::new(100);
        for _ in 0..18 {
            normal.observe(0, 30); // A
        }
        for _ in 0..2 {
            normal.observe(1, 30); // C (10% -> should fail max_normal_af)
        }

        let call = caller.call_position(&chrom, 100, b'A', Some(&tumor), Some(&normal));
        assert!(call.is_none());
    }
}
