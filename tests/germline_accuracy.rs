//! Simulated-truth germline ACCURACY harness.
//!
//! Inject a known set of het/hom SNVs into a reference, build diploid-sampled
//! coordinate-sorted BAMs directly (bypassing the toy aligner, so this measures
//! the CALLER not alignment), call with `call_germline_whole_genome`, and measure
//! detection precision/recall/F1 against the truth via the existing
//! `compare_callsets`. Deterministic (fixed-seed LCG, no `rand` dependency) so the
//! numbers are stable and the gate is non-flaky. Two regimes: a clean baseline
//! (40x / 0.5% error) and a low-coverage stress case (12x / 1.5% error). Also
//! exercises the unbiased depth-cap fix: a deep het site is recovered under a
//! depth cap below local depth.

use std::collections::BTreeMap;
use std::path::PathBuf;
use std::time::{SystemTime, UNIX_EPOCH};

use rust_htslib::bam;
use rust_htslib::bam::record::{Cigar, CigarString, Record};

use rosalind::call::Filter;
use rosalind::core::Locus;
use rosalind::genomics::{
    compare_callsets, sort_bam_deterministic, ComparisonReport, GenomeIndex, IndexReader,
    IndexWriter, VcfVariant,
};
use rosalind::{
    call_germline_whole_genome, GermlineCall, GermlineParams, PileupParams, StreamingBamSource,
};

const NUCS: [u8; 4] = [b'A', b'C', b'G', b'T'];
const N: usize = 20_000;
const READ_LEN: usize = 100;

/// Deterministic LCG (PCG-style multiplier) — reproducible, no `rand` dependency.
struct Lcg(u64);
impl Lcg {
    fn next_u32(&mut self) -> u32 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 32) as u32
    }
    fn below(&mut self, n: u32) -> u32 {
        self.next_u32() % n
    }
    fn frac(&mut self) -> f64 {
        self.next_u32() as f64 / (u32::MAX as f64 + 1.0)
    }
}

/// A deterministic alt allele distinct from the reference base.
fn alt_of(base: u8) -> u8 {
    match base {
        b'A' => b'C',
        b'C' => b'G',
        b'G' => b'T',
        _ => b'A',
    }
}

#[derive(Clone, Copy)]
struct TruthVar {
    pos: usize,
    refb: u8,
    altb: u8,
    het: bool,
}

fn tmp_dir() -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let d = std::env::temp_dir().join(format!("rosalind-acc-{nanos}"));
    std::fs::create_dir_all(&d).unwrap();
    d
}

struct AccuracyOutcome {
    /// Detection report over ALL emitted variant records (incl. LowQual/LowDepth).
    report_all: ComparisonReport,
    /// Detection report over PASS records only (the GIAB-standard metric).
    report_pass: ComparisonReport,
    deep_site: usize,
    /// Whether the deep-het probe site was recovered as a PASS call.
    deep_called: bool,
}

/// Run the full simulate -> call -> compare pipeline for one regime. Returns the
/// detection report plus whether the deep-het probe site was recovered.
fn run_accuracy(coverage: usize, error_rate: f64, cap: u32, seed: u64) -> AccuracyOutcome {
    let dir = tmp_dir();
    let mut rng = Lcg(seed);

    // 1. Deterministic reference.
    let reference: Vec<u8> = (0..N).map(|_| NUCS[rng.below(4) as usize]).collect();

    // 2. Truth variants: ~40 SNVs, spaced, away from edges; ~60% het / 40% hom-alt.
    let n_vars = 40usize;
    let margin = 200usize;
    let span = N - 2 * margin;
    let step = span / n_vars;
    let mut truth: Vec<TruthVar> = Vec::new();
    for k in 0..n_vars {
        let pos = margin + k * step + rng.below((step / 2) as u32) as usize;
        if pos >= N - margin {
            continue;
        }
        let refb = reference[pos];
        truth.push(TruthVar {
            pos,
            refb,
            altb: alt_of(refb),
            het: rng.frac() < 0.6,
        });
    }
    truth.sort_by_key(|v| v.pos);
    truth.dedup_by_key(|v| v.pos);

    // 3. Two haplotypes. hom-alt -> both haps; het -> hap2 only.
    let mut hap1 = reference.clone();
    let mut hap2 = reference.clone();
    for v in &truth {
        hap2[v.pos] = v.altb;
        if !v.het {
            hap1[v.pos] = v.altb;
        }
    }
    let deep_site = truth
        .iter()
        .find(|v| v.het)
        .map(|v| v.pos)
        .expect("at least one het site");

    // 4. Simulate reads -> BAM.
    let mut header = bam::Header::new();
    {
        let mut sq = bam::header::HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", &"chr1");
        sq.push_tag(b"LN", &(N as i64));
        header.push_record(&sq);
    }
    let raw_bam = dir.join("reads.bam");
    let sorted_bam = dir.join("sorted.bam");
    {
        let mut w = bam::Writer::from_path(&raw_bam, &header, bam::Format::Bam).unwrap();
        let mut rid = 0u64;
        let mut emit =
            |w: &mut bam::Writer, start: usize, hap: &[u8], rng: &mut Lcg, rid: &mut u64| {
                let end = (start + READ_LEN).min(N);
                if end <= start + 20 {
                    return;
                }
                let len = end - start;
                let mut seq: Vec<u8> = hap[start..end].to_vec();
                for b in seq.iter_mut() {
                    if rng.frac() < error_rate {
                        let mut nb = NUCS[rng.below(4) as usize];
                        while nb == *b {
                            nb = NUCS[rng.below(4) as usize];
                        }
                        *b = nb;
                    }
                }
                let cigar = CigarString(vec![Cigar::Match(len as u32)]);
                let qual = vec![40u8; len];
                let mut rec = Record::new();
                let name = format!("r{}", *rid);
                *rid += 1;
                rec.set(name.as_bytes(), Some(&cigar), &seq, &qual);
                rec.set_tid(0);
                rec.set_pos(start as i64);
                rec.set_mapq(60);
                w.write(&rec).unwrap();
            };

        // 4a. Genome-wide ~coverage x at random starts, random haplotype.
        let n_reads = coverage * N / READ_LEN;
        for _ in 0..n_reads {
            let start = rng.below((N - 20) as u32) as usize;
            let pick_hap2 = rng.frac() < 0.5;
            let hap: &[u8] = if pick_hap2 { &hap2 } else { &hap1 };
            emit(&mut w, start, hap, &mut rng, &mut rid);
        }
        // 4b. Deep-het probe: ref reads starting upstream + alt reads starting AT
        // the site, so the cap engages at a deep site (the unbiased reservoir must
        // keep the at-variant alt reads — the old leftmost-arrival cap dropped them).
        let up = deep_site.saturating_sub(READ_LEN - 10);
        for _ in 0..30 {
            emit(&mut w, up, &hap1, &mut rng, &mut rid);
        }
        for _ in 0..30 {
            emit(&mut w, deep_site, &hap2, &mut rng, &mut rid);
        }
    }
    sort_bam_deterministic(&raw_bam, &sorted_bam, 64 << 20).unwrap();

    // 5. Index + call.
    let idx = dir.join("ref.idx");
    let index =
        GenomeIndex::from_named_sequences(&[("chr1".to_string(), reference.clone())]).unwrap();
    IndexWriter::create(&idx)
        .unwrap()
        .write_genome_index(&index)
        .unwrap();
    let reader = IndexReader::open(&idx).unwrap();
    let rv = reader.reference_view().unwrap();
    let contigs = reader.contigs();

    let pp = PileupParams {
        max_depth: Some(cap),
        ..PileupParams::default()
    };
    let gp = GermlineParams::default();
    let src = StreamingBamSource::new(&sorted_bam, contigs).unwrap();
    let mut calls_all: Vec<VcfVariant> = Vec::new();
    let mut calls_pass: Vec<VcfVariant> = Vec::new();
    call_germline_whole_genome(src, &rv, contigs, pp, &gp, &mut |(locus, refb, call): (
        Locus,
        u8,
        GermlineCall,
    )| {
        let v = VcfVariant {
            chrom: "chr1".to_string(),
            pos0: locus.pos.0,
            reference: vec![refb],
            alternate: vec![call.alt_base],
            qual: Some(call.qual as f32),
            filter: format!("{:?}", call.filter),
            info: BTreeMap::new(),
        };
        if matches!(call.filter, Filter::Pass) {
            calls_pass.push(v.clone());
        }
        calls_all.push(v);
        Ok(())
    })
    .unwrap();

    // 6. Compare against truth (detection: pos + ref + alt).
    let truth_vcf: Vec<VcfVariant> = truth
        .iter()
        .map(|v| VcfVariant {
            chrom: "chr1".to_string(),
            pos0: v.pos as u32,
            reference: vec![v.refb],
            alternate: vec![v.altb],
            qual: None,
            filter: ".".to_string(),
            info: BTreeMap::new(),
        })
        .collect();

    let report_all = compare_callsets(&reference, &calls_all, &truth_vcf, None).unwrap();
    let report_pass = compare_callsets(&reference, &calls_pass, &truth_vcf, None).unwrap();
    let deep_called = calls_pass.iter().any(|c| c.pos0 as usize == deep_site);

    std::fs::remove_dir_all(&dir).ok();
    AccuracyOutcome {
        report_all,
        report_pass,
        deep_site,
        deep_called,
    }
}

fn f1(report: &ComparisonReport) -> f64 {
    let (p, r) = (report.precision(), report.recall());
    if p + r == 0.0 {
        0.0
    } else {
        2.0 * p * r / (p + r)
    }
}

fn log_report(label: &str, o: &AccuracyOutcome) {
    for (which, r) in [("PASS", &o.report_pass), ("all", &o.report_all)] {
        eprintln!(
            "germline accuracy [{label}/{which}]: truth={} calls={} tp={} fp={} fn={} precision={:.4} recall={:.4} f1={:.4}",
            r.total_truth,
            r.total_calls,
            r.true_positive,
            r.false_positive,
            r.false_negative,
            r.precision(),
            r.recall(),
            f1(r),
        );
    }
}

#[test]
fn detection_accuracy_clean_baseline_40x() {
    // Clean regime: 40x, 0.5% error, cap 50. Expect near-perfect PASS detection.
    let o = run_accuracy(40, 0.005, 50, 0x9E37_79B9_7F4A_7C15);
    log_report("40x/0.5%", &o);
    assert!(
        o.deep_called,
        "deep het site at {} must be recovered as a PASS call under the cap",
        o.deep_site
    );
    // Measured PASS: precision=1.0000, recall=1.0000 (2026-06-02). Gate a margin below.
    let r = &o.report_pass;
    assert!(
        r.recall() >= 0.97,
        "PASS recall {:.4} below 0.97",
        r.recall()
    );
    assert!(
        r.precision() >= 0.97,
        "PASS precision {:.4} below 0.97",
        r.precision()
    );
}

#[test]
fn detection_accuracy_low_coverage_stress_12x() {
    // Harder regime: 12x, 1.5% error, cap 50 — fewer reads per site and more noise.
    // Shows the harness produces graded, honest numbers and that the PASS filter
    // (min_qual 30, min_depth 8) controls the error-driven false positives that
    // flood the unfiltered (`all`) callset.
    let o = run_accuracy(12, 0.015, 50, 0xD1B5_4A32_D192_ED03);
    log_report("12x/1.5%", &o);
    // Calibrated to the measured PASS run (margin below — see the findings doc).
    let r = &o.report_pass;
    assert!(
        r.recall() >= LOW_COV_R_MIN,
        "PASS recall {:.4} below {LOW_COV_R_MIN}",
        r.recall()
    );
    assert!(
        r.precision() >= LOW_COV_P_MIN,
        "PASS precision {:.4} below {LOW_COV_P_MIN}",
        r.precision()
    );
}

// Calibrated to the measured 12x/1.5% PASS run (recall 0.90, precision 0.68 on
// 2026-06-02) — a margin below, so the gate catches regression without flaking.
const LOW_COV_R_MIN: f64 = 0.85;
const LOW_COV_P_MIN: f64 = 0.60;
