use std::env;
use std::path::PathBuf;
use std::sync::Arc;
use std::time::{SystemTime, UNIX_EPOCH};

use rosalind::genomics::sort_bam_deterministic;
use rust_htslib::bam;
use rust_htslib::bam::record::{Cigar, CigarString, Record};

fn temp_path(name: &str) -> PathBuf {
    let timestamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("time went backwards")
        .as_nanos();
    env::temp_dir().join(format!("rosalind-{name}-{timestamp}.bam"))
}

#[test]
fn somatic_snv_truth_mixture_simple() {
    let tumor_bam = temp_path("tumor_truth");
    let tumor_sorted = temp_path("tumor_truth_sorted");
    let normal_bam = temp_path("normal_truth");
    let normal_sorted = temp_path("normal_truth_sorted");

    // Single contig.
    let mut header = bam::Header::new();
    header.push_record(
        bam::header::HeaderRecord::new(b"SQ")
            .push_tag(b"SN", "chr1")
            .push_tag(b"LN", 1000),
    );

    // Reference is all 'A'.
    let reference: Arc<[u8]> = Arc::from(vec![b'A'; 200].into_boxed_slice());
    let pos0 = 50i64;

    // Tumor: 12 alt reads (C), 8 ref reads (A) at pos 50 => AF 0.6.
    {
        let mut w = bam::Writer::from_path(&tumor_bam, &header, bam::Format::Bam).unwrap();
        for i in 0..12 {
            let mut r = Record::new();
            let cigar = CigarString::from(vec![Cigar::Match(1)]);
            r.set(format!("t_alt_{i}").as_bytes(), Some(&cigar), b"C", b"I");
            r.set_tid(0);
            r.set_pos(pos0);
            r.set_flags(0);
            r.set_mapq(60);
            w.write(&r).unwrap();
        }
        for i in 0..8 {
            let mut r = Record::new();
            let cigar = CigarString::from(vec![Cigar::Match(1)]);
            r.set(format!("t_ref_{i}").as_bytes(), Some(&cigar), b"A", b"I");
            r.set_tid(0);
            r.set_pos(pos0);
            r.set_flags(0);
            r.set_mapq(60);
            w.write(&r).unwrap();
        }
    }

    // Normal: all ref reads (A).
    {
        let mut w = bam::Writer::from_path(&normal_bam, &header, bam::Format::Bam).unwrap();
        for i in 0..20 {
            let mut r = Record::new();
            let cigar = CigarString::from(vec![Cigar::Match(1)]);
            r.set(format!("n_ref_{i}").as_bytes(), Some(&cigar), b"A", b"I");
            r.set_tid(0);
            r.set_pos(pos0);
            r.set_flags(0);
            r.set_mapq(60);
            w.write(&r).unwrap();
        }
    }

    sort_bam_deterministic(&tumor_bam, &tumor_sorted, 1 << 20).unwrap();
    sort_bam_deterministic(&normal_bam, &normal_sorted, 1 << 20).unwrap();

    use rosalind::call::{call_somatic_region, SomaticParams};
    use rosalind::core::{ContigSet, Position};
    use rosalind::io::bam::BamSource;
    use rosalind::pileup::PileupParams;

    let mut contigs = ContigSet::new();
    contigs.push("chr1", 1000);
    let calls = call_somatic_region(
        BamSource::new(&tumor_sorted, &contigs).unwrap(),
        BamSource::new(&normal_sorted, &contigs).unwrap(),
        Arc::clone(&reference),
        0,
        0..200,
        PileupParams::default(),
        &SomaticParams {
            min_tumor_depth: 10,
            min_normal_depth: 10,
            min_tumor_af: 0.2,
            max_normal_af: 0.01,
            min_quality: 0.0,
            seq_error_rate: 1e-3,
        },
    )
    .unwrap();
    assert_eq!(calls.len(), 1);
    assert_eq!(calls[0].0.pos, Position(50));
    assert_eq!(calls[0].1.ref_base, b'A');
    assert_eq!(calls[0].1.alt_base, b'C');

    let _ = std::fs::remove_file(tumor_bam);
    let _ = std::fs::remove_file(tumor_sorted);
    let _ = std::fs::remove_file(normal_bam);
    let _ = std::fs::remove_file(normal_sorted);
}

#[test]
fn alignment_truth_exact_matches_map_to_expected_position() {
    use rosalind::genomics::BWTAligner;

    let reference = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
    let mut aligner = BWTAligner::new(reference).unwrap();

    for start in 0..(reference.len() - 8) {
        let read = &reference[start..start + 8];
        let res = aligner.align_read(read).unwrap();
        assert!(
            res.primary_position.is_some(),
            "expected mapping for start={start}"
        );
        let mapped = res.primary_position.unwrap() as usize;

        // This reference is highly repetitive; accept any exact-match location.
        let mut matches = Vec::new();
        for i in 0..=(reference.len() - read.len()) {
            if &reference[i..i + read.len()] == read {
                matches.push(i);
            }
        }
        assert!(
            matches.contains(&mapped),
            "mapped position {mapped} not in exact-match set {matches:?} for start={start}"
        );
    }
}
