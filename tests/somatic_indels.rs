use std::env;
use std::path::PathBuf;
use std::sync::Arc;
use std::time::{SystemTime, UNIX_EPOCH};

use rosalind::genomics::{SomaticCaller, SomaticCallerConfig};
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
fn somatic_indel_emits_simple_insertion() {
    let tumor_bam = temp_path("tumor");
    let normal_bam = temp_path("normal");

    let mut header = bam::Header::new();
    header.push_record(
        bam::header::HeaderRecord::new(b"SQ")
            .push_tag(b"SN", &"chr1")
            .push_tag(b"LN", &1000),
    );

    // Tumor: one read with 1bp insertion after position 10 (anchor at 10).
    {
        let mut writer = bam::Writer::from_path(&tumor_bam, &header, bam::Format::Bam).unwrap();
        let mut r = Record::new();
        let cigar = CigarString::from(vec![Cigar::Match(1), Cigar::Ins(1), Cigar::Match(3)]);
        r.set(b"t1", Some(&cigar), b"AACGT", b"IIIII");
        r.set_tid(0);
        r.set_pos(11); // after consuming 1M, insertion occurs at ref_pos 12; anchor at 11
        r.set_flags(0);
        r.set_mapq(60);
        writer.write(&r).unwrap();
    }

    // Normal: no indels.
    {
        let mut writer = bam::Writer::from_path(&normal_bam, &header, bam::Format::Bam).unwrap();
        let mut r = Record::new();
        let cigar = CigarString::from(vec![Cigar::Match(4)]);
        r.set(b"n1", Some(&cigar), b"ACGT", b"IIII");
        r.set_tid(0);
        r.set_pos(11);
        r.set_flags(0);
        r.set_mapq(60);
        writer.write(&r).unwrap();
    }

    let cfg = SomaticCallerConfig {
        min_indel_support: 1,
        ..SomaticCallerConfig::default()
    };
    let caller = SomaticCaller::new(cfg);
    let chrom: Arc<str> = Arc::from("chr1");
    let reference: Arc<[u8]> = Arc::from(vec![b'A'; 200].into_boxed_slice());
    let indels = caller
        .call_indels_from_sorted_bams(chrom, reference, 0, &tumor_bam, &normal_bam)
        .unwrap();
    assert_eq!(indels.len(), 1);
    assert_eq!(indels[0].tumor_support, 1);
    assert_eq!(indels[0].normal_support, 0);

    let _ = std::fs::remove_file(tumor_bam);
    let _ = std::fs::remove_file(normal_bam);
}


