use std::env;
use std::path::PathBuf;
use std::sync::Arc;
use std::time::{SystemTime, UNIX_EPOCH};

use rosalind::genomics::{sort_bam_deterministic, BamPileupStream};
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
fn pileup_stream_counts_depth() {
    let input = temp_path("pileup_unsorted");
    let sorted = temp_path("pileup_sorted");

    // Header with one contig.
    let mut header = bam::Header::new();
    header.push_record(
        bam::header::HeaderRecord::new(b"SQ")
            .push_tag(b"SN", &"chr1")
            .push_tag(b"LN", &1000),
    );

    // Two reads overlapping position 5.
    {
        let mut writer = bam::Writer::from_path(&input, &header, bam::Format::Bam).unwrap();

        let mut r1 = Record::new();
        r1.set(b"read1", Some(&CigarString::from(vec![Cigar::Match(4)])), b"ACGT", b"IIII");
        r1.set_tid(0);
        r1.set_pos(4); // covers 4..8
        r1.set_flags(0);
        r1.set_mapq(60);
        writer.write(&r1).unwrap();

        let mut r2 = Record::new();
        r2.set(b"read2", Some(&CigarString::from(vec![Cigar::Match(4)])), b"TTTT", b"IIII");
        r2.set_tid(0);
        r2.set_pos(5); // covers 5..9
        r2.set_flags(0);
        r2.set_mapq(60);
        writer.write(&r2).unwrap();
    }

    sort_bam_deterministic(&input, &sorted, 1 << 20).unwrap();

    let chrom: Arc<str> = Arc::from("chr1");
    let region = 4..7;
    let mut stream = BamPileupStream::new(&sorted, Arc::clone(&chrom), region).unwrap();
    let mut depths = Vec::new();
    while let Some(node) = stream.next() {
        let node = node.unwrap();
        depths.push((node.position, node.depth));
    }

    // Position 4: depth 1 (read1)
    // Position 5: depth 2 (read1 + read2)
    // Position 6: depth 2 (read1 + read2)
    assert!(depths.contains(&(4, 1)));
    assert!(depths.contains(&(5, 2)));
    assert!(depths.contains(&(6, 2)));

    let _ = std::fs::remove_file(input);
    let _ = std::fs::remove_file(sorted);
}


