use std::env;
use std::path::PathBuf;
use std::time::{SystemTime, UNIX_EPOCH};

use rosalind::genomics::sort_bam_deterministic;
use rust_htslib::bam;
use rust_htslib::bam::record::{Cigar, CigarString, Record};
use rust_htslib::bam::Read as BamRead;

fn temp_path(name: &str) -> PathBuf {
    let timestamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("time went backwards")
        .as_nanos();
    env::temp_dir().join(format!("rosalind-{name}-{timestamp}.bam"))
}

#[test]
fn deterministic_bam_sort_orders_by_tid_pos_strand_qname() {
    let input = temp_path("unsorted");
    let output = temp_path("sorted");

    // Minimal header with one contig.
    let mut header = bam::Header::new();
    header.push_record(
        bam::header::HeaderRecord::new(b"SQ")
            .push_tag(b"SN", "chr1")
            .push_tag(b"LN", 1000),
    );

    // Write unsorted records: positions 10, 5, 10 (different qnames).
    {
        let mut writer = bam::Writer::from_path(&input, &header, bam::Format::Bam).unwrap();

        let mut r1 = Record::new();
        r1.set(
            b"readB",
            Some(&CigarString::from(vec![Cigar::Match(4)])),
            b"ACGT",
            b"IIII",
        );
        r1.set_tid(0);
        r1.set_pos(10);
        r1.set_flags(0);
        r1.set_mapq(60);
        writer.write(&r1).unwrap();

        let mut r2 = Record::new();
        r2.set(
            b"readA",
            Some(&CigarString::from(vec![Cigar::Match(4)])),
            b"ACGT",
            b"IIII",
        );
        r2.set_tid(0);
        r2.set_pos(5);
        r2.set_flags(0);
        r2.set_mapq(60);
        writer.write(&r2).unwrap();

        let mut unplaced = Record::new();
        unplaced.set(b"unplaced", None, b"ACGT", &[30; 4]);
        unplaced.set_tid(-1);
        unplaced.set_pos(-1);
        unplaced.set_flags(4);
        writer.write(&unplaced).unwrap();

        let mut r3 = Record::new();
        r3.set(
            b"readC",
            Some(&CigarString::from(vec![Cigar::Match(4)])),
            b"ACGT",
            b"IIII",
        );
        r3.set_tid(0);
        r3.set_pos(10);
        r3.set_flags(0x10); // reverse
        r3.set_mapq(60);
        writer.write(&r3).unwrap();
    }

    sort_bam_deterministic(&input, &output, 1 << 20).unwrap();

    // Read sorted output and verify order.
    let mut reader = bam::Reader::from_path(&output).unwrap();
    let header_text = String::from_utf8_lossy(reader.header().as_bytes());
    assert!(
        header_text
            .lines()
            .next()
            .unwrap()
            .contains("SO:coordinate"),
        "the sorted output must declare coordinate order: {header_text}"
    );
    let mut qnames = Vec::new();
    let mut keys = Vec::new();
    for rec in reader.records() {
        let rec = rec.unwrap();
        qnames.push(String::from_utf8_lossy(rec.qname()).to_string());
        keys.push((rec.tid(), rec.pos(), rec.is_reverse()));
    }

    // Expected: pos=5 readA first, then pos=10 forward readB, then pos=10 reverse readC.
    assert_eq!(qnames, vec!["readA", "readB", "readC", "unplaced"]);
    assert_eq!(keys[0].1, 5);
    assert_eq!(keys[1].1, 10);
    assert!(!keys[1].2);
    assert_eq!(keys[2].1, 10);
    assert!(keys[2].2);
    assert_eq!(keys[3].0, -1);
    // An incorrect unmapped-first order declares SO:coordinate but cannot be
    // indexed, breaking the packaged indexed-evidence quickstart.
    bam::index::build(&output, None, bam::index::Type::Bai, 1).unwrap();

    let _ = std::fs::remove_file(input);
    let _ = std::fs::remove_file(format!("{}.bai", output.display()));
    let _ = std::fs::remove_file(output);
}
