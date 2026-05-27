//! Phase B1 DoD: the public reader API reads gzipped, multi-record FASTA and
//! produces a multi-contig ContigSet.

use std::io::{Cursor, Write};

use flate2::write::GzEncoder;
use flate2::Compression;
use rosalind::core::ContigSet;
use rosalind::io::decompress::maybe_decompress;
use rosalind::io::fasta::FastaReader;

fn gzip(bytes: &[u8]) -> Vec<u8> {
    let mut enc = GzEncoder::new(Vec::new(), Compression::default());
    enc.write_all(bytes).unwrap();
    enc.finish().unwrap()
}

#[test]
fn gzipped_multi_record_fasta_builds_a_contig_set() {
    let fasta = b">chr1 first\nACGTACGT\n>chr2\nAACC\n>chr3\nGGGGTTTT\n";
    let gz = gzip(fasta);

    let reader = maybe_decompress(Cursor::new(gz)).unwrap();
    let records: Vec<_> = FastaReader::new(reader)
        .collect::<Result<Vec<_>, _>>()
        .expect("gz FASTA should parse");

    assert_eq!(records.len(), 3);

    let mut contigs = ContigSet::new();
    for r in &records {
        contigs.push(r.name.clone(), r.sequence.len() as u32);
    }
    assert_eq!(contigs.len(), 3);
    assert_eq!(contigs.by_name("chr1").unwrap().global_offset, 0);
    assert_eq!(contigs.by_name("chr2").unwrap().global_offset, 8);
    assert_eq!(contigs.by_name("chr3").unwrap().global_offset, 12);
    assert_eq!(contigs.total_length(), 20);
}
