//! Phase B3b gates: the persisted, memory-mapped FM-index is byte-identical to
//! the in-RAM index, deterministic, self-contained, loads without rebuilding,
//! and keeps the index out of resident memory.

use std::env;
use std::time::{SystemTime, UNIX_EPOCH};

use rosalind::core::Locus;
use rosalind::genomics::{GenomeIndex, IndexReader, IndexWriter};

fn temp_path(suffix: &str) -> std::path::PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    env::temp_dir().join(format!("rosalind-b3b-it-{suffix}-{nanos}.idx"))
}

fn multi_contig_index() -> GenomeIndex {
    GenomeIndex::from_named_sequences(&[
        (
            "chr1".to_string(),
            b"ACGTACGTNNACGTACGTACGTAAGGCCTT".to_vec(),
        ),
        (
            "chr2".to_string(),
            b"TTTTGGGGCCCCAAAANNNNACGTACGTAC".to_vec(),
        ),
        (
            "chr3".to_string(),
            b"GATTACATTTTGATTACAGGGGGCCCCAAA".to_vec(),
        ),
    ])
    .expect("build")
}

#[test]
fn view_equals_in_ram_over_a_pattern_battery() {
    let idx = multi_contig_index();
    let path = temp_path("equiv");
    IndexWriter::create(&path)
        .unwrap()
        .write_genome_index(&idx)
        .unwrap();
    let loaded = IndexReader::open(&path).unwrap();
    let gv = loaded.genome_view().unwrap();

    let patterns: &[&[u8]] = &[
        b"A",
        b"C",
        b"G",
        b"T",
        b"N",
        b"GATTACA",
        b"ACGT",
        b"NNNN",
        b"GGGGG",
        b"TTTTGGGG",
        b"AAGGCCTT",
        b"CGTACGTAC",
        b"GATTACAGGGGG",
        b"acgt",
        b"NnAc",
    ];
    for &p in patterns {
        let expected: Vec<Locus> = idx.locate_exact(p, 4096);
        let got: Vec<Locus> = gv.locate_exact(p, 4096);
        assert_eq!(got, expected, "locate_exact mismatch for {p:?}");
    }
    let _ = std::fs::remove_file(path);
}

#[test]
fn build_is_deterministic() {
    // Two independent builds of the same reference (SA-IS build + serialize) must
    // produce byte-identical files — the end-to-end determinism gate.
    let p1 = temp_path("det1");
    let p2 = temp_path("det2");
    IndexWriter::create(&p1)
        .unwrap()
        .write_genome_index(&multi_contig_index())
        .unwrap();
    IndexWriter::create(&p2)
        .unwrap()
        .write_genome_index(&multi_contig_index())
        .unwrap();
    assert_eq!(std::fs::read(&p1).unwrap(), std::fs::read(&p2).unwrap());
    let _ = std::fs::remove_file(p1);
    let _ = std::fs::remove_file(p2);
}

#[test]
fn index_is_self_contained() {
    // Build, drop the in-RAM index, then query using only the file.
    let path = temp_path("selfcontained");
    {
        let idx = multi_contig_index();
        IndexWriter::create(&path)
            .unwrap()
            .write_genome_index(&idx)
            .unwrap();
    } // `idx` dropped here — the file is the only source.

    let loaded = IndexReader::open(&path).unwrap();
    let gv = loaded.genome_view().unwrap();
    // "GATTACA" occurs in chr3 at local pos 0 and pos 11 (both within chr3).
    let loci = gv.locate_exact(b"GATTACA", 4096);
    assert_eq!(loci.len(), 2, "expected two GATTACA hits in chr3");
    assert!(
        loci.iter().all(|l| l.contig == 2),
        "both hits must be in chr3"
    );
    let _ = std::fs::remove_file(path);
}

#[test]
fn view_is_a_small_borrow_not_an_owned_copy() {
    // Bounded residency (B3b witness; the enforced RSS gate is Phase C): the view
    // is a handful of slices + scalars, not an owned copy of the index. Its size
    // is independent of genome size.
    use rosalind::genomics::FmIndexView;
    assert!(
        std::mem::size_of::<FmIndexView<'_>>() <= 256,
        "FmIndexView must be a small borrow, got {} bytes",
        std::mem::size_of::<FmIndexView<'_>>()
    );
}
