use std::sync::Arc;

use rosalind::call::{call_germline_region, Genotype, GermlineParams};
use rosalind::core::{AlignedRead, CigarOp, CigarOpKind, Position, SamFlags};
use rosalind::pileup::{PileupParams, SliceSource};

fn read(pos: u32, seq: &[u8]) -> AlignedRead {
    AlignedRead {
        contig: 0,
        pos: Position(pos),
        mapq: 60,
        flags: SamFlags::default(),
        cigar: vec![CigarOp::new(CigarOpKind::Match, seq.len() as u32)],
        seq: Arc::from(seq.to_vec().into_boxed_slice()),
        qual: Arc::from(vec![35u8; seq.len()].into_boxed_slice()),
    }
}

#[test]
fn germline_pipeline_detects_a_het_snv() {
    let reference: Arc<[u8]> = Arc::from(b"ACGTACGT".to_vec().into_boxed_slice());
    // Position 3 (ref T): half the reads carry A.
    let reads = vec![
        read(0, b"ACGTACGT"),
        read(0, b"ACGAACGT"),
        read(0, b"ACGTACGT"),
        read(0, b"ACGAACGT"),
    ];
    let sites = call_germline_region(
        SliceSource::new(reads),
        reference,
        0,
        0..8,
        PileupParams::default(),
        &GermlineParams::default(),
    )
    .unwrap();
    let site = sites
        .iter()
        .find(|(l, _, _)| l.pos == Position(3))
        .expect("variant at pos 3");
    assert_eq!(site.1, b'T'); // ref base
    assert_eq!(site.2.alt_base, b'A');
    assert_eq!(site.2.genotype, Genotype::Het);
}
