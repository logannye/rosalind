use std::collections::HashSet;
use std::sync::Arc;

use blake3::hash;
use rosalind::call::{call_germline_region, GermlineParams};
use rosalind::core::{AlignedRead, CigarOp, CigarOpKind, ContigSet, Position, SamFlags};
use rosalind::io::vcf::{render_germline_vcf, GermlineRow};
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
fn germline_calling_and_vcf_are_deterministic() {
    let reference: Arc<[u8]> = Arc::from(b"ACGTACGTACGTACGT".to_vec().into_boxed_slice());
    let reads = vec![
        read(0, b"ACGTACGT"),
        read(2, b"GTAATCGT"),
        read(0, b"ACAATCGT"),
    ];
    let mut contigs = ContigSet::new();
    contigs.push("chrDet", reference.len() as u32);

    let mut fingerprints = HashSet::new();
    for _ in 0..5 {
        let sites = call_germline_region(
            SliceSource::new(reads.clone()),
            Arc::clone(&reference),
            0,
            0..reference.len() as u32,
            PileupParams::default(),
            &GermlineParams::default(),
        )
        .unwrap();
        let rows: Vec<GermlineRow> = sites
            .into_iter()
            .map(|(locus, ref_base, call)| GermlineRow {
                locus,
                ref_base,
                call,
            })
            .collect();
        let vcf = render_germline_vcf(&contigs, "S", &rows).unwrap();
        fingerprints.insert(hash(vcf.as_bytes()));
    }
    assert_eq!(fingerprints.len(), 1, "germline VCF diverged across runs");
}
