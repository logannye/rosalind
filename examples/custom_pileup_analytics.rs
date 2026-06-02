//! Cookbook: build your own bounded, deterministic per-locus analytics over the
//! `PileupColumn` substrate — no variant calling. Run with:
//!   cargo run --example custom_pileup_analytics

use std::sync::Arc;

use rosalind::core::{AlignedRead, CigarOp, CigarOpKind, Position, SamFlags};
use rosalind::{PileupColumn, PileupEngine, PileupParams, SliceSource};

fn read(pos: u32, seq: &[u8]) -> AlignedRead {
    AlignedRead {
        contig: 0,
        pos: Position(pos),
        mapq: 60,
        flags: SamFlags(0),
        cigar: vec![CigarOp::new(CigarOpKind::Match, seq.len() as u32)],
        seq: Arc::from(seq.to_vec().into_boxed_slice()),
        qual: Arc::from(vec![40u8; seq.len()].into_boxed_slice()),
    }
}

fn main() {
    let reference: Arc<[u8]> = Arc::from(b"ACGTACGTACGT".to_vec().into_boxed_slice());
    let reads = vec![read(0, b"ACGT"), read(2, b"GTAC"), read(4, b"ACGT")];

    // The substrate: one PileupColumn per covered position, bounded by coverage —
    // not by input size. PileupEngine<S: ReadSource> is an Iterator.
    let engine = PileupEngine::new(
        SliceSource::new(reads),
        Arc::clone(&reference),
        0,
        0..reference.len() as u32,
        PileupParams::default(),
    );

    // A custom per-locus metric — here, depth — computed without any calling.
    println!("pos\tref\tdepth");
    let mut total_depth = 0u64;
    for column in engine {
        let col: PileupColumn = column.expect("pileup column");
        total_depth += col.depth() as u64;
        println!(
            "{}\t{}\t{}",
            col.locus.pos.0,
            col.ref_base as char,
            col.depth()
        );
    }
    println!("# total observed depth across covered positions: {total_depth}");
}
