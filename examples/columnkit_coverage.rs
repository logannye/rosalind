//! Legacy ColumnKit cookbook: a custom per-locus coverage consumer.
//! The low-level driver provides ordered capacity-limited columns and writes to a
//! caller-owned sink. This example does not create a receipt or enforce a whole-
//! process budget. Use the public contract runner/scaffold for that lifecycle.
//!
//! Run with: `cargo run --example columnkit_coverage`

use std::collections::BTreeMap;
use std::io::{self, Write};
use std::sync::Arc;

use rosalind::call::{run_bounded_whole_genome, ColumnAnalyzer};
use rosalind::core::{AlignedRead, CigarOp, CigarOpKind, Position, SamFlags};
use rosalind::genomics::{GenomeIndex, IndexReader, IndexWriter};
use rosalind::{PileupColumn, PileupParams, SliceSource};

/// Emit a per-locus coverage track (contig, 1-based pos, depth). This is the
/// per-column domain logic; orchestration and receipts are separate concerns.
struct CoverageTrack;

impl ColumnAnalyzer for CoverageTrack {
    fn header(&self) -> Option<String> {
        Some("#contig\tpos\tdepth\n".to_string())
    }

    fn params(&self) -> BTreeMap<String, String> {
        BTreeMap::from([("analyzer".to_string(), "coverage".to_string())])
    }

    fn on_column(
        &mut self,
        col: &PileupColumn,
        contig: &str,
        out: &mut dyn Write,
    ) -> io::Result<()> {
        writeln!(out, "{contig}\t{}\t{}", col.locus.pos.0 + 1, col.depth())
    }
}

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

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // A tiny in-memory genome + a few reads. A real run would pass a persisted
    // index and a coordinate-sorted BAM (`StreamingBamSource`) — the driver and
    // the contract are identical either way.
    let dir = std::env::temp_dir().join("rosalind-columnkit-example");
    std::fs::create_dir_all(&dir)?;
    let idx_path = dir.join("ref.idx");
    let index =
        GenomeIndex::from_named_sequences(&[("chr1".to_string(), b"ACGTACGTACGTACGT".to_vec())])?;
    IndexWriter::create(&idx_path)?.write_genome_index(&index)?;
    let loaded = IndexReader::open(&idx_path)?;
    let ref_view = loaded.reference_view()?;
    let contigs = loaded.contigs();

    let reads = vec![
        read(0, b"ACGTACGT"),
        read(0, b"ACGTACGT"),
        read(4, b"ACGTACGT"),
    ];

    let mut analyzer = CoverageTrack;
    let mut out = io::stdout().lock();
    // emit_all_positions makes the coverage track REFERENCE-COMPLETE: every locus
    // emits, so uncovered positions (here ref 12–15, which no read reaches) report
    // depth 0 instead of being silently dropped. The memory bound is unchanged.
    let (ws, _skips) = run_bounded_whole_genome(
        &mut analyzer,
        SliceSource::new(reads),
        &ref_view,
        contigs,
        PileupParams {
            emit_all_positions: true,
            ..PileupParams::default()
        },
        &mut out,
    )?;
    out.flush()?;

    eprintln!(
        "\n# reported pileup working set: {} bytes (not a whole-process RSS measurement)",
        ws.bytes
    );

    std::fs::remove_dir_all(&dir).ok();
    Ok(())
}
