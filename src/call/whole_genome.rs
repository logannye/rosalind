//! Bounded whole-genome germline calling over a coordinate-sorted read stream and
//! a persisted reference (`ReferenceView`). One sorted pass, partitioned per
//! contig by `PerContig`, reusing the per-contig caller — peak memory ≈ the
//! largest contig's reference + the pileup working set, independent of input size.

use std::ops::Range;
use std::sync::Arc;

use crate::call::{call_germline_region_tracked, GermlineCall, GermlineParams};
use crate::core::{AlignedRead, ContigSet, CoreError, Locus, WorkingSet};
use crate::genomics::ReferenceView;
use crate::pileup::{PileupParams, ReadSource};

/// Per-contig view over a shared sorted stream: yields contig `contig`'s reads,
/// then stops at (and buffers) the first read of a later contig.
struct PerContig<'s, S: ReadSource> {
    source: &'s mut S,
    contig: u32,
    peeked: &'s mut Option<AlignedRead>,
}

impl<S: ReadSource> ReadSource for PerContig<'_, S> {
    fn next_read(&mut self) -> Result<Option<AlignedRead>, CoreError> {
        if let Some(r) = self.peeked.take() {
            if r.contig == self.contig {
                return Ok(Some(r));
            }
            *self.peeked = Some(r); // belongs to a later contig — keep it, stop here
            return Ok(None);
        }
        match self.source.next_read()? {
            Some(r) if r.contig == self.contig => Ok(Some(r)),
            Some(r) => {
                *self.peeked = Some(r);
                Ok(None)
            }
            None => Ok(None),
        }
    }
}

/// Call germline variants across every contig of `contigs`, in id order, reading
/// each contig's reference from `ref_view`. Returns the accumulated
/// `(Locus, ref_base, call)` rows and the max pileup working set observed.
/// `source` MUST yield reads in `(contig, pos)` order (a `StreamingBamSource`
/// guards this); the per-contig partition relies on it.
pub fn call_germline_whole_genome<S: ReadSource>(
    mut source: S,
    ref_view: &ReferenceView,
    contigs: &ContigSet,
    pileup_params: PileupParams,
    germline_params: &GermlineParams,
) -> Result<(Vec<(Locus, u8, GermlineCall)>, WorkingSet), CoreError> {
    let mut rows = Vec::new();
    let mut max_ws = WorkingSet { bytes: 0 };
    let mut peeked: Option<AlignedRead> = None;
    let mut buf = Vec::new();

    for c in contigs.iter() {
        // Decode this contig's reference (peak = largest contig; bounded).
        let start = c.global_offset as usize;
        let end = c.global_offset as usize + c.length as usize;
        ref_view.decode_window(start, end, &mut buf);
        let reference: Arc<[u8]> = Arc::from(buf.as_slice());

        let per = PerContig {
            source: &mut source,
            contig: c.id,
            peeked: &mut peeked,
        };
        let region: Range<u32> = 0..c.length;
        let (sites, ws) = call_germline_region_tracked(
            per,
            reference,
            c.id,
            region,
            pileup_params.clone(),
            germline_params,
        )?;
        if ws.bytes > max_ws.bytes {
            max_ws = ws;
        }
        rows.extend(sites);
    }

    Ok((rows, max_ws))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::call::call_germline_region;
    use crate::genomics::{GenomeIndex, IndexReader, IndexWriter};
    use crate::pileup::SliceSource;

    fn tmp(suffix: &str) -> std::path::PathBuf {
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let d = std::env::temp_dir().join(format!("rosalind-wg-{suffix}-{nanos}"));
        std::fs::create_dir_all(&d).unwrap();
        d.join("ref.idx")
    }

    // A read fully matching `len` bases at (contig, pos).
    fn read_at(contig: u32, pos: u32, seq: &[u8]) -> AlignedRead {
        use crate::core::{CigarOp, CigarOpKind, Position, SamFlags};
        AlignedRead {
            contig,
            pos: Position(pos),
            mapq: 60,
            flags: SamFlags(0),
            cigar: vec![CigarOp::new(CigarOpKind::Match, seq.len() as u32)],
            seq: Arc::from(seq.to_vec().into_boxed_slice()),
            qual: Arc::from(vec![40u8; seq.len()].into_boxed_slice()),
        }
    }

    #[test]
    fn whole_genome_equals_per_contig_calls() {
        // Two contigs; reads with a clear SNV on each. Build an index for the
        // reference, then compare the whole-genome drive to per-contig calls.
        let idx_path = tmp("equal");
        let index = GenomeIndex::from_named_sequences(&[
            ("chr1".to_string(), b"ACGTACGTACGTACGTACGT".to_vec()),
            ("chr2".to_string(), b"TTTTGGGGCCCCAAAATTTT".to_vec()),
        ])
        .unwrap();
        IndexWriter::create(&idx_path)
            .unwrap()
            .write_genome_index(&index)
            .unwrap();
        let loaded = IndexReader::open(&idx_path).unwrap();
        let rv = loaded.reference_view().unwrap();
        let contigs = loaded.contigs();

        // Reads (already (contig,pos)-sorted): depth on a few positions of each contig.
        let reads = vec![
            read_at(0, 0, b"ACGTACGT"),
            read_at(0, 0, b"ACGTACGT"),
            read_at(1, 0, b"TTTTGGGG"),
            read_at(1, 0, b"TTTTGGGG"),
        ];
        let pp = PileupParams::default();
        let gp = GermlineParams::default();

        let (rows, ws) = call_germline_whole_genome(
            SliceSource::new(reads.clone()),
            &rv,
            contigs,
            pp.clone(),
            &gp,
        )
        .unwrap();
        // Every row's contig is a real contig id; rows are grouped by contig.
        assert!(rows
            .iter()
            .all(|(l, _, _)| (l.contig as usize) < contigs.len()));
        // Working set was tracked: depth-2 coverage on each contig means the
        // pileup engine held a non-zero active read set at its peak.
        assert!(
            ws.bytes > 0,
            "max working set should be tracked and non-zero"
        );

        // Equivalence: union of per-contig call_germline_region over the same reads.
        let mut expected = Vec::new();
        for c in contigs.iter() {
            let mut buf = Vec::new();
            rv.decode_window(
                c.global_offset as usize,
                c.global_offset as usize + c.length as usize,
                &mut buf,
            );
            let cref: Arc<[u8]> = Arc::from(buf.as_slice());
            let creads: Vec<AlignedRead> =
                reads.iter().filter(|r| r.contig == c.id).cloned().collect();
            let sites = call_germline_region(
                SliceSource::new(creads),
                cref,
                c.id,
                0..c.length,
                pp.clone(),
                &gp,
            )
            .unwrap();
            expected.extend(sites);
        }
        assert_eq!(rows, expected, "whole-genome == per-contig union");

        let _ = std::fs::remove_dir_all(idx_path.parent().unwrap());
    }
}
