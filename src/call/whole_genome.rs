//! Bounded whole-genome germline calling over a coordinate-sorted read stream and
//! a persisted reference (`ReferenceView`). One sorted pass, partitioned per
//! contig by `PerContig`, reusing the per-contig caller — peak memory ≈ the
//! largest contig's reference + the pileup working set, independent of input size.

use std::ops::Range;
use std::sync::Arc;

use crate::call::{call_germline_region_streaming, GermlineCall, GermlineParams};
use crate::core::{governor, AlignedRead, ContigSet, CoreError, Locus, WorkingSet};
use crate::genomics::{ReferenceProvider, ReferenceSequence};
use crate::io::bam::IndexedBamRegionSource;
use crate::pileup::{PileupParams, ReadSource, SkipCounts};
use crate::selection::AnalysisSelection;

/// Per-contig view over a shared sorted stream: yields contig `contig`'s reads,
/// then stops at (and buffers) the first read of a later contig. Shared by the
/// germline-calling and feature-egress whole-genome drivers.
pub(crate) struct PerContig<'s, S: ReadSource> {
    pub(crate) source: &'s mut S,
    pub(crate) contig: u32,
    pub(crate) peeked: &'s mut Option<AlignedRead>,
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
/// each contig's reference from `ref_view` and streaming each emitted
/// `(Locus, ref_base, call)` to `on_row`. Returns the max pileup working set
/// observed; no genome-wide row buffer accumulates. `source` MUST yield reads in
/// `(contig, pos)` order (a `StreamingBamSource` guards this); the per-contig
/// partition relies on it.
pub fn call_germline_whole_genome<S: ReadSource>(
    mut source: S,
    ref_view: &dyn ReferenceSequence,
    contigs: &ContigSet,
    pileup_params: PileupParams,
    germline_params: &GermlineParams,
    on_row: &mut dyn FnMut((Locus, u8, GermlineCall)) -> Result<(), CoreError>,
) -> Result<(WorkingSet, SkipCounts), CoreError> {
    let mut max_ws = WorkingSet { bytes: 0 };
    let mut skips = SkipCounts::default();
    let mut peeked: Option<AlignedRead> = None;

    for c in contigs.iter() {
        // Cooperative budget check (no-op unless an --enforce governor is armed):
        // catch a breach before decoding the next contig's reference.
        governor::checkpoint()?;
        // Decode this contig's reference straight into the Arc — no intermediate
        // Vec, so peak resident reference memory is ONE copy, not the transient
        // two that `Arc::from(decoded)` (a reallocation) would briefly hold. Peak
        // = the largest contig; bounded, and the predicted-peak bound holds.
        let start = c.global_offset as usize;
        let end = c.global_offset as usize + c.length as usize;
        let reference: Arc<[u8]> = ref_view.decode_window_arc(start, end);

        let per = PerContig {
            source: &mut source,
            contig: c.id,
            peeked: &mut peeked,
        };
        let region: Range<u32> = 0..c.length;
        let (ws, contig_skips) = call_germline_region_streaming(
            per,
            reference,
            c.id,
            region,
            pileup_params.clone(),
            germline_params,
            on_row,
        )?;
        if ws.bytes > max_ws.bytes {
            max_ws = ws;
        }
        skips.accumulate(&contig_skips);
    }

    Ok((max_ws, skips))
}

/// Call germline variants over an indexed interval or shard selection.
pub fn call_germline_selected_bam(
    bam_path: &std::path::Path,
    reference: &dyn ReferenceProvider,
    selection: &AnalysisSelection,
    pileup_params: PileupParams,
    germline_params: &GermlineParams,
    on_row: &mut dyn FnMut((Locus, u8, GermlineCall)) -> Result<(), CoreError>,
) -> Result<(WorkingSet, SkipCounts), CoreError> {
    let intervals = selection.intervals().ok_or_else(|| {
        CoreError::MalformedRecord("selected germline driver requires intervals".into())
    })?;
    let mut max_working_set = WorkingSet { bytes: 0 };
    let mut skips = SkipCounts::default();
    for interval in intervals.intervals() {
        crate::core::governor::checkpoint()?;
        let contig = reference.contigs().by_id(interval.contig).ok_or_else(|| {
            CoreError::MalformedRecord(format!("unknown interval contig {}", interval.contig))
        })?;
        let start = contig.global_offset as usize + interval.start as usize;
        let end = contig.global_offset as usize + interval.end as usize;
        let sequence = reference.decode_window_arc(start, end);
        let source = IndexedBamRegionSource::new(bam_path, reference.contigs(), *interval)?;
        let (working_set, interval_skips) = call_germline_region_streaming(
            source,
            sequence,
            interval.contig,
            interval.start..interval.end,
            pileup_params.clone(),
            germline_params,
            on_row,
        )?;
        max_working_set.bytes = max_working_set.bytes.max(working_set.bytes);
        skips.accumulate(&interval_skips);
    }
    Ok((max_working_set, skips))
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

        let mut rows: Vec<(Locus, u8, GermlineCall)> = Vec::new();
        let (ws, _skips) = call_germline_whole_genome(
            SliceSource::new(reads.clone()),
            &rv,
            contigs,
            pp.clone(),
            &gp,
            &mut |row| {
                rows.push(row);
                Ok(())
            },
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

    #[test]
    fn working_set_is_bounded_by_reference_and_local_pressure_not_read_count() {
        // One small contig; extend the stream across nonoverlapping blocks at
        // fixed pressure of eight reads. The exact working set includes the
        // reference and retained observations, independent of total read count.
        let idx_path = tmp("bounded");
        let index =
            GenomeIndex::from_named_sequences(&[("chr1".to_string(), vec![b'A'; 2000])]).unwrap();
        IndexWriter::create(&idx_path)
            .unwrap()
            .write_genome_index(&index)
            .unwrap();
        let loaded = IndexReader::open(&idx_path).unwrap();
        let rv = loaded.reference_view().unwrap();
        let contigs = loaded.contigs();

        let params = PileupParams {
            max_depth: Some(8),
            ..PileupParams::default()
        };
        let gp = GermlineParams::default();

        let run = |n: usize| -> u64 {
            let reads: Vec<AlignedRead> = (0..n)
                .map(|read| read_at(0, ((read / 8) * 100) as u32, &[b'C'; 100]))
                .collect();
            call_germline_whole_genome(
                SliceSource::new(reads),
                &rv,
                contigs,
                params.clone(),
                &gp,
                &mut |_row| Ok(()),
            )
            .unwrap()
            .0
            .bytes
        };

        let ws_small = run(8);
        let ws_large = run(160);
        // (a) reference is counted: bound exceeds the 2000-byte reference.
        assert!(
            ws_small > 2000,
            "working set must include the reference bytes"
        );
        // (c) Twenty times more input, same exact local pressure and working set.
        assert_eq!(
            ws_small, ws_large,
            "working set must not grow with the number of input reads"
        );
        // (b) bounded by reference + capped active: ref(2000) + 8 reads *
        // (map 100*16 + seq 100 + qual 100 + 64) + 256, generously bounded.
        let bound = 2000 + 8 * (100 * 16 + 100 + 100 + 64) + 256;
        assert!(
            ws_small <= bound,
            "working set {ws_small} exceeded the analytic bound {bound}"
        );

        let _ = std::fs::remove_dir_all(idx_path.parent().unwrap());
    }
}
