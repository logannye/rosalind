//! Per-locus FEATURE egress: the bounded, deterministic `PileupColumn` stream
//! exposed as a tabular ML feature source. Same streaming engine, same memory
//! contract as germline calling — but every callable locus is emitted as a
//! feature row (not just variant sites). The output is byte-identical run-to-run
//! (canonical obs order + integer-mean formatting), so a BLAKE3 receipt over it
//! gives **bit-reproducible training inputs**.

use std::io::{self, Write};
use std::ops::Range;
use std::sync::Arc;

use crate::call::whole_genome::PerContig;
use crate::core::{AlignedRead, ContigSet, CoreError, WorkingSet};
use crate::genomics::ReferenceView;
use crate::pileup::{PileupColumn, PileupEngine, PileupParams, ReadSource, SkipCounts};

/// The tab-separated header for the per-locus feature table (written once).
pub const FEATURE_HEADER: &str = "#contig\tpos\tref\tdepth\traw_depth\ta\tc\tg\tt\t\
a_fwd\ta_rev\tc_fwd\tc_rev\tg_fwd\tg_rev\tt_fwd\tt_rev\tmean_bq\tmean_mapq";

/// Write the feature header line.
pub fn write_feature_header<W: Write>(out: &mut W) -> io::Result<()> {
    writeln!(out, "{FEATURE_HEADER}")
}

/// Write one feature row for a callable pileup column. `pos` is 1-based (VCF POS).
/// Means are integer sums over observations divided by callable depth, formatted
/// to 2 decimals so the row is byte-stable across runs.
pub fn write_feature_row<W: Write>(
    out: &mut W,
    contig_name: &str,
    col: &PileupColumn,
) -> io::Result<()> {
    let depth = col.depth();
    let ac = col.allele_counts();
    let sc = col.strand_counts();
    let (sum_bq, sum_mapq) = col.obs.iter().fold((0u64, 0u64), |(b, m), o| {
        (b + o.base_qual as u64, m + o.mapq as u64)
    });
    let (mean_bq, mean_mapq) = if depth > 0 {
        (sum_bq as f64 / depth as f64, sum_mapq as f64 / depth as f64)
    } else {
        (0.0, 0.0)
    };
    writeln!(
        out,
        "{contig}\t{pos}\t{refb}\t{depth}\t{raw}\t\
{a}\t{c}\t{g}\t{t}\t\
{afwd}\t{arev}\t{cfwd}\t{crev}\t{gfwd}\t{grev}\t{tfwd}\t{trev}\t\
{mean_bq:.2}\t{mean_mapq:.2}",
        contig = contig_name,
        pos = col.locus.pos.0 as u64 + 1,
        refb = col.ref_base as char,
        depth = depth,
        raw = col.raw_depth,
        a = ac[0],
        c = ac[1],
        g = ac[2],
        t = ac[3],
        afwd = sc[0][0],
        arev = sc[0][1],
        cfwd = sc[1][0],
        crev = sc[1][1],
        gfwd = sc[2][0],
        grev = sc[2][1],
        tfwd = sc[3][0],
        trev = sc[3][1],
    )
}

/// Stream feature columns over `region` of `contig` to a sink, returning the max
/// pileup working set + skip counts. The sink receives each callable column in
/// ascending position order; no genome-wide buffer accumulates.
pub fn stream_features_region<S: ReadSource>(
    source: S,
    reference: Arc<[u8]>,
    contig: u32,
    region: Range<u32>,
    pileup_params: PileupParams,
    on_row: &mut dyn FnMut(&PileupColumn) -> Result<(), CoreError>,
) -> Result<(WorkingSet, SkipCounts), CoreError> {
    let mut engine = PileupEngine::new(source, reference, contig, region, pileup_params);
    let mut max_ws = WorkingSet { bytes: 0 };
    while let Some(column) = engine.next() {
        let column = column?;
        let ws = engine.current_working_set();
        if ws.bytes > max_ws.bytes {
            max_ws = ws;
        }
        on_row(&column)?;
    }
    Ok((max_ws, engine.skip_counts()))
}

/// Stream feature columns across every contig (in id order), reading each contig's
/// reference from `ref_view`. The sink receives `(column, contig_name)`. Peak
/// memory ≈ the largest contig's reference + the pileup working set, independent
/// of input size. `source` MUST be `(contig, pos)`-ordered (a `StreamingBamSource`
/// guards this). Mirrors `call_germline_whole_genome`.
pub fn stream_features_whole_genome<S: ReadSource>(
    mut source: S,
    ref_view: &ReferenceView,
    contigs: &ContigSet,
    pileup_params: PileupParams,
    on_row: &mut dyn FnMut(&PileupColumn, &str) -> Result<(), CoreError>,
) -> Result<(WorkingSet, SkipCounts), CoreError> {
    let mut max_ws = WorkingSet { bytes: 0 };
    let mut skips = SkipCounts::default();
    let mut peeked: Option<AlignedRead> = None;

    for c in contigs.iter() {
        // Decode straight into the Arc (no intermediate Vec) — one resident copy,
        // not the transient two, so the predicted-peak bound holds (same as the
        // germline whole-genome driver).
        let start = c.global_offset as usize;
        let end = c.global_offset as usize + c.length as usize;
        let reference: Arc<[u8]> = ref_view.decode_window_arc(start, end);

        let per = PerContig {
            source: &mut source,
            contig: c.id,
            peeked: &mut peeked,
        };
        let region: Range<u32> = 0..c.length;
        let (ws, contig_skips) = stream_features_region(
            per,
            reference,
            c.id,
            region,
            pileup_params.clone(),
            &mut |col| on_row(col, &c.name),
        )?;
        if ws.bytes > max_ws.bytes {
            max_ws = ws;
        }
        skips.accumulate(&contig_skips);
    }

    Ok((max_ws, skips))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::{CigarOp, CigarOpKind, Position, SamFlags};
    use crate::pileup::SliceSource;

    fn mread(pos: u32, seq: &[u8], reverse: bool) -> AlignedRead {
        let flags = if reverse {
            SamFlags(SamFlags::REVERSE)
        } else {
            SamFlags::default()
        };
        AlignedRead {
            contig: 0,
            pos: Position(pos),
            mapq: 60,
            flags,
            cigar: vec![CigarOp::new(CigarOpKind::Match, seq.len() as u32)],
            seq: Arc::from(seq.to_vec().into_boxed_slice()),
            qual: Arc::from(vec![30u8; seq.len()].into_boxed_slice()),
        }
    }

    fn rows_for(reads: Vec<AlignedRead>, reference: &[u8]) -> String {
        let mut out: Vec<u8> = Vec::new();
        write_feature_header(&mut out).unwrap();
        stream_features_region(
            SliceSource::new(reads),
            Arc::from(reference.to_vec().into_boxed_slice()),
            0,
            0..reference.len() as u32,
            PileupParams::default(),
            &mut |col| write_feature_row(&mut out, "chr1", col).map_err(CoreError::from),
        )
        .unwrap();
        String::from_utf8(out).unwrap()
    }

    #[test]
    fn feature_row_has_the_expected_fields() {
        // Reference AAAA; at pos 0, two reads: A (ref, fwd) and C (alt, rev).
        let reference = b"AAAA";
        let reads = vec![mread(0, b"AAAA", false), mread(0, b"CAAA", true)];
        let text = rows_for(reads, reference);
        let line0 = text.lines().nth(1).expect("a data row"); // after header
                                                              // contig pos ref depth raw a c g t a_fwd a_rev c_fwd c_rev g_fwd g_rev t_fwd t_rev mean_bq mean_mapq
                                                              // pos 1 (1-based), ref A, depth 2, raw 2, a=1 c=1, a_fwd=1 c_rev=1, bq30 mapq60.
        let f: Vec<&str> = line0.split('\t').collect();
        assert_eq!(f[0], "chr1");
        assert_eq!(f[1], "1"); // 1-based
        assert_eq!(f[2], "A");
        assert_eq!(f[3], "2"); // depth
        assert_eq!(f[5], "1"); // a count
        assert_eq!(f[6], "1"); // c count
        assert_eq!(f[9], "1"); // a_fwd
        assert_eq!(f[12], "1"); // c_rev
        assert_eq!(f[17], "30.00"); // mean_bq
        assert_eq!(f[18], "60.00"); // mean_mapq
        assert_eq!(f.len(), 19);
    }

    #[test]
    fn feature_stream_is_deterministic_and_bounded() {
        use crate::core::MemoryBudget;
        // 50k single-base reads tiled at depth ~1: identical output twice, tiny WS.
        let reference = vec![b'A'; 50_000];
        let reads: Vec<AlignedRead> = (0..50_000u32).map(|p| mread(p, b"C", false)).collect();
        let a = rows_for(reads.clone(), &reference);
        let b = rows_for(reads, &reference);
        assert_eq!(a, b, "feature output must be byte-identical run-to-run");

        let reference2 = vec![b'A'; 50_000];
        let reads2: Vec<AlignedRead> = (0..50_000u32).map(|p| mread(p, b"C", false)).collect();
        let mut max_ws = 0u64;
        let (ws, _skips) = stream_features_region(
            SliceSource::new(reads2),
            Arc::from(reference2.into_boxed_slice()),
            0,
            0..50_000,
            PileupParams::default(),
            &mut |_col| Ok(()),
        )
        .unwrap();
        max_ws = max_ws.max(ws.bytes);
        // Bounded by coverage, not the 50k read count (reference ~50k + a tiny active set).
        assert!(
            ws.fits(MemoryBudget::from_mb(1)),
            "feature working set {max_ws} should fit 1 MiB"
        );
    }
}
