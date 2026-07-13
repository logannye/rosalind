//! ColumnKit: implement one trait, inherit the bounded memory contract.
//!
//! A builder who wants their OWN per-locus metric — methylation, a coverage/QC
//! track, custom ML features, a star-allele genotyper — would otherwise fork into
//! the raw pileup engine and re-implement, by hand, every surface that makes
//! Rosalind worth choosing: the whole-genome contig walk, the working-set bound
//! that `plan`/`--enforce` admit, and the canonical-JSON + BLAKE3 receipt.
//!
//! [`ColumnAnalyzer`] + [`run_bounded_whole_genome`] turn "extend the substrate"
//! into a first-class SDK: implement one trait, run it through the driver, and
//! inherit the SAME bounded per-contig stream, the SAME working-set bound, and
//! the SAME verifiable receipt the shipped `variants`/`features` subcommands
//! enjoy — for free. The trait is *welded* to the bounded kernel: the driver runs
//! the exact same column stream as `stream_features_whole_genome`, so the
//! estimator that admits a run upper-bounds the realized working set of
//! the builder's analyzer too. It is the contract made composable, not a feature
//! bolted beside it. (`FeatureAnalyzer` is the first impl — proof the trait
//! carries the real shipped analyzer, not a toy.)

use std::collections::BTreeMap;
use std::io::{self, Write};

use crate::call::features::{
    stream_features_region, stream_features_whole_genome, write_feature_row, FEATURE_HEADER,
};
use crate::core::{ContigSet, CoreError, WorkingSet};
use crate::genomics::{ReferenceProvider, ReferenceSequence};
use crate::io::bam::IndexedBamRegionSource;
use crate::pileup::{PileupColumn, PileupParams, ReadSource, SkipCounts};
use crate::selection::AnalysisSelection;

/// A per-locus analyzer over the bounded pileup-column stream. Implement this and
/// run it with [`run_bounded_whole_genome`] to inherit bounded memory,
/// byte-identical determinism, and a verifiable receipt.
pub trait ColumnAnalyzer {
    /// Optional header written once, before any column (e.g. a TSV header line,
    /// including its trailing newline). Default: none.
    fn header(&self) -> Option<String> {
        None
    }

    /// Parameters to fold into the run receipt (provenance). Default: none.
    fn params(&self) -> BTreeMap<String, String> {
        BTreeMap::new()
    }

    /// Process one callable pileup column, writing any output to `out`. Called in
    /// `(contig, position)` order; `contig` is the column's contig name. Keep
    /// per-call state O(1) (or bounded) to preserve the memory contract.
    fn on_column(
        &mut self,
        col: &PileupColumn,
        contig: &str,
        out: &mut dyn Write,
    ) -> io::Result<()>;

    /// Flush any bounded encoder state after the final column.
    fn finish(&mut self, _out: &mut dyn Write) -> io::Result<()> {
        Ok(())
    }
}

/// Drive `analyzer` over every callable column of a coordinate-sorted read stream
/// and a persisted reference, under the bounded memory contract: peak ≈ the
/// largest contig's reference + the pileup working set, independent of input
/// size. Writes the analyzer's header (if any), then streams columns to it in
/// `(contig, pos)` order. Returns the max pileup working set observed — the value
/// `plan`/`--enforce` reason about — and the skip counts, so the analyzer
/// inherits the receipt.
pub fn run_bounded_whole_genome<S: ReadSource>(
    analyzer: &mut dyn ColumnAnalyzer,
    source: S,
    ref_view: &dyn ReferenceSequence,
    contigs: &ContigSet,
    pileup_params: PileupParams,
    out: &mut dyn Write,
) -> Result<(WorkingSet, SkipCounts), CoreError> {
    if let Some(header) = analyzer.header() {
        out.write_all(header.as_bytes())?;
    }
    let result = stream_features_whole_genome(
        source,
        ref_view,
        contigs,
        pileup_params,
        &mut |col, contig| {
            analyzer
                .on_column(col, contig, &mut *out)
                .map_err(CoreError::from)
        },
    )?;
    analyzer.finish(out)?;
    Ok(result)
}

/// Drive an analyzer over one indexed interval selection. The header is written
/// once and every interval is fetched independently through the BAM's BAI.
pub fn run_bounded_selected_bam(
    analyzer: &mut dyn ColumnAnalyzer,
    bam_path: &std::path::Path,
    reference: &dyn ReferenceProvider,
    selection: &AnalysisSelection,
    pileup_params: PileupParams,
    out: &mut dyn Write,
) -> Result<(WorkingSet, SkipCounts), CoreError> {
    let intervals = selection.intervals().ok_or_else(|| {
        CoreError::MalformedRecord("selected BAM driver requires interval selection".into())
    })?;
    if let Some(header) = analyzer.header() {
        out.write_all(header.as_bytes())?;
    }
    let contigs = reference.contigs();
    let mut max_working_set = WorkingSet { bytes: 0 };
    let mut skips = SkipCounts::default();
    for interval in intervals.intervals() {
        crate::core::governor::checkpoint()?;
        let contig = contigs.by_id(interval.contig).ok_or_else(|| {
            CoreError::MalformedRecord(format!("unknown interval contig {}", interval.contig))
        })?;
        let global_start = contig.global_offset as usize + interval.start as usize;
        let global_end = contig.global_offset as usize + interval.end as usize;
        let sequence = reference.decode_window_arc(global_start, global_end);
        let source = IndexedBamRegionSource::new(bam_path, contigs, *interval)?;
        let (working_set, interval_skips) = stream_features_region(
            source,
            sequence,
            interval.contig,
            interval.start..interval.end,
            pileup_params.clone(),
            &mut |column| {
                analyzer
                    .on_column(column, &contig.name, &mut *out)
                    .map_err(CoreError::from)
            },
        )?;
        max_working_set.bytes = max_working_set.bytes.max(working_set.bytes);
        skips.accumulate(&interval_skips);
    }
    analyzer.finish(out)?;
    Ok((max_working_set, skips))
}

/// The shipped `features` egress, expressed as a [`ColumnAnalyzer`] — proof the
/// trait carries the real production analyzer. Counts the rows it emits (folded
/// into the receipt via [`ColumnAnalyzer::params`]).
#[derive(Debug, Default)]
pub struct FeatureAnalyzer {
    rows: u64,
}

impl FeatureAnalyzer {
    /// Number of feature rows written so far.
    pub fn rows(&self) -> u64 {
        self.rows
    }
}

impl ColumnAnalyzer for FeatureAnalyzer {
    fn header(&self) -> Option<String> {
        Some(format!("{FEATURE_HEADER}\n"))
    }

    fn params(&self) -> BTreeMap<String, String> {
        let mut p = BTreeMap::new();
        p.insert("feature_rows".to_string(), self.rows.to_string());
        p.insert("feature.schema".to_string(), "1".to_string());
        p.insert("artifact.format".to_string(), "tsv".to_string());
        p
    }

    fn on_column(
        &mut self,
        col: &PileupColumn,
        contig: &str,
        out: &mut dyn Write,
    ) -> io::Result<()> {
        write_feature_row(out, contig, col)?;
        self.rows += 1;
        Ok(())
    }
}

/// A per-locus coverage track: `(contig, 1-based pos, depth)`. The second first-party
/// [`ColumnAnalyzer`] — proof the SDK carries more than the feature egress. Its
/// `params()` names the analyzer so the receipt records which metric produced the run.
#[derive(Debug, Default)]
pub struct CoverageTrack;

impl ColumnAnalyzer for CoverageTrack {
    fn header(&self) -> Option<String> {
        Some("#contig\tpos\tdepth\n".to_string())
    }

    fn params(&self) -> BTreeMap<String, String> {
        BTreeMap::from([
            ("analyzer".to_string(), "coverage".to_string()),
            ("artifact.format".to_string(), "coverage-tsv".to_string()),
        ])
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::call::features::write_feature_header;

    #[test]
    fn coverage_track_header_and_params() {
        let t = CoverageTrack;
        assert_eq!(t.header().as_deref(), Some("#contig\tpos\tdepth\n"));
        assert_eq!(
            t.params().get("analyzer").map(String::as_str),
            Some("coverage")
        );
    }

    use crate::genomics::{GenomeIndex, IndexReader, IndexWriter};
    use crate::pileup::SliceSource;
    use std::sync::Arc;

    fn tmp(suffix: &str) -> std::path::PathBuf {
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let d = std::env::temp_dir().join(format!("rosalind-ck-{suffix}-{nanos}"));
        std::fs::create_dir_all(&d).unwrap();
        d.join("ref.idx")
    }

    fn read_at(contig: u32, pos: u32, seq: &[u8]) -> crate::core::AlignedRead {
        use crate::core::{CigarOp, CigarOpKind, Position, SamFlags};
        crate::core::AlignedRead {
            contig,
            pos: Position(pos),
            mapq: 60,
            flags: SamFlags(0),
            cigar: vec![CigarOp::new(CigarOpKind::Match, seq.len() as u32)],
            seq: Arc::from(seq.to_vec().into_boxed_slice()),
            qual: Arc::from(vec![40u8; seq.len()].into_boxed_slice()),
        }
    }

    // Running FeatureAnalyzer through the driver must be byte-identical to the
    // hand-wired features path (header + write_feature_row per column) — i.e. the
    // trait IS the production path, not a parallel one.
    #[test]
    fn feature_analyzer_via_driver_equals_the_direct_features_path() {
        let idx = tmp("equiv");
        let index = GenomeIndex::from_named_sequences(&[
            ("chr1".to_string(), b"ACGTACGTACGTACGTACGT".to_vec()),
            ("chr2".to_string(), b"TTTTGGGGCCCCAAAATTTT".to_vec()),
        ])
        .unwrap();
        IndexWriter::create(&idx)
            .unwrap()
            .write_genome_index(&index)
            .unwrap();
        let loaded = IndexReader::open(&idx).unwrap();
        let rv = loaded.reference_view().unwrap();
        let contigs = loaded.contigs();
        let reads = vec![
            read_at(0, 0, b"ACGTACGT"),
            read_at(0, 0, b"ACGTACGT"),
            read_at(1, 0, b"TTTTGGGG"),
        ];
        let pp = PileupParams::default();

        // Direct path: header + stream_features_whole_genome + write_feature_row.
        let mut direct: Vec<u8> = Vec::new();
        write_feature_header(&mut direct).unwrap();
        let (ws_direct, _) = stream_features_whole_genome(
            SliceSource::new(reads.clone()),
            &rv,
            contigs,
            pp.clone(),
            &mut |col, name| write_feature_row(&mut direct, name, col).map_err(CoreError::from),
        )
        .unwrap();

        // Driver path: FeatureAnalyzer through run_bounded_whole_genome.
        let mut analyzer = FeatureAnalyzer::default();
        let mut via_kit: Vec<u8> = Vec::new();
        let (ws_kit, _) = run_bounded_whole_genome(
            &mut analyzer,
            SliceSource::new(reads),
            &rv,
            contigs,
            pp,
            &mut via_kit,
        )
        .unwrap();

        assert_eq!(direct, via_kit, "driver output must be byte-identical");
        assert_eq!(
            ws_direct.bytes, ws_kit.bytes,
            "driver must report the same bounded working set"
        );
        assert!(analyzer.rows() > 0, "analyzer should count its rows");
        assert_eq!(
            analyzer.params().get("feature_rows").map(String::as_str),
            Some(analyzer.rows().to_string().as_str())
        );

        let _ = std::fs::remove_dir_all(idx.parent().unwrap());
    }

    // A tiny custom analyzer inherits the bounded walk: its working set is the
    // same coverage-bounded value, independent of how many reads stream in.
    #[test]
    fn a_custom_analyzer_inherits_the_bounded_working_set() {
        struct DepthTrack {
            max_depth_seen: u32,
        }
        impl ColumnAnalyzer for DepthTrack {
            fn on_column(
                &mut self,
                col: &PileupColumn,
                _contig: &str,
                _out: &mut dyn Write,
            ) -> io::Result<()> {
                self.max_depth_seen = self.max_depth_seen.max(col.depth());
                Ok(())
            }
        }

        let idx = tmp("custom");
        let index =
            GenomeIndex::from_named_sequences(&[("chr1".to_string(), vec![b'A'; 500])]).unwrap();
        IndexWriter::create(&idx)
            .unwrap()
            .write_genome_index(&index)
            .unwrap();
        let loaded = IndexReader::open(&idx).unwrap();
        let rv = loaded.reference_view().unwrap();
        let contigs = loaded.contigs();
        let pp = PileupParams {
            max_depth: Some(8),
            ..PileupParams::default()
        };

        let run = |n: usize| -> u64 {
            let reads: Vec<_> = (0..n).map(|_| read_at(0, 0, &[b'C'; 50])).collect();
            let mut a = DepthTrack { max_depth_seen: 0 };
            let mut sink = io::sink();
            run_bounded_whole_genome(
                &mut a,
                SliceSource::new(reads),
                &rv,
                contigs,
                pp.clone(),
                &mut sink,
            )
            .unwrap()
            .0
            .bytes
        };
        // Bounded: 100× more reads at one locus (capped at depth 8) → same WS.
        assert_eq!(
            run(20),
            run(2000),
            "analyzer working set must not grow with reads"
        );

        let _ = std::fs::remove_dir_all(idx.parent().unwrap());
    }
}
