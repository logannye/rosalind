//! Bounded gVCF: a banded reference-block + variant stream.
//!
//! A sites-only VCF lacks the reference confidence needed by cohort joiners —
//! they need to tell "reference here" from "not covered here". A gVCF says both:
//! every callable locus is either a variant record or part of a `<NON_REF>`
//! reference block (`END=` span) banded by genotype quality.
//!
//! gVCF is the ONE major caller output whose memory profile is intrinsically
//! streaming-friendly: reference blocks are locus-ordered and append-only, so the
//! whole-genome gVCF stays inside the same per-contig streaming envelope as
//! `variants` — the bander holds O(1) state (one open block). Combined with the
//! BLAKE3 receipt, each sample's gVCF is byte-reproducible — a property production
//! gVCF pipelines generally don't provide — and that reproducibility compounds at
//! cohort scale.

use std::io::{self, Write};
use std::ops::Range;
use std::sync::Arc;

use crate::call::germline::{genotype_column_gvcf, GvcfGenotype};
use crate::call::types::{Filter, Genotype, GermlineCall, GermlineParams};
use crate::call::whole_genome::PerContig;
use crate::core::{governor, AlignedRead, ContigSet, CoreError, Locus, WorkingSet};
use crate::genomics::{ReferenceProvider, ReferenceSequence};
use crate::io::bam::IndexedBamRegionSource;
use crate::pileup::{PileupEngine, PileupParams, ReadSource, SkipCounts};
use crate::selection::AnalysisSelection;

const GVCF_FORMAT: &str = "GT:DP:GQ:MIN_DP";

/// GQ band bucket (GATK-style): consecutive hom-ref columns coalesce into one
/// reference block while their GQ stays within the same 10-wide bucket.
fn gq_band(gq: u8) -> u8 {
    (gq / 10).min(9)
}

/// Accumulates consecutive hom-ref columns into banded `<NON_REF>` reference
/// blocks, writing variant records inline. Holds O(1) state (one open block), so
/// a whole-genome gVCF stays inside the bounded streaming envelope.
#[derive(Debug)]
pub struct GvcfBander<'c> {
    contigs: &'c ContigSet,
    open: Option<RefBlock>,
}

#[derive(Debug)]
struct RefBlock {
    contig: u32,
    ref_base: u8,
    start_pos0: u32,
    end_pos0: u32,
    band: u8,
    dp_start: u32,
    min_gq: u8,
    min_dp: u32,
}

impl<'c> GvcfBander<'c> {
    /// A bander resolving contig names through `contigs`.
    pub fn new(contigs: &'c ContigSet) -> Self {
        Self {
            contigs,
            open: None,
        }
    }

    fn chrom(&self, contig: u32) -> &str {
        self.contigs
            .by_id(contig)
            .map(|c| c.name.as_ref())
            .unwrap_or(".")
    }

    /// Feed one column's gVCF genotype. A variant flushes the open block then
    /// writes the variant record; a hom-ref column extends or opens a band; a
    /// no-call breaks the band (uncovered/non-callable positions are gaps).
    pub fn push<W: Write>(
        &mut self,
        locus: Locus,
        ref_base: u8,
        gt: GvcfGenotype,
        out: &mut W,
    ) -> io::Result<()> {
        match gt {
            GvcfGenotype::Variant(call) => {
                self.flush(out)?;
                self.write_variant(locus, ref_base, &call, out)?;
            }
            GvcfGenotype::HomRef { gq, dp } => {
                let band = gq_band(gq);
                let pos0 = locus.pos.0;
                let extend = matches!(&self.open, Some(b)
                    if b.contig == locus.contig && b.end_pos0 + 1 == pos0 && b.band == band);
                if extend {
                    let b = self.open.as_mut().expect("checked by matches!");
                    b.end_pos0 = pos0;
                    b.min_gq = b.min_gq.min(gq);
                    b.min_dp = b.min_dp.min(dp);
                } else {
                    self.flush(out)?;
                    self.open = Some(RefBlock {
                        contig: locus.contig,
                        ref_base,
                        start_pos0: pos0,
                        end_pos0: pos0,
                        band,
                        dp_start: dp,
                        min_gq: gq,
                        min_dp: dp,
                    });
                }
            }
            GvcfGenotype::NoCall { .. } => {
                self.flush(out)?;
            }
        }
        Ok(())
    }

    /// Write any open reference block. Call at contig boundaries and at the end of
    /// the stream so no block is left unflushed.
    pub fn flush<W: Write>(&mut self, out: &mut W) -> io::Result<()> {
        if let Some(b) = self.open.take() {
            let chrom = self.chrom(b.contig);
            writeln!(
                out,
                "{chrom}\t{}\t.\t{}\t<NON_REF>\t.\t.\tEND={}\t{GVCF_FORMAT}\t0/0:{}:{}:{}",
                b.start_pos0 as u64 + 1,
                b.ref_base as char,
                b.end_pos0 as u64 + 1,
                b.dp_start,
                b.min_gq,
                b.min_dp,
            )?;
        }
        Ok(())
    }

    fn write_variant<W: Write>(
        &self,
        locus: Locus,
        ref_base: u8,
        call: &GermlineCall,
        out: &mut W,
    ) -> io::Result<()> {
        let chrom = self.chrom(locus.contig);
        let gt = match call.genotype {
            Genotype::HomRef => "0/0",
            Genotype::Het => "0/1",
            Genotype::HomAlt => "1/1",
        };
        let filt = match call.filter {
            Filter::Pass => "PASS",
            Filter::LowQual => "LowQual",
            Filter::LowDepth => "LowDepth",
        };
        // `<NON_REF>` as a second ALT (the GATK gVCF convention) so the record
        // joins cleanly in cohort genotyping.
        writeln!(
            out,
            "{chrom}\t{}\t.\t{}\t{},<NON_REF>\t{:.1}\t{filt}\t.\t{GVCF_FORMAT}\t{gt}:{}:{}:{}",
            locus.pos.0 as u64 + 1,
            ref_base as char,
            call.alt_base as char,
            call.qual,
            call.dp,
            call.gq,
            call.dp,
        )
    }
}

/// `##fileformat` + contigs + the gVCF-specific `<NON_REF>` ALT / `END` / FORMAT
/// declarations + the `#CHROM` line.
pub fn write_gvcf_header<W: Write>(
    out: &mut W,
    contigs: &ContigSet,
    sample: &str,
) -> io::Result<()> {
    writeln!(out, "##fileformat=VCFv4.2")?;
    for c in contigs.iter() {
        writeln!(out, "##contig=<ID={},length={}>", c.name, c.length)?;
    }
    writeln!(
        out,
        r#"##ALT=<ID=NON_REF,Description="Represents any possible alternative allele">"#
    )?;
    writeln!(
        out,
        r#"##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the reference block">"#
    )?;
    writeln!(
        out,
        r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#
    )?;
    writeln!(
        out,
        r#"##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">"#
    )?;
    writeln!(
        out,
        r#"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">"#
    )?;
    writeln!(
        out,
        r#"##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP over the reference block">"#
    )?;
    writeln!(
        out,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}"
    )
}

/// Stream a banded gVCF for one contig region into `out`, feeding every callable
/// column to `bander`. Returns the max pileup working set observed.
#[allow(clippy::too_many_arguments)] // a streaming driver: source + region + params + sinks
fn stream_gvcf_region<S: ReadSource, W: Write>(
    bander: &mut GvcfBander,
    source: S,
    reference: Arc<[u8]>,
    contig: u32,
    region: Range<u32>,
    pileup_params: PileupParams,
    germline_params: &GermlineParams,
    out: &mut W,
) -> Result<(WorkingSet, SkipCounts), CoreError> {
    let mut engine = PileupEngine::new(source, reference, contig, region, pileup_params);
    let mut max_ws = WorkingSet { bytes: 0 };
    while let Some(column) = engine.next() {
        let column = column?;
        governor::checkpoint()?;
        let ws = engine.current_working_set();
        if ws.bytes > max_ws.bytes {
            max_ws = ws;
        }
        let gt = genotype_column_gvcf(&column, germline_params);
        bander
            .push(column.locus, column.ref_base, gt, out)
            .map_err(CoreError::from)?;
    }
    // Flush the contig's last reference block before the next contig.
    bander.flush(out).map_err(CoreError::from)?;
    Ok((max_ws, engine.skip_counts()))
}

/// Bounded whole-genome gVCF over a coordinate-sorted read stream and a persisted
/// reference. One sorted pass, partitioned per contig — peak ≈ the largest
/// contig's reference + the pileup working set, independent of input size (the
/// banding state is O(1)). The header must already be written to `out`.
pub fn stream_gvcf_whole_genome<S: ReadSource, W: Write>(
    mut source: S,
    ref_view: &dyn ReferenceSequence,
    contigs: &ContigSet,
    pileup_params: PileupParams,
    germline_params: &GermlineParams,
    out: &mut W,
) -> Result<(WorkingSet, SkipCounts), CoreError> {
    let mut bander = GvcfBander::new(contigs);
    let mut max_ws = WorkingSet { bytes: 0 };
    let mut skips = SkipCounts::default();
    let mut peeked: Option<AlignedRead> = None;

    for c in contigs.iter() {
        // Cooperative budget check (no-op unless an --enforce governor is armed):
        // catch a breach before decoding the next contig's reference.
        governor::checkpoint()?;
        let start = c.global_offset as usize;
        let end = c.global_offset as usize + c.length as usize;
        let reference: Arc<[u8]> = ref_view.decode_window_arc(start, end);
        let per = PerContig {
            source: &mut source,
            contig: c.id,
            peeked: &mut peeked,
        };
        let (ws, sk) = stream_gvcf_region(
            &mut bander,
            per,
            reference,
            c.id,
            0..c.length,
            pileup_params.clone(),
            germline_params,
            out,
        )?;
        if ws.bytes > max_ws.bytes {
            max_ws = ws;
        }
        skips.accumulate(&sk);
    }
    Ok((max_ws, skips))
}

/// Stream a banded gVCF over an indexed interval or shard selection.
pub fn stream_gvcf_selected_bam<W: Write>(
    bam_path: &std::path::Path,
    reference: &dyn ReferenceProvider,
    selection: &AnalysisSelection,
    pileup_params: PileupParams,
    germline_params: &GermlineParams,
    out: &mut W,
) -> Result<(WorkingSet, SkipCounts), CoreError> {
    let intervals = selection.intervals().ok_or_else(|| {
        CoreError::MalformedRecord("selected gVCF driver requires intervals".into())
    })?;
    let mut bander = GvcfBander::new(reference.contigs());
    let mut max_working_set = WorkingSet { bytes: 0 };
    let mut skips = SkipCounts::default();
    for interval in intervals.intervals() {
        governor::checkpoint()?;
        let contig = reference.contigs().by_id(interval.contig).ok_or_else(|| {
            CoreError::MalformedRecord(format!("unknown interval contig {}", interval.contig))
        })?;
        let start = contig.global_offset as usize + interval.start as usize;
        let end = contig.global_offset as usize + interval.end as usize;
        let sequence = reference.decode_window_arc(start, end);
        let source = IndexedBamRegionSource::new(bam_path, reference.contigs(), *interval)?;
        let (working_set, interval_skips) = stream_gvcf_region(
            &mut bander,
            source,
            sequence,
            interval.contig,
            interval.start..interval.end,
            pileup_params.clone(),
            germline_params,
            out,
        )?;
        max_working_set.bytes = max_working_set.bytes.max(working_set.bytes);
        skips.accumulate(&interval_skips);
    }
    Ok((max_working_set, skips))
}
