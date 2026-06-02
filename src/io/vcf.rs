//! One spec-valid VCFv4.2 writer for the calling layer: a typed header built
//! from a `ContigSet`, deterministic record ordering, germline single-sample
//! and somatic TUMOR/NORMAL output. Replaces the legacy non-conformant string
//! writers (no `##contig`/`##FORMAT`/sample columns).

use std::io::{self, Write};

use crate::call::{Filter, Genotype, GermlineCall, SomaticCall};
use crate::core::{ContigSet, Locus};

/// One germline row: the locus, its reference base, and the call there. The
/// call is locus-free, so the writer pairs it with coordinate + ref context.
#[derive(Debug, Clone)]
pub struct GermlineRow {
    /// Genomic coordinate.
    pub locus: Locus,
    /// Reference base (uppercase ASCII).
    pub ref_base: u8,
    /// The germline call at this locus.
    pub call: GermlineCall,
}

fn genotype_str(g: Genotype) -> &'static str {
    match g {
        Genotype::HomRef => "0/0",
        Genotype::Het => "0/1",
        Genotype::HomAlt => "1/1",
    }
}

fn filter_str(f: Filter) -> &'static str {
    match f {
        Filter::Pass => "PASS",
        Filter::LowQual => "LowQual",
        Filter::LowDepth => "LowDepth",
    }
}

/// `##fileformat` + one `##contig` line per contig (shared by both writers).
fn write_fileformat_and_contigs<W: Write>(out: &mut W, contigs: &ContigSet) -> io::Result<()> {
    writeln!(out, "##fileformat=VCFv4.2")?;
    for c in contigs.iter() {
        writeln!(out, "##contig=<ID={},length={}>", c.name, c.length)?;
    }
    Ok(())
}

/// Write the germline VCFv4.2 header (everything up to and including the
/// `#CHROM` line). Pair with [`write_germline_row`] to stream records.
pub fn write_germline_header<W: Write>(
    out: &mut W,
    contigs: &ContigSet,
    sample: &str,
) -> io::Result<()> {
    write_fileformat_and_contigs(out, contigs)?;
    writeln!(
        out,
        r#"##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth">"#
    )?;
    writeln!(
        out,
        r#"##INFO=<ID=AF,Number=A,Type=Float,Description="Alt allele fraction">"#
    )?;
    writeln!(
        out,
        r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#
    )?;
    writeln!(
        out,
        r#"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">"#
    )?;
    writeln!(
        out,
        r#"##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">"#
    )?;
    writeln!(
        out,
        r#"##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths (ref,alt)">"#
    )?;
    writeln!(
        out,
        r#"##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred genotype likelihoods">"#
    )?;
    writeln!(
        out,
        r#"##FILTER=<ID=LowQual,Description="QUAL below threshold">"#
    )?;
    writeln!(
        out,
        r#"##FILTER=<ID=LowDepth,Description="Depth below threshold">"#
    )?;
    writeln!(
        out,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}"
    )
}

/// Write one germline record line. The caller supplies rows in canonical
/// (contig, pos, ref, alt) order — this does not sort (use [`write_germline_vcf`]
/// for an unordered batch).
pub fn write_germline_row<W: Write>(
    out: &mut W,
    contigs: &ContigSet,
    r: &GermlineRow,
) -> io::Result<()> {
    let chrom = contigs
        .by_id(r.locus.contig)
        .map(|c| c.name.as_ref())
        .unwrap_or(".");
    let pos = r.locus.pos.0 as u64 + 1;
    let total = (r.call.ad[0] + r.call.ad[1]).max(1);
    let af = r.call.ad[1] as f64 / total as f64;
    writeln!(
        out,
        "{chrom}\t{pos}\t.\t{ref_b}\t{alt}\t{qual:.1}\t{filt}\tDP={dp};AF={af:.3}\tGT:GQ:DP:AD:PL\t{gt}:{gq}:{dp}:{ad0},{ad1}:{pl0},{pl1},{pl2}",
        ref_b = r.ref_base as char,
        alt = r.call.alt_base as char,
        qual = r.call.qual,
        filt = filter_str(r.call.filter),
        dp = r.call.dp,
        gt = genotype_str(r.call.genotype),
        gq = r.call.gq,
        ad0 = r.call.ad[0],
        ad1 = r.call.ad[1],
        pl0 = r.call.pl[0],
        pl1 = r.call.pl[1],
        pl2 = r.call.pl[2],
    )
}

/// Write a spec-valid germline (single-sample) VCFv4.2 to `out`. Records are
/// emitted in canonical (contig, pos, ref, alt) order regardless of input order.
pub fn write_germline_vcf<W: Write>(
    out: &mut W,
    contigs: &ContigSet,
    sample: &str,
    rows: &[GermlineRow],
) -> io::Result<()> {
    write_germline_header(out, contigs, sample)?;
    let mut ordered: Vec<&GermlineRow> = rows.iter().collect();
    ordered.sort_by(|a, b| {
        a.locus
            .cmp(&b.locus)
            .then_with(|| a.ref_base.cmp(&b.ref_base))
            .then_with(|| a.call.alt_base.cmp(&b.call.alt_base))
    });
    for r in ordered {
        write_germline_row(out, contigs, r)?;
    }
    out.flush()
}

/// Render a germline VCF into a `String` (tests/snapshots).
pub fn render_germline_vcf(
    contigs: &ContigSet,
    sample: &str,
    rows: &[GermlineRow],
) -> io::Result<String> {
    let mut buf = Vec::new();
    write_germline_vcf(&mut buf, contigs, sample, rows)?;
    Ok(String::from_utf8(buf).expect("VCF is valid UTF-8"))
}

/// Write a spec-valid somatic VCFv4.2 with `TUMOR` and `NORMAL` sample columns.
/// Records are emitted in canonical (contig, pos, ref, alt) order.
pub fn write_somatic_vcf<W: Write>(
    out: &mut W,
    contigs: &ContigSet,
    calls: &[(Locus, SomaticCall)],
) -> io::Result<()> {
    write_fileformat_and_contigs(out, contigs)?;
    writeln!(
        out,
        r#"##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic mutation">"#
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
        r#"##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths (ref,alt)">"#
    )?;
    writeln!(
        out,
        r#"##FORMAT=<ID=AF,Number=A,Type=Float,Description="Alt allele fraction">"#
    )?;
    writeln!(
        out,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\tNORMAL"
    )?;

    let mut ordered: Vec<&(Locus, SomaticCall)> = calls.iter().collect();
    ordered.sort_by(|a, b| {
        a.0.cmp(&b.0)
            .then_with(|| a.1.ref_base.cmp(&b.1.ref_base))
            .then_with(|| a.1.alt_base.cmp(&b.1.alt_base))
    });

    for (locus, c) in ordered {
        let chrom = contigs
            .by_id(locus.contig)
            .map(|ct| ct.name.as_ref())
            .unwrap_or(".");
        let pos = locus.pos.0 as u64 + 1;
        let t_ref = c.tumor_depth.saturating_sub(c.tumor_alt);
        let n_ref = c.normal_depth.saturating_sub(c.normal_alt);
        writeln!(
            out,
            "{chrom}\t{pos}\t.\t{ref_b}\t{alt}\t{qual:.1}\tPASS\tSOMATIC\tGT:DP:AD:AF\t0/1:{t_dp}:{t_ref},{t_alt}:{t_af:.3}\t0/0:{n_dp}:{n_ref},{n_alt}:{n_af:.3}",
            ref_b = c.ref_base as char,
            alt = c.alt_base as char,
            qual = c.quality,
            t_dp = c.tumor_depth,
            t_alt = c.tumor_alt,
            t_af = c.tumor_af,
            n_dp = c.normal_depth,
            n_alt = c.normal_alt,
            n_af = c.normal_af,
        )?;
    }
    out.flush()
}

/// Render a somatic VCF into a `String` (tests/snapshots).
pub fn render_somatic_vcf(
    contigs: &ContigSet,
    calls: &[(Locus, SomaticCall)],
) -> io::Result<String> {
    let mut buf = Vec::new();
    write_somatic_vcf(&mut buf, contigs, calls)?;
    Ok(String::from_utf8(buf).expect("VCF is valid UTF-8"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::Position;

    fn contigs() -> ContigSet {
        let mut c = ContigSet::new();
        c.push("chr1", 248_956_422);
        c.push("chr2", 100);
        c
    }

    fn row(contig: u32, pos0: u32, ref_base: u8, call: GermlineCall) -> GermlineRow {
        GermlineRow {
            locus: Locus {
                contig,
                pos: Position(pos0),
            },
            ref_base,
            call,
        }
    }

    fn het_call() -> GermlineCall {
        GermlineCall {
            genotype: Genotype::Het,
            alt_base: b'G',
            qual: 48.0,
            gq: 48,
            pl: [48, 0, 49],
            ad: [15, 15],
            dp: 30,
            filter: Filter::Pass,
        }
    }

    #[test]
    fn header_has_required_tags() {
        let vcf = render_germline_vcf(&contigs(), "SAMPLE", &[]).unwrap();
        assert!(vcf.starts_with("##fileformat=VCFv4.2\n"));
        assert!(vcf.contains("##contig=<ID=chr1,length=248956422>\n"));
        assert!(vcf.contains("##contig=<ID=chr2,length=100>\n"));
        assert!(vcf.contains("##INFO=<ID=DP,Number=1,Type=Integer,"));
        assert!(vcf.contains("##INFO=<ID=AF,Number=A,Type=Float,"));
        assert!(vcf.contains("##FORMAT=<ID=GT,Number=1,Type=String,"));
        assert!(vcf.contains("##FORMAT=<ID=AD,Number=R,Type=Integer,"));
        assert!(vcf.contains("##FORMAT=<ID=PL,Number=G,Type=Integer,"));
        assert!(vcf.contains("##FILTER=<ID=LowQual,"));
        assert!(vcf.contains("##FILTER=<ID=LowDepth,"));
        assert!(vcf.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"));
    }

    #[test]
    fn germline_record_is_exact() {
        let vcf =
            render_germline_vcf(&contigs(), "SAMPLE", &[row(0, 100, b'A', het_call())]).unwrap();
        let last = vcf.lines().last().unwrap();
        // POS is 1-based (100 -> 101); AF = 15/30 = 0.500.
        assert_eq!(
            last,
            "chr1\t101\t.\tA\tG\t48.0\tPASS\tDP=30;AF=0.500\tGT:GQ:DP:AD:PL\t0/1:48:30:15,15:48,0,49"
        );
    }

    #[test]
    fn records_are_sorted_deterministically() {
        let r1 = row(0, 100, b'A', het_call());
        let r2 = row(0, 50, b'A', het_call());
        let r3 = row(1, 10, b'A', het_call());
        let forward =
            render_germline_vcf(&contigs(), "S", &[r1.clone(), r2.clone(), r3.clone()]).unwrap();
        let shuffled = render_germline_vcf(&contigs(), "S", &[r3, r1, r2]).unwrap();
        assert_eq!(forward, shuffled);
        // chr1:51 sorts before chr1:101 before chr2:11.
        let body: Vec<&str> = forward.lines().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(body[0].split('\t').next().unwrap(), "chr1");
        assert!(body[0].contains("\t51\t"));
        assert!(body[1].contains("\t101\t"));
        assert_eq!(body[2].split('\t').next().unwrap(), "chr2");
    }

    #[test]
    fn filter_and_genotype_strings_render() {
        let mut c = het_call();
        c.genotype = Genotype::HomAlt;
        c.filter = Filter::LowDepth;
        let vcf = render_germline_vcf(&contigs(), "S", &[row(0, 0, b'C', c)]).unwrap();
        let last = vcf.lines().last().unwrap();
        assert!(last.contains("\tLowDepth\t"));
        assert!(last.contains("\t1/1:"));
    }

    fn somatic_call() -> SomaticCall {
        SomaticCall {
            ref_base: b'A',
            alt_base: b'T',
            tumor_alt: 12,
            tumor_depth: 40,
            normal_alt: 0,
            normal_depth: 38,
            tumor_af: 0.3,
            normal_af: 0.0,
            quality: 55.0,
        }
    }

    #[test]
    fn somatic_header_has_two_samples() {
        let vcf = render_somatic_vcf(&contigs(), &[]).unwrap();
        assert!(vcf.starts_with("##fileformat=VCFv4.2\n"));
        assert!(vcf.contains("##INFO=<ID=SOMATIC,"));
        assert!(vcf.contains("##FORMAT=<ID=AF,Number=A,Type=Float,"));
        assert!(
            vcf.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\tNORMAL\n")
        );
    }

    #[test]
    fn somatic_record_is_exact() {
        let vcf = render_somatic_vcf(
            &contigs(),
            &[(
                Locus {
                    contig: 0,
                    pos: Position(200),
                },
                somatic_call(),
            )],
        )
        .unwrap();
        let last = vcf.lines().last().unwrap();
        // TUMOR AD = [40-12, 12] = 28,12; NORMAL AD = [38, 0].
        assert_eq!(
            last,
            "chr1\t201\t.\tA\tT\t55.0\tPASS\tSOMATIC\tGT:DP:AD:AF\t0/1:40:28,12:0.300\t0/0:38:38,0:0.000"
        );
    }

    #[test]
    fn header_then_streamed_rows_equals_batch_write() {
        let r1 = row(0, 100, b'A', het_call());
        let r2 = row(0, 50, b'A', het_call());
        let r3 = row(1, 10, b'A', het_call());
        // Batch writer (sorts internally).
        let batch =
            render_germline_vcf(&contigs(), "S", &[r1.clone(), r2.clone(), r3.clone()]).unwrap();
        // Streaming: header once, then rows in already-sorted (contig,pos) order.
        let mut buf = Vec::new();
        write_germline_header(&mut buf, &contigs(), "S").unwrap();
        for r in [&r2, &r1, &r3] {
            write_germline_row(&mut buf, &contigs(), r).unwrap();
        }
        let streamed = String::from_utf8(buf).unwrap();
        assert_eq!(
            streamed, batch,
            "streamed header+rows must equal the batch write byte-for-byte"
        );
    }
}
