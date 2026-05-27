#[path = "common/mod.rs"]
mod common;
use common::assert_snapshot;
use rosalind::call::{Filter, Genotype, GermlineCall};
use rosalind::core::{ContigSet, Locus, Position};
use rosalind::io::vcf::{render_germline_vcf, GermlineRow};

#[test]
fn germline_vcf_matches_golden() {
    let mut contigs = ContigSet::new();
    contigs.push("chr1", 100_000);
    let rows = vec![
        GermlineRow {
            locus: Locus {
                contig: 0,
                pos: Position(99),
            },
            ref_base: b'T',
            call: GermlineCall {
                genotype: Genotype::Het,
                alt_base: b'A',
                qual: 42.0,
                gq: 40,
                pl: [42, 0, 60],
                ad: [6, 6],
                dp: 12,
                filter: Filter::Pass,
            },
        },
        GermlineRow {
            locus: Locus {
                contig: 0,
                pos: Position(199),
            },
            ref_base: b'G',
            call: GermlineCall {
                genotype: Genotype::HomAlt,
                alt_base: b'C',
                qual: 88.0,
                gq: 60,
                pl: [120, 60, 0],
                ad: [0, 8],
                dp: 8,
                filter: Filter::Pass,
            },
        },
    ];
    let actual = render_germline_vcf(&contigs, "SAMPLE", &rows).unwrap();
    assert_snapshot("variants/simple.vcf", &actual);
}
