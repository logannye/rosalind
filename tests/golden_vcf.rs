#[path = "common/mod.rs"]
mod common;
use common::assert_snapshot;
use rosalind::genomics::{render_vcf, Variant};
use std::sync::Arc;

#[test]
fn render_vcf_matches_golden() {
    let variants = vec![
        Variant {
            chrom: Arc::from("chr1"),
            position: 3,
            reference: b'T',
            alternate: b'A',
            depth: 12,
            quality: 42.0,
            allele_fraction: 0.75,
        },
        Variant {
            chrom: Arc::from("chr2"),
            position: 9,
            reference: b'G',
            alternate: b'C',
            depth: 8,
            quality: 18.5,
            allele_fraction: 0.625,
        },
    ];

    let actual = render_vcf(&variants).expect("VCF rendering should succeed");
    assert_snapshot("variants/simple.vcf", &actual);
}
