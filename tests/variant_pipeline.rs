use std::sync::Arc;

use sqrt_space_sim::genomics::{AlignedRead, CigarOp, CigarOpKind, StreamingVariantCaller};

#[test]
fn streaming_variant_caller_detects_simple_variant() {
    let chrom = Arc::from("chrTest");
    let reference = Arc::from(b"ACGTACGT".to_vec().into_boxed_slice());

    let reads = vec![
        AlignedRead::new(
            Arc::clone(&chrom),
            0,
            vec![CigarOp::new(CigarOpKind::Match, 4)],
            b"ACGT".to_vec(),
            vec![30; 4],
            false,
        ),
        AlignedRead::new(
            Arc::clone(&chrom),
            3,
            vec![CigarOp::new(CigarOpKind::Match, 4)],
            b"TAAA".to_vec(),
            vec![35; 4],
            false,
        ),
    ];

    let mut caller =
        StreamingVariantCaller::new(Arc::clone(&chrom), Arc::clone(&reference), 0, 4, 5.0, 1e-6)
            .unwrap();

    let variants = caller.call_variants(reads).unwrap();
    assert!(!variants.is_empty());
    let alt = &variants[0];
    assert_eq!(alt.chrom.as_ref(), "chrTest");
    assert!(alt.allele_fraction > 0.0);
    assert!(alt.quality >= 5.0);
}
