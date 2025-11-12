use std::collections::HashSet;
use std::sync::Arc;

use blake3::hash;
use rosalind::genomics::{render_vcf, AlignedRead, CigarOp, CigarOpKind, StreamingVariantCaller};

#[test]
fn streaming_variant_caller_is_deterministic() {
    let chrom = Arc::from("chrDeterministic");
    let reference = Arc::from(b"ACGTACGTACGTACGT".to_vec().into_boxed_slice());
    let read = AlignedRead::new(
        Arc::clone(&chrom),
        4,
        60,
        vec![CigarOp::new(CigarOpKind::Match, 6)],
        b"ACGTAC".to_vec(),
        vec![30; 6],
        false,
    );

    let mut fingerprints = HashSet::new();
    for _ in 0..5 {
        let mut caller = StreamingVariantCaller::new(
            Arc::clone(&chrom),
            Arc::clone(&reference),
            0,
            16,
            5.0,
            1e-6,
        )
        .expect("caller initialises");

        let variants = caller
            .call_variants(vec![read.clone()])
            .expect("variant calling succeeds");
        let vcf = render_vcf(&variants).expect("rendering succeeds");
        fingerprints.insert(hash(vcf.as_bytes()));
    }

    assert_eq!(fingerprints.len(), 1, "outputs diverged across runs");
}
