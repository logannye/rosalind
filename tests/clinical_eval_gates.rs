use rosalind::genomics::{compare_callsets, read_vcf_variants, BedIndex};
use std::collections::BTreeMap;

#[test]
fn eval_compare_smoke_and_thresholds() {
    // Small deterministic reference (single contig 'chr1').
    let references = BTreeMap::from([("chr1".to_string(), vec![b'A'; 100])]);

    let calls_vcf = "\
##fileformat=VCFv4.3
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t11\t.\tA\tC\t50\tPASS\t.
chr1\t21\t.\tA\tG\t50\tPASS\t.
";
    let truth_vcf = "\
##fileformat=VCFv4.3
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t11\t.\tA\tC\t50\tPASS\t.
chr1\t31\t.\tA\tT\t50\tPASS\t.
";

    let calls = read_vcf_variants(calls_vcf).unwrap();
    let truth = read_vcf_variants(truth_vcf).unwrap();

    // Mask to only evaluate positions 0..40.
    let bed = BedIndex::from_str("chr1\t0\t40\n").unwrap();
    let rep = compare_callsets(&references, &calls, &truth, Some(&bed)).unwrap();

    // TP: pos11 A>C
    // FP: pos21 A>G
    // FN: pos31 A>T
    assert_eq!(rep.true_positive, 1);
    assert_eq!(rep.false_positive, 1);
    assert_eq!(rep.false_negative, 1);

    // Example CI gate: require non-zero precision/recall on this smoke test.
    assert!(rep.precision() > 0.0);
    assert!(rep.recall() > 0.0);
}
