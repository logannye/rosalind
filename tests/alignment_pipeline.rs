use sqrt_space_sim::genomics::BWTAligner;

#[test]
fn align_batch_produces_results() {
    let reference = b"ACGTACGTACGT";
    let mut aligner = BWTAligner::new(reference).expect("aligner should initialize");
    let reads = vec![b"ACGT".as_slice(), b"CGTA".as_slice(), b"GTAC".as_slice()];
    let results = aligner
        .align_batch(reads.iter().copied())
        .expect("batch alignment should succeed");
    assert_eq!(results.len(), reads.len());
    assert!(results.iter().all(|res| res.interval.width() > 0));
}

#[test]
fn alignment_interval_is_within_bounds() {
    let reference = b"TTTACGTAAA";
    let mut aligner = BWTAligner::new(reference).expect("aligner should initialize");
    let result = aligner
        .align_read(b"ACGT")
        .expect("alignment should succeed");
    assert!(!result.interval.is_empty());
    assert!(result.interval.lower < result.interval.upper);
    assert!(result.score > 0.0);
}
