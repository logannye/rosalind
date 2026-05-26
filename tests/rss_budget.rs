use rosalind::genomics::BWTAligner;
use rosalind::util::rss::peak_rss_bytes;

#[test]
fn rss_peak_is_bounded_for_small_workloads() {
    // A real RSS guardrail. CI can tighten this via env var if desired.
    let cap_mb: u64 = std::env::var("ROSALIND_RSS_CAP_MB")
        .ok()
        .and_then(|s| s.parse::<u64>().ok())
        .unwrap_or(2048);
    let cap_bytes = cap_mb.saturating_mul(1024 * 1024);

    // Small-ish workload: build FM-index and align a handful of reads.
    let reference = vec![b'A'; 50_000];
    let mut aligner = BWTAligner::new(&reference).expect("aligner build");
    for _ in 0..100 {
        let _ = aligner.align_read(b"AAAAAAAAAAAAAAAAAAAA").expect("align");
    }

    let peak = peak_rss_bytes();
    assert!(
        peak > 0,
        "peak_rss_bytes() returned 0; RSS measurement unavailable"
    );
    assert!(
        peak <= cap_bytes,
        "peak RSS {peak} bytes exceeded cap {cap_bytes} bytes ({} MiB). Set ROSALIND_RSS_CAP_MB to adjust.",
        cap_mb
    );
}


