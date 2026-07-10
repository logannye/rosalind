//! `rosalind analyze <kind>`: a registered ColumnAnalyzer with a verifiable receipt
//! whose claim records the analyzer's own params.

use std::env;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{SystemTime, UNIX_EPOCH};

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_rosalind")
}

fn tmpdir() -> PathBuf {
    static COUNTER: AtomicU64 = AtomicU64::new(0);
    let n = COUNTER.fetch_add(1, Ordering::Relaxed);
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let d = env::temp_dir().join(format!("rosalind-analyze-{nanos}-{n}"));
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn run(args: &[&str]) -> std::process::Output {
    Command::new(bin())
        .args(args)
        .output()
        .expect("spawn rosalind")
}

fn build_index_and_sorted_bam(dir: &Path) -> (PathBuf, PathBuf) {
    let seq = "ACGTACGTACGTACGTACGTACGTACGTACGT";
    let fa = dir.join("ref.fa");
    std::fs::write(&fa, format!(">chr1\n{seq}\n")).unwrap();
    let fq = dir.join("reads.fq");
    let mut s = String::new();
    for (i, &start) in [0usize, 0, 8].iter().enumerate() {
        let read = &seq[start..start + 16];
        let qual: String = std::iter::repeat_n('I', 16).collect();
        s.push_str(&format!("@r{i}\n{read}\n+\n{qual}\n"));
    }
    std::fs::write(&fq, s).unwrap();
    let idx = dir.join("ref.idx");
    let raw = dir.join("raw.bam");
    let bam = dir.join("sorted.bam");
    assert!(run(&[
        "index",
        "--reference",
        fa.to_str().unwrap(),
        "--output",
        idx.to_str().unwrap()
    ])
    .status
    .success());
    assert!(run(&[
        "align",
        "--reference",
        fa.to_str().unwrap(),
        "--reads",
        fq.to_str().unwrap(),
        "--format",
        "bam",
        "--output",
        raw.to_str().unwrap()
    ])
    .status
    .success());
    assert!(run(&[
        "sort",
        "--input",
        raw.to_str().unwrap(),
        "--output",
        bam.to_str().unwrap()
    ])
    .status
    .success());
    (idx, bam)
}

#[test]
fn analyze_coverage_writes_a_verifiable_receipt_with_analyzer_params() {
    let d = tmpdir();
    let (idx, bam) = build_index_and_sorted_bam(&d);
    let tsv = d.join("cov.tsv");

    let out = run(&[
        "analyze",
        "coverage",
        "--index",
        idx.to_str().unwrap(),
        "--alignments",
        bam.to_str().unwrap(),
        "-o",
        tsv.to_str().unwrap(),
    ]);
    assert!(
        out.status.success(),
        "analyze coverage failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    // The output TSV carries the coverage header.
    let tsv_text = std::fs::read_to_string(&tsv).unwrap();
    assert!(
        tsv_text.starts_with("#contig\tpos\tdepth"),
        "tsv header: {tsv_text}"
    );

    // The receipt's CLAIM records the analyzer's params under the analyzer. prefix.
    let manifest = format!("{}.manifest.json", tsv.display());
    let json = std::fs::read_to_string(&manifest).unwrap();
    assert!(
        json.contains("\"analyzer.analyzer\":\"coverage\""),
        "receipt missing analyzer.analyzer=coverage: {json}"
    );
    for identity in [
        "\"producer.name\":\"rosalind\"",
        "\"producer.binary\":\"rosalind\"",
        "\"analyzer.id\":\"coverage\"",
        "\"analyzer.version\":",
        "\"replay_schema\":\"3\"",
        "\"replay.kind\":\"rosalind\"",
        "\"command_argv\":",
    ] {
        assert!(
            json.contains(identity),
            "receipt missing {identity}: {json}"
        );
    }

    // And it verifies.
    let v = run(&["verify", "--manifest", &manifest]);
    assert!(
        v.status.success(),
        "verify must pass: {}",
        String::from_utf8_lossy(&v.stderr)
    );

    std::fs::remove_dir_all(&d).ok();
}
