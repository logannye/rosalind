//! `rosalind diff <a> <b>`: claim-level divergence localization between two receipts.

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
    let d = env::temp_dir().join(format!("rosalind-diff-{nanos}-{n}"));
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

fn features(dir: &Path, idx: &Path, bam: &Path, out: &str, max_depth: &str) -> PathBuf {
    let tsv = dir.join(out);
    assert!(run(&[
        "features",
        "--index",
        idx.to_str().unwrap(),
        "--alignments",
        bam.to_str().unwrap(),
        "--max-depth",
        max_depth,
        "-o",
        tsv.to_str().unwrap(),
    ])
    .status
    .success());
    PathBuf::from(format!("{}.manifest.json", tsv.display()))
}

#[test]
fn diff_localizes_a_science_param_change_and_exits_one() {
    let d = tmpdir();
    let (idx, bam) = build_index_and_sorted_bam(&d);
    let a = features(&d, &idx, &bam, "a.tsv", "1000");
    let b = features(&d, &idx, &bam, "b.tsv", "500");

    let out = run(&["diff", a.to_str().unwrap(), b.to_str().unwrap()]);
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert_eq!(
        out.status.code(),
        Some(1),
        "differing claims must exit 1. stdout:\n{stdout}"
    );
    assert!(
        stdout.contains("max_depth"),
        "must localize the param change: {stdout}"
    );
    assert!(stdout.contains("DIFFER"), "{stdout}");

    std::fs::remove_dir_all(&d).ok();
}

#[test]
fn diff_of_a_receipt_against_itself_is_identical_and_exits_zero() {
    let d = tmpdir();
    let (idx, bam) = build_index_and_sorted_bam(&d);
    let a = features(&d, &idx, &bam, "a.tsv", "1000");

    let out = run(&["diff", a.to_str().unwrap(), a.to_str().unwrap()]);
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert_eq!(
        out.status.code(),
        Some(0),
        "identical claims must exit 0: {stdout}"
    );
    assert!(stdout.contains("IDENTICAL"), "{stdout}");

    std::fs::remove_dir_all(&d).ok();
}
