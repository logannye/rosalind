//! CLI: `rosalind features` streams a bounded, byte-identical per-locus feature
//! TSV under the memory contract, with a verifiable receipt.

use std::path::PathBuf;
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_rosalind")
}

fn unique_dir(prefix: &str) -> PathBuf {
    static COUNTER: AtomicU64 = AtomicU64::new(0);
    let n = COUNTER.fetch_add(1, Ordering::Relaxed);
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let d = std::env::temp_dir().join(format!("{prefix}-{nanos}-{n}"));
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn run(args: &[&str]) -> std::process::Output {
    Command::new(bin())
        .args(args)
        .output()
        .expect("spawn rosalind")
}

// index -> align (bam) -> sort, returning (dir, index, sorted.bam). Single-contig.
fn build_sorted_bam_fixture() -> (PathBuf, PathBuf, PathBuf) {
    let dir = unique_dir("rosalind-features");
    let seq = "ACGTACGTACGTACGTACGTACGTACGTACGT"; // 32 bp
    let fa = dir.join("ref.fa");
    std::fs::write(&fa, format!(">chr1\n{seq}\n")).unwrap();
    let fq = dir.join("reads.fq");
    let mut s = String::new();
    for (i, &start) in [0usize, 0, 4, 4, 8].iter().enumerate() {
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
    (dir, idx, bam)
}

#[test]
fn features_emit_a_schema_correct_byte_identical_table_with_a_receipt() {
    let (dir, idx, bam) = build_sorted_bam_fixture();
    let tsv = dir.join("features.tsv");
    let tsv2 = dir.join("features2.tsv");
    let manifest = dir.join("features.tsv.manifest.json");

    let out = Command::new(bin())
        .args(["features", "--index"])
        .arg(&idx)
        .arg("--alignments")
        .arg(&bam)
        .arg("-o")
        .arg(&tsv)
        .output()
        .unwrap();
    assert!(out.status.success(), "features failed: {out:?}");

    let text = std::fs::read_to_string(&tsv).unwrap();
    let mut lines = text.lines();
    let header = lines.next().expect("a header line");
    assert!(
        header.starts_with("#contig\tpos\tref\tdepth\traw_depth\t"),
        "header: {header}"
    );
    let first = lines.next().expect("at least one data row");
    let cols: Vec<&str> = first.split('\t').collect();
    assert_eq!(cols.len(), 19, "expected 19 columns, got: {first}");
    assert_eq!(cols[0], "chr1");

    // Byte-identical across a second run = bit-reproducible features.
    let out2 = Command::new(bin())
        .args(["features", "--index"])
        .arg(&idx)
        .arg("--alignments")
        .arg(&bam)
        .arg("-o")
        .arg(&tsv2)
        .output()
        .unwrap();
    assert!(out2.status.success());
    assert_eq!(
        std::fs::read(&tsv).unwrap(),
        std::fs::read(&tsv2).unwrap(),
        "feature TSV must be byte-identical run-to-run"
    );

    // Receipt records the row count; verify re-checks it without re-running.
    let json = std::fs::read_to_string(&manifest).unwrap();
    assert!(
        json.contains("\"feature_rows\":"),
        "receipt missing feature_rows: {json}"
    );
    let m = rosalind::provenance::RunManifest::from_canonical_json(&json).unwrap();
    let rows: u64 = m.params.get("feature_rows").unwrap().parse().unwrap();
    assert!(rows > 0, "expected feature rows, got {rows}");

    let verify = Command::new(bin())
        .args(["verify", "--manifest"])
        .arg(&manifest)
        .output()
        .unwrap();
    assert!(verify.status.success(), "verify failed: {verify:?}");
    assert!(String::from_utf8_lossy(&verify.stdout).contains("verify: OK"));

    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn features_honor_the_memory_contract() {
    // A 1 MiB budget is below the process baseline -> --enforce refuses up front.
    let (dir, idx, bam) = build_sorted_bam_fixture();
    let out = Command::new(bin())
        .args(["features", "--index"])
        .arg(&idx)
        .arg("--alignments")
        .arg(&bam)
        .args(["--memory-budget-mb", "1", "--enforce"])
        .output()
        .unwrap();
    assert_eq!(
        out.status.code(),
        Some(3),
        "expected refuse exit 3: {out:?}"
    );
    assert!(String::from_utf8_lossy(&out.stderr).contains("REFUSE"));
    std::fs::remove_dir_all(&dir).ok();
}
