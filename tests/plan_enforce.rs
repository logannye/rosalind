//! CLI contract surface: `rosalind plan` predicts feasibility, and
//! `variants --index --enforce` refuses up front / passes within budget.

use std::path::PathBuf;
use std::process::Command;

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_rosalind")
}

// Build a tiny 2-contig index in a fresh temp dir; return (dir, index_path).
fn build_index() -> (PathBuf, PathBuf) {
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!("rosalind-plan-{nanos}"));
    std::fs::create_dir_all(&dir).unwrap();
    let fa = dir.join("ref.fa");
    std::fs::write(
        &fa,
        b">chr1\nACGTACGTACGTACGTACGT\n>chr2\nTTTTGGGGCCCCAAAATTTT\n",
    )
    .unwrap();
    let idx = dir.join("ref.idx");
    let out = Command::new(bin())
        .args(["index", "--reference"])
        .arg(&fa)
        .arg("--output")
        .arg(&idx)
        .output()
        .unwrap();
    assert!(out.status.success(), "index build failed: {out:?}");
    (dir, idx)
}

#[test]
fn plan_index_reports_a_breakdown_and_fits_a_generous_budget() {
    let (dir, idx) = build_index();
    let out = Command::new(bin())
        .args(["plan", "--index"])
        .arg(&idx)
        .args(["--budget-mb", "4096"])
        .output()
        .unwrap();
    assert!(out.status.success());
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(
        stdout.contains("predicted peak"),
        "missing breakdown: {stdout}"
    );
    assert!(stdout.contains("[FITS]"), "generous budget should FIT: {stdout}");
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn plan_reference_reports_build_estimate() {
    let (dir, _idx) = build_index();
    let fa = dir.join("ref.fa");
    let out = Command::new(bin())
        .args(["plan", "--reference"])
        .arg(&fa)
        .args(["--budget-mb", "4096"])
        .output()
        .unwrap();
    assert!(out.status.success());
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(stdout.contains("plan:"), "missing build plan line: {stdout}");
    std::fs::remove_dir_all(&dir).ok();
}
