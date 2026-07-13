use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicU64, Ordering};

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_rosalind")
}

fn run(args: &[&str]) -> Output {
    Command::new(bin()).args(args).output().unwrap()
}

fn must_run(args: &[&str]) {
    let output = run(args);
    assert!(
        output.status.success(),
        "command failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

fn unique_dir() -> PathBuf {
    static COUNTER: AtomicU64 = AtomicU64::new(0);
    let id = COUNTER.fetch_add(1, Ordering::Relaxed);
    let path = std::env::temp_dir().join(format!("rosalind-doctor-{}-{}", std::process::id(), id));
    std::fs::create_dir_all(&path).unwrap();
    path
}

fn fixture(dir: &Path) -> (PathBuf, PathBuf) {
    let reference = dir.join("reference.fa");
    let reads = dir.join("reads.fastq");
    let index = dir.join("reference.idx");
    let raw = dir.join("raw.bam");
    let sorted = dir.join("sorted.bam");
    std::fs::write(&reference, ">chr1\nACGTACGTACGTACGTACGTACGTACGT\n").unwrap();
    std::fs::write(&reads, "@r1\nACGTACGT\n+\nIIIIIIII\n").unwrap();
    must_run(&[
        "index",
        "--reference",
        reference.to_str().unwrap(),
        "--output",
        index.to_str().unwrap(),
    ]);
    must_run(&[
        "align",
        "--reference",
        reference.to_str().unwrap(),
        "--reads",
        reads.to_str().unwrap(),
        "--format",
        "bam",
        "--output",
        raw.to_str().unwrap(),
    ]);
    must_run(&[
        "sort",
        "--input",
        raw.to_str().unwrap(),
        "--output",
        sorted.to_str().unwrap(),
    ]);
    (index, sorted)
}

#[test]
fn doctor_proves_a_ready_fixture_and_reports_output_safety() {
    let dir = unique_dir();
    let (index, bam) = fixture(&dir);
    let output_path = dir.join("calls.vcf");
    let output = run(&[
        "doctor",
        "--index",
        index.to_str().unwrap(),
        "--alignments",
        bam.to_str().unwrap(),
        "--output",
        output_path.to_str().unwrap(),
        "--budget-mb",
        "128",
        "--deep",
        "--json",
    ]);
    assert!(
        output.status.success(),
        "stdout: {}\nstderr: {}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );
    let json = String::from_utf8(output.stdout).unwrap();
    for expected in [
        "\"ok\":true",
        "\"declared_sort_order\":\"coordinate\"",
        "\"coordinate_order_proven\":true",
        "\"output_safe\":true",
        "\"budget_feasible\":true",
    ] {
        assert!(json.contains(expected), "missing {expected}: {json}");
    }
    std::fs::remove_dir_all(dir).ok();
}

#[test]
fn doctor_reports_an_existing_output_as_actionable() {
    let dir = unique_dir();
    let (index, bam) = fixture(&dir);
    let output_path = dir.join("calls.vcf");
    std::fs::write(&output_path, b"existing\n").unwrap();
    let output = run(&[
        "doctor",
        "--index",
        index.to_str().unwrap(),
        "--alignments",
        bam.to_str().unwrap(),
        "--output",
        output_path.to_str().unwrap(),
        "--json",
    ]);
    assert_eq!(output.status.code(), Some(2));
    let json = String::from_utf8(output.stdout).unwrap();
    assert!(json.contains("\"output_safe\":false"), "{json}");
    assert!(json.contains("output reservation collides"), "{json}");
    std::fs::remove_dir_all(dir).ok();
}
