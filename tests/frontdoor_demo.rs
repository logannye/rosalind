//! Proves the README's in-house contract demo actually runs end-to-end on the
//! bundled single-contig fixture: index → sort → plan → variants --enforce → verify.

use std::path::PathBuf;
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_rosalind")
}

fn unique_dir() -> PathBuf {
    static C: AtomicU64 = AtomicU64::new(0);
    let n = C.fetch_add(1, Ordering::Relaxed);
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let d = std::env::temp_dir().join(format!("rosalind-frontdoor-{nanos}-{n}"));
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn run(args: &[&str]) -> std::process::Output {
    Command::new(bin())
        .args(args)
        .output()
        .expect("spawn rosalind")
}

#[test]
fn readme_inhouse_contract_demo_runs_end_to_end() {
    // CARGO_MANIFEST_DIR points at the crate root; the fixture is bundled there.
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let fa = root.join("examples/data/illumina_toy/reference.fa");
    let bam = root.join("examples/data/illumina_toy/alignments.bam");
    assert!(fa.exists(), "bundled reference missing: {}", fa.display());
    assert!(
        bam.exists(),
        "bundled alignments missing: {}",
        bam.display()
    );

    let dir = unique_dir();
    let idx = dir.join("toy.idx");
    let sorted = dir.join("toy.sorted.bam");
    let vcf = dir.join("toy.vcf");
    let manifest = dir.join("toy.vcf.manifest.json");

    let out = run(&[
        "index",
        "--reference",
        fa.to_str().unwrap(),
        "--output",
        idx.to_str().unwrap(),
    ]);
    assert!(
        out.status.success(),
        "index: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    let out = run(&[
        "sort",
        "--input",
        bam.to_str().unwrap(),
        "--output",
        sorted.to_str().unwrap(),
    ]);
    assert!(
        out.status.success(),
        "sort: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    let out = run(&[
        "plan",
        "--index",
        idx.to_str().unwrap(),
        "--budget-mb",
        "512",
    ]);
    assert!(
        out.status.success(),
        "plan: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert!(String::from_utf8_lossy(&out.stdout).contains("predicted peak"));

    let out = run(&[
        "variants",
        "--index",
        idx.to_str().unwrap(),
        "--alignments",
        sorted.to_str().unwrap(),
        "--memory-budget-mb",
        "512",
        "--enforce",
        "-o",
        vcf.to_str().unwrap(),
    ]);
    assert!(
        out.status.success(),
        "variants --enforce: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert!(manifest.exists(), "receipt sidecar must be written");

    let out = run(&["verify", "--manifest", manifest.to_str().unwrap()]);
    assert!(
        out.status.success(),
        "verify: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert!(String::from_utf8_lossy(&out.stdout).contains("verify: OK"));

    std::fs::remove_dir_all(&dir).ok();
}
