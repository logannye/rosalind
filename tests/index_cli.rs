//! Phase B3c CLI gates: `rosalind index` + `rosalind locate` (build-once → query).

use std::env;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{SystemTime, UNIX_EPOCH};

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_rosalind")
}

fn tmpdir() -> PathBuf {
    // Cargo runs these tests as parallel threads in one process, so a timestamp
    // alone can collide; a per-call atomic counter guarantees uniqueness.
    static COUNTER: AtomicU64 = AtomicU64::new(0);
    let n = COUNTER.fetch_add(1, Ordering::Relaxed);
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let d = env::temp_dir().join(format!("rosalind-b3c-{nanos}-{n}"));
    std::fs::create_dir_all(&d).unwrap();
    d
}

// chr1/chr2/chr3 each 30 bp (90 bp total). GATTACA occurs in chr3 at local
// positions 0 and 11 (2 hits, both in chr3).
fn write_fasta(dir: &Path) -> PathBuf {
    let p = dir.join("ref.fa");
    std::fs::write(
        &p,
        ">chr1\nACGTACGTNNACGTACGTACGTAAGGCCTT\n\
         >chr2\nTTTTGGGGCCCCAAAANNNNACGTACGTAC\n\
         >chr3\nGATTACATTTTGATTACAGGGGGCCCCAAA\n",
    )
    .unwrap();
    p
}

#[test]
fn index_builds_a_loadable_artifact() {
    let dir = tmpdir();
    let fa = write_fasta(&dir);
    let idx = dir.join("ref.idx");

    let out = Command::new(bin())
        .args([
            "index",
            "--reference",
            fa.to_str().unwrap(),
            "--output",
            idx.to_str().unwrap(),
        ])
        .output()
        .expect("run index");
    assert!(
        out.status.success(),
        "index failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert!(idx.exists(), "index file not written");

    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(
        stdout.contains("contigs: 3 (90 bp total)"),
        "receipt: {stdout}"
    );
    assert!(stdout.contains("  chr1\t30"), "receipt: {stdout}");
    assert!(stdout.contains("index_bytes: "), "receipt: {stdout}");

    let loaded = rosalind::genomics::IndexReader::open(&idx).expect("open");
    assert_eq!(loaded.contigs().len(), 3);

    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn index_memory_budget_prints_plan_line_and_never_refuses() {
    let dir = tmpdir();
    let fa = write_fasta(&dir);
    let idx = dir.join("ref.idx");

    let out = Command::new(bin())
        .args([
            "index",
            "--reference",
            fa.to_str().unwrap(),
            "--output",
            idx.to_str().unwrap(),
            "--memory-budget-mb",
            "0",
        ])
        .output()
        .expect("run index");
    assert!(
        out.status.success(),
        "index must succeed even when over budget"
    );
    assert!(idx.exists(), "index written despite over-budget");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("plan:") && stderr.contains("[OVER]"),
        "expected an [OVER] plan line on stderr, got: {stderr}"
    );

    std::fs::remove_dir_all(&dir).ok();
}
