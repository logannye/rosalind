//! `rosalind pack` — the contract's prediction turned into a placement decision.
//! Each job's peak is read from its index header (no run); the schedule proves
//! every node fits before any job launches, or refuses (exit 3).

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

fn build_index(dir: &std::path::Path, name: &str) -> PathBuf {
    let fa = dir.join(format!("{name}.fa"));
    // A small single contig — enough for a real index + a real predicted peak.
    std::fs::write(&fa, format!(">{name}\n{}\n", "ACGTACGTAC".repeat(40))).unwrap();
    let idx = dir.join(format!("{name}.idx"));
    let out = Command::new(bin())
        .args(["index", "--reference"])
        .arg(&fa)
        .arg("--output")
        .arg(&idx)
        .output()
        .unwrap();
    assert!(out.status.success(), "index build failed: {out:?}");
    idx
}

#[test]
fn pack_proves_a_co_location_fits_and_refuses_when_it_cannot() {
    let dir = unique_dir("rosalind-pack");
    let a = build_index(&dir, "a");
    let b = build_index(&dir, "b");
    let jobs = dir.join("jobs.tsv");
    std::fs::write(
        &jobs,
        format!(
            "{}\n# a comment line\n\n{}\t500\t150\n",
            a.display(),
            b.display()
        ),
    )
    .unwrap();

    // A generous node fits both jobs on one node, proven up front.
    let out = Command::new(bin())
        .args(["pack", "--jobs"])
        .arg(&jobs)
        .args(["--node-mb", "512"])
        .output()
        .unwrap();
    assert!(out.status.success(), "generous node should pack: {out:?}");
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(
        stdout.contains("every node proven within capacity"),
        "missing proof line: {stdout}"
    );

    // JSON form is machine-readable.
    let outj = Command::new(bin())
        .args(["pack", "--jobs"])
        .arg(&jobs)
        .args(["--node-mb", "512", "--json"])
        .output()
        .unwrap();
    assert!(outj.status.success());
    let j = String::from_utf8_lossy(&outj.stdout);
    assert!(
        j.contains("\"node_mb\":512") && j.contains("\"nodes\":["),
        "unexpected json: {j}"
    );

    // A tiny node cannot hold even one job (predicted peak ≫ 4 MiB) → exit 3.
    let refuse = Command::new(bin())
        .args(["pack", "--jobs"])
        .arg(&jobs)
        .args(["--node-mb", "4"])
        .output()
        .unwrap();
    assert_eq!(
        refuse.status.code(),
        Some(3),
        "an impossible packing must exit 3: {refuse:?}"
    );
    assert!(
        String::from_utf8_lossy(&refuse.stderr).contains("REFUSE"),
        "missing REFUSE message"
    );

    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn plan_json_emits_machine_readable_prediction() {
    let dir = unique_dir("rosalind-planjson");
    let idx = build_index(&dir, "p");
    let out = Command::new(bin())
        .args(["plan", "--index"])
        .arg(&idx)
        .args(["--budget-mb", "512", "--json"])
        .output()
        .unwrap();
    assert!(out.status.success());
    let j = String::from_utf8_lossy(&out.stdout);
    for needle in [
        "\"predicted_peak_rss_bytes\":",
        "\"largest_contig_len\":",
        "\"verdict\":\"fits\"",
    ] {
        assert!(j.contains(needle), "plan --json missing {needle}: {j}");
    }
    std::fs::remove_dir_all(&dir).ok();
}
