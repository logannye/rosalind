//! CLI contract surface: `rosalind plan` predicts feasibility, and
//! `variants --index --enforce` refuses up front / passes within budget.

use std::path::PathBuf;
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_rosalind")
}

// A fresh, collision-free temp dir (atomic counter + nanos — concurrent tests
// must not share a directory).
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

// Build a tiny 2-contig index in a fresh temp dir; return (dir, index_path).
fn build_index() -> (PathBuf, PathBuf) {
    let dir = unique_dir("rosalind-plan");
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
    assert!(
        stdout.contains("[FITS]"),
        "generous budget should FIT: {stdout}"
    );
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
    assert!(
        stdout.contains("plan:"),
        "missing build plan line: {stdout}"
    );
    std::fs::remove_dir_all(&dir).ok();
}

// ---- enforce tests: need a real coordinate-sorted BAM via the CLI pipeline ----

fn run(args: &[&str]) -> std::process::Output {
    Command::new(bin())
        .args(args)
        .output()
        .expect("spawn rosalind")
}

// `index` -> `align --format bam` -> `sort`, mirroring tests/variants_index.rs.
// Returns (dir, index_path, sorted_bam_path). Single-contig (aligner is single-contig).
fn build_sorted_bam_fixture() -> (PathBuf, PathBuf, PathBuf) {
    let dir = unique_dir("rosalind-enforce");
    let seq = "ACGTACGTACGTACGTACGTACGTACGTACGT"; // 32 bp
    let fa = dir.join("ref.fa");
    std::fs::write(&fa, format!(">chr1\n{seq}\n")).unwrap();
    let fq = dir.join("reads.fq");
    let mut s = String::new();
    for (i, &start) in [0usize, 0, 4, 4, 8].iter().enumerate() {
        let read = &seq[start..start + 16];
        let qual: String = std::iter::repeat('I').take(16).collect();
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
fn enforce_refuses_up_front_when_budget_below_predicted() {
    // A 1 MiB budget is below the process baseline alone, so the pre-run check
    // refuses with exit 3 before doing any calling — and writes no VCF.
    let (dir, idx, bam) = build_sorted_bam_fixture();
    let out = Command::new(bin())
        .args(["variants", "--index"])
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
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("REFUSE"),
        "missing refuse message: {stderr}"
    );
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn enforce_passes_within_a_generous_budget() {
    let (dir, idx, bam) = build_sorted_bam_fixture();
    // --manifest into the temp dir (stdout run would otherwise drop the cwd-default
    // sidecar into the repo root).
    let manifest = dir.join("run.manifest.json");
    let out = Command::new(bin())
        .args(["variants", "--index"])
        .arg(&idx)
        .arg("--alignments")
        .arg(&bam)
        .args(["--memory-budget-mb", "4096", "--enforce", "--manifest"])
        .arg(&manifest)
        .output()
        .unwrap();
    assert!(out.status.success(), "generous budget should pass: {out:?}");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("contract: OK"), "missing OK line: {stderr}");
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn stdout_run_persists_a_self_describing_receipt() {
    let (dir, idx, bam) = build_sorted_bam_fixture();
    let manifest = dir.join("run.manifest.json");
    // stdout output (no -o), explicit --manifest so we know where to look.
    let out = Command::new(bin())
        .args(["variants", "--index"])
        .arg(&idx)
        .arg("--alignments")
        .arg(&bam)
        .args(["--memory-budget-mb", "4096", "--enforce", "--manifest"])
        .arg(&manifest)
        .output()
        .unwrap();
    assert!(out.status.success(), "run failed: {out:?}");
    let json = std::fs::read_to_string(&manifest).expect("manifest written");
    for needle in [
        "\"contract_verdict\":\"within\"",
        "\"enforced\":\"true\"",
        "\"max_depth\":\"1000\"",
        "\"memory_budget_mb\":\"4096\"",
        "\"peak_rss_bytes\":",
        "\"max_working_set_bytes\":",
    ] {
        assert!(json.contains(needle), "manifest missing {needle}: {json}");
    }
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn verify_passes_on_an_untampered_run_and_fails_on_a_tampered_output() {
    let (dir, idx, bam) = build_sorted_bam_fixture();
    let vcf = dir.join("calls.vcf");
    let manifest = dir.join("calls.vcf.manifest.json");
    let out = Command::new(bin())
        .args(["variants", "--index"])
        .arg(&idx)
        .arg("--alignments")
        .arg(&bam)
        .args(["--memory-budget-mb", "4096", "--enforce", "-o"])
        .arg(&vcf)
        .output()
        .unwrap();
    assert!(out.status.success(), "run failed: {out:?}");

    // Untampered → verify OK (exit 0).
    let ok = Command::new(bin())
        .args(["verify", "--manifest"])
        .arg(&manifest)
        .output()
        .unwrap();
    assert!(
        ok.status.success(),
        "verify should pass: {}",
        String::from_utf8_lossy(&ok.stderr)
    );
    assert!(String::from_utf8_lossy(&ok.stdout).contains("verify: OK"));

    // Tamper with the output VCF → verify FAILS (exit 5).
    std::fs::write(&vcf, b"##tampered\n").unwrap();
    let bad = Command::new(bin())
        .args(["verify", "--manifest"])
        .arg(&manifest)
        .output()
        .unwrap();
    assert_eq!(
        bad.status.code(),
        Some(5),
        "tampered output must fail verify: {bad:?}"
    );
    assert!(String::from_utf8_lossy(&bad.stderr).contains("hash mismatch"));
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn enforce_breach_exits_4_after_writing_output_and_receipt() {
    // Budget 4096 MiB passes the pre-run exit-3 gate (tiny fixture), but a forced
    // realized peak of 8 GiB trips the post-run breach -> exit 4, with the VCF +
    // receipt still written (the documented "output + receipt written" semantics).
    let (dir, idx, bam) = build_sorted_bam_fixture();
    let vcf = dir.join("calls.vcf");
    let manifest = dir.join("calls.vcf.manifest.json");
    let out = Command::new(bin())
        .args(["variants", "--index"])
        .arg(&idx)
        .arg("--alignments")
        .arg(&bam)
        .args(["--memory-budget-mb", "4096", "--enforce", "-o"])
        .arg(&vcf)
        .env("ROSALIND_FORCE_PEAK_RSS_BYTES", "8589934592") // 8 GiB
        .output()
        .unwrap();
    assert_eq!(
        out.status.code(),
        Some(4),
        "expected breach exit 4: {out:?}"
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("VIOLATED"),
        "missing VIOLATED line: {stderr}"
    );
    // Output + receipt were still written before the breach exit.
    assert!(vcf.exists(), "VCF must be written before exit 4");
    let json = std::fs::read_to_string(&manifest).expect("receipt written");
    assert!(
        json.contains("\"contract_verdict\":\"over\""),
        "verdict should be over: {json}"
    );
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn enforce_aborts_on_a_read_longer_than_declared_max_read_len() {
    // Standard fixture has 16 bp reads. Declare --max-read-len 8 under --enforce
    // with a generous budget: the pre-run estimate (using 8) fits, so the run
    // proceeds to ingest and hits the over-long-read abort.
    let (dir, idx, bam) = build_sorted_bam_fixture();
    let out = Command::new(bin())
        .args(["variants", "--index"])
        .arg(&idx)
        .arg("--alignments")
        .arg(&bam)
        .args([
            "--memory-budget-mb",
            "4096",
            "--max-read-len",
            "8",
            "--max-depth",
            "1000",
            "--enforce",
        ])
        .output()
        .unwrap();
    assert!(!out.status.success(), "over-long read must abort: {out:?}");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("exceeds declared --max-read-len"),
        "missing clear message: {stderr}"
    );
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn receipt_records_skip_counts() {
    // The standard fixture (5 reads, well below the default max-depth 1000) drops
    // nothing → over_max_depth 0. A tight --max-depth 1 forces drops → > 0.
    let (dir, idx, bam) = build_sorted_bam_fixture();
    let manifest = dir.join("run.manifest.json");
    let out = Command::new(bin())
        .args(["variants", "--index"])
        .arg(&idx)
        .arg("--alignments")
        .arg(&bam)
        .args(["--max-depth", "1", "--manifest"])
        .arg(&manifest)
        .output()
        .unwrap();
    assert!(out.status.success(), "run failed: {out:?}");
    let json = std::fs::read_to_string(&manifest).unwrap();
    assert!(
        json.contains("\"over_max_depth\":") && json.contains("\"reads_skipped_total\":"),
        "manifest missing skip fields: {json}"
    );
    // At pos 0 the fixture stacks >1 read; --max-depth 1 must drop at least one.
    let m = rosalind::provenance::RunManifest::from_canonical_json(&json).unwrap();
    let over: u64 = m.params.get("over_max_depth").unwrap().parse().unwrap();
    assert!(over > 0, "tight cap should drop reads, got {over}");
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn estimator_upper_bounds_the_realized_working_set() {
    // Run the real pipeline, read max_working_set_bytes from the receipt, and
    // assert the pure estimator (same shared constants) is a true upper bound for
    // the declared --max-depth / --max-read-len. Deterministic (working-set
    // numbers, not process RSS).
    let (dir, idx, bam) = build_sorted_bam_fixture();
    let vcf = dir.join("calls.vcf");
    let manifest = dir.join("calls.vcf.manifest.json");
    let out = Command::new(bin())
        .args(["variants", "--index"])
        .arg(&idx)
        .arg("--alignments")
        .arg(&bam)
        .args(["--max-depth", "1000", "--max-read-len", "250", "-o"])
        .arg(&vcf)
        .output()
        .unwrap();
    assert!(out.status.success(), "run failed: {out:?}");

    let text = std::fs::read_to_string(&manifest).unwrap();
    let m = rosalind::provenance::RunManifest::from_canonical_json(&text).unwrap();
    let realized: u64 = m
        .params
        .get("max_working_set_bytes")
        .unwrap()
        .parse()
        .unwrap();

    // The fixture's single contig is 32 bp; the estimator's bound at the declared
    // cap must dominate the realized working set.
    let predicted = rosalind::call::estimate_variants_working_set(32, 1000, 250).bytes;
    assert!(
        predicted >= realized,
        "estimator bound {predicted} must be >= realized working set {realized}"
    );
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn verify_rejects_an_internally_inconsistent_manifest() {
    // A generous budget → a real run records contract_verdict 'within'. Flip the
    // recorded verdict to 'over' (a tampered/corrupt receipt) and `verify` must
    // reject it on internal consistency — even though the file hashes still match
    // (the manifest has no self-hash yet; these cross-checks are the first line).
    let (dir, idx, bam) = build_sorted_bam_fixture();
    let vcf = dir.join("calls.vcf");
    let manifest = dir.join("calls.vcf.manifest.json");
    let out = Command::new(bin())
        .args(["variants", "--index"])
        .arg(&idx)
        .arg("--alignments")
        .arg(&bam)
        .args(["--memory-budget-mb", "512", "--enforce", "-o"])
        .arg(&vcf)
        .output()
        .unwrap();
    assert!(out.status.success(), "run failed: {out:?}");

    let text = std::fs::read_to_string(&manifest).unwrap();
    assert!(
        text.contains("\"contract_verdict\":\"within\""),
        "expected a 'within' verdict to flip: {text}"
    );
    let tampered = text.replace(
        "\"contract_verdict\":\"within\"",
        "\"contract_verdict\":\"over\"",
    );
    std::fs::write(&manifest, tampered).unwrap();

    let v = Command::new(bin())
        .args(["verify", "--manifest"])
        .arg(&manifest)
        .args(["--budget-mb", "512"])
        .output()
        .unwrap();
    assert_eq!(
        v.status.code(),
        Some(5),
        "an inconsistent manifest must fail verify: {v:?}"
    );
    let stderr = String::from_utf8_lossy(&v.stderr);
    assert!(
        stderr.contains("inconsistent"),
        "expected an internal-inconsistency error: {stderr}"
    );
    std::fs::remove_dir_all(&dir).ok();
}
