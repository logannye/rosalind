//! Provenance-DAG gates: the `index` receipt (PR1) and `chain verify` (PR2).
//! Driven through the real CLI binary (no rust-htslib dev-dependency), reusing the
//! `index -> align -> sort -> variants` toy pipeline from `tests/reproduce.rs`.

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
    let d = env::temp_dir().join(format!("rosalind-chain-{nanos}-{n}"));
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn run(args: &[&str]) -> std::process::Output {
    Command::new(bin())
        .args(args)
        .output()
        .expect("spawn rosalind")
}

fn write_fasta(dir: &Path, name: &str, seq: &str) -> PathBuf {
    let p = dir.join("ref.fa");
    std::fs::write(&p, format!(">{name}\n{seq}\n")).unwrap();
    p
}

fn write_fastq(dir: &Path, seq: &str, starts: &[usize], len: usize) -> PathBuf {
    let p = dir.join("reads.fq");
    let mut s = String::new();
    for (i, &start) in starts.iter().enumerate() {
        let read = &seq[start..start + len];
        let qual: String = std::iter::repeat_n('I', len).collect();
        s.push_str(&format!("@r{i}\n{read}\n+\n{qual}\n"));
    }
    std::fs::write(&p, s).unwrap();
    p
}

/// `index` -> `align --format bam` -> `sort` -> `(idx, sorted.bam)`, all in `dir`.
fn build_index_and_sorted_bam(dir: &Path, fa: &Path, fq: &Path) -> (PathBuf, PathBuf) {
    let idx = dir.join("ref.idx");
    let raw = dir.join("raw.bam");
    let sorted = dir.join("sorted.bam");
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
        raw.to_str().unwrap(),
    ])
    .status
    .success());
    assert!(run(&[
        "sort",
        "--input",
        raw.to_str().unwrap(),
        "--output",
        sorted.to_str().unwrap()
    ])
    .status
    .success());
    (idx, sorted)
}

#[test]
fn index_writes_a_self_verifying_receipt() {
    let d = tmpdir();
    let fa = write_fasta(&d, "chr1", "ACGTACGTACGTACGTACGTACGTACGTACGT");
    let idx = d.join("ref.idx");

    let out = run(&[
        "index",
        "--reference",
        fa.to_str().unwrap(),
        "--output",
        idx.to_str().unwrap(),
    ]);
    assert!(
        out.status.success(),
        "index failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    let manifest = d.join("ref.idx.manifest.json");
    assert!(
        manifest.exists(),
        "index must write a receipt sidecar at {}",
        manifest.display()
    );

    let v = run(&["verify", "--manifest", manifest.to_str().unwrap()]);
    assert!(
        v.status.success(),
        "verify must pass on the index receipt. stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&v.stdout),
        String::from_utf8_lossy(&v.stderr)
    );

    std::fs::remove_dir_all(&d).ok();
}

#[test]
fn index_output_hash_equals_variants_index_input_hash() {
    use rosalind::provenance::RunManifest;

    let d = tmpdir();
    let seq = "ACGTACGTACGTACGTACGTACGTACGTACGT";
    let fa = write_fasta(&d, "chr1", seq);
    let fq = write_fastq(&d, seq, &[0, 0, 8], 16);
    let (idx, bam) = build_index_and_sorted_bam(&d, &fa, &fq);
    let vcf = d.join("calls.vcf");

    let made = run(&[
        "variants",
        "--index",
        idx.to_str().unwrap(),
        "--alignments",
        bam.to_str().unwrap(),
        "-o",
        vcf.to_str().unwrap(),
    ]);
    assert!(
        made.status.success(),
        "variants: {}",
        String::from_utf8_lossy(&made.stderr)
    );

    let index_m = RunManifest::from_canonical_json(
        &std::fs::read_to_string(d.join("ref.idx.manifest.json")).unwrap(),
    )
    .unwrap();
    let variants_m = RunManifest::from_canonical_json(
        &std::fs::read_to_string(format!("{}.manifest.json", vcf.display())).unwrap(),
    )
    .unwrap();

    // The index's recorded .idx output hash == the variants' recorded --index input hash.
    let index_out = &index_m.outputs[0].blake3;
    let variants_index_in = variants_m
        .inputs
        .iter()
        .find(|f| f.path.ends_with("ref.idx"))
        .expect("variants receipt records the index as an input")
        .blake3
        .clone();
    assert_eq!(
        *index_out, variants_index_in,
        "the provenance edge must resolve: index outputs[0] == variants inputs[--index]"
    );

    std::fs::remove_dir_all(&d).ok();
}

#[test]
fn chain_verify_is_intact_on_a_real_index_variants_chain() {
    let d = tmpdir();
    let seq = "ACGTACGTACGTACGTACGTACGTACGTACGT";
    let fa = write_fasta(&d, "chr1", seq);
    let fq = write_fastq(&d, seq, &[0, 0, 8], 16);
    let (idx, bam) = build_index_and_sorted_bam(&d, &fa, &fq);
    let vcf = d.join("calls.vcf");
    assert!(run(&[
        "variants",
        "--index",
        idx.to_str().unwrap(),
        "--alignments",
        bam.to_str().unwrap(),
        "-o",
        vcf.to_str().unwrap(),
    ])
    .status
    .success());

    // d now holds ref.idx.manifest.json + calls.vcf.manifest.json.
    let out = run(&["chain", "verify", d.to_str().unwrap()]);
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(
        out.status.success(),
        "chain verify must exit 0 (INTACT). stdout:\n{stdout}\nstderr:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert!(stdout.contains("CHAIN INTACT"), "verdict line: {stdout}");

    let j = run(&["chain", "verify", d.to_str().unwrap(), "--json"]);
    let js = String::from_utf8_lossy(&j.stdout);
    assert!(js.contains("\"intact\":true"), "json: {js}");
    assert!(js.contains("\"edges_resolved\":1"), "json: {js}");

    std::fs::remove_dir_all(&d).ok();
}

#[test]
fn chain_verify_breaks_when_the_index_receipt_is_tampered() {
    let d = tmpdir();
    let seq = "ACGTACGTACGTACGTACGTACGTACGTACGT";
    let fa = write_fasta(&d, "chr1", seq);
    let fq = write_fastq(&d, seq, &[0, 0, 8], 16);
    let (idx, bam) = build_index_and_sorted_bam(&d, &fa, &fq);
    let vcf = d.join("calls.vcf");
    assert!(run(&[
        "variants",
        "--index",
        idx.to_str().unwrap(),
        "--alignments",
        bam.to_str().unwrap(),
        "-o",
        vcf.to_str().unwrap(),
    ])
    .status
    .success());

    // Tamper with the index receipt's claim: flip a digit of total_bp without
    // recomputing manifest_blake3 → the node self-hash must fail.
    let mpath = d.join("ref.idx.manifest.json");
    let text = std::fs::read_to_string(&mpath).unwrap();
    let tampered = text.replacen("\"total_bp\":\"", "\"total_bp\":\"9", 1);
    assert_ne!(
        text, tampered,
        "the receipt must contain a total_bp field to tamper"
    );
    std::fs::write(&mpath, tampered).unwrap();

    let out = run(&["chain", "verify", d.to_str().unwrap()]);
    assert_eq!(out.status.code(), Some(5), "a tampered node must exit 5");
    assert!(
        String::from_utf8_lossy(&out.stdout).contains("CHAIN BROKEN"),
        "stdout: {}",
        String::from_utf8_lossy(&out.stdout)
    );

    std::fs::remove_dir_all(&d).ok();
}

#[test]
fn chain_verify_reports_the_bam_input_as_external_not_a_failure() {
    let d = tmpdir();
    let seq = "ACGTACGTACGTACGTACGTACGTACGTACGT";
    let fa = write_fasta(&d, "chr1", seq);
    let fq = write_fastq(&d, seq, &[0, 0, 8], 16);
    let (idx, bam) = build_index_and_sorted_bam(&d, &fa, &fq);
    let vcf = d.join("calls.vcf");
    assert!(run(&[
        "variants",
        "--index",
        idx.to_str().unwrap(),
        "--alignments",
        bam.to_str().unwrap(),
        "-o",
        vcf.to_str().unwrap(),
    ])
    .status
    .success());

    let out = run(&["chain", "verify", d.to_str().unwrap()]);
    let stdout = String::from_utf8_lossy(&out.stdout);
    // The BAM alignments input has no producing receipt → external/integrity-only,
    // and it must NOT fail the chain.
    assert!(
        out.status.success(),
        "external inputs must not break the chain: {stdout}"
    );
    assert!(
        stdout.contains("--alignments-->  (external)  [integrity-only]"),
        "the BAM edge must be reported external: {stdout}"
    );

    std::fs::remove_dir_all(&d).ok();
}
