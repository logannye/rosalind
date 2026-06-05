//! Track D gates: `rosalind reproduce` — third-party byte re-derivation from a receipt.
//! Exercised through the real `index -> align -> sort -> variants -> reproduce` CLI
//! pipeline (no rust-htslib dev-dependency). The DIVERGED classification is unit-tested
//! in `reproduce::tests` (a same-binary same-input run cannot non-deterministically
//! diverge by design); here we cover REPRODUCED + the INCONCLUSIVE preconditions.

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
    let d = env::temp_dir().join(format!("rosalind-repro-{nanos}-{n}"));
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

/// `index` -> `align --format bam` -> `sort` -> `(idx, sorted.bam)`.
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
fn reproduces_a_variants_index_vcf_byte_for_byte() {
    let dir = tmpdir();
    let seq = "ACGTACGTACGTACGTACGTACGTACGTACGT";
    let fa = write_fasta(&dir, "chr1", seq);
    let fq = write_fastq(&dir, seq, &[0, 0, 8], 16);
    let (idx, bam) = build_index_and_sorted_bam(&dir, &fa, &fq);
    let vcf = dir.join("calls.vcf");

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
        "make: {}",
        String::from_utf8_lossy(&made.stderr)
    );

    let manifest = format!("{}.manifest.json", vcf.display());
    let out = run(&[
        "reproduce",
        "--manifest",
        &manifest,
        "--inputs",
        dir.to_str().unwrap(),
    ]);
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(
        out.status.success(),
        "reproduce should exit 0 (REPRODUCED). stdout:\n{stdout}\nstderr:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert!(stdout.contains("REPRODUCED"), "verdict line: {stdout}");
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn inconclusive_when_an_input_is_missing() {
    let dir = tmpdir();
    let seq = "ACGTACGTACGTACGTACGTACGTACGTACGT";
    let fa = write_fasta(&dir, "chr1", seq);
    let fq = write_fastq(&dir, seq, &[0, 0, 8], 16);
    let (idx, bam) = build_index_and_sorted_bam(&dir, &fa, &fq);
    let vcf = dir.join("calls.vcf");
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

    // Point --inputs at an EMPTY dir → the recorded inputs cannot be content-located.
    let empty = tmpdir();
    let manifest = format!("{}.manifest.json", vcf.display());
    let out = run(&[
        "reproduce",
        "--manifest",
        &manifest,
        "--inputs",
        empty.to_str().unwrap(),
    ]);
    assert_eq!(
        out.status.code(),
        Some(7),
        "missing input must be INCONCLUSIVE (exit 7)"
    );
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(stdout.contains("INCONCLUSIVE"), "{stdout}");
    assert!(
        stdout.contains("not located"),
        "names the unlocatable input: {stdout}"
    );
    std::fs::remove_dir_all(&dir).ok();
    std::fs::remove_dir_all(&empty).ok();
}

#[test]
fn inconclusive_on_a_pre_schema_5_receipt() {
    let dir = tmpdir();
    // A minimal, well-formed receipt with no recorded `command` (pre-schema-5 shape) and
    // no manifest_blake3 (so the integrity gate notes-and-skips rather than failing).
    let receipt = dir.join("old.manifest.json");
    std::fs::write(
        &receipt,
        r#"{"inputs":[],"outputs":[{"blake3":"aa","path":"x.vcf"}],"params":{},"subcommand":"variants","tool_version":"0.1.0"}"#,
    )
    .unwrap();
    let out = run(&[
        "reproduce",
        "--manifest",
        receipt.to_str().unwrap(),
        "--inputs",
        dir.to_str().unwrap(),
    ]);
    assert_eq!(
        out.status.code(),
        Some(7),
        "pre-schema-5 receipt must be INCONCLUSIVE"
    );
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(stdout.contains("pre-schema-5"), "explains why: {stdout}");
    std::fs::remove_dir_all(&dir).ok();
}
