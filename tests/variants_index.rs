//! Phase B4 gates: `rosalind variants --index` — bounded, self-contained,
//! germline calling, parity with `--reference` on the same sorted BAM. Exercised
//! through the real `index -> align -> sort -> variants` CLI pipeline (no
//! rust-htslib dev-dependency; `--index` is BAM-only). Multi-contig calling is
//! gated by the `call::whole_genome` library test.

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
    let d = env::temp_dir().join(format!("rosalind-b4v-{nanos}-{n}"));
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

// A FASTQ of reads that are substrings of `seq` (so they align), giving depth.
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

/// `index` -> `align --format bam` -> `sort` -> returns `(idx, sorted.bam)`.
/// Single-contig (the aligner is single-contig); produces a real sorted BAM via
/// the repo's own pipeline so the streaming `--index` path has valid input.
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
    let a = run(&[
        "align",
        "--reference",
        fa.to_str().unwrap(),
        "--reads",
        fq.to_str().unwrap(),
        "--format",
        "bam",
        "--output",
        raw.to_str().unwrap(),
    ]);
    assert!(
        a.status.success(),
        "align: {}",
        String::from_utf8_lossy(&a.stderr)
    );
    let s = run(&[
        "sort",
        "--input",
        raw.to_str().unwrap(),
        "--output",
        sorted.to_str().unwrap(),
    ]);
    assert!(
        s.status.success(),
        "sort: {}",
        String::from_utf8_lossy(&s.stderr)
    );
    (idx, sorted)
}

fn vcf_records(out: &[u8]) -> Vec<String> {
    String::from_utf8_lossy(out)
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.to_string())
        .collect()
}

#[test]
fn variants_index_matches_reference_on_the_same_sorted_bam() {
    let dir = tmpdir();
    let seq = "ACGTACGTACGTACGTACGTACGTACGTACGT"; // 32 bp
    let fa = write_fasta(&dir, "chr1", seq);
    let fq = write_fastq(&dir, seq, &[0, 0, 4, 4, 8], 16);
    let (idx, bam) = build_index_and_sorted_bam(&dir, &fa, &fq);

    // Same sorted BAM through both paths → identical reads → calls must match.
    let by_ref = run(&[
        "variants",
        "--reference",
        fa.to_str().unwrap(),
        "--alignments",
        bam.to_str().unwrap(),
    ]);
    assert!(
        by_ref.status.success(),
        "--reference: {}",
        String::from_utf8_lossy(&by_ref.stderr)
    );
    let by_idx = run(&[
        "variants",
        "--index",
        idx.to_str().unwrap(),
        "--alignments",
        bam.to_str().unwrap(),
    ]);
    assert!(
        by_idx.status.success(),
        "--index: {}",
        String::from_utf8_lossy(&by_idx.stderr)
    );

    // Record lines only (headers legitimately differ: sample name / ##contig set).
    assert_eq!(
        vcf_records(&by_ref.stdout),
        vcf_records(&by_idx.stdout),
        "variants --index records must match --reference on the same BAM"
    );
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn variants_index_is_self_contained() {
    let dir = tmpdir();
    let seq = "ACGTACGTACGTACGTACGTACGTACGTACGT";
    let fa = write_fasta(&dir, "chr1", seq);
    let fq = write_fastq(&dir, seq, &[0, 0, 8], 16);
    let (idx, bam) = build_index_and_sorted_bam(&dir, &fa, &fq);
    std::fs::remove_file(&fa).unwrap(); // reference FASTA gone

    let out = run(&[
        "variants",
        "--index",
        idx.to_str().unwrap(),
        "--alignments",
        bam.to_str().unwrap(),
    ]);
    assert!(
        out.status.success(),
        "must call from the index alone: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn variants_requires_exactly_one_reference_source() {
    let dir = tmpdir();
    let seq = "ACGTACGTACGTACGT";
    let fa = write_fasta(&dir, "chr1", seq);
    let fq = write_fastq(&dir, seq, &[0], 8);
    let (_idx, bam) = build_index_and_sorted_bam(&dir, &fa, &fq);
    // Neither --index nor --reference → error.
    let neither = run(&["variants", "--alignments", bam.to_str().unwrap()]);
    assert!(
        !neither.status.success(),
        "must require one of --index/--reference"
    );
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn variants_index_rejects_sam() {
    // --index requires a BAM; a .sam alignments file errors with guidance.
    let dir = tmpdir();
    let fa = write_fasta(&dir, "chr1", "ACGTACGTACGTACGT");
    let idx = dir.join("ref.idx");
    assert!(run(&[
        "index",
        "--reference",
        fa.to_str().unwrap(),
        "--output",
        idx.to_str().unwrap()
    ])
    .status
    .success());
    let sam = dir.join("reads.sam");
    std::fs::write(&sam, "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:16\n").unwrap();
    let out = run(&[
        "variants",
        "--index",
        idx.to_str().unwrap(),
        "--alignments",
        sam.to_str().unwrap(),
    ]);
    assert!(!out.status.success(), "--index with SAM must error");
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn variants_index_memory_budget_reports_and_never_refuses() {
    let dir = tmpdir();
    let seq = "ACGTACGTACGTACGTACGTACGTACGTACGT";
    let fa = write_fasta(&dir, "chr1", seq);
    let fq = write_fastq(&dir, seq, &[0, 0, 8], 16);
    let (idx, bam) = build_index_and_sorted_bam(&dir, &fa, &fq);

    // Budget 0 -> reported EXCEEDED, but the call still completes (record-only).
    let out = run(&[
        "variants",
        "--index",
        idx.to_str().unwrap(),
        "--alignments",
        bam.to_str().unwrap(),
        "--memory-budget-mb",
        "0",
    ]);
    assert!(
        out.status.success(),
        "must complete even when over budget: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("peak RSS"),
        "receipt must report realized peak: {stderr}"
    );
    assert!(
        stderr.contains("EXCEEDED") || stderr.contains("budget"),
        "budget verdict: {stderr}"
    );
    std::fs::remove_dir_all(&dir).ok();
}
