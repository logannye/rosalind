//! Regression test: `align --format bam` must write raw Phred (0–93) into the
//! BAM QUAL field, not ASCII Phred+33 (33–126).
//!
//! Bug: the FASTQ parser stored quality bytes as raw ASCII (e.g. b'I' = 73 for
//! Q40).  The BAM write path previously passed those ASCII bytes unchanged into
//! `Record::set(...)`, so BAM files stored 73 instead of 40.  The germline
//! variant caller then read those inflated values and produced over-confident
//! calls.  The SAM path was unaffected because SAM is ASCII text and the SAM
//! reader explicitly subtracts 33 on read-back.
//!
//! This test:
//!   1. Writes a tiny FASTA reference and a matching FASTQ with known qualities.
//!   2. Runs `rosalind align --format bam` (via the compiled binary).
//!   3. Reads the BAM back via `rosalind::io::bam::read_bam_as_core_reads`.
//!   4. Asserts every quality byte equals the expected *raw* Phred value, NOT
//!      the ASCII-encoded value.

use std::io::Write;
use std::path::PathBuf;
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

use rosalind::core::ContigSet;
use rosalind::io::bam::read_bam_as_core_reads;

fn ts() -> u128 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos()
}

fn tmp(label: &str, ext: &str) -> PathBuf {
    std::env::temp_dir().join(format!("rosalind-bam-qual-{label}-{}.{ext}", ts()))
}

fn write_text(path: &PathBuf, contents: &str) {
    let mut f = std::fs::File::create(path).expect("create temp file");
    f.write_all(contents.as_bytes())
        .expect("write temp file contents");
}

/// ASCII-encoded FASTQ quality byte for Phred Q `p` is `p + 33`.
fn ascii_for_phred(p: u8) -> u8 {
    p + 33
}

#[test]
fn bam_qual_field_stores_raw_phred_not_ascii() {
    // ── reference ────────────────────────────────────────────────────────────
    // A short but unambiguous reference so our reads align uniquely.
    let reference_seq = "ACGTACGTACGTACGT";
    let fasta_path = tmp("ref", "fa");
    write_text(&fasta_path, &format!(">chr1\n{reference_seq}\n"));

    // ── reads with known qualities ────────────────────────────────────────────
    // Two reads that exactly match a prefix of the reference.
    //   read1: ACGTACGT  qualities: IIIIIIII  (Q40 each, ASCII 73)
    //   read2: CGTACGTA  qualities: !!!!!!!!  (Q0  each, ASCII 33)
    let expected_q_read1: Vec<u8> = vec![40; 8]; // raw Phred
    let expected_q_read2: Vec<u8> = vec![0; 8]; // raw Phred

    let qual_str_read1: String = expected_q_read1
        .iter()
        .map(|q| ascii_for_phred(*q) as char)
        .collect();
    let qual_str_read2: String = expected_q_read2
        .iter()
        .map(|q| ascii_for_phred(*q) as char)
        .collect();

    let fastq_contents = format!(
        "@read1\nACGTACGT\n+\n{}\n@read2\nCGTACGTA\n+\n{}\n",
        qual_str_read1, qual_str_read2
    );
    let fastq_path = tmp("reads", "fastq");
    write_text(&fastq_path, &fastq_contents);

    // ── run `rosalind align --format bam` ────────────────────────────────────
    let bam_path = tmp("out", "bam");

    let bin = env!("CARGO_BIN_EXE_rosalind");
    let status = Command::new(bin)
        .args([
            "align",
            "--reference",
            fasta_path.to_str().unwrap(),
            "--reads",
            fastq_path.to_str().unwrap(),
            "--format",
            "bam",
            "--output",
            bam_path.to_str().unwrap(),
        ])
        .status()
        .expect("failed to run rosalind align");

    assert!(
        status.success(),
        "rosalind align exited with non-zero status: {status}"
    );
    assert!(
        bam_path.exists(),
        "BAM output file was not created at {bam_path:?}"
    );

    // ── read the BAM back and verify qualities are raw Phred ──────────────────
    let mut contigs = ContigSet::new();
    contigs.push("chr1", reference_seq.len() as u32);

    let reads =
        read_bam_as_core_reads(&bam_path, &contigs).expect("should be able to read BAM output");

    // Both reads should have mapped (reference is long enough).
    assert!(
        !reads.is_empty(),
        "expected at least one mapped read in BAM; got none"
    );

    for read in &reads {
        let name = std::str::from_utf8(read.seq.as_ref()).unwrap_or("<non-utf8>");
        let _ = name;

        for (i, &qual) in read.qual.iter().enumerate() {
            assert!(
                qual <= 93,
                "BAM quality byte[{i}] = {qual} is ≥ 94, which is not a valid raw Phred score \
                 (it looks like an ASCII Phred+33 value was stored unmodified)"
            );
        }
    }

    // Specifically verify read1 (Q40 = raw 40, not 73).
    let read1 = reads.iter().find(|r| {
        // read1 starts at pos 0, read2 at pos 1 (CGTACGTA offset by 1).
        r.pos.0 == 0
    });
    if let Some(r) = read1 {
        for (i, &qual) in r.qual.iter().enumerate() {
            assert_eq!(
                qual, 40,
                "read1 quality byte[{i}]: expected raw Phred 40 (Q40), got {qual} \
                 (73 would indicate ASCII 'I' was written verbatim instead of raw Phred)"
            );
        }
    }

    // ── cleanup ──────────────────────────────────────────────────────────────
    let _ = std::fs::remove_file(&fasta_path);
    let _ = std::fs::remove_file(&fastq_path);
    let _ = std::fs::remove_file(&bam_path);
}

/// Guard against double-correction: if the FASTQ reader were ever changed to
/// subtract 33 itself, `fastq_quals_to_phred` would produce values ≤ 60 even
/// for Q40, but specifically 7 (i.e. 40 - 33) rather than 40.  This test
/// confirms the BAM does NOT store 7.
#[test]
fn bam_qual_field_is_not_double_decoded() {
    let reference_seq = "ACGTACGTACGTACGT";
    let fasta_path = tmp("ref2", "fa");
    write_text(&fasta_path, &format!(">chr1\n{reference_seq}\n"));

    // Q40 reads; raw Phred 40, ASCII 'I' (73).
    let fastq_path = tmp("reads2", "fastq");
    write_text(&fastq_path, "@read1\nACGTACGT\n+\nIIIIIIII\n");

    let bam_path = tmp("out2", "bam");
    let bin = env!("CARGO_BIN_EXE_rosalind");
    let status = Command::new(bin)
        .args([
            "align",
            "--reference",
            fasta_path.to_str().unwrap(),
            "--reads",
            fastq_path.to_str().unwrap(),
            "--format",
            "bam",
            "--output",
            bam_path.to_str().unwrap(),
        ])
        .status()
        .expect("failed to run rosalind align");
    assert!(status.success());

    let mut contigs = ContigSet::new();
    contigs.push("chr1", reference_seq.len() as u32);
    let reads = read_bam_as_core_reads(&bam_path, &contigs).unwrap();

    for read in &reads {
        for (i, &qual) in read.qual.iter().enumerate() {
            assert_ne!(
                qual, 7,
                "BAM quality byte[{i}] = 7 suggests double-decoding (40 - 33 = 7); \
                 each byte should be exactly 40"
            );
            // Not ASCII either.
            assert_ne!(
                qual, 73,
                "BAM quality byte[{i}] = 73 means ASCII 'I' was stored verbatim (the original bug)"
            );
        }
    }

    let _ = std::fs::remove_file(&fasta_path);
    let _ = std::fs::remove_file(&fastq_path);
    let _ = std::fs::remove_file(&bam_path);
}
