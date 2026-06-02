//! `eval-germline` must normalize each variant against its OWN contig. Loading
//! only the first FASTA record silently miscompared (or crashed with an opaque
//! out-of-bounds error) on any multi-contig benchmark — the headline truth-
//! comparison surface, and exactly what a real GIAB run looks like.

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

// chr1 = 20 bp, chr2 = 40 bp.
const REFERENCE: &str =
    ">chr1\nACGTACGTACGTACGTACGT\n>chr2\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";

#[test]
fn eval_germline_compares_each_variant_against_its_own_contig() {
    let dir = unique_dir("rosalind-eval-mc");
    let fa = dir.join("ref.fa");
    std::fs::write(&fa, REFERENCE).unwrap();

    // A chr1 SNV (pos 5) and a chr2 SNV at pos 30 — past chr1's length (20), so
    // the old first-contig-only path crashed out-of-bounds on the chr2 variant.
    let header = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    let calls = dir.join("calls.vcf");
    std::fs::write(
        &calls,
        format!("{header}chr1\t5\t.\tA\tC\t50\tPASS\t.\nchr2\t30\t.\tG\tT\t50\tPASS\t.\n"),
    )
    .unwrap();
    let truth = dir.join("truth.vcf");
    std::fs::write(
        &truth,
        format!("{header}chr1\t5\t.\tA\tC\t50\tPASS\t.\nchr2\t30\t.\tG\tT\t50\tPASS\t.\n"),
    )
    .unwrap();

    let out = Command::new(bin())
        .args(["eval-germline", "--reference"])
        .arg(&fa)
        .arg("--calls")
        .arg(&calls)
        .arg("--truth")
        .arg(&truth)
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "multi-contig eval must succeed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8_lossy(&out.stdout);
    // Both contigs' variants match truth → 2 TP, 0 FP, 0 FN.
    assert!(
        stdout.contains("tp=2"),
        "expected 2 true positives: {stdout}"
    );
    assert!(
        stdout.contains("fp=0"),
        "expected 0 false positives: {stdout}"
    );
    assert!(
        stdout.contains("fn=0"),
        "expected 0 false negatives: {stdout}"
    );
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn eval_germline_errors_clearly_on_a_contig_naming_mismatch() {
    let dir = unique_dir("rosalind-eval-naming");
    let fa = dir.join("ref.fa");
    std::fs::write(&fa, REFERENCE).unwrap();

    // A variant on a contig the FASTA does not contain (Ensembl '1' vs UCSC
    // 'chr1') must error clearly, not silently miscompare.
    let header = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    let calls = dir.join("calls.vcf");
    std::fs::write(&calls, format!("{header}1\t5\t.\tA\tC\t50\tPASS\t.\n")).unwrap();
    let truth = dir.join("truth.vcf");
    std::fs::write(&truth, header).unwrap();

    let out = Command::new(bin())
        .args(["eval-germline", "--reference"])
        .arg(&fa)
        .arg("--calls")
        .arg(&calls)
        .arg("--truth")
        .arg(&truth)
        .output()
        .unwrap();
    assert!(!out.status.success(), "naming mismatch must fail");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("naming scheme") || stderr.contains("no matching sequence"),
        "expected a clear naming-mismatch error: {stderr}"
    );
    std::fs::remove_dir_all(&dir).ok();
}
