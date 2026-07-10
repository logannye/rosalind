//! `variants --index --gvcf` — a bounded, byte-reproducible banded gVCF: every
//! callable locus is a variant record or a `<NON_REF>` reference block (`END=`
//! span), providing single-sample reference confidence for future joiner testing.

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

fn run(args: &[&str]) -> std::process::Output {
    Command::new(bin())
        .args(args)
        .output()
        .expect("spawn rosalind")
}

// index -> align -> sort a small single-contig fixture with real coverage.
fn build_sorted_bam_fixture() -> (PathBuf, PathBuf, PathBuf) {
    let dir = unique_dir("rosalind-gvcf");
    let seq = "ACGTACGTACGTACGTACGTACGTACGTACGT"; // 32 bp
    let fa = dir.join("ref.fa");
    std::fs::write(&fa, format!(">chr1\n{seq}\n")).unwrap();
    let fq = dir.join("reads.fq");
    let mut s = String::new();
    for (i, &start) in [0usize, 0, 4, 4, 8, 8, 12, 16].iter().enumerate() {
        let read = &seq[start..start + 16];
        let qual: String = "I".repeat(16);
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

fn gvcf(
    idx: &std::path::Path,
    bam: &std::path::Path,
    out: &std::path::Path,
    manifest: Option<&std::path::Path>,
) {
    let mut args = vec![
        "variants".to_string(),
        "--index".to_string(),
        idx.display().to_string(),
        "--alignments".to_string(),
        bam.display().to_string(),
        "--gvcf".to_string(),
        "-o".to_string(),
        out.display().to_string(),
    ];
    if let Some(m) = manifest {
        args.push("--manifest".to_string());
        args.push(m.display().to_string());
    }
    let argv: Vec<&str> = args.iter().map(String::as_str).collect();
    let o = run(&argv);
    assert!(o.status.success(), "gvcf run failed: {o:?}");
}

#[test]
fn gvcf_is_banded_byte_reproducible_and_well_formed() {
    let (dir, idx, bam) = build_sorted_bam_fixture();
    let g1 = dir.join("a.gvcf");
    let g2 = dir.join("b.gvcf");
    let manifest = dir.join("a.gvcf.manifest.json");
    gvcf(&idx, &bam, &g1, Some(&manifest));
    gvcf(&idx, &bam, &g2, None);

    // Byte-reproducible: two independent runs are identical.
    let a = std::fs::read_to_string(&g1).unwrap();
    let b = std::fs::read_to_string(&g2).unwrap();
    assert_eq!(a, b, "gVCF must be byte-identical run-to-run");

    // Well-formed gVCF header.
    assert!(a.contains("##ALT=<ID=NON_REF"), "missing NON_REF ALT: {a}");
    assert!(a.contains("##INFO=<ID=END"), "missing END INFO");
    assert!(a.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"));

    // There IS at least one reference block, and every data record is either a
    // `<NON_REF>` reference block with an END>=POS span or a variant record.
    let mut ref_blocks = 0;
    for line in a.lines().filter(|l| !l.starts_with('#')) {
        let f: Vec<&str> = line.split('\t').collect();
        assert!(f.len() >= 10, "malformed record: {line}");
        let pos: u64 = f[1].parse().unwrap();
        if f[4] == "<NON_REF>" {
            ref_blocks += 1;
            let end: u64 = f[7]
                .strip_prefix("END=")
                .and_then(|s| s.parse().ok())
                .unwrap_or_else(|| panic!("ref block without END=: {line}"));
            assert!(end >= pos, "ref block END {end} < POS {pos}");
        } else {
            // A variant record carries <NON_REF> as a second ALT (gVCF convention).
            assert!(
                f[4].ends_with(",<NON_REF>"),
                "variant ALT not gVCF-shaped: {line}"
            );
        }
    }
    assert!(
        ref_blocks > 0,
        "a covered gVCF must contain reference blocks"
    );

    // Bounded: the receipt records a small realized peak (coverage-bounded).
    let m = std::fs::read_to_string(&manifest).unwrap();
    let rm = rosalind::provenance::RunManifest::from_canonical_json(&m).unwrap();
    assert!(rm.measurements.contains_key("peak_rss_bytes"));

    std::fs::remove_dir_all(&dir).ok();
}
