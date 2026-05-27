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

#[test]
fn locate_matches_in_ram_ground_truth() {
    use rosalind::genomics::GenomeIndex;

    let dir = tmpdir();
    let fa = write_fasta(&dir);
    let idx = dir.join("ref.idx");
    let build = Command::new(bin())
        .args([
            "index",
            "--reference",
            fa.to_str().unwrap(),
            "--output",
            idx.to_str().unwrap(),
        ])
        .output()
        .expect("index");
    assert!(build.status.success());

    // In-RAM ground truth over the same sequences.
    let gi = GenomeIndex::from_named_sequences(&[
        (
            "chr1".to_string(),
            b"ACGTACGTNNACGTACGTACGTAAGGCCTT".to_vec(),
        ),
        (
            "chr2".to_string(),
            b"TTTTGGGGCCCCAAAANNNNACGTACGTAC".to_vec(),
        ),
        (
            "chr3".to_string(),
            b"GATTACATTTTGATTACAGGGGGCCCCAAA".to_vec(),
        ),
    ])
    .unwrap();

    // Battery: multi-hit, lowercase (case-insensitive), cross-contig boundary
    // straddle ("CTTTTTT" spans chr1's tail into chr2 -> rejected, 0 hits),
    // N-bearing, and an invalid base ("ZZZZ" -> 0 hits).
    for pat in [
        "GATTACA", "gattaca", "ACGT", "GGGGG", "NNNN", "TTTTGGGG", "CTTTTTT", "ZZZZ",
    ] {
        let out = Command::new(bin())
            .args(["locate", "--index", idx.to_str().unwrap(), "--pattern", pat])
            .output()
            .expect("locate");
        assert!(
            out.status.success(),
            "locate {pat} failed: {}",
            String::from_utf8_lossy(&out.stderr)
        );
        let stdout = String::from_utf8_lossy(&out.stdout);

        let mut expected: Vec<String> = gi
            .locate_exact(pat.as_bytes(), 1024)
            .into_iter()
            .map(|l| {
                let name = gi.contigs().by_id(l.contig).unwrap().name.to_string();
                format!("{name}\t{}", l.pos.0)
            })
            .collect();
        expected.sort();
        let mut got: Vec<String> = stdout.lines().map(|s| s.to_string()).collect();
        got.sort();
        assert_eq!(got, expected, "locate mismatch for {pat}");
    }
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn index_build_is_deterministic_via_cli() {
    let dir = tmpdir();
    let fa = write_fasta(&dir);
    let idx1 = dir.join("a.idx");
    let idx2 = dir.join("b.idx");
    for out in [&idx1, &idx2] {
        let r = Command::new(bin())
            .args([
                "index",
                "--reference",
                fa.to_str().unwrap(),
                "--output",
                out.to_str().unwrap(),
            ])
            .output()
            .expect("index");
        assert!(r.status.success());
    }
    assert_eq!(
        std::fs::read(&idx1).unwrap(),
        std::fs::read(&idx2).unwrap(),
        "two CLI builds of the same reference must be byte-identical"
    );
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn locate_works_from_the_artifact_alone() {
    // Build, delete the source FASTA, then locate from the index file alone —
    // proving the load path is self-contained and never rebuilds.
    let dir = tmpdir();
    let fa = write_fasta(&dir);
    let idx = dir.join("ref.idx");
    assert!(Command::new(bin())
        .args([
            "index",
            "--reference",
            fa.to_str().unwrap(),
            "--output",
            idx.to_str().unwrap(),
        ])
        .output()
        .expect("index")
        .status
        .success());
    std::fs::remove_file(&fa).unwrap(); // the only source of the sequence is now gone

    let out = Command::new(bin())
        .args([
            "locate",
            "--index",
            idx.to_str().unwrap(),
            "--pattern",
            "GATTACA",
        ])
        .output()
        .expect("locate");
    assert!(out.status.success());
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert_eq!(
        stdout.lines().count(),
        2,
        "GATTACA: 2 hits in chr3, served from the artifact alone"
    );
    std::fs::remove_dir_all(&dir).ok();
}
