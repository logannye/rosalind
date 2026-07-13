use std::path::PathBuf;
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

use rosalind::genomics::{ReferencePackReader, ReferenceProvider, ReferenceSequence};

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_rosalind")
}

fn tempdir(name: &str) -> PathBuf {
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let path = std::env::temp_dir().join(format!("rosalind-reference-cli-{name}-{nonce}"));
    std::fs::create_dir_all(&path).unwrap();
    path
}

#[test]
fn build_and_inspect_json_are_deterministic() {
    let root = tempdir("build");
    let fasta = root.join("reference with spaces.fa");
    let pack = root.join("reference with spaces.rref");
    std::fs::write(&fasta, b">chr1\nACGTN\n>chr2 unusual\nttuu\n").unwrap();

    let build = Command::new(bin())
        .args(["reference", "build", "--fasta"])
        .arg(&fasta)
        .arg("--output")
        .arg(&pack)
        .output()
        .unwrap();
    assert!(
        build.status.success(),
        "{}",
        String::from_utf8_lossy(&build.stderr)
    );

    let inspect = || {
        Command::new(bin())
            .args(["reference", "inspect", "--reference-pack"])
            .arg(&pack)
            .arg("--json")
            .output()
            .unwrap()
    };
    let first = inspect();
    let second = inspect();
    assert!(first.status.success());
    assert_eq!(first.stdout, second.stdout);
    let json = String::from_utf8(first.stdout).unwrap();
    assert!(json.contains("\"format\":\"rref\""));
    assert!(json.contains("\"total_bases\":9"));

    let reference = ReferencePackReader::open(&pack).unwrap();
    assert_eq!(reference.contigs().len(), 2);
    let mut sequence = Vec::new();
    reference.decode_window(0, reference.len(), &mut sequence);
    assert_eq!(sequence, b"ACGTNTTTT");
    std::fs::remove_dir_all(root).unwrap();
}

#[test]
fn legacy_index_conversion_preserves_reference_identity_and_bytes() {
    let root = tempdir("convert");
    let fasta = root.join("reference.fa");
    let index = root.join("reference.idx");
    let pack = root.join("reference.rref");
    std::fs::write(&fasta, b">chr1\nACGTNNACGT\n").unwrap();

    let built = Command::new(bin())
        .args(["index", "--reference"])
        .arg(&fasta)
        .arg("--output")
        .arg(&index)
        .output()
        .unwrap();
    assert!(built.status.success());
    let converted = Command::new(bin())
        .args(["reference", "convert", "--index"])
        .arg(&index)
        .arg("--output")
        .arg(&pack)
        .output()
        .unwrap();
    assert!(
        converted.status.success(),
        "{}",
        String::from_utf8_lossy(&converted.stderr)
    );

    let legacy = rosalind::IndexReader::open(&index).unwrap();
    let converted = ReferencePackReader::open(&pack).unwrap();
    assert_eq!(
        legacy.header.reference_blake3,
        converted.source_reference_blake3()
    );
    let mut legacy_bytes = Vec::new();
    legacy
        .reference_view()
        .unwrap()
        .decode_window(0, converted.len(), &mut legacy_bytes);
    let mut pack_bytes = Vec::new();
    converted.decode_window(0, converted.len(), &mut pack_bytes);
    assert_eq!(legacy_bytes, pack_bytes);
    std::fs::remove_dir_all(root).unwrap();
}

#[test]
fn reference_pack_and_legacy_index_produce_identical_feature_bytes() {
    let root = tempdir("analysis-equivalence");
    let fixture = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("examples/data/illumina_toy");
    let fasta = fixture.join("reference.fa");
    let bam = fixture.join("alignments.bam");
    let index = root.join("reference.idx");
    let pack = root.join("reference.rref");
    let sorted = root.join("sorted.bam");
    let legacy_output = root.join("legacy.tsv");
    let pack_output = root.join("pack.tsv");

    let setup_commands: Vec<Vec<std::ffi::OsString>> = vec![
        vec![
            "index".into(),
            "--reference".into(),
            fasta.as_os_str().into(),
            "--output".into(),
            index.as_os_str().into(),
        ],
        vec![
            "reference".into(),
            "build".into(),
            "--fasta".into(),
            fasta.as_os_str().into(),
            "--output".into(),
            pack.as_os_str().into(),
        ],
        vec![
            "sort".into(),
            "--input".into(),
            bam.as_os_str().into(),
            "--output".into(),
            sorted.as_os_str().into(),
        ],
    ];
    for arguments in setup_commands {
        let output = Command::new(bin()).args(arguments).output().unwrap();
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
    }

    let legacy = Command::new(bin())
        .args(["features", "--index"])
        .arg(&index)
        .arg("--alignments")
        .arg(&sorted)
        .arg("--output")
        .arg(&legacy_output)
        .output()
        .unwrap();
    assert!(
        legacy.status.success(),
        "{}",
        String::from_utf8_lossy(&legacy.stderr)
    );
    let packed = Command::new(bin())
        .args(["features", "--reference-pack"])
        .arg(&pack)
        .arg("--alignments")
        .arg(&sorted)
        .arg("--output")
        .arg(&pack_output)
        .output()
        .unwrap();
    assert!(
        packed.status.success(),
        "{}",
        String::from_utf8_lossy(&packed.stderr)
    );
    assert_eq!(
        std::fs::read(&legacy_output).unwrap(),
        std::fs::read(&pack_output).unwrap()
    );
    let receipt = std::fs::read_to_string(pack_output.with_extension("tsv.manifest.json")).unwrap();
    assert!(receipt.contains("--reference-pack"));
    std::fs::remove_dir_all(root).unwrap();
}

#[test]
fn reference_pack_receipt_verifies_and_replays() {
    let work = tempdir("receipt-replay");
    let fasta = work.join("reference.fa");
    let pack = work.join("reference.rref");
    std::fs::write(&fasta, ">chr1\nACGTNNACGT\n").unwrap();
    let build = Command::new(bin())
        .args([
            "reference",
            "build",
            "--fasta",
            fasta.to_str().unwrap(),
            "--output",
            pack.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(build.success());
    let receipt = PathBuf::from(format!("{}.manifest.json", pack.display()));
    let verify = Command::new(bin())
        .args(["verify", "--manifest", receipt.to_str().unwrap()])
        .status()
        .unwrap();
    assert!(verify.success());
    let replay = Command::new(bin())
        .args([
            "reproduce",
            "--manifest",
            receipt.to_str().unwrap(),
            "--inputs",
            work.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        replay.status.success(),
        "{}",
        String::from_utf8_lossy(&replay.stderr)
    );
    assert!(String::from_utf8_lossy(&replay.stdout).contains("REPRODUCED"));
    std::fs::remove_dir_all(work).unwrap();
}
