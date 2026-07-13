use rust_htslib::bam::Read;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_rosalind")
}

fn tempdir() -> PathBuf {
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let path = std::env::temp_dir().join(format!("rosalind-partition-{nonce}"));
    std::fs::create_dir_all(&path).unwrap();
    path
}

fn run(arguments: &[&str], paths: &[&Path]) -> std::process::Output {
    let mut command = Command::new(bin());
    command.args(arguments);
    for path in paths {
        command.arg(path);
    }
    let output = command.output().unwrap();
    assert!(
        output.status.success(),
        "command failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

#[test]
fn arrow_shards_merge_byte_identically_to_whole_genome() {
    let root = tempdir();
    let fixture = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("examples/data/illumina_toy");
    let fasta = fixture.join("reference.fa");
    let bam = fixture.join("alignments.bam");
    let reference = root.join("reference.rref");
    let sorted_all = root.join("sorted-all.bam");
    let sorted = root.join("sorted.bam");
    run(
        &["reference", "build", "--fasta"],
        &[&fasta, Path::new("--output"), &reference],
    );
    run(
        &["sort", "--input"],
        &[&bam, Path::new("--output"), &sorted_all],
    );
    let mut reader = rust_htslib::bam::Reader::from_path(&sorted_all).unwrap();
    let header = rust_htslib::bam::Header::from_template(reader.header());
    let mut writer =
        rust_htslib::bam::Writer::from_path(&sorted, &header, rust_htslib::bam::Format::Bam)
            .unwrap();
    for record in reader.records() {
        let record = record.unwrap();
        if !record.is_unmapped() {
            writer.write(&record).unwrap();
        }
    }
    drop(writer);
    rust_htslib::bam::index::build(
        &sorted,
        None::<&PathBuf>,
        rust_htslib::bam::index::Type::Bai,
        1,
    )
    .unwrap();

    let whole = root.join("whole.arrow");
    let output = Command::new(bin())
        .args(["features", "--reference-pack"])
        .arg(&reference)
        .arg("--alignments")
        .arg(&sorted)
        .args(["--format", "arrow-ipc", "--output"])
        .arg(&whole)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );

    let mut manifests = Vec::new();
    for index in 0..3 {
        let shard = root.join(format!("shard-{index}.arrow"));
        let output = Command::new(bin())
            .args(["features", "--reference-pack"])
            .arg(&reference)
            .arg("--alignments")
            .arg(&sorted)
            .args([
                "--format",
                "arrow-ipc",
                "--shard-count",
                "3",
                "--shard-index",
            ])
            .arg(index.to_string())
            .arg("--output")
            .arg(&shard)
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        manifests.push(PathBuf::from(format!("{}.manifest.json", shard.display())));
    }
    let merged = root.join("merged.arrow");
    let mut command = Command::new(bin());
    command.arg("merge");
    for manifest in &manifests {
        command.arg("--manifest").arg(manifest);
    }
    let output = command
        .arg("--inputs")
        .arg(&root)
        .arg("--output")
        .arg(&merged)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        std::fs::read(whole).unwrap(),
        std::fs::read(merged).unwrap()
    );

    let collision_output = root.join("collision.arrow");
    let collision_receipt = root.join("existing-merge.manifest.json");
    std::fs::write(&collision_receipt, "keep me").unwrap();
    let mut collision = Command::new(bin());
    collision.arg("merge");
    for manifest in &manifests {
        collision.arg("--manifest").arg(manifest);
    }
    let collision = collision
        .arg("--inputs")
        .arg(&root)
        .arg("--output")
        .arg(&collision_output)
        .arg("--output-manifest")
        .arg(&collision_receipt)
        .status()
        .unwrap();
    assert_eq!(collision.code(), Some(2));
    assert!(!collision_output.exists());
    assert_eq!(
        std::fs::read_to_string(collision_receipt).unwrap(),
        "keep me"
    );

    let missing = Command::new(bin())
        .arg("merge")
        .arg("--manifest")
        .arg(&manifests[0])
        .arg("--manifest")
        .arg(&manifests[1])
        .arg("--inputs")
        .arg(&root)
        .arg("--output")
        .arg(root.join("missing.arrow"))
        .status()
        .unwrap();
    assert_eq!(missing.code(), Some(3));

    let tampered = root.join("tampered.manifest.json");
    let original = std::fs::read_to_string(&manifests[2]).unwrap();
    let marker = "\"manifest_blake3\":\"";
    let offset = original.find(marker).unwrap() + marker.len();
    let mut bytes = original.into_bytes();
    bytes[offset] = if bytes[offset] == b'0' { b'1' } else { b'0' };
    std::fs::write(&tampered, bytes).unwrap();
    let tampered_result = Command::new(bin())
        .arg("merge")
        .arg("--manifest")
        .arg(&manifests[0])
        .arg("--manifest")
        .arg(&manifests[1])
        .arg("--manifest")
        .arg(&tampered)
        .arg("--inputs")
        .arg(&root)
        .arg("--output")
        .arg(root.join("tampered.arrow"))
        .status()
        .unwrap();
    assert_eq!(tampered_result.code(), Some(5));

    let empty_bed = root.join("empty.bed");
    std::fs::write(&empty_bed, "# deliberately empty\nchrToy\t10\t10\n").unwrap();
    for (command, format, extension) in [
        ("features", Some("tsv"), "tsv"),
        ("features", Some("arrow-ipc"), "arrow"),
        ("variants", None, "vcf"),
    ] {
        let destination = root.join(format!("empty.{extension}"));
        let mut invocation = Command::new(bin());
        invocation
            .arg(command)
            .arg("--reference-pack")
            .arg(&reference)
            .arg("--alignments")
            .arg(&sorted)
            .arg("--regions")
            .arg(&empty_bed);
        if let Some(format) = format {
            invocation.arg("--format").arg(format);
        }
        let result = invocation
            .arg("--output")
            .arg(&destination)
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        assert!(destination.exists());
        assert!(PathBuf::from(format!("{}.manifest.json", destination.display())).exists());
    }

    for (label, prefix, extension) in [
        ("features-tsv", vec!["features"], "tsv"),
        ("coverage", vec!["analyze", "coverage"], "tsv"),
        ("sites", vec!["variants"], "vcf"),
        ("gvcf", vec!["variants", "--gvcf"], "g.vcf"),
    ] {
        let whole = root.join(format!("{label}-whole.{extension}"));
        let mut invocation = Command::new(bin());
        let result = invocation
            .args(&prefix)
            .arg("--reference-pack")
            .arg(&reference)
            .arg("--alignments")
            .arg(&sorted)
            .arg("--output")
            .arg(&whole)
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{label} whole: {}",
            String::from_utf8_lossy(&result.stderr)
        );
        let mut codec_manifests = Vec::new();
        for index in 0..3 {
            let shard = root.join(format!("{label}-shard-{index}.{extension}"));
            let mut invocation = Command::new(bin());
            let result = invocation
                .args(&prefix)
                .arg("--reference-pack")
                .arg(&reference)
                .arg("--alignments")
                .arg(&sorted)
                .args(["--shard-count", "3", "--shard-index"])
                .arg(index.to_string())
                .arg("--output")
                .arg(&shard)
                .output()
                .unwrap();
            assert!(
                result.status.success(),
                "{label} shard {index}: {}",
                String::from_utf8_lossy(&result.stderr)
            );
            codec_manifests.push(PathBuf::from(format!("{}.manifest.json", shard.display())));
        }
        let merged = root.join(format!("{label}-merged.{extension}"));
        let mut merge = Command::new(bin());
        merge.arg("merge");
        for manifest in codec_manifests {
            merge.arg("--manifest").arg(manifest);
        }
        let result = merge
            .arg("--inputs")
            .arg(&root)
            .arg("--output")
            .arg(&merged)
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{label} merge: {}",
            String::from_utf8_lossy(&result.stderr)
        );
        assert_eq!(
            std::fs::read(&whole).unwrap(),
            std::fs::read(&merged).unwrap(),
            "{label} canonical merge must equal the unsharded bytes"
        );
    }
    std::fs::remove_dir_all(root).unwrap();
}

#[test]
fn sparse_selection_without_bai_refuses_before_output() {
    let root = tempdir();
    let fixture = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("examples/data/illumina_toy");
    let reference = root.join("reference.rref");
    let output = root.join("must-not-exist.tsv");
    let build = Command::new(bin())
        .args(["reference", "build", "--fasta"])
        .arg(fixture.join("reference.fa"))
        .arg("--output")
        .arg(&reference)
        .output()
        .unwrap();
    assert!(build.status.success());
    let result = Command::new(bin())
        .args(["features", "--reference-pack"])
        .arg(&reference)
        .arg("--alignments")
        .arg(fixture.join("alignments.bam"))
        .args(["--region", "chrToy:1-10", "--output"])
        .arg(&output)
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert!(!output.exists());
    std::fs::remove_dir_all(root).unwrap();
}
