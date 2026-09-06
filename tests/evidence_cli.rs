use arrow_ipc::reader::StreamReader;
use rosalind::provenance::RunManifest;
use rust_htslib::bam::{
    self,
    header::HeaderRecord,
    record::{Cigar, CigarString},
    Header, Record,
};
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

struct Fixture {
    root: PathBuf,
    reference: PathBuf,
    bam: PathBuf,
    bed: PathBuf,
}
impl Fixture {
    fn new() -> Self {
        static NEXT: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
        let root = std::env::temp_dir().join(format!(
            "rosalind-evidence-cli-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let root = root.with_extension(
            NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed)
                .to_string(),
        );
        std::fs::create_dir(&root).unwrap();
        let reference = root.join("ref.fa");
        std::fs::write(&reference, format!(">chr1\n{}\n", "A".repeat(4099))).unwrap();
        std::fs::write(root.join("ref.fa.fai"), "chr1\t4099\t6\t4099\t4100\n").unwrap();
        let bam = root.join("reads.bam");
        let mut header = Header::new();
        header.push_record(
            HeaderRecord::new(b"HD")
                .push_tag(b"VN", "1.6")
                .push_tag(b"SO", "coordinate"),
        );
        header.push_record(
            HeaderRecord::new(b"SQ")
                .push_tag(b"SN", "chr1")
                .push_tag(b"LN", 4099),
        );
        let mut writer = bam::Writer::from_path(&bam, &header, bam::Format::Bam).unwrap();
        for (i, (pos, length, flags, mq, bq)) in [
            (0, 2, 0, 60, 40),
            (0, 8, 16, 60, 40),
            (0, 8, 1024, 60, 40),
            (0, 8, 0, 255, 40),
            (0, 8, 0, 60, 255),
            (4090, 9, 0, 60, 40),
        ]
        .into_iter()
        .enumerate()
        {
            let mut record = Record::new();
            record.set(
                format!("r{i}").as_bytes(),
                Some(&CigarString(vec![Cigar::Match(length)])),
                &vec![b'A'; length as usize],
                &vec![bq; length as usize],
            );
            record.set_tid(0);
            record.set_pos(pos);
            record.set_flags(flags);
            record.set_mapq(mq);
            writer.write(&record).unwrap();
        }
        drop(writer);
        bam::index::build(&bam, None::<&PathBuf>, bam::index::Type::Bai, 1).unwrap();
        let bed = root.join("panel.bed");
        std::fs::write(&bed, "chr1\t0\t4099\twhole\nchr1\t2\t8\toverlap\n").unwrap();
        Self {
            root,
            reference,
            bam,
            bed,
        }
    }
    fn command(&self, kind: &str) -> Command {
        let mut cmd = Command::new(env!("CARGO_BIN_EXE_rosalind"));
        cmd.args(["analyze", kind, "--reference"])
            .arg(&self.reference)
            .arg("--alignments")
            .arg(&self.bam)
            .arg("--regions")
            .arg(&self.bed);
        cmd
    }
}
impl Drop for Fixture {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.root);
    }
}
fn success(output: Output) -> Output {
    assert!(
        output.status.success(),
        "exit {:?}; stdout: {}; stderr: {}",
        output.status.code(),
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    output
}
fn receipt(path: &Path) -> RunManifest {
    RunManifest::from_canonical_json(
        &std::fs::read_to_string(format!("{}.manifest.json", path.display())).unwrap(),
    )
    .unwrap()
}

#[test]
fn evidence_cli_arrow_is_invariant_and_receipts_replay() {
    let f = Fixture::new();
    let mut outputs = Vec::new();
    let mut identities = Vec::new();
    for (i, (budget, tile)) in [(64, 1), (96, 137), (128, 16384)].into_iter().enumerate() {
        let path = f.root.join(format!("run-{i}.arrow"));
        success(
            f.command("evidence")
                .args([
                    "--format",
                    "arrow-ipc",
                    "--memory-budget-mb",
                    &budget.to_string(),
                    "--tile-bases",
                    &tile.to_string(),
                    "-o",
                ])
                .arg(&path)
                .output()
                .unwrap(),
        );
        let manifest = receipt(&path);
        assert_eq!(manifest.self_hash_ok(), Some(true));
        assert_eq!(manifest.measurement_hash_ok(), Some(true));
        assert_eq!(manifest.params["run_status"], "completed");
        assert_eq!(
            manifest.params["evidence.profile"],
            "shortread-dna-readcount-v1"
        );
        identities.push(manifest.params["evidence.science_blake3"].clone());
        outputs.push(std::fs::read(path).unwrap());
    }
    assert!(outputs.windows(2).all(|w| w[0] == w[1]));
    assert!(identities.windows(2).all(|w| w[0] == w[1]));
    let batches: Vec<_> = StreamReader::try_new(outputs[0].as_slice(), None)
        .unwrap()
        .map(|b| b.unwrap().num_rows())
        .collect();
    assert_eq!(batches, [1024, 1024, 1024, 1024, 3]);
    let manifest = f.root.join("run-2.arrow.manifest.json");
    success(
        Command::new(env!("CARGO_BIN_EXE_rosalind"))
            .args(["verify", "--manifest"])
            .arg(&manifest)
            .output()
            .unwrap(),
    );
    success(
        Command::new(env!("CARGO_BIN_EXE_rosalind"))
            .args(["reproduce", "--manifest"])
            .arg(&manifest)
            .arg("--inputs")
            .arg(&f.root)
            .arg("--no-attest")
            .output()
            .unwrap(),
    );
}

#[test]
fn panel_fuses_position_evidence_and_uses_full_target_denominators() {
    let f = Fixture::new();
    let panel = f.root.join("z-panel.tsv");
    let positions = f.root.join("a-positions.arrow");
    success(
        f.command("panel-qc")
            .args(["--min-callable-depth", "2", "-o"])
            .arg(&panel)
            .arg("--position-output")
            .arg(&positions)
            .output()
            .unwrap(),
    );
    let text = std::fs::read_to_string(&panel).unwrap();
    let mut lines = text.lines();
    let header: Vec<_> = lines.next().unwrap().split('\t').collect();
    let index = |name| header.iter().position(|field| *field == name).unwrap();
    let whole: Vec<_> = lines.next().unwrap().split('\t').collect();
    assert_eq!(whole[index("length")], "4099");
    assert_eq!(whole[index("callable_depth_sum")], "19");
    assert_eq!(whole[index("callable_positions")], "2");
    assert_eq!(whole[index("min_callable_depth")], "0");
    let overlap: Vec<_> = lines.next().unwrap().split('\t').collect();
    assert_eq!(overlap[index("length")], "6");
    assert_eq!(overlap[index("callable_depth_sum")], "6");
    let standalone = f.root.join("standalone.arrow");
    success(
        f.command("evidence")
            .args(["--format", "arrow-ipc", "-o"])
            .arg(&standalone)
            .output()
            .unwrap(),
    );
    assert_eq!(
        std::fs::read(positions).unwrap(),
        std::fs::read(standalone).unwrap()
    );
    assert_eq!(receipt(&panel).outputs.len(), 2);
    let manifest = receipt(&panel);
    assert!(manifest.outputs[0].path.ends_with("a-positions.arrow"));
    assert_eq!(manifest.params["artifact.output.0.role"], "exact-evidence");
    assert_eq!(manifest.params["artifact.output.0.format"], "arrow-ipc");
    assert_eq!(manifest.params["artifact.output.1.role"], "panel-summary");
    assert_eq!(manifest.params["artifact.output.1.format"], "tsv");
}

#[test]
fn evidence_resource_and_input_failures_never_publish_success() {
    let f = Fixture::new();
    let path = f.root.join("refused.arrow");
    let output = f
        .command("evidence")
        .args(["--memory-budget-mb", "20", "--format", "arrow-ipc", "-o"])
        .arg(&path)
        .output()
        .unwrap();
    assert_eq!(
        output.status.code(),
        Some(3),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(!path.exists());
    let path = f.root.join("capacity.arrow");
    let output = f
        .command("evidence")
        .args(["--max-read-len", "4", "--format", "arrow-ipc", "-o"])
        .arg(&path)
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(4));
    assert!(!path.exists());
    assert!(f.root.join("capacity.arrow.partial").exists());
    let manifest = receipt(&path);
    assert_eq!(manifest.params["failure.kind"], "declared-capacity");
    assert_eq!(manifest.self_hash_ok(), Some(true));
    let output = f
        .command("evidence")
        .arg("--sites")
        .arg(f.root.join("unused.vcf"))
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(2));
}

#[test]
fn arrow_legacy_plan_includes_encoder_flush_transients() {
    let f = Fixture::new();
    let reference = f.root.join("ref.rref");
    success(
        Command::new(env!("CARGO_BIN_EXE_rosalind"))
            .args(["reference", "build", "--fasta"])
            .arg(&f.reference)
            .arg("-o")
            .arg(&reference)
            .output()
            .unwrap(),
    );
    let plan = success(
        Command::new(env!("CARGO_BIN_EXE_rosalind"))
            .args(["plan", "--reference-pack"])
            .arg(reference)
            .args([
                "--max-depth",
                "1",
                "--format",
                "arrow-ipc",
                "--budget-mb",
                "32",
                "--json",
            ])
            .output()
            .unwrap(),
    );
    let text = String::from_utf8(plan.stdout).unwrap();
    assert!(text.contains("\"encoder_additional_bytes\":"));
    assert!(text.contains("\"verdict\":\"refuse\""), "{text}");
}

#[test]
fn cache_reuse_and_worker_cli_preserve_canonical_output() {
    let f = Fixture::new();
    let cache = f.root.join("cache");
    let first = f.root.join("first.arrow");
    success(
        f.command("evidence")
            .args(["--format", "arrow-ipc", "--workers", "2", "--cache-dir"])
            .arg(&cache)
            .arg("-o")
            .arg(&first)
            .output()
            .unwrap(),
    );
    assert_eq!(
        receipt(&first).measurements["execution.computed_partitions"],
        "1"
    );
    let second = f.root.join("second.arrow");
    success(
        f.command("evidence")
            .args([
                "--format",
                "arrow-ipc",
                "--workers",
                "1",
                "--tile-bases",
                "17",
                "--cache-dir",
            ])
            .arg(&cache)
            .arg("-o")
            .arg(&second)
            .output()
            .unwrap(),
    );
    let manifest = receipt(&second);
    assert_eq!(manifest.measurements["execution.reused_partitions"], "1");
    assert_eq!(manifest.measurements["execution.record_visits"], "0");
    assert_eq!(
        std::fs::read(first).unwrap(),
        std::fs::read(second).unwrap()
    );
}

#[test]
fn finalization_breach_demotes_both_panel_artifacts() {
    let f = Fixture::new();
    let panel = f.root.join("final-panel.tsv");
    let positions = f.root.join("final-positions.arrow");
    let output = f
        .command("panel-qc")
        .args(["--memory-budget-mb", "128", "-o"])
        .arg(&panel)
        .arg("--position-output")
        .arg(&positions)
        .env("ROSALIND_FORCE_FINAL_RSS_BYTES", (256u64 << 20).to_string())
        .output()
        .unwrap();
    assert_eq!(
        output.status.code(),
        Some(4),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(!panel.exists());
    assert!(!positions.exists());
    assert!(f.root.join("final-panel.tsv.partial").exists());
    assert!(f.root.join("final-positions.arrow.partial").exists());
    let manifest = receipt(&panel);
    assert_eq!(manifest.params["run_status"], "resource-failed");
    assert_eq!(manifest.params["failure.kind"], "runtime-memory");
    assert!(manifest
        .outputs
        .iter()
        .all(|output| output.path.ends_with(".partial")));
    assert_eq!(manifest.self_hash_ok(), Some(true));
    assert_eq!(manifest.measurement_hash_ok(), Some(true));
}

#[test]
fn explicitly_supplied_inapplicable_defaults_are_rejected() {
    let f = Fixture::new();
    for (kind, flag, value) in [
        ("features", "--base-quality-threshold", "20"),
        ("evidence", "--max-depth", "1000"),
        ("evidence", "--min-callable-depth", "10"),
    ] {
        let output = f.command(kind).args([flag, value]).output().unwrap();
        assert_eq!(output.status.code(), Some(2));
    }
}

#[test]
fn legacy_arrow_full_batch_refuses_small_budget_and_verifies_admitted_run() {
    let f = Fixture::new();
    std::fs::write(&f.reference, format!(">chr1\n{}\n", "A".repeat(70_000))).unwrap();
    let pack = f.root.join("large.rref");
    success(
        Command::new(env!("CARGO_BIN_EXE_rosalind"))
            .args(["reference", "build", "--fasta"])
            .arg(&f.reference)
            .arg("--output")
            .arg(&pack)
            .output()
            .unwrap(),
    );
    let mut header = Header::new();
    header.push_record(HeaderRecord::new(b"HD").push_tag(b"SO", "coordinate"));
    header.push_record(
        HeaderRecord::new(b"SQ")
            .push_tag(b"SN", "chr1")
            .push_tag(b"LN", 70_000),
    );
    let mut writer = bam::Writer::from_path(&f.bam, &header, bam::Format::Bam).unwrap();
    for pos in (0..70_000).step_by(250) {
        let mut record = Record::new();
        record.set(
            format!("r{pos}").as_bytes(),
            Some(&CigarString(vec![Cigar::Match(250)])),
            &[b'A'; 250],
            &[40; 250],
        );
        record.set_tid(0);
        record.set_pos(pos);
        record.set_flags(0);
        record.set_mapq(60);
        writer.write(&record).unwrap();
    }
    drop(writer);
    for budget in [32, 256] {
        let path = f.root.join(format!("legacy-{budget}.arrow"));
        let output = Command::new(env!("CARGO_BIN_EXE_rosalind"))
            .args(["features", "--reference-pack"])
            .arg(&pack)
            .arg("--alignments")
            .arg(&f.bam)
            .args(["--format", "arrow-ipc", "--enforce", "--memory-budget-mb"])
            .arg(budget.to_string())
            .arg("-o")
            .arg(&path)
            .output()
            .unwrap();
        if budget == 32 {
            assert_eq!(output.status.code(), Some(3));
            assert!(!path.exists());
        } else {
            success(output);
            let sizes: Vec<_> = StreamReader::try_new(std::fs::File::open(&path).unwrap(), None)
                .unwrap()
                .map(|batch| batch.unwrap().num_rows())
                .collect();
            assert_eq!(sizes, [65_536, 4_464]);
            success(
                Command::new(env!("CARGO_BIN_EXE_rosalind"))
                    .args(["verify", "--manifest"])
                    .arg(format!("{}.manifest.json", path.display()))
                    .output()
                    .unwrap(),
            );
        }
    }
}
