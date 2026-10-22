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
        Self::with_length(4099)
    }
    fn with_length(length: usize) -> Self {
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
        std::fs::write(&reference, format!(">chr1\n{}\n", "A".repeat(length))).unwrap();
        std::fs::write(
            root.join("ref.fa.fai"),
            format!("chr1\t{length}\t6\t{length}\t{}\n", length + 1),
        )
        .unwrap();
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
                .push_tag(b"LN", length),
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
        std::fs::write(
            &bed,
            format!("chr1\t0\t{length}\twhole\nchr1\t2\t8\toverlap\n"),
        )
        .unwrap();
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
fn dataset_coverage_is_bounded_metadata_only_and_leaves_saved_queries_strict() {
    let f = Fixture::new();
    std::fs::write(&f.bed, "chr1\t0\t2\nchr1\t20\t21\n").unwrap();
    let artifact = f.root.join("sparse.arrow");
    success(
        f.command("evidence")
            .args([
                "--format",
                "arrow-ipc",
                "--fields",
                "depths,alleles",
                "--cache-dir",
            ])
            .arg(f.root.join("cache"))
            .arg("-o")
            .arg(&artifact)
            .output()
            .unwrap(),
    );
    let parent =
        PathBuf::from(&receipt(&artifact).measurements["execution.evidence_dataset_manifest"]);
    let sites = f.root.join("candidates.vcf");
    let header = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=4099>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    std::fs::write(
        &sites,
        format!("{header}chr1\t21\t.\tA\tC\t.\t.\t.\nchr1\t10\t.\tA\tG\t.\t.\t.\nchr1\t1\t.\tA\tG\t.\t.\t.\nchr1\t1\t.\tA\tC\t.\t.\t.\n"),
    )
    .unwrap();
    std::fs::remove_file(&f.reference).unwrap();
    std::fs::remove_file(&f.bam).unwrap();
    let command = |action: &str| {
        let mut command = Command::new(env!("CARGO_BIN_EXE_rosalind"));
        command
            .args(["dataset", action, "--dataset"])
            .arg(&parent)
            .arg("--sites")
            .arg(&sites);
        command
    };
    let response = success(command("coverage").output().unwrap());
    let coverage: serde_json::Value = serde_json::from_slice(&response.stdout).unwrap();
    assert_eq!(coverage["schema"], "rosalind-dataset-coverage-v1");
    assert_eq!(coverage["verification_level"], "metadata_only");
    assert_eq!(coverage["partition_bodies_verified"], false);
    assert_eq!(coverage["original_sources_rehashed"], false);
    assert_eq!(coverage["reference_validation"], "not_performed");
    assert_eq!(coverage["requested_loci"], 3);
    assert_eq!(coverage["covered_loci"], 2);
    assert_eq!(coverage["missing_loci"], 1);
    assert_eq!(coverage["query"]["fields"], 3);
    assert_eq!(coverage["covered"]["sites"][0]["position"], 0);
    assert_eq!(
        coverage["covered"]["sites"][0]["alternates"],
        serde_json::json!([67, 71])
    );
    assert_eq!(coverage["covered"]["sites"][1]["position"], 20);
    assert_eq!(coverage["missing"]["sites"][0]["position"], 9);
    assert_eq!(
        coverage["parent"]["blake3"],
        rosalind::provenance::blake3_file(&parent).unwrap()
    );
    assert_eq!(
        coverage["descriptor"]["blake3"],
        rosalind::provenance::blake3_file(
            &parent.parent().unwrap().join("dataset.descriptor.json")
        )
        .unwrap()
    );
    // The coverage adapter must not relax the existing reader's missingness rule.
    let strict = command("extract").arg("--plan").output().unwrap();
    assert!(!strict.status.success());
    assert!(String::from_utf8_lossy(&strict.stderr).contains("absent from the dataset"));
    let bounded = command("coverage")
        .args(["--max-output-bytes", "64"])
        .output()
        .unwrap();
    assert!(!bounded.status.success());
    assert!(bounded.stdout.is_empty());
    assert!(String::from_utf8_lossy(&bounded.stderr).contains("max-output-bytes"));
    let exact_bound = success(
        command("coverage")
            .arg("--max-output-bytes")
            .arg(response.stdout.len().to_string())
            .output()
            .unwrap(),
    );
    assert_eq!(exact_bound.stdout, response.stdout);
    let below_bound = command("coverage")
        .arg("--max-output-bytes")
        .arg((response.stdout.len() - 1).to_string())
        .output()
        .unwrap();
    assert!(!below_bound.status.success());
    assert!(below_bound.stdout.is_empty());
    let absent_fields = command("coverage")
        .args(["--fields", "strands"])
        .output()
        .unwrap();
    assert!(!absent_fields.status.success());
    assert!(String::from_utf8_lossy(&absent_fields.stderr).contains("fields absent"));
    // Reordering and merging duplicate ALT annotations preserves native identity.
    std::fs::write(
        &sites,
        format!("{header}chr1\t1\t.\tA\tC,G\t.\t.\t.\nchr1\t10\t.\tA\tG\t.\t.\t.\nchr1\t21\t.\tA\tC\t.\t.\t.\n"),
    )
    .unwrap();
    let repeat: serde_json::Value =
        serde_json::from_slice(&success(command("coverage").output().unwrap()).stdout).unwrap();
    assert_eq!(coverage["query_blake3"], repeat["query_blake3"]);
    assert_eq!(coverage["covered"], repeat["covered"]);
    std::fs::write(&sites, format!("{header}chr1\t10\t.\tA\tG\t.\t.\t.\n")).unwrap();
    let missing: serde_json::Value =
        serde_json::from_slice(&success(command("coverage").output().unwrap()).stdout).unwrap();
    assert_eq!(missing["covered_loci"], 0);
    assert_eq!(missing["missing_loci"], 1);
    assert_eq!(missing["covered"]["sites"], serde_json::json!([]));
    std::fs::write(&sites, header).unwrap();
    let empty: serde_json::Value =
        serde_json::from_slice(&success(command("coverage").output().unwrap()).stdout).unwrap();
    assert_eq!(empty["requested_loci"], 0);
    assert_eq!(empty["covered_loci"], 0);
    assert_eq!(empty["missing_loci"], 0);
    // Native BED conversion and normalization are shared with strict extraction.
    std::fs::write(&f.bed, "chr1\t0\t10\nchr1\t1\t3\n").unwrap();
    let bed_result = success(
        Command::new(env!("CARGO_BIN_EXE_rosalind"))
            .args(["dataset", "coverage", "--dataset"])
            .arg(&parent)
            .arg("--regions")
            .arg(&f.bed)
            .output()
            .unwrap(),
    );
    let bed: serde_json::Value = serde_json::from_slice(&bed_result.stdout).unwrap();
    assert_eq!(bed["covered_loci"], 2);
    assert_eq!(bed["missing_loci"], 8);
    assert!(bed["covered"]["sites"].is_null());
    assert_eq!(
        bed["missing"]["intervals"],
        serde_json::json!([{"contig":0,"start":2,"end":10}])
    );
    // REF checks require observed saved rows, and are deliberately deferred.
    std::fs::write(&sites, format!("{header}chr1\t1\t.\tT\tC\t.\t.\t.\n")).unwrap();
    success(command("coverage").output().unwrap());
    let wrong_ref = command("extract").output().unwrap();
    assert!(!wrong_ref.status.success());
    // Missing Arrow bodies do not undermine the honest metadata-only statement.
    std::fs::write(&sites, format!("{header}chr1\t1\t.\tA\tC\t.\t.\t.\n")).unwrap();
    let descriptor: serde_json::Value = serde_json::from_slice(
        &std::fs::read(parent.parent().unwrap().join("dataset.descriptor.json")).unwrap(),
    )
    .unwrap();
    let arrow = parent.parent().unwrap().join(
        descriptor["partitions"][0]["arrow"]["path"]
            .as_str()
            .unwrap(),
    );
    std::fs::remove_file(arrow).unwrap();
    success(command("coverage").output().unwrap());
    assert!(!command("extract").output().unwrap().status.success());
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

#[test]
fn panel_boundary_buffers_reuse_capacity_under_three_process_budgets() {
    // A full tile followed by smaller boundary windows formerly retained both
    // allocations and breached an admitted 64 MiB budget on macOS.
    let fixture = Fixture::with_length(200_000);
    let mut expected = None;
    for fields in ["all", "depths,quality-sums"] {
        for budget in [64, 96, 128] {
            let output = fixture.root.join(format!("pressure-{fields}-{budget}.tsv"));
            success(
                fixture
                    .command("panel-qc")
                    .args([
                        "--fields",
                        fields,
                        "--memory-budget-mb",
                        &budget.to_string(),
                        "--output",
                    ])
                    .arg(&output)
                    .output()
                    .unwrap(),
            );
            let bytes = std::fs::read(&output).unwrap();
            if let Some(before) = &expected {
                assert_eq!(&bytes, before);
            } else {
                expected = Some(bytes);
            }
            let manifest = receipt(&output);
            assert!(
                manifest.measurements["peak_rss_bytes"]
                    .parse::<u64>()
                    .unwrap()
                    <= budget * (1 << 20)
            );
        }
    }
}
