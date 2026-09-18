//! Exercise the preview through the public binary from a directory outside the
//! checkout, with separately generated samples and independently known counts.
use rosalind::dataset::{
    publish_evidence_dataset, run_dataset_with_snapshot, DatasetOptions, DescriptorLimits,
    VerifiedInputSession,
};
use rosalind::evidence::{
    EvidenceBatch, EvidenceCallback, EvidenceEngine, EvidenceFields, EvidenceRequest,
    EvidenceSelection,
};
use rosalind::selection::GenomicInterval;
use rust_htslib::bam::{
    self,
    header::HeaderRecord,
    record::{Aux, Cigar, CigarString},
};
use rust_htslib::bcf::{self, Read as BcfRead};
use serde_json::Value;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

struct Fixture {
    root: PathBuf,
    cohort: PathBuf,
    snapshot: String,
    sources: PathBuf,
}
impl Fixture {
    fn new() -> Self {
        Self::with_length(32)
    }
    fn with_length(length: u32) -> Self {
        Self::with_length_and_ref_reads(length, false)
    }
    fn with_length_and_ref_reads(length: u32, paired: bool) -> Self {
        static NEXT: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
        let root = std::env::temp_dir().join(format!(
            "rosalind-cohort-cli-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed)
        ));
        fs::create_dir(&root).unwrap();
        let root = fs::canonicalize(root).unwrap();
        let reference = root.join("reference.fa");
        fs::write(
            &reference,
            format!(">chr1\n{}\n", "A".repeat(length as usize)),
        )
        .unwrap();
        let fai = root.join("reference.fa.fai");
        fs::write(
            &fai,
            format!("chr1\t{length}\t6\t{length}\t{}\n", length + 1),
        )
        .unwrap();
        let fields = EvidenceFields::DEPTHS.union(EvidenceFields::ALLELES);
        let mut members = String::from("id\tmanifest\tgroup\n");
        let mut mappings = String::from("id\trole\tpath\n");
        for (sample, count) in [("A", 1), ("B", 2)] {
            let alignment = root.join(format!("{sample}.bam"));
            let index = root.join(format!("{sample}.bam.bai"));
            let mut header = bam::Header::new();
            header.push_record(HeaderRecord::new(b"HD").push_tag(b"SO", "coordinate"));
            header.push_record(
                HeaderRecord::new(b"SQ")
                    .push_tag(b"SN", "chr1")
                    .push_tag(b"LN", length),
            );
            header.push_record(
                HeaderRecord::new(b"RG")
                    .push_tag(b"ID", "rg1")
                    .push_tag(b"SM", sample),
            );
            let mut writer = bam::Writer::from_path(&alignment, &header, bam::Format::Bam).unwrap();
            for position in [1, 8] {
                for read in 0..count {
                    let mut record = bam::Record::new();
                    record.set(
                        format!("read-{position}-{read}").as_bytes(),
                        Some(&CigarString(vec![Cigar::Match(1)])),
                        if paired && sample == "B" && read == 1 {
                            b"A"
                        } else {
                            b"C"
                        },
                        &[35],
                    );
                    record.set_tid(0);
                    record.set_pos(position);
                    record.set_flags(0);
                    record.set_mapq(60);
                    record.push_aux(b"RG", Aux::String("rg1")).unwrap();
                    writer.write(&record).unwrap();
                }
            }
            drop(writer);
            bam::index::build(&alignment, Some(&index), bam::index::Type::Bai, 1).unwrap();
            let mut request = EvidenceRequest::new(&alignment, &reference);
            request.alignment_index = Some(index.clone());
            request.reference_fai = Some(fai.clone());
            request.fields = fields;
            request.selection = EvidenceSelection::Intervals(vec![GenomicInterval {
                contig: 0,
                start: 0,
                end: 4,
            }]);
            let roles = vec![
                ("alignments".into(), alignment),
                ("alignment-index".into(), index),
                ("reference".into(), reference.clone()),
                ("reference-fai".into(), fai.clone()),
            ];
            for (role, path) in &roles {
                mappings.push_str(&format!(
                    "{sample}\t{role}\t{}\n",
                    path.file_name().unwrap().to_str().unwrap()
                ));
            }
            let session = VerifiedInputSession::open(roles).unwrap();
            let mut engine = EvidenceEngine::open(request).unwrap();
            let namespace = session.dataset_namespace(&engine).unwrap();
            let mut callback = EvidenceCallback::with_fields(|_: &EvidenceBatch| Ok(()), 0, fields);
            let outcome = run_dataset_with_snapshot(
                &mut engine,
                &namespace,
                &DatasetOptions {
                    cache_dir: root.join(format!("cache-{sample}")),
                    workers: 1,
                    resume: false,
                },
                &mut callback,
                session.snapshot(),
            )
            .unwrap();
            let manifest =
                publish_evidence_dataset(&engine, &outcome, &session, DescriptorLimits::default())
                    .unwrap();
            members.push_str(&format!(
                "{sample}\t{}\tstudy\n",
                manifest.strip_prefix(&root).unwrap().display()
            ));
        }
        fs::write(root.join("members.tsv"), members).unwrap();
        let sources = root.join("sources.tsv");
        fs::write(&sources, mappings).unwrap();
        let cohort = root.join("cohort");
        let result = run(
            &root,
            &[
                "cohort",
                "create",
                "--cohort",
                "cohort",
                "--members",
                "members.tsv",
            ],
        );
        let snapshot = json_ok(result)["snapshot_id"].as_str().unwrap().to_owned();
        fs::write(root.join("candidates.vcf"), "##fileformat=VCFv4.3\n##contig=<ID=chr1,length=32>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t2\t.\tA\tC,G\t.\t.\t.\nchr1\t3\t.\tA\tC\t.\t.\t.\nchr1\t9\t.\tA\tC\t.\t.\t.\n").unwrap();
        Self {
            root,
            cohort,
            snapshot,
            sources,
        }
    }
    fn query(&self, command: &str, extra: &[&str]) -> Output {
        let mut args = vec![
            "cohort",
            command,
            "--cohort",
            self.cohort.to_str().unwrap(),
            "--snapshot",
            &self.snapshot,
            "--sites",
            "candidates.vcf",
        ];
        args.extend_from_slice(extra);
        run(&self.root, &args)
    }
}
impl Drop for Fixture {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.root);
    }
}
fn run(root: &Path, args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_rosalind"))
        .args(args)
        .current_dir(root)
        .output()
        .unwrap()
}
fn json_ok(output: Output) -> Value {
    assert!(
        output.status.success(),
        "stdout={} stderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    serde_json::from_slice(&output.stdout).unwrap()
}
fn rows(path: &Path) -> Vec<Vec<String>> {
    fs::read_to_string(path)
        .unwrap()
        .lines()
        .map(|line| line.split('\t').map(str::to_owned).collect())
        .collect()
}

#[test]
fn saved_only_cli_plans_missingness_and_produces_exact_relocated_reports() {
    let mut fixture = Fixture::new();
    let planned = fixture.query("extract", &["--plan"]);
    let plan = json_ok(planned);
    assert_eq!(plan["status"], "blocked");
    assert_eq!(plan["candidate_rows"], 4);
    assert_eq!(plan["members"][0]["missing_loci"], 1);
    let strict = fixture.query("extract", &["--output", "strict.tsv"]);
    assert!(!strict.status.success());
    assert!(!fixture.root.join("strict.tsv").exists());
    for name in [
        "A.bam",
        "A.bam.bai",
        "B.bam",
        "B.bam.bai",
        "reference.fa",
        "reference.fa.fai",
    ] {
        fs::remove_file(fixture.root.join(name)).unwrap();
    }
    let relocated = fixture.root.join("relocated");
    fs::rename(&fixture.cohort, &relocated).unwrap();
    fixture.cohort = relocated;
    let outcome = json_ok(fixture.query(
        "extract",
        &[
            "--missing",
            "partial",
            "--output",
            "evidence.tsv",
            "--min-callable-depth",
            "1",
        ],
    ));
    assert_eq!(outcome["sample_candidate_rows"], 8);
    assert_eq!(outcome["unmeasured_rows"], 2);
    assert_eq!(outcome["original_alignment_records_decoded"], 0);
    let result = rows(&fixture.root.join("evidence.tsv"));
    let col = |name: &str| result[0].iter().position(|field| field == name).unwrap();
    assert_eq!(result[1][col("alt_count")], "1");
    assert_eq!(
        result[2][col("alt_count")],
        "0",
        "new ALT at stored locus is observed"
    );
    assert_eq!(
        result[3][col("callable_depth")],
        "0",
        "stored zero is observed"
    );
    assert_eq!(
        result[4][col("callable_depth")],
        ".",
        "absent locus is unmeasured"
    );
    json_ok(fixture.query(
        "summarize",
        &[
            "--missing",
            "partial",
            "--output",
            "summary.tsv",
            "--min-callable-depth",
            "1",
        ],
    ));
    let summary = rows(&fixture.root.join("summary.tsv"));
    let col = |name: &str| summary[0].iter().position(|field| field == name).unwrap();
    assert_eq!(summary[1][col("n_requested")], "2");
    assert_eq!(summary[1][col("n_alt_supported")], "2");
    assert_eq!(summary[1][col("alt_total")], "3");
    assert_eq!(summary[4][col("n_observed")], "0");
    assert_eq!(summary[4][col("observed_alt_fraction_denominator")], ".");
    let existing = fs::read(fixture.root.join("evidence.tsv")).unwrap();
    assert!(!fixture
        .query(
            "extract",
            &["--missing", "partial", "--output", "evidence.tsv"]
        )
        .status
        .success());
    assert_eq!(
        existing,
        fs::read(fixture.root.join("evidence.tsv")).unwrap()
    );
    let verified = json_ok(run(
        &fixture.root,
        &[
            "cohort",
            "verify",
            "--cohort-root",
            fixture.cohort.to_str().unwrap(),
            "--snapshot",
            &fixture.snapshot,
        ],
    ));
    assert_eq!(verified["current_payloads_verified"], true);
    assert_eq!(verified["original_sources_opened"], false);
}

#[test]
fn extension_is_explicit_and_preserves_the_original_snapshot() {
    let fixture = Fixture::new();
    let old = fs::read(
        fixture
            .cohort
            .join("snapshots")
            .join(&fixture.snapshot)
            .join("snapshot.json"),
    )
    .unwrap();
    let plan = json_ok(fixture.query("extend", &["--sources", "sources.tsv", "--plan"]));
    assert_eq!(plan["raw_sources_opened"], false);
    assert_eq!(plan["affected_members"], serde_json::json!(["A", "B"]));
    let outcome =
        json_ok(fixture.query("extend", &["--sources", fixture.sources.to_str().unwrap()]));
    assert_eq!(outcome["changed"], true);
    let new = outcome["snapshot_id"].as_str().unwrap();
    assert_ne!(new, fixture.snapshot);
    assert_eq!(outcome["members"][0]["computed_loci"], 1);
    assert_eq!(
        old,
        fs::read(
            fixture
                .cohort
                .join("snapshots")
                .join(&fixture.snapshot)
                .join("snapshot.json")
        )
        .unwrap()
    );
    json_ok(run(
        &fixture.root,
        &[
            "cohort",
            "extract",
            "--cohort",
            "cohort",
            "--snapshot",
            new,
            "--sites",
            "candidates.vcf",
            "--output",
            "extended.tsv",
        ],
    ));
    let result = rows(&fixture.root.join("extended.tsv"));
    let alt = result[0]
        .iter()
        .position(|field| field == "alt_count")
        .unwrap();
    assert_eq!(result[4][alt], "1");
    assert_eq!(result[8][alt], "2");
    assert_eq!(
        json_ok(fixture.query("extract", &["--plan"]))["status"],
        "blocked"
    );
}

#[test]
fn compressed_candidates_replay_operands_and_invalid_tables_are_checked() {
    let fixture = Fixture::new();
    for (name, format) in [
        ("candidates.vcf.gz", bcf::Format::Vcf),
        ("candidates.bcf", bcf::Format::Bcf),
    ] {
        let mut reader = bcf::Reader::from_path(fixture.root.join("candidates.vcf")).unwrap();
        let header = bcf::Header::from_template(reader.header());
        let mut writer =
            bcf::Writer::from_path(fixture.root.join(name), &header, false, format).unwrap();
        for record in reader.records() {
            writer.write(&record.unwrap()).unwrap();
        }
    }
    for (name, output) in [
        ("candidates.vcf", "a.tsv"),
        ("candidates.vcf.gz", "b.tsv"),
        ("candidates.bcf", "c.tsv"),
    ] {
        json_ok(run(
            &fixture.root,
            &[
                "cohort",
                "extract",
                "--cohort",
                "cohort",
                "--snapshot",
                &fixture.snapshot,
                "--sites",
                name,
                "--missing",
                "partial",
                "--output",
                output,
            ],
        ));
    }
    assert_eq!(
        fs::read(fixture.root.join("a.tsv")).unwrap(),
        fs::read(fixture.root.join("b.tsv")).unwrap()
    );
    assert_eq!(
        fs::read(fixture.root.join("a.tsv")).unwrap(),
        fs::read(fixture.root.join("c.tsv")).unwrap()
    );
    let manifest = fixture
        .cohort
        .join("snapshots")
        .join(&fixture.snapshot)
        .join("snapshot.json");
    let plan = json_ok(run(
        &fixture.root,
        &[
            "cohort",
            "extract",
            "--cohort-snapshot-manifest",
            manifest.to_str().unwrap(),
            "--snapshot",
            &fixture.snapshot,
            "--sites",
            "candidates.vcf",
            "--cohort-members",
            "[]",
            "--plan",
            "--memory-budget-bytes",
            "536870912",
            "--enforce",
            "--max-microtile-bases",
            "7",
        ],
    ));
    assert_eq!(plan["output_rows"], 0);
    fs::write(
        fixture.root.join("bad.tsv"),
        "id\tmanifest\tgroup\nA\tx\tfirst\nA\ty\tsecond\n",
    )
    .unwrap();
    let result = run(
        &fixture.root,
        &[
            "cohort",
            "create",
            "--cohort",
            "bad-store",
            "--members",
            "bad.tsv",
        ],
    );
    assert!(!result.status.success());
    assert!(!fixture.root.join("bad-store").exists());
    assert!(!run(
        &fixture.root,
        &[
            "cohort",
            "create",
            "--cohort",
            "small-store",
            "--members",
            "members.tsv",
            "--max-table-bytes",
            "1"
        ]
    )
    .status
    .success());
    assert!(!fixture.root.join("small-store").exists());
}

#[test]
fn cohort_receipts_replay_all_output_modes_after_relocation_and_reject_changed_inputs() {
    let fixture = Fixture::new();
    for command in ["extract", "summarize"] {
        for format in ["tsv", "arrow-ipc"] {
            let output = format!("{command}-{format}");
            json_ok(fixture.query(
                command,
                &[
                    "--missing",
                    "partial",
                    "--format",
                    format,
                    "--output",
                    &output,
                ],
            ));
        }
    }
    for name in [
        "A.bam",
        "A.bam.bai",
        "B.bam",
        "B.bam.bai",
        "reference.fa",
        "reference.fa.fai",
    ] {
        fs::remove_file(fixture.root.join(name)).unwrap();
    }
    let replay = fixture.root.join("relocated-inputs");
    fs::create_dir(&replay).unwrap();
    fs::rename(&fixture.cohort, replay.join("cohort")).unwrap();
    fs::rename(
        fixture.root.join("candidates.vcf"),
        replay.join("candidates.vcf"),
    )
    .unwrap();
    for command in ["extract", "summarize"] {
        for format in ["tsv", "arrow-ipc"] {
            let manifest = format!("{command}-{format}.manifest.json");
            let result = run(
                &fixture.root,
                &[
                    "reproduce",
                    "--manifest",
                    &manifest,
                    "--inputs",
                    "relocated-inputs",
                    "--binary",
                    env!("CARGO_BIN_EXE_rosalind"),
                    "--no-attest",
                ],
            );
            assert!(
                result.status.success(),
                "{command}/{format}: stdout={} stderr={}",
                String::from_utf8_lossy(&result.stdout),
                String::from_utf8_lossy(&result.stderr)
            );
            assert!(String::from_utf8_lossy(&result.stdout).contains("REPRODUCED"));
        }
    }
    fn first_arrow(path: &Path) -> Option<PathBuf> {
        for entry in fs::read_dir(path).unwrap() {
            let path = entry.unwrap().path();
            if path.is_dir() {
                if let Some(found) = first_arrow(&path) {
                    return Some(found);
                }
            } else if path
                .extension()
                .is_some_and(|extension| extension == "arrow")
            {
                return Some(path);
            }
        }
        None
    }
    let arrow = first_arrow(&replay.join("cohort/objects")).unwrap();
    let mut bytes = fs::read(&arrow).unwrap();
    bytes[0] ^= 1;
    fs::write(&arrow, bytes).unwrap();
    let changed = run(
        &fixture.root,
        &[
            "reproduce",
            "--manifest",
            "extract-tsv.manifest.json",
            "--inputs",
            "relocated-inputs",
            "--binary",
            env!("CARGO_BIN_EXE_rosalind"),
            "--no-attest",
        ],
    );
    assert!(
        !changed.status.success(),
        "changed consumed evidence must refuse replay"
    );
}

#[test]
fn replay_transport_only_accepts_bounded_matching_cohort_requests() {
    let fixture = Fixture::new();
    let path = fixture.root.join("request.json");
    let valid = vec![
        "cohort",
        "extract",
        "--cohort",
        "cohort",
        "--snapshot",
        &fixture.snapshot,
        "--sites",
        "candidates.vcf",
        "--missing",
        "partial",
        "--plan",
    ];
    fs::write(&path, serde_json::to_vec(&valid).unwrap()).unwrap();
    let plan = json_ok(run(
        &fixture.root,
        &["cohort", "replay", "--request", "request.json"],
    ));
    assert_eq!(plan["status"], "ready");
    let mismatched = run(
        &fixture.root,
        &[
            "cohort",
            "replay",
            "--request",
            "request.json",
            "--memory-budget-bytes",
            "536870912",
            "--enforce",
        ],
    );
    assert!(!mismatched.status.success());
    assert!(String::from_utf8_lossy(&mismatched.stderr).contains("exactly match"));
    let tiny = run(
        &fixture.root,
        &[
            "cohort",
            "replay",
            "--request",
            "request.json",
            "--memory-budget-bytes",
            "1",
            "--enforce",
        ],
    );
    assert_eq!(
        tiny.status.code(),
        Some(3),
        "transport must admit before expanding argv"
    );
    for forbidden in [
        vec!["cohort", "replay", "--request", "request.json"],
        vec!["cohort", "create"],
        vec!["variants"],
        vec!["sh", "-c", "echo unsafe"],
    ] {
        fs::write(&path, serde_json::to_vec(&forbidden).unwrap()).unwrap();
        let result = run(
            &fixture.root,
            &["cohort", "replay", "--request", "request.json"],
        );
        assert!(!result.status.success());
        assert!(String::from_utf8_lossy(&result.stderr).contains("extract or cohort summarize"));
    }
    fs::File::create(&path)
        .unwrap()
        .set_len((32 << 20) + 1)
        .unwrap();
    let oversized = run(
        &fixture.root,
        &["cohort", "replay", "--request", "request.json"],
    );
    assert!(!oversized.status.success());
    assert!(String::from_utf8_lossy(&oversized.stderr).contains("exceeds 32 MiB"));
}

#[test]
fn file_based_replay_supports_normalized_candidate_queries_larger_than_inline_envelope() {
    let fixture = Fixture::with_length(4096);
    let mut candidates = String::from("##fileformat=VCFv4.3\n##contig=<ID=chr1,length=4096>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    for position in 1..=2048 {
        candidates.push_str(&format!("chr1\t{position}\t.\tA\tC\t.\t.\t.\n"));
    }
    fs::write(fixture.root.join("candidates.vcf"), candidates).unwrap();
    json_ok(fixture.query(
        "extract",
        &["--missing", "partial", "--output", "large.tsv"],
    ));
    let receipt: Value =
        serde_json::from_slice(&fs::read(fixture.root.join("large.tsv.manifest.json")).unwrap())
            .unwrap();
    assert!(
        receipt["params"]["cohort.normalized_query"]
            .as_str()
            .unwrap()
            .len()
            > 32 << 10
    );
    let argv: Vec<String> =
        serde_json::from_str(receipt["params"]["command_argv"].as_str().unwrap()).unwrap();
    assert!(argv.iter().any(|token| token == "--sites"));
    assert!(!argv.iter().any(|token| token == "--cohort-query"));
    let replay = fixture.root.join("relocated-inputs");
    fs::create_dir(&replay).unwrap();
    fs::rename(&fixture.cohort, replay.join("cohort")).unwrap();
    fs::rename(
        fixture.root.join("candidates.vcf"),
        replay.join("candidates.vcf"),
    )
    .unwrap();
    for name in [
        "A.bam",
        "A.bam.bai",
        "B.bam",
        "B.bam.bai",
        "reference.fa",
        "reference.fa.fai",
    ] {
        fs::remove_file(fixture.root.join(name)).unwrap();
    }
    let result = run(
        &fixture.root,
        &[
            "reproduce",
            "--manifest",
            "large.tsv.manifest.json",
            "--inputs",
            "relocated-inputs",
            "--binary",
            env!("CARGO_BIN_EXE_rosalind"),
            "--no-attest",
        ],
    );
    assert!(
        result.status.success(),
        "stdout={} stderr={}",
        String::from_utf8_lossy(&result.stdout),
        String::from_utf8_lossy(&result.stderr)
    );
    assert!(String::from_utf8_lossy(&result.stdout).contains("REPRODUCED"));
}

#[test]
fn paired_reports_preserve_direction_missingness_and_budget_invariance() {
    let fixture = Fixture::with_length_and_ref_reads(32, true);
    fs::write(
        fixture.root.join("pairs.tsv"),
        "id\tleft\tright\nreverse-first\tB\tA\nforward-second\tA\tB\n",
    )
    .unwrap();
    let plan = json_ok(fixture.query("compare-pairs", &["--pairs", "pairs.tsv", "--plan"]));
    assert_eq!(plan["status"], "blocked");
    assert_eq!(plan["paired_candidate_rows"], 8);
    assert_eq!(plan["pairs"][0]["id"], "reverse-first");
    assert_eq!(plan["pairs"][0]["left"], "B");
    assert_eq!(plan["pair_direction"], "right-minus-left");
    assert!(!fixture
        .query(
            "compare-pairs",
            &["--pairs", "pairs.tsv", "-o", "strict.tsv"]
        )
        .status
        .success());
    assert!(!fixture.root.join("strict.tsv").exists());
    for format in ["tsv", "arrow-ipc"] {
        let mut expected = None;
        for (index, budget, tile) in [(0, "256", "1"), (1, "512", "2"), (2, "1024", "16384")] {
            let filename = format!("pairs-{format}-{index}");
            let result = json_ok(fixture.query(
                "compare-pairs",
                &[
                    "--pairs",
                    "pairs.tsv",
                    "--missing",
                    "partial",
                    "--format",
                    format,
                    "--memory-budget-mb",
                    budget,
                    "--enforce",
                    "--tile-bases",
                    tile,
                    "-o",
                    &filename,
                ],
            ));
            assert_eq!(result["paired_candidate_rows"], 8);
            let bytes = fs::read(fixture.root.join(&filename)).unwrap();
            if let Some(expected) = &expected {
                assert_eq!(&bytes, expected);
            } else {
                expected = Some(bytes);
            }
        }
    }
    let parsed = rows(&fixture.root.join("pairs-tsv-0"));
    let column = |name: &str| parsed[0].iter().position(|value| value == name).unwrap();
    let cell = |row: usize, name: &str| parsed[row][column(name)].as_str();
    assert_eq!(cell(1, "pair_id"), "reverse-first");
    assert_eq!(cell(1, "left_callable_depth"), "2");
    assert_eq!(cell(1, "left_alt_count"), "1");
    assert_eq!(cell(1, "right_callable_depth"), "1");
    assert_eq!(cell(1, "difference_negative"), "false");
    assert_eq!(cell(1, "difference_numerator"), "1");
    assert_eq!(cell(1, "difference_denominator"), "2");
    assert_eq!(cell(1, "both_depth_eligible"), "false");
    assert_eq!(cell(3, "left_status"), "observed");
    assert_eq!(cell(3, "left_callable_depth"), "0");
    assert_eq!(cell(3, "difference_numerator"), ".");
    assert_eq!(cell(4, "left_status"), "unmeasured");
    assert_eq!(cell(4, "left_callable_depth"), ".");
    assert_eq!(cell(4, "both_depth_eligible"), ".");
    assert_eq!(cell(5, "pair_id"), "forward-second");
    assert_eq!(cell(5, "difference_negative"), "true");
    assert_eq!(cell(5, "difference_numerator"), "1");
    let receipt: Value =
        serde_json::from_slice(&fs::read(fixture.root.join("pairs-tsv-0.manifest.json")).unwrap())
            .unwrap();
    assert_eq!(
        receipt["params"]["cohort.result_semantics"],
        "cohort-pairs-v1"
    );
    assert_eq!(
        receipt["params"]["cohort.pair_direction"],
        "right-minus-left"
    );
    assert!(receipt["params"]["command_argv"]
        .as_str()
        .unwrap()
        .contains("--pairs"));
}

#[test]
fn pair_tables_refuse_ambiguity_and_replay_after_relocation() {
    let fixture = Fixture::with_length_and_ref_reads(32, true);
    for invalid in [
        "id\tleft\tright\nx\tA\tA\n",
        "id\tleft\tright\nx\tA\tB\nx\tB\tA\n",
        "id\tleft\tright\nx\tA\tunknown\n",
        "left\tright\nA\tB\n",
        "id\tleft\tright\nx\tA\tB\textra\n",
        "id\tleft\tright\n\tA\tB\n",
    ] {
        fs::write(fixture.root.join("bad.tsv"), invalid).unwrap();
        assert!(!fixture
            .query(
                "compare-pairs",
                &[
                    "--pairs",
                    "bad.tsv",
                    "--missing",
                    "partial",
                    "-o",
                    "bad-result.tsv"
                ]
            )
            .status
            .success());
        assert!(!fixture.root.join("bad-result.tsv").exists());
    }
    fs::write(
        fixture.root.join("pairs.tsv"),
        "id\tleft\tright\nexplicit\tA\tB\n",
    )
    .unwrap();
    assert!(!fixture
        .query(
            "compare-pairs",
            &["--pairs", "pairs.tsv", "--member", "A", "--plan"]
        )
        .status
        .success());
    assert!(!fixture
        .query(
            "compare-pairs",
            &[
                "--pairs",
                "pairs.tsv",
                "--max-pair-table-bytes",
                "1",
                "--plan"
            ]
        )
        .status
        .success());
    for format in ["tsv", "arrow-ipc"] {
        json_ok(fixture.query(
            "compare-pairs",
            &[
                "--pairs",
                "pairs.tsv",
                "--missing",
                "partial",
                "--format",
                format,
                "-o",
                &format!("paired-{format}"),
            ],
        ));
    }
    let relocated = fixture.root.join("relocated-inputs");
    fs::create_dir(&relocated).unwrap();
    fs::rename(&fixture.cohort, relocated.join("cohort")).unwrap();
    for name in ["pairs.tsv", "candidates.vcf"] {
        fs::rename(fixture.root.join(name), relocated.join(name)).unwrap();
    }
    for name in [
        "A.bam",
        "A.bam.bai",
        "B.bam",
        "B.bam.bai",
        "reference.fa",
        "reference.fa.fai",
    ] {
        fs::remove_file(fixture.root.join(name)).unwrap();
    }
    for format in ["tsv", "arrow-ipc"] {
        let output = run(
            &fixture.root,
            &[
                "reproduce",
                "--manifest",
                &format!("paired-{format}.manifest.json"),
                "--inputs",
                "relocated-inputs",
                "--binary",
                env!("CARGO_BIN_EXE_rosalind"),
                "--no-attest",
            ],
        );
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(String::from_utf8_lossy(&output.stdout).contains("REPRODUCED"));
    }
    fs::write(
        relocated.join("pairs.tsv"),
        "id\tleft\tright\nexplicit\tB\tA\n",
    )
    .unwrap();
    let changed = run(
        &fixture.root,
        &[
            "reproduce",
            "--manifest",
            "paired-tsv.manifest.json",
            "--inputs",
            "relocated-inputs",
            "--binary",
            env!("CARGO_BIN_EXE_rosalind"),
            "--no-attest",
        ],
    );
    assert!(
        !changed.status.success(),
        "changed pair direction must not replay under the original identity"
    );
}
