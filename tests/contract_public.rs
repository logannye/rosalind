//! A downstream crate can call the public orchestration layer directly and
//! receives typed outcomes rather than process termination.

use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use rosalind::contract::{
    run_column_analysis, AnalyzerIdentity, AnalyzerMemoryModel, ContractRunError, ContractRunSpec,
    ContractVerdict, EnforcementMode, OutputPolicy, OutputTarget, ProducerIdentity,
    ReplayInvocation,
};
use rosalind::{ColumnAnalyzer, PileupColumn};

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_rosalind")
}

fn unique_dir() -> PathBuf {
    static COUNTER: AtomicU64 = AtomicU64::new(0);
    let n = COUNTER.fetch_add(1, Ordering::Relaxed);
    let dir = std::env::temp_dir().join(format!(
        "rosalind-public-contract-{}-{n}",
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn cli(args: &[&str]) {
    let output = Command::new(bin()).args(args).output().unwrap();
    assert!(
        output.status.success(),
        "CLI fixture failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

fn fixture(dir: &Path) -> (PathBuf, PathBuf) {
    let sequence = "ACGTACGTACGTACGTACGTACGTACGTACGT";
    let fasta = dir.join("ref.fa");
    let fastq = dir.join("reads.fastq");
    std::fs::write(&fasta, format!(">chr1\n{sequence}\n")).unwrap();
    std::fs::write(
        &fastq,
        format!("@r1\n{}\n+\n{}\n", &sequence[..16], "I".repeat(16)),
    )
    .unwrap();
    let index = dir.join("ref.idx");
    let raw = dir.join("raw.bam");
    let sorted = dir.join("sorted.bam");
    cli(&[
        "index",
        "--reference",
        fasta.to_str().unwrap(),
        "--output",
        index.to_str().unwrap(),
    ]);
    cli(&[
        "align",
        "--reference",
        fasta.to_str().unwrap(),
        "--reads",
        fastq.to_str().unwrap(),
        "--format",
        "bam",
        "--output",
        raw.to_str().unwrap(),
    ]);
    cli(&[
        "sort",
        "--input",
        raw.to_str().unwrap(),
        "--output",
        sorted.to_str().unwrap(),
    ]);
    (index, sorted)
}

#[derive(Default)]
struct DepthAnalyzer;

impl ColumnAnalyzer for DepthAnalyzer {
    fn header(&self) -> Option<String> {
        Some("#contig\tpos\tdepth\n".to_string())
    }

    fn on_column(
        &mut self,
        column: &PileupColumn,
        contig: &str,
        output: &mut dyn Write,
    ) -> std::io::Result<()> {
        writeln!(
            output,
            "{contig}\t{}\t{}",
            column.locus.pos.0 + 1,
            column.depth()
        )
    }
}

fn spec(index: PathBuf, alignments: PathBuf, output: PathBuf) -> ContractRunSpec {
    ContractRunSpec {
        producer: ProducerIdentity {
            name: "external-test".to_string(),
            version: "1.2.3".to_string(),
            repository: Some("https://example.test/external".to_string()),
            binary: "external-test".to_string(),
        },
        analyzer: AnalyzerIdentity::new("depth", "4.5.6"),
        analyzer_memory: AnalyzerMemoryModel::Unknown,
        invocation: ReplayInvocation::new(["run"]).option("--label", "value with spaces"),
        index,
        alignments,
        output: OutputTarget::File(output),
        output_policy: OutputPolicy::CreateNewAtomic,
        manifest: None,
        mapq_threshold: 0,
        max_depth: 1000,
        max_read_len: 250,
        memory_budget_mb: None,
        enforcement: EnforcementMode::RecordOnly,
    }
}

#[test]
fn direct_runner_completes_and_seals_external_identity() {
    let dir = unique_dir();
    let (index, bam) = fixture(&dir);
    let output = dir.join("depth.tsv");
    let outcome =
        run_column_analysis(&mut DepthAnalyzer, spec(index, bam, output.clone())).unwrap();
    assert_eq!(outcome.verdict, ContractVerdict::Unset);
    assert!(outcome.claim_hash.is_some());
    let receipt = std::fs::read_to_string(format!("{}.manifest.json", output.display())).unwrap();
    for expected in [
        "\"tool_version\":\"1.2.3\"",
        "\"producer.name\":\"external-test\"",
        "\"analyzer.id\":\"depth\"",
        "value with spaces",
    ] {
        assert!(receipt.contains(expected), "missing {expected}: {receipt}");
    }
    std::fs::remove_dir_all(dir).ok();
}

#[test]
fn direct_runner_refuses_before_primary_output_creation() {
    // Refusal starts the process-global governor before reference validation.
    // Its deliberately tiny budget must not trip concurrent record-only runs
    // in this test binary. Keep the real public API call in a child process.
    const ISOLATED: &str = "ROSALIND_TEST_CONTRACT_REFUSAL_ISOLATED";
    if std::env::var_os(ISOLATED).is_none() {
        let output = Command::new(std::env::current_exe().unwrap())
            .args([
                "--exact",
                "direct_runner_refuses_before_primary_output_creation",
                "--nocapture",
            ])
            .env(ISOLATED, "1")
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "isolated contract refusal failed:\n{}\n{}",
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        );
        return;
    }
    let dir = unique_dir();
    let (index, bam) = fixture(&dir);
    let output = dir.join("must-not-exist.tsv");
    let mut configured = spec(index, bam, output.clone());
    configured.memory_budget_mb = Some(1);
    configured.analyzer_memory = AnalyzerMemoryModel::Fixed {
        model_id: "test-v1".to_string(),
        max_additional_bytes: 0,
    };
    configured.enforcement = EnforcementMode::Cooperative;
    let error = run_column_analysis(&mut DepthAnalyzer, configured).unwrap_err();
    assert!(matches!(error, ContractRunError::Refused(_)));
    assert!(!output.exists());
    std::fs::remove_dir_all(dir).ok();
}

#[test]
fn enforced_unknown_analyzer_is_rejected_before_output_creation() {
    let dir = unique_dir();
    let (index, bam) = fixture(&dir);
    let output = dir.join("unknown-must-not-exist.tsv");
    let mut configured = spec(index, bam, output.clone());
    configured.memory_budget_mb = Some(128);
    configured.enforcement = EnforcementMode::Cooperative;
    let error = run_column_analysis(&mut DepthAnalyzer, configured).unwrap_err();
    assert!(matches!(error, ContractRunError::UnknownAnalyzerBound));
    assert!(!output.exists());
    assert!(!PathBuf::from(format!("{}.partial", output.display())).exists());
    std::fs::remove_dir_all(dir).ok();
}

#[test]
fn fixed_analyzer_contribution_is_reported_exactly_once() {
    let dir = unique_dir();
    let (index, bam) = fixture(&dir);
    let output = dir.join("bounded.tsv");
    let mut configured = spec(index, bam, output);
    configured.analyzer_memory = AnalyzerMemoryModel::Fixed {
        model_id: "retained-v1".to_string(),
        max_additional_bytes: 12_345,
    };
    let outcome = run_column_analysis(&mut DepthAnalyzer, configured).unwrap();
    assert_eq!(outcome.analyzer_predicted_bytes, Some(12_345));
    let receipt = std::fs::read_to_string(outcome.manifest_path.unwrap()).unwrap();
    assert!(receipt.contains("\"analyzer.max_additional_bytes\":\"12345\""));
    std::fs::remove_dir_all(dir).ok();
}

#[test]
fn create_new_policy_preserves_existing_output() {
    let dir = unique_dir();
    let (index, bam) = fixture(&dir);
    let output = dir.join("existing.tsv");
    std::fs::write(&output, b"keep me\n").unwrap();
    let error =
        run_column_analysis(&mut DepthAnalyzer, spec(index, bam, output.clone())).unwrap_err();
    assert!(matches!(error, ContractRunError::OutputExists(path) if path == output));
    assert_eq!(std::fs::read(&output).unwrap(), b"keep me\n");
    std::fs::remove_dir_all(dir).ok();
}

#[test]
fn replace_policy_atomically_replaces_existing_output() {
    let dir = unique_dir();
    let (index, bam) = fixture(&dir);
    let output = dir.join("replace.tsv");
    std::fs::write(&output, b"old bytes\n").unwrap();
    let mut configured = spec(index, bam, output.clone());
    configured.output_policy = OutputPolicy::ReplaceAtomic;
    let outcome = run_column_analysis(&mut DepthAnalyzer, configured).unwrap();
    assert_eq!(outcome.output_path.as_deref(), Some(output.as_path()));
    assert_ne!(std::fs::read(&output).unwrap(), b"old bytes\n");
    std::fs::remove_dir_all(dir).ok();
}

#[test]
fn analyzer_replay_options_cannot_shadow_contract_flags() {
    let dir = unique_dir();
    let mut configured = spec(
        dir.join("missing.idx"),
        dir.join("missing.bam"),
        dir.join("missing.tsv"),
    );
    configured.invocation = ReplayInvocation::new(["run"]).option("--max-depth", 7);
    let error = run_column_analysis(&mut DepthAnalyzer, configured).unwrap_err();
    assert!(matches!(
        error,
        ContractRunError::InvalidConfiguration(message) if message.contains("collides")
    ));
    std::fs::remove_dir_all(dir).ok();
}

#[test]
fn required_os_limit_accepts_a_finite_cgroup_ceiling() {
    let dir = unique_dir();
    let (index, bam) = fixture(&dir);
    let output = dir.join("cgroup.tsv");
    let result = Command::new(bin())
        .args([
            "features",
            "--index",
            index.to_str().unwrap(),
            "--alignments",
            bam.to_str().unwrap(),
            "--memory-budget-mb",
            "128",
            "--enforce",
            "--require-os-limit",
            "--output",
            output.to_str().unwrap(),
        ])
        .env("ROSALIND_TEST_CGROUP_MEMORY_MAX", "134217728")
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let receipt = std::fs::read_to_string(format!("{}.manifest.json", output.display())).unwrap();
    assert!(receipt.contains("\"contract.assurance\":\"cgroup-v2\""));
    assert!(receipt.contains("\"os.memory_limit_bytes\":\"134217728\""));
    std::fs::remove_dir_all(dir).ok();
}

#[test]
fn required_os_limit_refuses_unavailable_or_broad_limits_before_output() {
    let dir = unique_dir();
    let (index, bam) = fixture(&dir);
    for (label, limit) in [("unlimited", "max"), ("too-broad", "268435456")] {
        let output = dir.join(format!("{label}.tsv"));
        let result = Command::new(bin())
            .args([
                "features",
                "--index",
                index.to_str().unwrap(),
                "--alignments",
                bam.to_str().unwrap(),
                "--memory-budget-mb",
                "128",
                "--enforce",
                "--require-os-limit",
                "--output",
                output.to_str().unwrap(),
            ])
            .env("ROSALIND_TEST_CGROUP_MEMORY_MAX", limit)
            .output()
            .unwrap();
        assert_eq!(result.status.code(), Some(3), "{label}");
        assert!(!output.exists(), "{label}");
        assert!(
            String::from_utf8_lossy(&result.stderr).contains("OS enforcement unavailable"),
            "{label}: {}",
            String::from_utf8_lossy(&result.stderr)
        );
    }
    std::fs::remove_dir_all(dir).ok();
}

#[test]
fn exact_capacity_failure_preserves_partial_and_distinct_receipt_without_rss_budget() {
    use rust_htslib::bam;
    use rust_htslib::bam::Read;
    let dir = unique_dir();
    let (index, alignments) = fixture(&dir);
    let reader = bam::Reader::from_path(&alignments).unwrap();
    let header = bam::Header::from_template(reader.header());
    drop(reader);
    let mut writer = bam::Writer::from_path(&alignments, &header, bam::Format::Bam).unwrap();
    for (name, sequence) in [
        (b"short".as_slice(), b"AA".as_slice()),
        (b"long".as_slice(), b"CCCCCCCC".as_slice()),
    ] {
        let mut record = bam::Record::new();
        let cigar =
            bam::record::CigarString(vec![bam::record::Cigar::Match(sequence.len() as u32)]);
        record.set(name, Some(&cigar), sequence, &vec![40; sequence.len()]);
        record.set_tid(0);
        record.set_pos(0);
        record.set_flags(0);
        record.set_mapq(60);
        writer.write(&record).unwrap();
    }
    drop(writer);
    let output = dir.join("capacity.tsv");
    let mut configured = spec(index, alignments, output.clone());
    configured.max_depth = 1;
    let Err(ContractRunError::Breached(outcome)) =
        run_column_analysis(&mut DepthAnalyzer, configured)
    else {
        panic!("exact capacity must produce a typed breached outcome");
    };
    let capacity = outcome.capacity_exceeded.unwrap();
    assert_eq!(
        (
            capacity.contig,
            capacity.position,
            capacity.capacity,
            capacity.required
        ),
        (0, 0, 1, 2)
    );
    assert_eq!(
        outcome.verdict,
        ContractVerdict::Unset,
        "capacity must not be reported as excess RSS"
    );
    assert!(!output.exists());
    assert!(outcome.partial_output_path.as_ref().unwrap().is_file());
    let receipt = rosalind::provenance::RunManifest::from_canonical_json(
        &std::fs::read_to_string(outcome.manifest_path.unwrap()).unwrap(),
    )
    .unwrap();
    assert_eq!(receipt.params["pileup.semantics"], "exact-or-fail-v1");
    assert_eq!(receipt.params["run_status"], "capacity-exceeded");
    assert_eq!(receipt.params["failure.kind"], "capacity-exceeded");
    assert_eq!(receipt.params["over_max_depth"], "0");
    assert_eq!(receipt.self_hash_ok(), Some(true));
    assert_eq!(
        receipt.outputs[0].blake3,
        rosalind::provenance::blake3_file(&outcome.partial_output_path.unwrap()).unwrap()
    );
    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn invalid_receipt_claim_does_not_publish_or_replace_output() {
    struct InvalidClaim;
    impl ColumnAnalyzer for InvalidClaim {
        fn params(&self) -> std::collections::BTreeMap<String, String> {
            std::collections::BTreeMap::from([("pileup.semantics".into(), "invalid".into())])
        }
        fn on_column(
            &mut self,
            _: &PileupColumn,
            _: &str,
            output: &mut dyn Write,
        ) -> std::io::Result<()> {
            writeln!(output, "would have been published")
        }
    }
    let dir = unique_dir();
    let (index, bam) = fixture(&dir);
    let output = dir.join("keep.tsv");
    std::fs::write(&output, b"previous valid artifact").unwrap();
    let mut configured = spec(index, bam, output.clone());
    configured.output_policy = OutputPolicy::ReplaceAtomic;
    configured.analyzer = AnalyzerIdentity::new("invalid", "1").with_param_prefix("");
    let result = run_column_analysis(&mut InvalidClaim, configured);
    assert!(matches!(
        result,
        Err(ContractRunError::InvalidConfiguration(_))
    ));
    assert_eq!(std::fs::read(&output).unwrap(), b"previous valid artifact");
    assert!(!PathBuf::from(format!("{}.manifest.json", output.display())).exists());
    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn receipt_staging_failure_preserves_existing_output() {
    let dir = unique_dir();
    let (index, bam) = fixture(&dir);
    let output = dir.join("keep.tsv");
    std::fs::write(&output, b"previous valid artifact").unwrap();
    let mut configured = spec(index, bam, output.clone());
    configured.output_policy = OutputPolicy::ReplaceAtomic;
    configured.manifest = Some(dir.join("missing-parent/receipt.json"));
    assert!(matches!(
        run_column_analysis(&mut DepthAnalyzer, configured),
        Err(ContractRunError::Io(_))
    ));
    assert_eq!(std::fs::read(&output).unwrap(), b"previous valid artifact");
    std::fs::remove_dir_all(dir).unwrap();
}
