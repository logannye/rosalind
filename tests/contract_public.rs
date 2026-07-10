//! A downstream crate can call the public orchestration layer directly and
//! receives typed outcomes rather than process termination.

use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use rosalind::contract::{
    run_column_analysis, AnalyzerIdentity, ContractRunError, ContractRunSpec, ContractVerdict,
    OutputTarget, ProducerIdentity, ReplayInvocation,
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
        invocation: ReplayInvocation::new(["run"]).option("--label", "value with spaces"),
        index,
        alignments,
        output: OutputTarget::File(output),
        manifest: None,
        mapq_threshold: 0,
        max_depth: 1000,
        max_read_len: 250,
        memory_budget_mb: None,
        enforce: false,
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
    let dir = unique_dir();
    let (index, bam) = fixture(&dir);
    let output = dir.join("must-not-exist.tsv");
    let mut configured = spec(index, bam, output.clone());
    configured.memory_budget_mb = Some(1);
    configured.enforce = true;
    let error = run_column_analysis(&mut DepthAnalyzer, configured).unwrap_err();
    assert!(matches!(error, ContractRunError::Refused(_)));
    assert!(!output.exists());
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
