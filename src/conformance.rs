//! Self-contained conformance harness for external analyzer binaries.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicU64, Ordering};

use anyhow::{anyhow, Context, Result};

use crate::provenance::{diff_receipts, RunManifest};

/// Machine-readable analyzer conformance result.
#[derive(Debug, Clone)]
pub struct ConformanceReport {
    /// Whether every required check passed.
    pub passed: bool,
    /// Stable check name to boolean result.
    pub checks: BTreeMap<String, bool>,
    /// Diagnostic details for failed checks.
    pub failures: Vec<String>,
}

impl ConformanceReport {
    /// JSON suitable for CI storage or a badge endpoint.
    pub fn to_json(&self) -> String {
        let checks = self
            .checks
            .iter()
            .map(|(name, passed)| format!("\"{}\":{passed}", escape(name)))
            .collect::<Vec<_>>()
            .join(",");
        let failures = self
            .failures
            .iter()
            .map(|failure| format!("\"{}\"", escape(failure)))
            .collect::<Vec<_>>()
            .join(",");
        format!(
            "{{\"schema\":1,\"kind\":\"rosalind-analyzer-conformance\",\"passed\":{},\"checks\":{{{checks}}},\"failures\":[{failures}]}}",
            self.passed
        )
    }
}

/// Run the complete embedded analyzer conformance suite without networking.
pub fn conform_analyzer(binary: &Path) -> Result<ConformanceReport> {
    let binary = std::fs::canonicalize(binary)
        .with_context(|| format!("cannot resolve analyzer binary {}", binary.display()))?;
    let rosalind = std::env::current_exe().context("cannot locate the rosalind executable")?;
    let root = temporary_dir()?;
    let demo = root.join("demo");
    let index = demo.join("ref.idx");
    let alignments = demo.join("sorted.bam");
    let mut checks = BTreeMap::new();
    let mut failures = Vec::new();

    let demo_result = Command::new(&rosalind)
        .args(["demo", "--output-dir"])
        .arg(&demo)
        .arg("--json")
        .output()?;
    check(
        &mut checks,
        &mut failures,
        "embedded_fixture",
        demo_result.status.success(),
        output_detail(&demo_result),
    );
    if !demo_result.status.success() {
        let report = ConformanceReport {
            passed: false,
            checks,
            failures,
        };
        std::fs::remove_dir_all(root).ok();
        return Ok(report);
    }

    let first = root.join("first.tsv");
    let second = root.join("second.tsv");
    let first_run = analyzer_run(
        &binary,
        &index,
        &alignments,
        &first,
        1,
        Some(128),
        true,
        false,
    );
    let second_run = analyzer_run(
        &binary,
        &index,
        &alignments,
        &second,
        1,
        Some(128),
        true,
        false,
    );
    check(
        &mut checks,
        &mut failures,
        "known_bound_enforced_fit",
        first_run.status.success() && second_run.status.success(),
        format!(
            "first: {}; second: {}",
            output_detail(&first_run),
            output_detail(&second_run)
        ),
    );
    let repeat_equal = std::fs::read(&first).ok() == std::fs::read(&second).ok()
        && first.exists()
        && second.exists();
    check(
        &mut checks,
        &mut failures,
        "repeat_run_byte_equality",
        repeat_equal,
        "repeat outputs differ or are absent".to_string(),
    );

    let first_receipt_path = sidecar(&first, ".manifest.json");
    let first_receipt = read_manifest(&first_receipt_path);
    let receipt_ok = first_receipt
        .as_ref()
        .is_ok_and(|manifest| manifest.self_hash_ok() == Some(true));
    check(
        &mut checks,
        &mut failures,
        "receipt_integrity",
        receipt_ok,
        first_receipt
            .as_ref()
            .err()
            .map(ToString::to_string)
            .unwrap_or_else(|| "receipt self-hash failed".to_string()),
    );
    let relocation_ok = first_receipt.as_ref().is_ok_and(|manifest| {
        let before = manifest.content_hash();
        let mut relocated = (*manifest).clone();
        for (index, file) in relocated.inputs.iter_mut().enumerate() {
            file.path = format!("/relocated/input-{index}");
        }
        for (index, file) in relocated.outputs.iter_mut().enumerate() {
            file.path = format!("/relocated/output-{index}");
        }
        before == relocated.content_hash()
    });
    check(
        &mut checks,
        &mut failures,
        "path_relocation",
        relocation_ok,
        "claim changed after path relocation".to_string(),
    );

    let relocated = root.join("relocated-inputs");
    std::fs::create_dir_all(&relocated)?;
    std::fs::copy(&index, relocated.join("renamed index.idx"))?;
    std::fs::copy(&alignments, relocated.join("renamed alignments.bam"))?;
    let dry_run = Command::new(&rosalind)
        .args(["reproduce", "--manifest"])
        .arg(&first_receipt_path)
        .arg("--inputs")
        .arg(&relocated)
        .arg("--binary")
        .arg(&binary)
        .args(["--dry-run", "--json"])
        .output()?;
    check(
        &mut checks,
        &mut failures,
        "safe_replay_plan",
        dry_run.status.success()
            && String::from_utf8_lossy(&dry_run.stdout).contains("\"validated\":true"),
        output_detail(&dry_run),
    );
    let reproduction = Command::new(&rosalind)
        .args(["reproduce", "--manifest"])
        .arg(&first_receipt_path)
        .arg("--inputs")
        .arg(&relocated)
        .arg("--binary")
        .arg(&binary)
        .args(["--no-attest", "--json"])
        .output()?;
    check(
        &mut checks,
        &mut failures,
        "external_reproduction",
        reproduction.status.success()
            && String::from_utf8_lossy(&reproduction.stdout).contains("REPRODUCED"),
        output_detail(&reproduction),
    );

    let sanitized = root.join("sanitized.json");
    let sanitize_ok =
        crate::sanitize_receipt(&first_receipt_path, &sanitized).is_ok() && sanitized.exists();
    check(
        &mut checks,
        &mut failures,
        "sanitized_receipt",
        sanitize_ok,
        "receipt sanitization failed".to_string(),
    );

    let scaled = root.join("scaled.tsv");
    let scaled_run = analyzer_run(&binary, &index, &alignments, &scaled, 2, None, false, false);
    let diff_ok = scaled_run.status.success()
        && read_manifest(&sidecar(&scaled, ".manifest.json"))
            .ok()
            .zip(first_receipt.as_ref().ok())
            .is_some_and(|(scaled, first)| {
                diff_receipts(first, &scaled)
                    .science_params
                    .iter()
                    .any(|change| change.key.contains("scale"))
            });
    check(
        &mut checks,
        &mut failures,
        "causal_parameter_diff",
        diff_ok,
        output_detail(&scaled_run),
    );

    let refused = root.join("refused.tsv");
    let refusal = analyzer_run(
        &binary,
        &index,
        &alignments,
        &refused,
        1,
        Some(1),
        true,
        false,
    );
    check(
        &mut checks,
        &mut failures,
        "upfront_refusal",
        refusal.status.code() == Some(3) && !refused.exists(),
        output_detail(&refusal),
    );

    let breached = root.join("breached.tsv");
    let breach = analyzer_command(
        &binary,
        &index,
        &alignments,
        &breached,
        1,
        Some(128),
        true,
        false,
    )
    .env("ROSALIND_FORCE_PEAK_RSS_BYTES", "8589934592")
    .output()?;
    check(
        &mut checks,
        &mut failures,
        "forced_breach_partial",
        breach.status.code() == Some(4)
            && !breached.exists()
            && sidecar(&breached, ".partial").exists(),
        output_detail(&breach),
    );

    let collision = analyzer_run(&binary, &index, &alignments, &first, 1, None, false, false);
    check(
        &mut checks,
        &mut failures,
        "safe_output_collision",
        collision.status.code() == Some(2),
        output_detail(&collision),
    );

    let passed = checks.values().all(|passed| *passed);
    let report = ConformanceReport {
        passed,
        checks,
        failures,
    };
    std::fs::remove_dir_all(root).ok();
    Ok(report)
}

#[allow(clippy::too_many_arguments)]
fn analyzer_command(
    binary: &Path,
    index: &Path,
    alignments: &Path,
    output: &Path,
    scale: u32,
    budget_mb: Option<u64>,
    enforce: bool,
    force: bool,
) -> Command {
    let mut command = Command::new(binary);
    command
        .arg("run")
        .arg("--index")
        .arg(index)
        .arg("--alignments")
        .arg(alignments)
        .arg("--output")
        .arg(output)
        .arg("--scale")
        .arg(scale.to_string());
    if let Some(budget) = budget_mb {
        command.arg("--memory-budget-mb").arg(budget.to_string());
    }
    if enforce {
        command.arg("--enforce");
    }
    if force {
        command.arg("--force");
    }
    command
}

#[allow(clippy::too_many_arguments)]
fn analyzer_run(
    binary: &Path,
    index: &Path,
    alignments: &Path,
    output: &Path,
    scale: u32,
    budget_mb: Option<u64>,
    enforce: bool,
    force: bool,
) -> Output {
    analyzer_command(
        binary, index, alignments, output, scale, budget_mb, enforce, force,
    )
    .output()
    .unwrap_or_else(|error| panic!("cannot start analyzer {}: {error}", binary.display()))
}

fn read_manifest(path: &Path) -> Result<RunManifest> {
    let text = std::fs::read_to_string(path)
        .with_context(|| format!("cannot read receipt {}", path.display()))?;
    RunManifest::from_canonical_json(&text).map_err(|error| anyhow!(error))
}

fn check(
    checks: &mut BTreeMap<String, bool>,
    failures: &mut Vec<String>,
    name: &str,
    passed: bool,
    detail: String,
) {
    checks.insert(name.to_string(), passed);
    if !passed {
        failures.push(format!("{name}: {detail}"));
    }
}

fn temporary_dir() -> Result<PathBuf> {
    static COUNTER: AtomicU64 = AtomicU64::new(0);
    let path = std::env::temp_dir().join(format!(
        "rosalind-conformance-{}-{}",
        std::process::id(),
        COUNTER.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir(&path)?;
    Ok(path)
}

fn sidecar(path: &Path, suffix: &str) -> PathBuf {
    let mut value = path.as_os_str().to_os_string();
    value.push(suffix);
    PathBuf::from(value)
}

fn output_detail(output: &Output) -> String {
    format!(
        "exit {:?}; stdout={:?}; stderr={:?}",
        output.status.code(),
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    )
}

fn escape(value: &str) -> String {
    value
        .replace('\\', "\\\\")
        .replace('"', "\\\"")
        .replace('\n', "\\n")
        .replace('\r', "\\r")
        .replace('\t', "\\t")
}
