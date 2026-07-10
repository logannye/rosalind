//! Project scaffolding for downstream, statically linked analyzer binaries.

use std::path::{Path, PathBuf};

/// Files written by [`create_analyzer_project`].
#[derive(Debug, Clone)]
pub struct ScaffoldReport {
    /// Project root.
    pub root: PathBuf,
    /// Created paths, relative to `root`.
    pub files: Vec<PathBuf>,
}

/// Create a standalone analyzer crate. The destination must be absent or an
/// empty directory; existing files are never overwritten.
pub fn create_analyzer_project(name: &str, output: &Path) -> std::io::Result<ScaffoldReport> {
    validate_name(name)?;
    if output.exists() && std::fs::read_dir(output)?.next().is_some() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::AlreadyExists,
            format!("destination is not empty: {}", output.display()),
        ));
    }
    std::fs::create_dir_all(output.join("src"))?;
    std::fs::create_dir_all(output.join("tests"))?;
    std::fs::create_dir_all(output.join("scripts"))?;
    std::fs::create_dir_all(output.join(".github/workflows"))?;

    let crate_ident = name.replace('-', "_");
    let version = env!("CARGO_PKG_VERSION");
    let replacements = |template: &str| {
        template
            .replace("__PACKAGE_NAME__", name)
            .replace("__CRATE_IDENT__", &crate_ident)
            .replace("__ROSALIND_VERSION__", version)
    };
    let files = [
        ("Cargo.toml", replacements(CARGO_TEMPLATE)),
        ("build.rs", BUILD_TEMPLATE.to_string()),
        ("src/main.rs", replacements(MAIN_TEMPLATE)),
        ("tests/contract.rs", TEST_TEMPLATE.to_string()),
        ("scripts/contract-check.sh", replacements(CHECK_TEMPLATE)),
        (".github/workflows/ci.yml", CI_TEMPLATE.to_string()),
        ("README.md", replacements(README_TEMPLATE)),
    ];
    let mut created = Vec::new();
    for (relative, body) in files {
        let path = output.join(relative);
        std::fs::write(&path, body)?;
        created.push(PathBuf::from(relative));
    }
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        let path = output.join("scripts/contract-check.sh");
        let mut permissions = std::fs::metadata(&path)?.permissions();
        permissions.set_mode(0o755);
        std::fs::set_permissions(path, permissions)?;
    }
    Ok(ScaffoldReport {
        root: output.to_path_buf(),
        files: created,
    })
}

fn validate_name(name: &str) -> std::io::Result<()> {
    let valid = !name.is_empty()
        && name
            .bytes()
            .all(|b| b.is_ascii_lowercase() || b.is_ascii_digit() || b == b'-')
        && name.as_bytes()[0].is_ascii_lowercase()
        && !name.ends_with('-')
        && !name.contains("--");
    if valid {
        Ok(())
    } else {
        Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "analyzer name must be lowercase kebab-case and start with a letter",
        ))
    }
}

const CARGO_TEMPLATE: &str = r#"[package]
name = "__PACKAGE_NAME__"
version = "0.1.0"
edition = "2021"
rust-version = "1.83"
build = "build.rs"

[dependencies]
anyhow = "1"
clap = { version = "4.5", features = ["derive"] }
rosalind-bio = { version = "=__ROSALIND_VERSION__", features = ["contract-testkit"] }

[build-dependencies]
rosalind-build-info = "=0.1.0"
"#;

const BUILD_TEMPLATE: &str = "fn main() {\n    rosalind_build_info::emit();\n}\n";

const MAIN_TEMPLATE: &str = r##"use std::collections::BTreeMap;
use std::io::Write;
use std::path::PathBuf;

use anyhow::Result;
use clap::{Parser, Subcommand};
use rosalind::contract::{
    run_column_analysis, AnalyzerIdentity, AnalyzerMemoryModel, ContractRunError, ContractRunSpec,
    EnforcementMode, OutputPolicy, OutputTarget, ProducerIdentity, ReplayInvocation,
};
use rosalind::{ColumnAnalyzer, PileupColumn};

#[derive(Parser)]
#[command(name = "__PACKAGE_NAME__")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Run this analyzer under Rosalind's bounded, receipted contract.
    Run {
        #[arg(long)] index: PathBuf,
        #[arg(long)] alignments: PathBuf,
        #[arg(short, long)] output: PathBuf,
        #[arg(long)] manifest: Option<PathBuf>,
        #[arg(long)] memory_budget_mb: Option<u64>,
        #[arg(long, default_value_t = 1000)] max_depth: u32,
        #[arg(long, default_value_t = 250)] max_read_len: u32,
        #[arg(long, default_value_t = 0)] mapq_threshold: u8,
        #[arg(long)] enforce: bool,
        /// Analyzer-specific multiplier, recorded in replay argv and the claim.
        #[arg(long, default_value_t = 1)] scale: u32,
    },
}

struct CustomAnalyzer {
    scale: u32,
}

impl ColumnAnalyzer for CustomAnalyzer {
    fn header(&self) -> Option<String> {
        Some("#contig\tpos\tcustom_value\n".to_string())
    }

    fn params(&self) -> BTreeMap<String, String> {
        BTreeMap::from([
            ("metric".to_string(), "scaled-depth".to_string()),
            ("scale".to_string(), self.scale.to_string()),
        ])
    }

    fn on_column(
        &mut self,
        column: &PileupColumn,
        contig: &str,
        out: &mut dyn Write,
    ) -> std::io::Result<()> {
        writeln!(out, "{contig}\t{}\t{}", column.locus.pos.0 + 1, column.depth().saturating_mul(self.scale))
    }
}

fn main() -> Result<()> {
    rosalind::provenance::set_build_identity(rosalind::provenance::BuildIdentity {
        code_git_sha: env!("ROSALIND_GIT_SHA").to_string(),
        code_dirty: env!("ROSALIND_GIT_DIRTY").to_string(),
        rustc_version: env!("ROSALIND_RUSTC_VERSION").to_string(),
        target_triple: env!("ROSALIND_TARGET").to_string(),
        deps_lock_blake3: env!("ROSALIND_DEPS_LOCK_BLAKE3").to_string(),
    });
    let Cli { command } = Cli::parse();
    let Command::Run {
        index, alignments, output, manifest, memory_budget_mb,
        max_depth, max_read_len, mapq_threshold, enforce, scale,
    } = command;
    let mut analyzer = CustomAnalyzer { scale };
    let spec = ContractRunSpec {
        producer: ProducerIdentity {
            name: "__PACKAGE_NAME__".to_string(),
            version: env!("CARGO_PKG_VERSION").to_string(),
            repository: None,
            binary: "__PACKAGE_NAME__".to_string(),
        },
        analyzer: AnalyzerIdentity::new("__PACKAGE_NAME__", env!("CARGO_PKG_VERSION")),
        analyzer_memory: AnalyzerMemoryModel::Fixed {
            model_id: "fixed-additional-v1".to_string(),
            max_additional_bytes: 0,
        },
        invocation: ReplayInvocation::new(["run"]).option("--scale", scale),
        index,
        alignments,
        output: OutputTarget::File(output),
        output_policy: OutputPolicy::CreateNewAtomic,
        manifest,
        mapq_threshold,
        max_depth,
        max_read_len,
        memory_budget_mb,
        enforcement: if enforce {
            EnforcementMode::Cooperative
        } else {
            EnforcementMode::RecordOnly
        },
    };
    match run_column_analysis(&mut analyzer, spec) {
        Ok(outcome) => {
            eprintln!("claim: {}", outcome.claim_hash.as_deref().unwrap_or("(none)"));
            Ok(())
        }
        Err(ContractRunError::Refused(report)) => {
            eprintln!("contract: REFUSE — predicted {} MiB > budget {} MiB",
                report.predicted_peak_rss_bytes / (1 << 20), report.budget_mb);
            std::process::exit(3);
        }
        Err(ContractRunError::Breached(outcome)) => {
            eprintln!("contract: VIOLATED — peak {} MiB", outcome.peak_rss_bytes / (1 << 20));
            std::process::exit(4);
        }
        Err(error) => Err(error.into()),
    }
}
"##;

const TEST_TEMPLATE: &str = r#"use rosalind::contract::testkit::assert_paths_are_portable;
use rosalind::provenance::{FileHash, RunManifest};

#[test]
fn receipt_claim_is_path_portable() {
    let mut manifest = RunManifest::new("run");
    manifest.inputs.push(FileHash { path: "/one/input.bam".into(), blake3: "abc".into() });
    manifest.outputs.push(FileHash { path: "/one/output.tsv".into(), blake3: "def".into() });
    manifest.finalize();
    assert_paths_are_portable(&manifest);
}
"#;

const CHECK_TEMPLATE: &str = r#"#!/usr/bin/env bash
set -euo pipefail
ROSALIND_BIN="${ROSALIND_BIN:-rosalind}"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

"$ROSALIND_BIN" demo --output-dir "$WORK/demo" --json >/dev/null
cargo build --release
BIN="target/release/__PACKAGE_NAME__"

for n in 1 2; do
  "$BIN" run --index "$WORK/demo/ref.idx" --alignments "$WORK/demo/sorted.bam" \
    --memory-budget-mb 128 --enforce --output "$WORK/out$n.tsv"
done
cmp "$WORK/out1.tsv" "$WORK/out2.tsv"
"$ROSALIND_BIN" verify --manifest "$WORK/out1.tsv.manifest.json" --json >/dev/null
"$ROSALIND_BIN" reproduce --manifest "$WORK/out1.tsv.manifest.json" \
  --inputs "$WORK/demo" --binary "$BIN" --json >/dev/null

"$BIN" run --index "$WORK/demo/ref.idx" --alignments "$WORK/demo/sorted.bam" \
  --scale 2 --output "$WORK/scaled.tsv"
set +e
diff_report=$("$ROSALIND_BIN" diff "$WORK/out1.tsv.manifest.json" \
  "$WORK/scaled.tsv.manifest.json")
diff_code=$?
set -e
test "$diff_code" -eq 1
printf '%s' "$diff_report" | grep -q 'analyzer.scale'

set +e
"$BIN" run --index "$WORK/demo/ref.idx" --alignments "$WORK/demo/sorted.bam" \
  --memory-budget-mb 1 --enforce --output "$WORK/refused.tsv"
code=$?
set -e
test "$code" -eq 3
test ! -e "$WORK/refused.tsv"
"#;

const CI_TEMPLATE: &str = r#"name: CI
on: [push, pull_request]
jobs:
  contract:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: dtolnay/rust-toolchain@stable
      - run: cargo install rosalind-bio --locked
      - run: cargo test
      - run: bash scripts/contract-check.sh
"#;

const README_TEMPLATE: &str = r#"# __PACKAGE_NAME__

A standalone Rosalind `ColumnAnalyzer` that inherits bounded execution,
deterministic output, the live RSS governor, and a content-addressed receipt.

```sh
cargo build --release
rosalind demo --output-dir /tmp/rosalind-demo
target/release/__PACKAGE_NAME__ run \
  --index /tmp/rosalind-demo/ref.idx \
  --alignments /tmp/rosalind-demo/sorted.bam \
  --memory-budget-mb 128 --enforce -o output.tsv
rosalind verify --manifest output.tsv.manifest.json
rosalind reproduce --manifest output.tsv.manifest.json \
  --inputs /tmp/rosalind-demo --binary target/release/__PACKAGE_NAME__
```

Run `bash scripts/contract-check.sh` for determinism, verification, external
reproduction, and refusal checks.
"#;

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicU64, Ordering};

    fn temporary(name: &str) -> PathBuf {
        static COUNTER: AtomicU64 = AtomicU64::new(0);
        std::env::temp_dir().join(format!(
            "rosalind-scaffold-{name}-{}-{}",
            std::process::id(),
            COUNTER.fetch_add(1, Ordering::Relaxed)
        ))
    }

    #[test]
    fn creates_the_complete_analyzer_project_and_refuses_overwrite() {
        let destination = temporary("complete");
        let report = create_analyzer_project("depth-track", &destination).unwrap();
        for required in [
            "Cargo.toml",
            "build.rs",
            "src/main.rs",
            "tests/contract.rs",
            "scripts/contract-check.sh",
            ".github/workflows/ci.yml",
            "README.md",
        ] {
            assert!(report.files.contains(&PathBuf::from(required)));
            assert!(destination.join(required).is_file());
        }
        let cargo = std::fs::read_to_string(destination.join("Cargo.toml")).unwrap();
        assert!(cargo.contains(&format!("version = \"={}\"", env!("CARGO_PKG_VERSION"))));
        assert!(cargo.contains("contract-testkit"));
        let error = create_analyzer_project("depth-track", &destination).unwrap_err();
        assert_eq!(error.kind(), std::io::ErrorKind::AlreadyExists);
        std::fs::remove_dir_all(destination).ok();
    }

    #[test]
    fn rejects_names_that_are_not_cargo_kebab_case() {
        let destination = temporary("invalid");
        let error = create_analyzer_project("Not Valid", &destination).unwrap_err();
        assert_eq!(error.kind(), std::io::ErrorKind::InvalidInput);
        assert!(!destination.exists());
    }
}
