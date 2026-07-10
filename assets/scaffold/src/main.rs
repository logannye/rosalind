use std::collections::BTreeMap;
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
        #[arg(long)]
        index: PathBuf,
        #[arg(long)]
        alignments: PathBuf,
        #[arg(short, long)]
        output: PathBuf,
        #[arg(long)]
        manifest: Option<PathBuf>,
        #[arg(long)]
        memory_budget_mb: Option<u64>,
        #[arg(long, default_value_t = 1000)]
        max_depth: u32,
        #[arg(long, default_value_t = 250)]
        max_read_len: u32,
        #[arg(long, default_value_t = 0)]
        mapq_threshold: u8,
        #[arg(long)]
        enforce: bool,
        #[arg(long, requires = "enforce")]
        require_os_limit: bool,
        #[arg(long)]
        force: bool,
        /// Analyzer-specific multiplier, recorded in replay argv and the claim.
        #[arg(long, default_value_t = 1)]
        scale: u32,
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
        writeln!(
            out,
            "{contig}\t{}\t{}",
            column.locus.pos.0 + 1,
            column.depth().saturating_mul(self.scale)
        )
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
        index,
        alignments,
        output,
        manifest,
        memory_budget_mb,
        max_depth,
        max_read_len,
        mapq_threshold,
        enforce,
        require_os_limit,
        force,
        scale,
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
        output_policy: if force {
            OutputPolicy::ReplaceAtomic
        } else {
            OutputPolicy::CreateNewAtomic
        },
        manifest,
        mapq_threshold,
        max_depth,
        max_read_len,
        memory_budget_mb,
        enforcement: if require_os_limit {
            EnforcementMode::RequireOsLimit
        } else if enforce {
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
            eprintln!(
                "contract: REFUSE — predicted {} MiB > budget {} MiB",
                report.predicted_peak_rss_bytes / (1 << 20),
                report.budget_mb
            );
            std::process::exit(3);
        }
        Err(ContractRunError::Breached(outcome)) => {
            eprintln!(
                "contract: VIOLATED — peak {} MiB; partial {:?}",
                outcome.peak_rss_bytes / (1 << 20),
                outcome.partial_output_path
            );
            std::process::exit(4);
        }
        Err(ContractRunError::OutputExists(path)) => {
            eprintln!("output already exists: {}", path.display());
            std::process::exit(2);
        }
        Err(ContractRunError::UnknownAnalyzerBound) => {
            eprintln!("contract: REFUSE — analyzer memory bound is unknown");
            std::process::exit(3);
        }
        Err(ContractRunError::OsEnforcementUnavailable(message)) => {
            eprintln!("contract: REFUSE — {message}");
            std::process::exit(3);
        }
        Err(error) => Err(error.into()),
    }
}
