//! A complete external analyzer: customize the reducer, keep the artifact runner.
use std::collections::BTreeMap;
use std::io::Write;
use std::path::PathBuf;

use clap::{Args, Parser, Subcommand};
use rosalind::contract::{
    AnalyzerIdentity, EnforcementMode, OutputPolicy, ProducerIdentity, ReplayInvocation,
};
use rosalind::dataset::{DatasetQuery, DatasetReadLimits};
use rosalind::evidence::{
    parse_artifact_query, run_evidence_artifact, ArtifactSelection, EvidenceAnalyzer,
    EvidenceArtifactFactory, EvidenceArtifactSource, EvidenceArtifactSpec, EvidenceBatch,
    EvidenceError, EvidenceExecution, EvidenceFields, EvidenceProfile, EvidenceRequest,
    EvidenceRequirements, EvidenceSampleSelection, EvidenceSelection,
};
use rosalind::variant_io::VariantLimits;

#[derive(Parser)]
#[command(name = "rosalind-example-evidence-analyzer", version)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}
#[derive(Subcommand)]
enum Command {
    /// Write a bounded, atomic summary and a replayable receipt.
    Run(Box<Run>),
}

#[derive(Args)]
struct Run {
    #[arg(long, required_unless_present = "dataset", conflicts_with = "dataset")]
    alignments: Option<PathBuf>,
    #[arg(long)]
    reference: Option<PathBuf>,
    #[arg(long)]
    alignment_index: Option<PathBuf>,
    #[arg(long)]
    reference_fai: Option<PathBuf>,
    #[arg(long)]
    cram_reference: Option<PathBuf>,
    #[arg(long)]
    cram_reference_fai: Option<PathBuf>,
    /// Portable evidence-dataset.manifest.json; original alignments are unnecessary.
    #[arg(long)]
    dataset: Option<PathBuf>,
    /// Explicit lineage operands supplied by Rosalind's relocated replay.
    #[arg(long, requires = "dataset")]
    dataset_artifact: Vec<PathBuf>,
    #[arg(long, conflicts_with_all = ["sites", "whole_genome", "whole_dataset", "query_json"])]
    regions: Option<PathBuf>,
    #[arg(long, conflicts_with_all = ["whole_genome", "whole_dataset", "query_json"])]
    sites: Option<PathBuf>,
    #[arg(long, conflicts_with_all = ["whole_dataset", "query_json"])]
    whole_genome: bool,
    #[arg(long, requires = "dataset", conflicts_with = "query_json")]
    whole_dataset: bool,
    /// Bounded canonical query JSON, normally supplied by replay.
    #[arg(long)]
    query_json: Option<String>,
    /// Physical field mask; this reducer requires depth(1) and allele(2) groups.
    #[arg(long)]
    fields: Option<u32>,
    #[arg(long)]
    mapq_threshold: Option<u8>,
    #[arg(long)]
    base_quality_threshold: Option<u8>,
    #[arg(long)]
    include_secondary: bool,
    #[arg(long)]
    include_supplementary: bool,
    #[arg(long)]
    include_qc_fail: bool,
    #[arg(long)]
    include_duplicates: bool,
    #[arg(long, conflicts_with = "pool_samples")]
    sample: Option<String>,
    #[arg(long)]
    pool_samples: bool,
    #[arg(long, default_value_t = 16384)]
    tile_bases: u32,
    #[arg(long, default_value_t = 250)]
    max_read_len: usize,
    #[arg(long, default_value_t = 1048576)]
    max_record_bytes: usize,
    #[arg(long, default_value_t = 8388608)]
    max_variant_header_bytes: usize,
    #[arg(long, default_value_t = 1048576)]
    max_variant_record_bytes: usize,
    #[arg(long, default_value_t = 33554432)]
    max_receipt_bytes: usize,
    #[arg(long, default_value_t = 33554432)]
    max_dataset_metadata_bytes: usize,
    #[arg(long)]
    memory_budget_mb: Option<u64>,
    #[arg(long, conflicts_with = "memory_budget_mb")]
    memory_budget_bytes: Option<u64>,
    #[arg(long)]
    enforce: bool,
    #[arg(long, requires = "enforce")]
    require_os_limit: bool,
    #[arg(short, long)]
    output: PathBuf,
    #[arg(long)]
    manifest: Option<PathBuf>,
    #[arg(long)]
    force: bool,
    /// Example analyzer parameter included in scientific identity and replay.
    #[arg(long, default_value_t = 1)]
    scale: u64,
}

fn fields() -> EvidenceFields {
    EvidenceFields::DEPTHS.union(EvidenceFields::ALLELES)
}
fn requirements() -> EvidenceRequirements {
    EvidenceRequirements {
        fields: fields(),
        requires_reference: true,
        context_bases: 0,
        // Counters, scale, borrowed writer, and Box allocation/alignment allowance.
        // The runner separately reserves its output buffer and receipt machinery.
        retained_bytes: Some(128),
    }
}

struct SummaryFactory {
    scale: u64,
}
impl EvidenceArtifactFactory for SummaryFactory {
    fn requirements(&self) -> EvidenceRequirements {
        requirements()
    }
    fn params(&self) -> BTreeMap<String, String> {
        BTreeMap::from([
            ("metric".into(), "candidate-summary-v1".into()),
            ("scale".into(), self.scale.to_string()),
        ])
    }
    fn create<'a>(
        &'a mut self,
        out: &'a mut dyn Write,
    ) -> Result<Box<dyn EvidenceAnalyzer + 'a>, EvidenceError> {
        Ok(Box::new(CandidateSummary {
            out,
            scale: self.scale,
            counts: [0; 3],
        }))
    }
}

struct CandidateSummary<'a> {
    out: &'a mut dyn Write,
    scale: u64,
    counts: [u64; 3],
}
fn add(counter: &mut u64, value: u64) -> Result<(), EvidenceError> {
    *counter = counter
        .checked_add(value)
        .ok_or(EvidenceError::CounterOverflow)?;
    Ok(())
}
impl EvidenceAnalyzer for CandidateSummary<'_> {
    fn requirements(&self) -> EvidenceRequirements {
        requirements()
    }
    fn on_batch(&mut self, batch: &EvidenceBatch) -> Result<(), EvidenceError> {
        if !batch.fields().contains(fields()) {
            return Err(EvidenceError::Analyzer(
                "candidate summary requires depth and allele fields".into(),
            ));
        }
        for row in batch.rows() {
            let depths = row
                .depths
                .ok_or_else(|| EvidenceError::Analyzer("depth group is absent".into()))?;
            let alleles = row
                .alleles
                .ok_or_else(|| EvidenceError::Analyzer("allele group is absent".into()))?;
            add(&mut self.counts[0], 1)?;
            add(&mut self.counts[1], depths.callable_depth)?;
            for alternate in row.requested_alts {
                let index = b"ACGT"
                    .iter()
                    .position(|base| base == alternate)
                    .ok_or_else(|| EvidenceError::Analyzer("invalid SNV allele".into()))?;
                add(&mut self.counts[2], alleles.allele_counts[index])?;
            }
        }
        Ok(())
    }
    fn finish(&mut self) -> Result<(), EvidenceError> {
        let [loci, callable, alternates] = self.counts.map(|value| value.checked_mul(self.scale));
        let (loci, callable, alternates) = (
            loci.ok_or(EvidenceError::CounterOverflow)?,
            callable.ok_or(EvidenceError::CounterOverflow)?,
            alternates.ok_or(EvidenceError::CounterOverflow)?,
        );
        writeln!(
            self.out,
            "selected_loci\tcallable_read_observations\tcandidate_alt_read_observations"
        )?;
        writeln!(self.out, "{loci}\t{callable}\t{alternates}")?;
        Ok(())
    }
}

impl Run {
    fn spec(self) -> Result<EvidenceArtifactSpec, Box<dyn std::error::Error>> {
        let execution = EvidenceExecution {
            memory_budget_bytes: self.memory_budget_bytes.or(self
                .memory_budget_mb
                .map(|value| {
                    value
                        .checked_mul(1 << 20)
                        .ok_or("memory budget overflows bytes")
                })
                .transpose()?),
            max_microtile_bases: self.tile_bases,
            max_read_len: self.max_read_len,
            max_record_bytes: self.max_record_bytes,
            analyzer_bytes: 0,
        };
        let query = self
            .query_json
            .as_deref()
            .map(parse_artifact_query)
            .transpose()?;
        let fields = match (self.fields, query.as_ref()) {
            (Some(bits), Some(query)) if bits != query.fields.bits() => {
                return Err("--fields disagrees with --query-json".into())
            }
            (Some(bits), _) => EvidenceFields::from_bits(bits)?,
            (None, Some(query)) => query.fields,
            (None, None) => fields(),
        };
        let selection = if let Some(path) = self.regions {
            ArtifactSelection::Bed(path)
        } else if let Some(path) = self.sites {
            ArtifactSelection::Variants(
                path,
                VariantLimits {
                    max_header_bytes: self.max_variant_header_bytes,
                    max_record_bytes: self.max_variant_record_bytes,
                },
            )
        } else if self.dataset.is_some() && !self.whole_genome && query.is_none() {
            ArtifactSelection::Stored
        } else {
            ArtifactSelection::Request
        };
        let query = query.unwrap_or(DatasetQuery {
            selection: EvidenceSelection::WholeGenome,
            fields,
        });
        let source = if let Some(manifest) = self.dataset {
            if self.reference.is_some()
                || self.alignment_index.is_some()
                || self.reference_fai.is_some()
                || self.cram_reference.is_some()
                || self.cram_reference_fai.is_some()
                || self.mapq_threshold.is_some()
                || self.base_quality_threshold.is_some()
                || self.include_secondary
                || self.include_supplementary
                || self.include_qc_fail
                || self.include_duplicates
                || self.sample.is_some()
                || self.pool_samples
            {
                return Err("--dataset uses its stored reference/profile/sample identity; native input and filter flags are incompatible".into());
            }
            EvidenceArtifactSource::Dataset {
                manifest,
                query,
                selection,
                execution,
                artifacts: self.dataset_artifact,
                limits: DatasetReadLimits {
                    max_manifest_bytes: self.max_dataset_metadata_bytes,
                    max_descriptor_bytes: self.max_dataset_metadata_bytes,
                    ..DatasetReadLimits::default()
                },
            }
        } else {
            let mut request = EvidenceRequest::coverage(
                self.alignments
                    .ok_or("--alignments or --dataset is required")?,
            );
            request.reference = self.reference;
            request.alignment_index = self.alignment_index;
            request.reference_fai = self.reference_fai;
            request.cram_reference = self.cram_reference;
            request.cram_reference_fai = self.cram_reference_fai;
            request.fields = fields;
            request.selection = query.selection;
            request.execution = execution;
            request.profile = EvidenceProfile {
                min_mapq: self.mapq_threshold.unwrap_or(20),
                min_base_quality: self.base_quality_threshold.unwrap_or(20),
                exclude_secondary: !self.include_secondary,
                exclude_supplementary: !self.include_supplementary,
                exclude_qc_fail: !self.include_qc_fail,
                exclude_duplicates: !self.include_duplicates,
            };
            request.sample_selection = if self.pool_samples {
                EvidenceSampleSelection::Pool
            } else if let Some(name) = self.sample {
                EvidenceSampleSelection::Named(name)
            } else {
                EvidenceSampleSelection::Auto
            };
            EvidenceArtifactSource::Native { request, selection }
        };
        let mut spec = EvidenceArtifactSpec::new(
            source,
            self.output,
            ProducerIdentity {
                name: "rosalind-example-evidence-analyzer".into(),
                version: env!("CARGO_PKG_VERSION").into(),
                repository: None,
                binary: "rosalind-example-evidence-analyzer".into(),
            },
            AnalyzerIdentity::new("candidate-summary", env!("CARGO_PKG_VERSION")),
            ReplayInvocation::new(["run"]).option("--scale", self.scale),
        );
        spec.manifest = self.manifest;
        spec.max_receipt_bytes = self.max_receipt_bytes;
        spec.output_policy = if self.force {
            OutputPolicy::ReplaceAtomic
        } else {
            OutputPolicy::CreateNewAtomic
        };
        spec.enforcement = if self.require_os_limit {
            EnforcementMode::RequireOsLimit
        } else if self.enforce {
            EnforcementMode::Cooperative
        } else {
            EnforcementMode::RecordOnly
        };
        spec.handle_signals = true;
        Ok(spec)
    }
}

fn main() {
    rosalind::provenance::set_build_identity(rosalind::provenance::BuildIdentity {
        code_git_sha: env!("ROSALIND_GIT_SHA").into(),
        code_dirty: env!("ROSALIND_GIT_DIRTY").into(),
        rustc_version: env!("ROSALIND_RUSTC_VERSION").into(),
        target_triple: env!("ROSALIND_TARGET").into(),
        deps_lock_blake3: env!("ROSALIND_DEPS_LOCK_BLAKE3").into(),
    });
    let Command::Run(args) = Cli::parse().command;
    let mut factory = SummaryFactory { scale: args.scale };
    let spec = args.spec().unwrap_or_else(|error| {
        eprintln!("{error}");
        std::process::exit(2)
    });
    match run_evidence_artifact(&mut factory, spec) {
        Ok(_) => eprintln!("analysis completed; artifact and receipt published"),
        Err(error) => {
            eprintln!("{error}");
            std::process::exit(error.exit_code());
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rosalind::evidence::EvidenceRow;

    #[test]
    fn checked_reducer_includes_zero_loci_and_multiple_requested_alts() {
        let batch = EvidenceBatch::from_full_rows(
            0,
            "chr1",
            0,
            vec![
                EvidenceRow {
                    callable_depth: 10,
                    allele_counts: [6, 3, 1, 0],
                    requested_alts: vec![b'C', b'G'],
                    ..EvidenceRow::default()
                },
                EvidenceRow::default(),
            ],
        );
        let mut output = Vec::new();
        let mut factory = SummaryFactory { scale: 1 };
        let expected = factory.requirements();
        {
            let mut analyzer = factory.create(&mut output).unwrap();
            assert_eq!(analyzer.requirements(), expected);
            analyzer.on_batch(&batch).unwrap();
            analyzer.finish().unwrap();
        }
        assert_eq!(
            output,
            b"selected_loci\tcallable_read_observations\tcandidate_alt_read_observations\n2\t10\t4\n"
        );
    }
    #[test]
    fn omitted_fields_are_not_zero_and_overflow_is_an_error() {
        let mut output = Vec::new();
        let mut analyzer = CandidateSummary {
            out: &mut output,
            scale: 2,
            counts: [u64::MAX, 0, 0],
        };
        assert!(matches!(
            analyzer.finish(),
            Err(EvidenceError::CounterOverflow)
        ));
        let batch = EvidenceBatch::new(0, "chr1", 0, EvidenceFields::DEPTHS, vec![]);
        assert!(matches!(
            analyzer.on_batch(&batch),
            Err(EvidenceError::Analyzer(_))
        ));
        assert!(output.is_empty());
    }

    #[test]
    fn dataset_cannot_silently_change_stored_filters() {
        let Command::Run(args) = Cli::try_parse_from([
            "analyzer",
            "run",
            "--dataset",
            "dataset.json",
            "--output",
            "out.tsv",
            "--mapq-threshold",
            "30",
        ])
        .unwrap()
        .command;
        assert!(args
            .spec()
            .unwrap_err()
            .to_string()
            .contains("stored reference/profile/sample"));
    }

    #[test]
    fn replay_preserves_exact_byte_budget_and_rejects_mismatched_fields() {
        let Command::Run(args) = Cli::try_parse_from([
            "analyzer",
            "run",
            "--alignments",
            "sample.bam",
            "--reference",
            "ref.fa",
            "--output",
            "out.tsv",
            "--memory-budget-bytes",
            "134217729",
        ])
        .unwrap()
        .command;
        let EvidenceArtifactSource::Native { request, .. } = args.spec().unwrap().source else {
            panic!("expected native source")
        };
        assert_eq!(request.execution.memory_budget_bytes, Some(134217729));

        let query = r#"{"version":1,"fields":3,"selection":{"intervals":[],"sites":null}}"#;
        let Command::Run(args) = Cli::try_parse_from([
            "analyzer",
            "run",
            "--dataset",
            "dataset.json",
            "--output",
            "out.tsv",
            "--query-json",
            query,
            "--fields",
            "1",
        ])
        .unwrap()
        .command;
        assert!(args.spec().unwrap_err().to_string().contains("disagrees"));
    }
}
