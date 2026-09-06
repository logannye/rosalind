//! Offline scientific queries over verified portable evidence partitions.
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};

use crate::FeatureFormat;
use anyhow::{bail, Result};
use clap::{Args, Subcommand};
use rosalind::core::governor::{checkpoint, MemoryGovernor};
use rosalind::core::CoreError;
use rosalind::dataset::*;
use rosalind::evidence::*;
use rosalind::provenance::{CommandCapture, RunManifest};
use rosalind::util::atomic::{commit_group, ensure_destination, AtomicFile};
use rosalind::util::rss::peak_rss_bytes;

#[derive(Subcommand, Debug)]
pub(crate) enum DatasetAction {
    /// Inspect verified metadata; does not scan partition bodies or original sources.
    Inspect(DatasetInput),
    /// Verify every stored partition and evidence row without original alignments.
    Verify(DatasetInput),
    /// Query stored evidence as bounded Arrow batches or TSV.
    Extract {
        #[command(flatten)]
        query: QueryArgs,
        #[command(flatten)]
        output: OutputArgs,
        #[arg(long, value_enum, default_value_t = FeatureFormat::ArrowIpc)]
        format: FeatureFormat,
    },
    /// Calculate panel summaries from stored evidence without alignment decoding.
    PanelQc {
        #[command(flatten)]
        query: QueryArgs,
        #[command(flatten)]
        output: OutputArgs,
        #[arg(long, default_value_t = DEFAULT_MIN_CALLABLE_DEPTH)]
        min_callable_depth: u64,
    },
    /// Export bounded Parquet parts and a verified lineage receipt into a new directory.
    Export {
        #[command(flatten)]
        query: QueryArgs,
        #[arg(short, long)]
        output: PathBuf,
    },
}

#[derive(Args, Debug)]
pub(crate) struct DatasetInput {
    /// Portable evidence-dataset.manifest.json created by analyze --cache-dir.
    #[arg(long)]
    dataset: PathBuf,
    /// Explicit verified dependencies captured for strict byte replay.
    #[arg(long = "dataset-artifact", hide = true)]
    artifacts: Vec<PathBuf>,
    /// Largest portable descriptor or receipt accepted; parsing memory is planned separately.
    #[arg(long, default_value_t = 33_554_432)]
    max_dataset_metadata_bytes: usize,
    /// Whole native process budget in MiB; independent of downstream consumers.
    #[arg(long)]
    memory_budget_mb: Option<u64>,
}
#[derive(Args, Debug)]
pub(crate) struct QueryArgs {
    #[command(flatten)]
    input: DatasetInput,
    /// BED subset; omitted selection means every stored locus, including uncovered loci.
    #[arg(long, conflicts_with = "sites")]
    regions: Option<PathBuf>,
    /// SNV VCF, VCF.gz or BCF subset with validated stored reference alleles.
    #[arg(long, conflicts_with = "regions")]
    sites: Option<PathBuf>,
    /// Comma-separated evidence groups, all, all-supported, or none; defaults to stored fields.
    #[arg(long, value_delimiter = ',')]
    fields: Option<Vec<String>>,
    /// Print the complete memory plan without creating artifacts.
    #[arg(long)]
    plan: bool,
}
#[derive(Args, Debug)]
pub(crate) struct OutputArgs {
    /// File destination; stdout when omitted (stdout is not a verified artifact).
    #[arg(short, long)]
    output: Option<PathBuf>,
    /// Receipt destination; defaults to <output>.manifest.json.
    #[arg(long)]
    manifest: Option<PathBuf>,
    /// Atomically replace the file and receipt after successful completion.
    #[arg(long)]
    force: bool,
}

pub(crate) fn run(action: DatasetAction) -> Result<()> {
    if let Err(error) = run_inner(action) {
        let code = if let Some(error) = error.downcast_ref::<DatasetError>() {
            error.exit_code()
        } else {
            match error.downcast_ref::<EvidenceError>() {
                Some(EvidenceError::Refused { .. }) => 3,
                Some(EvidenceError::Core(CoreError::BudgetExceeded { .. }))
                | Some(EvidenceError::RecordLimit(_)) => 4,
                Some(EvidenceError::Analyzer(_)) | Some(EvidenceError::CounterOverflow) => 1,
                _ => 2,
            }
        };
        eprintln!("{error:#}");
        std::process::exit(code);
    }
    Ok(())
}
fn run_inner(action: DatasetAction) -> Result<()> {
    let input = match &action {
        DatasetAction::Inspect(input) | DatasetAction::Verify(input) => input,
        DatasetAction::Extract { query, .. }
        | DatasetAction::PanelQc { query, .. }
        | DatasetAction::Export { query, .. } => &query.input,
    };
    let started = Instant::now();
    let budget = input
        .memory_budget_mb
        .map(|mb| {
            mb.checked_mul(1 << 20)
                .filter(|b| *b > 0)
                .ok_or_else(|| anyhow::anyhow!("memory budget must be positive and fit u64"))
        })
        .transpose()?;
    let _governor = budget
        .map(|b| MemoryGovernor::start(b, Duration::from_millis(100), peak_rss_bytes))
        .transpose()?;
    checkpoint().map_err(EvidenceError::from)?;
    let mut dataset = VerifiedEvidenceDataset::open(
        &input.dataset,
        DatasetReadLimits {
            memory_budget_bytes: budget,
            max_manifest_bytes: input.max_dataset_metadata_bytes,
            max_descriptor_bytes: input.max_dataset_metadata_bytes,
            ..Default::default()
        },
    )?;
    let execution = EvidenceExecution {
        memory_budget_bytes: budget,
        ..Default::default()
    };
    match action {
        DatasetAction::Inspect(_) => {
            dataset.verify_unchanged()?;
            println!("{}", serde_json::to_string_pretty(dataset.descriptor())?);
        }
        DatasetAction::Verify(_) => {
            let query = DatasetQuery {
                selection: dataset.selection(),
                fields: EvidenceFields::from_bits(0)?,
            };
            let mut sink =
                EvidenceCallback::with_fields(|_: &EvidenceBatch| Ok(()), 0, query.fields);
            let stats = dataset.visit_batches(&query, &mut sink, &execution)?;
            dataset.verify_unchanged()?;
            println!(
                "{}",
                serde_json::json!({"status":"verified", "verified_loci":stats.emitted_loci, "verified_partition_files":dataset.verified_partition_hashes()?.len(), "original_sources_rehashed":false, "dataset_namespace":dataset.descriptor().dataset_namespace})
            );
        }
        DatasetAction::Extract {
            query,
            output,
            format,
        } => query_output(
            &mut dataset,
            query,
            output,
            format,
            None,
            &execution,
            started,
        )?,
        DatasetAction::PanelQc {
            query,
            output,
            min_callable_depth,
        } => query_output(
            &mut dataset,
            query,
            output,
            FeatureFormat::Tsv,
            Some(min_callable_depth),
            &execution,
            started,
        )?,
        DatasetAction::Export { query, output } => {
            let plan_only = query.plan;
            let max_receipt_bytes = query.input.max_dataset_metadata_bytes;
            let (query, _, selection_guard) = make_query(&dataset, &query, None)?;
            let options = ParquetExportOptions {
                max_receipt_bytes,
                ..Default::default()
            };
            let plan = plan_parquet_export_with_inputs(
                &dataset,
                &query,
                &options,
                &execution,
                selection_guard.as_ref(),
            )?;
            if plan_only {
                println!("{}", plan_json(&plan));
                return Ok(());
            }
            let destination = absolute_destination(&output)?;
            let source_root = Path::new(&dataset.source_hashes()[0].path)
                .parent()
                .unwrap();
            if destination.starts_with(source_root) {
                bail!("write exports outside the immutable dataset directory");
            }
            // Selection inputs must stay unchanged through export publication.
            if let Some(guard) = &selection_guard {
                guard.verify()?;
            }
            let result = export_parquet_dataset_with_inputs(
                &mut dataset,
                &query,
                output,
                &options,
                &execution,
                selection_guard.as_ref(),
            )?;
            eprintln!(
                "Parquet dataset: {} ({} loci, planned {} bytes)",
                result.manifest_path.display(),
                result.stats.emitted_loci,
                plan.predicted_peak_rss_bytes
            );
        }
    }
    Ok(())
}

fn parse_fields(names: &[String]) -> Result<EvidenceFields> {
    let mut fields = EvidenceFields::from_bits(0)?;
    for name in names {
        let group = match name.as_str() {
            "all" if names.len() == 1 => EvidenceFields::ALL,
            "all-supported" if names.len() == 1 => EvidenceFields::ALL_SUPPORTED,
            "none" if names.len() == 1 => EvidenceFields::from_bits(0)?,
            "depths" => EvidenceFields::DEPTHS,
            "alleles" => EvidenceFields::ALLELES,
            "strands" => EvidenceFields::STRANDS,
            "quality-sums" => EvidenceFields::QUALITY_SUMS,
            "quality-histograms" => EvidenceFields::QUALITY_HISTOGRAMS,
            "read-position" => EvidenceFields::READ_POSITION,
            "allele-quality" => EvidenceFields::ALLELE_QUALITY,
            _ => bail!("invalid evidence field group {name:?}"),
        };
        fields = fields.union(group);
    }
    Ok(fields)
}
fn make_query(
    dataset: &VerifiedEvidenceDataset,
    args: &QueryArgs,
    panel: Option<u64>,
) -> Result<(DatasetQuery, CommandCapture, Option<VerifiedInputSession>)> {
    let mut capture = CommandCapture::new(if panel.is_some() {
        "dataset panel-qc"
    } else {
        "dataset extract"
    });
    let parent = dataset
        .source_hashes()
        .first()
        .ok_or_else(|| anyhow::anyhow!("verified parent receipt identity absent"))?;
    capture.input_hashed("--dataset", &parent.path, &parent.blake3);
    capture.opt(
        "--max-dataset-metadata-bytes",
        args.input.max_dataset_metadata_bytes,
    );
    let selected_input = args
        .regions
        .as_ref()
        .map(|p| ("regions", p))
        .or_else(|| args.sites.as_ref().map(|p| ("sites", p)));
    let selection_guard = selected_input
        .map(|(role, path)| VerifiedInputSession::open([(role.to_string(), path.clone())]))
        .transpose()?;
    if let Some(guard) = &selection_guard {
        for source in guard.identities() {
            capture.input_hashed(&format!("--{}", source.role), &source.path, &source.blake3);
        }
    }
    let selection = if let Some(path) = &args.regions {
        EvidenceSelection::from_bed(path, dataset.contigs())?
    } else if let Some(path) = &args.sites {
        EvidenceSelection::from_vcf(path, dataset.contigs())?
    } else {
        dataset.selection()
    };
    if let Some(guard) = &selection_guard {
        guard.verify()?;
    }
    let panel_fields = EvidenceFields::DEPTHS.union(EvidenceFields::QUALITY_SUMS);
    let fields = args
        .fields
        .as_ref()
        .map(|names| parse_fields(names))
        .transpose()?
        .unwrap_or_else(|| {
            if panel.is_some() {
                panel_fields
            } else {
                dataset.fields()
            }
        });
    if panel.is_some() && !fields.contains(panel_fields) {
        bail!("panel-qc requires depths and quality-sums");
    }
    let names = fields.names().join(",");
    capture.opt("--fields", if names.is_empty() { "none" } else { &names });
    if let Some(mb) = args.input.memory_budget_mb {
        capture.opt("--memory-budget-mb", mb);
    }
    Ok((DatasetQuery { selection, fields }, capture, selection_guard))
}
fn plan_json(plan: &DatasetReadPlan) -> serde_json::Value {
    serde_json::json!({"model":plan.model_id,"baseline_rss_bytes":plan.baseline_rss_bytes,"metadata_bytes":plan.metadata_bytes,"source_decoder_bytes":plan.source_decoder_bytes,"projection_bytes":plan.projection_bytes,"analyzer_bytes":plan.analyzer_bytes,"predicted_peak_rss_bytes":plan.predicted_peak_rss_bytes,"selected_loci":plan.selected_loci,"source_fields":plan.source_fields.bits(),"fields":plan.output_fields.bits(),"output_batch_rows":plan.output_batch_rows})
}

#[allow(clippy::too_many_arguments)]
fn query_output(
    dataset: &mut VerifiedEvidenceDataset,
    args: QueryArgs,
    output: OutputArgs,
    format: FeatureFormat,
    panel_depth: Option<u64>,
    execution: &EvidenceExecution,
    started: Instant,
) -> Result<()> {
    let (query, mut capture, guard) = make_query(dataset, &args, panel_depth)?;
    let mut panel = if let Some(depth) = panel_depth {
        let bed = args
            .regions
            .as_ref()
            .ok_or_else(|| anyhow::anyhow!("dataset panel-qc requires --regions BED"))?;
        capture.opt("--min-callable-depth", depth);
        Some(
            PanelQcAnalyzer::new(PanelTarget::from_bed(bed, dataset.contigs())?)?
                .with_min_callable_depth(depth),
        )
    } else {
        None
    };
    if panel_depth.is_none() {
        capture.opt(
            "--format",
            if format == FeatureFormat::ArrowIpc {
                "arrow-ipc"
            } else {
                "tsv"
            },
        );
    }
    let encoder_bytes = if let Some(panel) = &panel {
        panel.additional_memory_bytes().unwrap()
    } else if format == FeatureFormat::ArrowIpc {
        EvidenceArrowWriter::with_fields(io::sink(), query.fields)
            .additional_memory_bytes()
            .unwrap()
    } else {
        EvidenceTsvWriter::with_fields(io::sink(), query.fields)
            .additional_memory_bytes()
            .unwrap()
    };
    let reserve = encoder_bytes
        .checked_add((dataset.descriptor().partitions.len() as u64).saturating_mul(24 << 10))
        .and_then(|n| n.checked_add(2 << 20))
        .ok_or_else(|| anyhow::anyhow!("receipt memory envelope overflow"))?;
    let planner = EvidenceCallback::with_fields(|_: &EvidenceBatch| Ok(()), reserve, query.fields);
    let plan = dataset.plan(&query, &planner, execution)?;
    if args.plan {
        println!("{}", plan_json(&plan));
        return Ok(());
    }
    let receipt_path = output.manifest.or_else(|| {
        output
            .output
            .as_ref()
            .map(|p| PathBuf::from(format!("{}.manifest.json", p.display())))
    });
    let mut destinations = Vec::new();
    for path in [output.output.as_ref(), receipt_path.as_ref()]
        .into_iter()
        .flatten()
    {
        ensure_destination(path, output.force)?;
        let absolute = absolute_destination(path)?;
        if destinations.contains(&absolute) {
            bail!("artifact and receipt destinations must differ");
        }
        // Guard all descriptor-owned files, even partitions the query does not touch.
        let root = Path::new(&dataset.source_hashes()[0].path)
            .parent()
            .unwrap();
        if absolute.starts_with(root) {
            bail!("write derived artifacts outside the immutable dataset directory");
        }
        if let Some(guard) = &guard {
            if guard
                .identities()
                .iter()
                .any(|s| Path::new(&s.path) == absolute)
            {
                bail!("output must not replace a selection input");
            }
        }
        destinations.push(absolute);
    }
    let mut pending = output
        .output
        .as_ref()
        .map(|p| AtomicFile::create(p))
        .transpose()?;
    let stats = {
        let stdout = io::stdout();
        let mut writer: Box<dyn Write + '_> = if let Some(file) = &mut pending {
            Box::new(io::BufWriter::new(file.file_mut()))
        } else {
            Box::new(stdout.lock())
        };
        let mut bounded_execution = execution.clone();
        bounded_execution.analyzer_bytes = reserve;
        let stats = if let Some(panel) = &mut panel {
            let stats = dataset.visit_batches(&query, panel, &bounded_execution)?;
            panel.write_tsv(&mut writer, dataset.contigs())?;
            stats
        } else if format == FeatureFormat::ArrowIpc {
            dataset.visit_batches(
                &query,
                &mut EvidenceArrowWriter::with_fields(&mut writer, query.fields),
                &bounded_execution,
            )?
        } else {
            dataset.visit_batches(
                &query,
                &mut EvidenceTsvWriter::with_fields(&mut writer, query.fields),
                &bounded_execution,
            )?
        };
        writer.flush()?;
        stats
    };
    if let Some(guard) = &guard {
        guard.verify()?;
    }
    dataset.verify_unchanged()?;
    if let (Some(file), Some(path)) = (&mut pending, &output.output) {
        file.file_mut().sync_all()?;
        let hash = rosalind::provenance::blake3_file(file.temporary_path())?;
        capture.output_hashed("-o", &path.display().to_string(), &hash);
    }
    let dependencies = dataset
        .source_hashes()
        .iter()
        .cloned()
        .chain(dataset.verified_partition_hashes()?)
        .collect::<Vec<_>>();
    for requested in &args.input.artifacts {
        let path = std::fs::canonicalize(requested)?;
        if !dependencies.iter().any(|s| Path::new(&s.path) == path) {
            let hash = rosalind::provenance::blake3_file(&path)?;
            if !dependencies.iter().any(|s| s.blake3 == hash) {
                bail!("--dataset-artifact must match a verified dependency of this query");
            }
        }
    }
    for source in dependencies.iter().skip(1) {
        capture.input_hashed("--dataset-artifact", &source.path, &source.blake3);
    }
    let mut receipt = RunManifest::new(if panel_depth.is_some() {
        "dataset panel-qc"
    } else {
        "dataset extract"
    });
    capture.record_into(&mut receipt);
    for (key, value) in [
        ("run_status", "completed".to_string()),
        ("producer.name", "rosalind".into()),
        ("producer.binary", "rosalind".into()),
        ("producer.version", env!("CARGO_PKG_VERSION").into()),
        (
            "dataset.namespace",
            dataset.descriptor().dataset_namespace.clone(),
        ),
        (
            "dataset.compatibility_blake3",
            dataset.descriptor().compatibility_blake3.clone(),
        ),
        ("evidence.fields", query.fields.bits().to_string()),
        (
            "evidence.fields_version",
            query.fields.mask_version().to_string(),
        ),
        ("evidence.schema", query.fields.schema_version().to_string()),
        ("evidence.profile", EvidenceProfile::ID.into()),
        (
            "evidence.sample_scope",
            dataset.descriptor().sample_scope.canonical_json()?,
        ),
        ("evidence.semantics", EVIDENCE_SEMANTICS_VERSION.into()),
        (
            "dataset.query_semantics",
            "persisted-evidence-query-v1".into(),
        ),
        (
            "artifact.output.0.role",
            if panel_depth.is_some() {
                "panel-qc"
            } else {
                "evidence"
            }
            .into(),
        ),
        (
            "artifact.output.0.format",
            if format == FeatureFormat::ArrowIpc {
                "arrow-ipc"
            } else {
                "tsv"
            }
            .into(),
        ),
        (
            "contract.assurance",
            if execution.memory_budget_bytes.is_some() {
                "declared-bound-cooperative"
            } else {
                "observed-only"
            }
            .into(),
        ),
    ] {
        receipt.params.insert(key.into(), value);
    }
    receipt.tool_version = env!("CARGO_PKG_VERSION").into();
    receipt.params.insert(
        "dataset.query".into(),
        canonical_dataset_query(&query, dataset.contigs())?,
    );
    let query_digest = dataset_query_digest(&query, dataset.contigs())?;
    let identity = format!(
        "persisted-analysis-v1\n{}\n{}\n{}\n{:?}\n",
        dataset.descriptor().compatibility_blake3,
        query_digest,
        receipt.subcommand,
        panel_depth
    );
    receipt.params.insert(
        "science.blake3".into(),
        blake3::hash(identity.as_bytes()).to_hex().to_string(),
    );
    if pending.is_none() {
        receipt
            .params
            .retain(|key, _| !key.starts_with("artifact.output."));
    }
    receipt
        .measurements
        .insert("execution.alignment_record_visits".into(), "0".into());
    receipt.measurements.insert(
        "execution.emitted_loci".into(),
        stats.emitted_loci.to_string(),
    );
    receipt.measurements.insert(
        "execution.elapsed_ms".into(),
        started.elapsed().as_millis().to_string(),
    );
    receipt
        .measurements
        .insert("execution.plan".into(), plan_json(&plan).to_string());
    receipt
        .measurements
        .insert("original_sources_rehashed".into(), "false".into());
    receipt
        .measurements
        .insert("peak_rss_bytes".into(), peak_rss_bytes().to_string());
    receipt.finalize();
    let mut receipt_file = receipt_path
        .as_ref()
        .map(|p| AtomicFile::create(p))
        .transpose()?;
    if let Some(file) = &mut receipt_file {
        file.file_mut()
            .write_all(receipt.to_canonical_json().as_bytes())?;
        file.file_mut().sync_all()?;
    }
    // Include initial receipt encoding/sync in the observation, then refresh its
    // measurement block before the final guarded publication.
    receipt
        .measurements
        .insert("peak_rss_bytes".into(), peak_rss_bytes().to_string());
    receipt.measurements.insert(
        "execution.elapsed_ms".into(),
        started.elapsed().as_millis().to_string(),
    );
    receipt.finalize();
    if let Some(file) = &mut receipt_file {
        use std::io::{Seek, SeekFrom};
        file.file_mut().set_len(0)?;
        file.file_mut().seek(SeekFrom::Start(0))?;
        file.file_mut()
            .write_all(receipt.to_canonical_json().as_bytes())?;
        file.file_mut().sync_all()?;
    }
    if let Some(guard) = &guard {
        guard.verify()?;
    }
    dataset.verify_unchanged()?;
    checkpoint().map_err(EvidenceError::from)?;
    if let Some(budget) = execution.memory_budget_bytes {
        let peak = peak_rss_bytes();
        if peak > budget {
            return Err(EvidenceError::Core(CoreError::BudgetExceeded {
                needed: peak,
                budget,
            })
            .into());
        }
    }
    let mut group = Vec::new();
    if let (Some(file), Some(path)) = (pending, output.output) {
        group.push((file, path));
    }
    if let (Some(file), Some(path)) = (receipt_file, receipt_path) {
        group.push((file, path));
    }
    commit_group(group, output.force)?;
    Ok(())
}
fn absolute_destination(path: &Path) -> Result<PathBuf> {
    if path.exists() {
        return Ok(std::fs::canonicalize(path)?);
    }
    let parent = path
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    Ok(std::fs::canonicalize(parent)?.join(
        path.file_name()
            .ok_or_else(|| anyhow::anyhow!("destination must name a file"))?,
    ))
}
