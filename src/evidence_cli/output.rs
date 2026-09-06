use super::*;
use rosalind::dataset::{plan_dataset, run_dataset_with_snapshot, DatasetError, DatasetOptions};

struct TemporaryCache(Option<PathBuf>);
impl Drop for TemporaryCache {
    fn drop(&mut self) {
        if let Some(path) = &self.0 {
            let _ = std::fs::remove_dir_all(path);
        }
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn execute_outputs(
    command: EvidenceCommand,
    engine: &mut EvidenceEngine,
    mut capture: CommandCapture,
    receipt_path: Option<PathBuf>,
    started: Instant,
    initial_rss: u64,
    hashing_ms: u128,
    os_limit: Option<u64>,
    science_digest: String,
    input_snapshot: InputSnapshot,
) -> Result<()> {
    input_snapshot.verify()?;
    let mut panel = if command.panel {
        let targets = PanelTarget::from_bed(
            command.selection.regions.as_ref().unwrap(),
            engine.contigs(),
        )?;
        let analyzer = PanelQcAnalyzer::new(targets)?
            .with_min_callable_depth(command.options.min_callable_depth);
        Some(analyzer)
    } else {
        None
    };
    let arrow_memory = EvidenceArrowWriter::new(io::sink())
        .additional_memory_bytes()
        .unwrap();
    let tsv_memory = EvidenceTsvWriter::new(io::sink())
        .additional_memory_bytes()
        .unwrap();
    let consumer_bytes = panel
        .as_ref()
        .and_then(EvidenceAnalyzer::additional_memory_bytes)
        .unwrap_or(0)
        + if command.options.position_output.is_some()
            || (!command.panel && command.options.format == FeatureFormat::ArrowIpc)
        {
            arrow_memory
        } else {
            tsv_memory
        };
    let planner = EvidenceCallback::new(|_: &EvidenceBatch| Ok(()), consumer_bytes);
    let plan = engine.plan_for_analyzer(&planner)?.clone();
    let temporary_cache = TemporaryCache(
        if command.options.cache_dir.is_none() && command.options.workers > 1 {
            Some(std::env::temp_dir().join(format!(
                    "rosalind-workers-{}-{}",
                    std::process::id(),
                    std::time::SystemTime::now()
                        .duration_since(std::time::UNIX_EPOCH)?
                        .as_nanos()
                )))
        } else {
            None
        },
    );
    let dataset_options = command
        .options
        .cache_dir
        .as_ref()
        .or(temporary_cache.0.as_ref())
        .map(|path| DatasetOptions {
            cache_dir: path.clone(),
            resume: command.options.resume,
            workers: command.options.workers,
        });
    let cache_digest = if dataset_options.is_some() {
        let binary_hash = hash_file(&std::env::current_exe()?)?;
        blake3::hash(
            format!("rosalind-cache-producer-v1\n{science_digest}\n{binary_hash}\n").as_bytes(),
        )
        .to_hex()
        .to_string()
    } else {
        science_digest.clone()
    };
    let dataset_plan = dataset_options
        .as_ref()
        .map(|options| plan_dataset(engine, options, &planner))
        .transpose()?;
    let predicted_peak = dataset_plan
        .as_ref()
        .map_or(plan.predicted_peak_rss_bytes, |plan| {
            plan.predicted_peak_rss_bytes
        });
    if command.options.plan {
        input_snapshot.verify()?;
        println!("{{\"model\":\"{}\",\"baseline_rss_bytes\":{},\"fixed_bytes\":{},\"bytes_per_locus\":{},\"microtile_bases\":{},\"canonical_tile_bases\":{},\"analyzer_bytes\":{},\"selected_loci\":{},\"predicted_peak_rss_bytes\":{},\"science_blake3\":\"{}\"}}",
            plan.model_id, plan.baseline_rss_bytes, plan.fixed_bytes, plan.bytes_per_locus,
            dataset_plan.as_ref().map_or(plan.microtile_bases, |plan| plan.microtile_bases), plan.canonical_tile_bases, plan.analyzer_bytes,
            plan.selected_loci, predicted_peak, science_digest);
        return Ok(());
    }
    let mut primary = command
        .output
        .as_ref()
        .map(|path| AtomicFile::create(path))
        .transpose()?;
    let mut positional = command
        .options
        .position_output
        .as_ref()
        .map(|path| AtomicFile::create(path))
        .transpose()?;
    let setup_ms = started.elapsed().as_millis();
    let mut dataset_outcome = None;
    let mut drive = |engine: &mut EvidenceEngine,
                     analyzer: &mut dyn EvidenceAnalyzer|
     -> Result<EvidenceRunStats> {
        if let Some(options) = &dataset_options {
            let outcome = run_dataset_with_snapshot(
                engine,
                &cache_digest,
                options,
                analyzer,
                &input_snapshot,
            )
            .map_err(|error| match error {
                DatasetError::Evidence(error) => anyhow::Error::from(error),
                other => anyhow::Error::from(other),
            })?;
            let stats = outcome.stats.clone();
            dataset_outcome = Some(outcome);
            Ok(stats)
        } else {
            Ok(engine.run(analyzer)?)
        }
    };
    let result = {
        let stdout = io::stdout();
        let mut output: Box<dyn Write + '_> = match &mut primary {
            Some(file) => Box::new(io::BufWriter::new(file.file_mut())),
            None => Box::new(stdout.lock()),
        };
        let result = if let Some(panel) = &mut panel {
            let result = if let Some(file) = &mut positional {
                let mut writer = EvidenceArrowWriter::new(io::BufWriter::new(file.file_mut()));
                let result = drive(engine, &mut FusedAnalyzers::new(vec![panel, &mut writer]));
                if result.is_err() {
                    let _ = writer.finish();
                }
                result
            } else {
                drive(engine, panel)
            };
            if result.is_ok()
                || result.as_ref().err().is_some_and(|error| {
                    error
                        .downcast_ref::<EvidenceError>()
                        .is_some_and(resource_failure)
                })
            {
                panel.write_tsv(&mut output, engine.contigs())?;
            }
            result
        } else if command.options.format == FeatureFormat::ArrowIpc {
            let mut writer = EvidenceArrowWriter::new(&mut output);
            let result = drive(engine, &mut writer);
            if result.is_err() {
                let _ = writer.finish();
            }
            result
        } else {
            let mut writer = EvidenceTsvWriter::new(&mut output);
            drive(engine, &mut writer)
        };
        output.flush()?;
        result
    };
    input_snapshot.verify()?;
    let analysis_ms = started.elapsed().as_millis().saturating_sub(setup_ms);
    let (stats, mut failure) = match result {
        Ok(stats) => (stats, None),
        Err(error)
            if error
                .downcast_ref::<EvidenceError>()
                .is_some_and(resource_failure) =>
        {
            (
                EvidenceRunStats::default(),
                Some(error.downcast::<EvidenceError>().unwrap()),
            )
        }
        Err(error) => return Err(error),
    };
    let mut artifacts = Vec::new();
    for (flag, pending, requested) in [
        ("-o", primary.as_ref(), command.output.as_ref()),
        (
            "--position-output",
            positional.as_ref(),
            command.options.position_output.as_ref(),
        ),
    ] {
        if let (Some(file), Some(path)) = (pending, requested) {
            let hash = match hash_file(file.temporary_path()) {
                Ok(hash) => hash,
                Err(error)
                    if error
                        .downcast_ref::<EvidenceError>()
                        .is_some_and(resource_failure) =>
                {
                    if failure.is_none() {
                        failure = error.downcast::<EvidenceError>().ok();
                    }
                    rosalind::provenance::blake3_file(file.temporary_path())?
                }
                Err(error) => return Err(error),
            };
            artifacts.push((flag, path.clone(), hash));
        }
    }
    if failure.is_none() {
        if let Err(error) = checkpoint() {
            failure = Some(error.into());
        }
        if let Some(budget) = engine.request().execution.memory_budget_bytes {
            let peak = peak_rss_bytes();
            if peak > budget {
                failure = Some(
                    CoreError::BudgetExceeded {
                        needed: peak,
                        budget,
                    }
                    .into(),
                );
            }
        }
    }
    let mut failed = failure.is_some();
    for (flag, requested, hash) in &artifacts {
        let path = if failed {
            partial_path(requested)
        } else {
            requested.clone()
        };
        capture.output_hashed(flag, &path.display().to_string(), hash);
    }
    capture.opt("--mapq-threshold", engine.request().profile.min_mapq);
    capture.opt(
        "--base-quality-threshold",
        engine.request().profile.min_base_quality,
    );
    capture.opt("--max-read-len", command.max_read_len);
    capture.opt("--max-record-bytes", command.options.max_record_bytes);
    capture.opt("--tile-bases", command.options.tile_bases);
    capture.opt("--workers", command.options.workers);
    capture.opt(
        "--format",
        if command.options.format == FeatureFormat::ArrowIpc {
            "arrow-ipc"
        } else {
            "tsv"
        },
    );
    if command.panel {
        capture.opt("--min-callable-depth", command.options.min_callable_depth);
    }
    if let Some(mb) = command.memory_budget_mb {
        capture.opt("--memory-budget-mb", mb);
    }
    capture.flag_if(command.require_os_limit, "--require-os-limit");
    capture.flag_if(command.require_os_limit || command.enforce, "--enforce");
    let mut manifest = RunManifest::new(if command.panel {
        "analyze panel-qc"
    } else {
        "analyze evidence"
    });
    manifest.tool_version = env!("CARGO_PKG_VERSION").to_string();
    capture.record_into(&mut manifest);
    let analysis_identity = format!(
        "evidence={science_digest};analyzer={};schema=1;callable_depth={};panel_selection={}",
        if command.panel {
            "panel-qc"
        } else {
            "evidence"
        },
        if command.panel {
            command.options.min_callable_depth
        } else {
            0
        },
        if command.panel {
            manifest
                .inputs
                .iter()
                .find(|input| {
                    command.selection.regions.as_ref().is_some_and(|p| {
                        std::fs::canonicalize(p).is_ok_and(|path| path == Path::new(&input.path))
                    })
                })
                .map_or("", |input| input.blake3.as_str())
        } else {
            ""
        }
    );
    manifest.params.insert(
        "science.blake3".into(),
        blake3::hash(analysis_identity.as_bytes())
            .to_hex()
            .to_string(),
    );
    manifest.params.insert(
        "evidence.fields".into(),
        engine.request().fields.bits().to_string(),
    );
    manifest
        .params
        .insert("evidence.sampling".into(), "none".into());
    manifest.params.insert(
        "evidence.excluded_flags".into(),
        "unmapped,secondary,supplementary,duplicate,qc-failed".into(),
    );
    for (key, value) in [
        ("producer.name", "rosalind".to_string()),
        ("producer.version", env!("CARGO_PKG_VERSION").to_string()),
        ("producer.binary", "rosalind".to_string()),
        ("replay.kind", "rosalind".to_string()),
        ("evidence.schema", EVIDENCE_SCHEMA_VERSION.to_string()),
        ("evidence.profile", EvidenceProfile::ID.to_string()),
        ("evidence.counting_unit", "read".to_string()),
        ("evidence.science_blake3", science_digest),
        (
            "analyzer.id",
            if command.panel {
                "panel-qc"
            } else {
                "evidence"
            }
            .to_string(),
        ),
        ("analyzer.version", "1".to_string()),
        (
            "run_status",
            if failed {
                "resource-failed"
            } else {
                "completed"
            }
            .to_string(),
        ),
        (
            "contract.assurance",
            if os_limit.is_some() {
                "cgroup-v2"
            } else if command.memory_budget_mb.is_some() {
                "declared-bound-cooperative"
            } else {
                "observed-only"
            }
            .to_string(),
        ),
    ] {
        manifest.params.insert(key.to_string(), value);
    }
    if let Some(error) = &failure {
        manifest.params.insert(
            "failure.kind".into(),
            if matches!(error, EvidenceError::RecordLimit(_)) {
                "declared-capacity"
            } else {
                "runtime-memory"
            }
            .into(),
        );
        manifest
            .params
            .insert("failure.message".into(), error.to_string());
    }
    assign_artifact_roles(&mut manifest, &artifacts, &command, failed);
    let peak = peak_rss_bytes();
    manifest.params.insert(
        "predicted_working_set_bytes".into(),
        predicted_peak
            .saturating_sub(
                dataset_plan
                    .as_ref()
                    .map_or(plan.baseline_rss_bytes, |plan| plan.baseline_rss_bytes),
            )
            .to_string(),
    );
    if let Some(outcome) = dataset_outcome {
        manifest
            .params
            .insert("dataset.science_blake3".into(), outcome.science_digest);
        manifest.measurements.insert(
            "execution.reused_partitions".into(),
            outcome.reused_partitions.to_string(),
        );
        manifest.measurements.insert(
            "execution.computed_partitions".into(),
            outcome.computed_partitions.to_string(),
        );
        manifest.measurements.insert(
            "execution.worker_count".into(),
            outcome.plan.worker_count.to_string(),
        );
        manifest.measurements.insert(
            "execution.worker_bytes".into(),
            outcome.plan.worker_bytes.to_string(),
        );
        manifest.measurements.insert(
            "execution.reducer_bytes".into(),
            outcome.plan.reducer_bytes.to_string(),
        );
        manifest.measurements.insert(
            "execution.metadata_queue_bytes".into(),
            outcome.plan.metadata_queue_bytes.to_string(),
        );
        if temporary_cache.0.is_none() {
            manifest.measurements.insert(
                "execution.dataset_manifest".into(),
                outcome.dataset_manifest.display().to_string(),
            );
        }
    }
    for (key, value) in [
        ("baseline_rss_bytes", initial_rss),
        ("peak_rss_bytes", peak),
        ("predicted_peak_rss_bytes", predicted_peak),
        ("execution.record_visits", stats.record_visits),
        ("execution.microtiles", stats.microtiles),
        ("execution.emitted_loci", stats.emitted_loci),
        ("execution.analyzer_bytes", plan.analyzer_bytes),
        ("execution.hashing_ms", hashing_ms as u64),
        ("execution.setup_ms", setup_ms as u64),
        ("execution.analysis_encoding_ms", analysis_ms as u64),
        ("execution.elapsed_ms", started.elapsed().as_millis() as u64),
    ] {
        manifest.measurements.insert(key.into(), value.to_string());
    }
    manifest.params.insert(
        "contract_verdict".into(),
        if failed {
            "failed"
        } else if command.memory_budget_mb.is_some() {
            "within"
        } else {
            "unset"
        }
        .into(),
    );
    if let Some(limit) = os_limit {
        manifest
            .params
            .insert("os.memory_limit_bytes".into(), limit.to_string());
    }
    manifest.finalize();
    let mut receipt = receipt_path
        .as_ref()
        .map(|path| AtomicFile::create(path))
        .transpose()?;
    if let Some(file) = &mut receipt {
        file.file_mut()
            .write_all(manifest.to_canonical_json().as_bytes())?;
        file.file_mut().sync_all()?;
    }
    // Include receipt encoding and sync in the final resource observation. Any
    // newly detected failure changes output paths before the group is published.
    let final_peak = std::env::var("ROSALIND_FORCE_FINAL_RSS_BYTES")
        .ok()
        .and_then(|value| value.parse::<u64>().ok())
        .unwrap_or_else(peak_rss_bytes)
        .max(peak);
    if !failed {
        failure = checkpoint().err().map(EvidenceError::from);
        if let Some(budget) = engine.request().execution.memory_budget_bytes {
            if final_peak > budget {
                failure = Some(
                    CoreError::BudgetExceeded {
                        needed: final_peak,
                        budget,
                    }
                    .into(),
                );
            }
        }
        if let Some(error) = &failure {
            failed = true;
            manifest
                .params
                .insert("run_status".into(), "resource-failed".into());
            manifest
                .params
                .insert("failure.kind".into(), "runtime-memory".into());
            manifest
                .params
                .insert("failure.message".into(), error.to_string());
            manifest
                .measurements
                .insert("contract_verdict".into(), "failed".into());
            for (index, output) in manifest.outputs.iter_mut().enumerate() {
                output.path = partial_path(Path::new(&output.path)).display().to_string();
                manifest.params.insert(
                    format!("artifact.output.{index}.role"),
                    "partial-evidence".into(),
                );
            }
        }
    }
    manifest
        .measurements
        .insert("peak_rss_bytes".into(), final_peak.to_string());
    manifest.measurements.insert(
        "execution.elapsed_ms".into(),
        started.elapsed().as_millis().to_string(),
    );
    manifest.measurements.insert(
        "execution.finalization_ms".into(),
        started
            .elapsed()
            .as_millis()
            .saturating_sub(setup_ms + analysis_ms)
            .to_string(),
    );
    assign_artifact_roles(&mut manifest, &artifacts, &command, failed);
    manifest.finalize();
    if let Some(file) = &mut receipt {
        use std::io::{Seek, SeekFrom};
        file.file_mut().set_len(0)?;
        file.file_mut().seek(SeekFrom::Start(0))?;
        file.file_mut()
            .write_all(manifest.to_canonical_json().as_bytes())?;
        file.file_mut().sync_all()?;
    }
    let mut group = Vec::new();
    for (file, requested) in [
        (primary.take(), command.output.as_ref()),
        (positional.take(), command.options.position_output.as_ref()),
    ] {
        if let (Some(file), Some(path)) = (file, requested) {
            group.push((
                file,
                if failed {
                    partial_path(path)
                } else {
                    path.clone()
                },
            ));
        }
    }
    if let (Some(file), Some(path)) = (receipt, receipt_path.as_ref()) {
        group.push((file, path.clone()));
    }
    input_snapshot.verify()?;
    rosalind::util::atomic::commit_group(group, command.force)?;
    if let Some(error) = failure {
        return Err(error.into());
    }
    eprintln!(
        "evidence: {} loci; peak RSS {} MiB; {} record visits; receipt {}",
        stats.emitted_loci,
        final_peak / (1 << 20),
        stats.record_visits,
        receipt_path.map_or_else(
            || "not requested (stdout)".to_string(),
            |p| p.display().to_string()
        )
    );
    Ok(())
}

fn resource_failure(error: &EvidenceError) -> bool {
    matches!(
        error,
        EvidenceError::RecordLimit(_) | EvidenceError::Core(CoreError::BudgetExceeded { .. })
    )
}

fn partial_path(path: &Path) -> PathBuf {
    PathBuf::from(format!("{}.partial", path.display()))
}

// Receipt serialization orders outputs by path. Bind indexed metadata to that
// order, including the changed names after a late resource failure.
fn assign_artifact_roles(
    manifest: &mut RunManifest,
    artifacts: &[(&str, PathBuf, String)],
    command: &EvidenceCommand,
    failed: bool,
) {
    let mut ordered: Vec<_> = artifacts
        .iter()
        .map(|(flag, path, _)| {
            let final_path = if failed {
                partial_path(path)
            } else {
                path.clone()
            };
            (final_path, *flag)
        })
        .collect();
    ordered.sort_by(|a, b| a.0.cmp(&b.0));
    for (index, (_, flag)) in ordered.iter().enumerate() {
        let role = if failed {
            "partial-evidence"
        } else if command.panel && *flag == "-o" {
            "panel-summary"
        } else {
            "exact-evidence"
        };
        let format =
            if *flag == "--position-output" || command.options.format == FeatureFormat::ArrowIpc {
                "arrow-ipc"
            } else {
                "tsv"
            };
        manifest
            .params
            .insert(format!("artifact.output.{index}.role"), role.into());
        manifest
            .params
            .insert(format!("artifact.output.{index}.format"), format.into());
    }
}
