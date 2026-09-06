use super::*;
use rosalind::dataset::{
    plan_dataset, publication_memory_bytes, publish_evidence_dataset, run_dataset_with_snapshot,
    DatasetError, DatasetOptions, DescriptorLimits,
};

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
    input_session: VerifiedInputSession,
) -> Result<()> {
    input_session.verify()?;
    let mut reuse_source = command
        .options
        .reuse_dataset
        .as_ref()
        .map(|path| {
            rosalind::dataset::VerifiedEvidenceDataset::open(
                path,
                rosalind::dataset::DatasetReadLimits {
                    memory_budget_bytes: engine.request().execution.memory_budget_bytes,
                    max_manifest_bytes: command.options.max_dataset_metadata_bytes,
                    max_descriptor_bytes: command.options.max_dataset_metadata_bytes,
                    ..Default::default()
                },
            )
        })
        .transpose()?;
    if let Some(source) = &reuse_source {
        let parent = &source.source_hashes()[0];
        capture.input_hashed("--reuse-dataset", &parent.path, &parent.blake3);
    }
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
    let fields = engine.request().fields;
    let arrow_memory = EvidenceArrowWriter::with_fields(io::sink(), fields)
        .additional_memory_bytes()
        .unwrap();
    let tsv_memory = EvidenceTsvWriter::with_fields(io::sink(), fields)
        .additional_memory_bytes()
        .unwrap();
    let annotation_bytes = if command.options.annotated_variants.is_some() {
        rosalind::variant_annotation::annotation_memory_bytes(
            engine,
            command.options.variant_limits(),
        )
    } else {
        0
    };
    let encoder_bytes = if command.options.position_output.is_some()
        || (!command.panel && command.options.format == FeatureFormat::ArrowIpc)
    {
        arrow_memory
    } else {
        tsv_memory
    };
    let publication_bytes = if command.options.cache_dir.is_some() {
        publication_memory_bytes(
            engine,
            DescriptorLimits {
                max_bytes: command.options.max_dataset_metadata_bytes,
            },
        )?
    } else {
        0
    };
    let consumer_bytes = annotation_bytes
        .checked_add(publication_bytes)
        .ok_or_else(|| {
            EvidenceError::InvalidRequest("dataset publication memory envelope is too large".into())
        })?
        .checked_add(
            panel
                .as_ref()
                .and_then(EvidenceAnalyzer::additional_memory_bytes)
                .unwrap_or(0),
        )
        .and_then(|bytes| bytes.checked_add(encoder_bytes))
        .ok_or_else(|| {
            EvidenceError::InvalidRequest("annotation/consumer memory envelope is too large".into())
        })?;
    let planner = EvidenceCallback::with_fields(|_: &EvidenceBatch| Ok(()), consumer_bytes, fields);
    let plan = engine.plan_for_analyzer(&planner)?.clone();
    let reuse_options = rosalind::dataset::ReuseOptions {
        workers: command.options.workers,
    };
    let reuse_plan = reuse_source
        .as_ref()
        .map(|source| {
            rosalind::dataset::plan_reuse(engine, &input_session, source, &planner, &reuse_options)
        })
        .transpose()?;
    let temporary_cache = TemporaryCache(
        if command.options.cache_dir.is_none()
            && (command.options.workers > 1 || command.options.annotated_variants.is_some())
        {
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
    let cache_digest = input_session.dataset_namespace(engine)?;
    let dataset_plan = dataset_options
        .as_ref()
        .map(|options| plan_dataset(engine, options, &planner))
        .transpose()?;
    let predicted_peak = dataset_plan
        .as_ref()
        .map_or(plan.predicted_peak_rss_bytes, |plan| {
            plan.predicted_peak_rss_bytes
        });
    let predicted_peak = reuse_plan
        .as_ref()
        .map_or(predicted_peak, |p| p.predicted_peak_rss_bytes);
    if command.options.plan {
        input_session.verify()?;
        if let Some(p) = &reuse_plan {
            println!(
                "{}",
                serde_json::json!({"model":p.model_id,"reused_loci":p.reused_loci,"computed_loci":p.computed_loci,"metadata_bytes":p.metadata_bytes,"source_decoder_bytes":p.source_decoder_bytes,"merge_bytes":p.merge_bytes,"analyzer_bytes":p.analyzer_bytes,"native_bytes":p.native_bytes,"fields":p.output_fields.bits(),"source_fields":p.source_fields.bits(),"predicted_peak_rss_bytes":p.predicted_peak_rss_bytes,"science_blake3":science_digest})
            );
            return Ok(());
        }
        println!("{{\"model\":\"{}\",\"baseline_rss_bytes\":{},\"fixed_bytes\":{},\"bytes_per_locus\":{},\"microtile_bases\":{},\"canonical_tile_bases\":{},\"analyzer_bytes\":{},\"selected_loci\":{},\"fields\":{},\"schema\":{},\"predicted_peak_rss_bytes\":{},\"science_blake3\":\"{}\"}}",
            plan.model_id, plan.baseline_rss_bytes, plan.fixed_bytes, plan.bytes_per_locus,
            dataset_plan.as_ref().map_or(plan.microtile_bases, |plan| plan.microtile_bases), plan.canonical_tile_bases, plan.analyzer_bytes,
            plan.selected_loci, fields.bits(), fields.schema_version(), predicted_peak, science_digest);
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
    let mut annotated = command
        .options
        .annotated_variants
        .as_ref()
        .map(|path| AtomicFile::create(path))
        .transpose()?;
    let setup_ms = started.elapsed().as_millis();
    let mut dataset_outcome = None;
    let mut reuse_outcome = None;
    let mut drive = |engine: &mut EvidenceEngine,
                     analyzer: &mut dyn EvidenceAnalyzer|
     -> Result<EvidenceRunStats> {
        if let Some(source) = &mut reuse_source {
            let outcome = rosalind::dataset::run_reusing_dataset(
                engine,
                &input_session,
                source,
                analyzer,
                &reuse_options,
            )
            .map_err(|error| match error {
                DatasetError::Evidence(error) => anyhow::Error::from(error),
                other => anyhow::Error::from(other),
            })?;
            let stats = outcome.stats.clone();
            reuse_outcome = Some(outcome);
            Ok(stats)
        } else if let Some(options) = &dataset_options {
            let outcome = run_dataset_with_snapshot(
                engine,
                &cache_digest,
                options,
                analyzer,
                input_session.snapshot(),
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
    let mut result = {
        let stdout = io::stdout();
        let mut output: Box<dyn Write + '_> = match &mut primary {
            Some(file) => Box::new(io::BufWriter::new(file.file_mut())),
            None => Box::new(stdout.lock()),
        };
        let result = if let Some(panel) = &mut panel {
            let result = if let Some(file) = &mut positional {
                let mut writer =
                    EvidenceArrowWriter::with_fields(io::BufWriter::new(file.file_mut()), fields);
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
            let mut writer = EvidenceArrowWriter::with_fields(&mut output, fields);
            let result = drive(engine, &mut writer);
            if result.is_err() {
                let _ = writer.finish();
            }
            result
        } else {
            let mut writer = EvidenceTsvWriter::with_fields(&mut output, fields);
            drive(engine, &mut writer)
        };
        output.flush()?;
        result
    };
    let mut persisted_manifest = None;
    if result.is_ok() && command.options.cache_dir.is_some() {
        let published = publish_evidence_dataset(
            engine,
            dataset_outcome.as_ref().unwrap(),
            &input_session,
            DescriptorLimits {
                max_bytes: command.options.max_dataset_metadata_bytes,
            },
        );
        match published {
            Ok(path) => persisted_manifest = Some(path),
            Err(DatasetError::Evidence(error)) => result = Err(error.into()),
            Err(error) => result = Err(error.into()),
        }
    }
    let mut annotated_records = None;
    if result.is_ok() {
        if let (Some(file), Some(requested)) = (&mut annotated, &command.options.annotated_variants)
        {
            file.close_for_path_writer();
            let annotation_result = rosalind::variant_annotation::annotate_variants(
                command.options.sites.as_ref().unwrap(),
                file.temporary_path(),
                rosalind::variant_annotation::annotation_format(requested)?,
                command.options.variant_limits(),
                engine,
                dataset_outcome.as_ref().unwrap(),
                &science_digest,
            )
            .map_err(|error| match error {
                DatasetError::Evidence(error) => anyhow::Error::from(error),
                other => anyhow::Error::from(other),
            });
            match annotation_result {
                Ok(count) => annotated_records = Some(count),
                Err(error) => result = Err(error),
            }
        }
    } else {
        // No annotation bytes exist if extraction failed. Drop its staging file.
        annotated = None;
    }
    input_session.verify()?;
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
        (
            "--annotated-variants",
            annotated.as_ref(),
            command.options.annotated_variants.as_ref(),
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
    if let Some(source) = &reuse_source {
        source.verify_unchanged()?;
        let dependencies = source
            .source_hashes()
            .iter()
            .cloned()
            .chain(source.verified_partition_hashes()?)
            .collect::<Vec<_>>();
        for path in &command.options.reuse_artifacts {
            let path = std::fs::canonicalize(path)?;
            if !dependencies.iter().any(|s| Path::new(&s.path) == path) {
                let hash = hash_file(&path)?;
                if !dependencies.iter().any(|s| s.blake3 == hash) {
                    bail!("--reuse-artifact must match a verified source dependency");
                }
            }
        }
        for file in dependencies.iter().skip(1) {
            capture.input_hashed("--reuse-artifact", &file.path, &file.blake3);
        }
    }
    capture.opt("--mapq-threshold", engine.request().profile.min_mapq);
    capture.opt(
        "--base-quality-threshold",
        engine.request().profile.min_base_quality,
    );
    capture.opt("--max-read-len", command.max_read_len);
    capture.opt("--max-record-bytes", command.options.max_record_bytes);
    if command.options.sites.is_some() {
        capture.opt(
            "--max-variant-record-bytes",
            command.options.max_variant_record_bytes,
        );
        capture.opt(
            "--max-variant-header-bytes",
            command.options.max_variant_header_bytes,
        );
    }
    capture.opt("--tile-bases", command.options.tile_bases);
    capture.opt("--workers", command.options.workers);
    if command.options.cache_dir.is_some() || command.options.reuse_dataset.is_some() {
        capture.opt(
            "--max-dataset-metadata-bytes",
            command.options.max_dataset_metadata_bytes,
        );
    }
    let field_names = fields.names().join(",");
    capture.opt(
        "--fields",
        if field_names.is_empty() {
            "none"
        } else {
            &field_names
        },
    );
    if let Some(sample) = &command.options.sample {
        capture.opt("--sample", sample);
    }
    capture.flag_if(command.options.pool_samples, "--pool-samples");
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
    let mut analysis_identity = format!(
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
    if command.options.annotated_variants.is_some() {
        let sites_path = std::fs::canonicalize(command.options.sites.as_ref().unwrap())?;
        let sites_hash = &manifest
            .inputs
            .iter()
            .find(|input| Path::new(&input.path) == sites_path)
            .ok_or_else(|| anyhow::anyhow!("annotation input identity missing"))?
            .blake3;
        analysis_identity.push_str(&format!(
            ";annotation={};variant_bytes={sites_hash}",
            rosalind::variant_annotation::ANNOTATION_SEMANTICS
        ));
        manifest.params.insert(
            "annotation.semantics".into(),
            rosalind::variant_annotation::ANNOTATION_SEMANTICS.into(),
        );
        manifest
            .params
            .insert("annotation.variant_blake3".into(), sites_hash.clone());
        manifest.measurements.insert(
            "execution.annotation_bytes".into(),
            annotation_bytes.to_string(),
        );
        if let Some(count) = annotated_records {
            manifest
                .measurements
                .insert("execution.annotated_records".into(), count.to_string());
        }
    }
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
        ("evidence.schema", fields.schema_version().to_string()),
        ("evidence.fields_version", fields.mask_version().to_string()),
        ("evidence.profile", EvidenceProfile::ID.to_string()),
        ("evidence.counting_unit", "read".to_string()),
        (
            "evidence.sample_scope",
            engine.sample_scope().canonical_json(),
        ),
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
    if let Some(path) = persisted_manifest {
        manifest.measurements.insert(
            "execution.evidence_dataset_manifest".into(),
            path.display().to_string(),
        );
        manifest.measurements.insert(
            "execution.dataset_publication_bytes".into(),
            publication_bytes.to_string(),
        );
    }
    if let Some(outcome) = reuse_outcome {
        manifest.measurements.insert(
            "execution.reused_loci".into(),
            outcome.reused_loci.to_string(),
        );
        manifest.measurements.insert(
            "execution.computed_loci".into(),
            outcome.computed_loci.to_string(),
        );
        manifest.measurements.insert(
            "execution.reuse_model".into(),
            reuse_plan.as_ref().unwrap().model_id.into(),
        );
    }
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
        (
            "execution.sample_filtered_record_visits",
            stats.sample_filtered_record_visits,
        ),
        ("execution.sample_scope_bytes", plan.sample_scope_bytes),
        ("execution.microtiles", stats.microtiles),
        (
            "execution.microtile_bases",
            u64::from(
                dataset_plan
                    .as_ref()
                    .map_or(plan.microtile_bases, |value| value.microtile_bases),
            ),
        ),
        ("execution.bytes_per_locus", plan.bytes_per_locus),
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
    if let Some(source) = &reuse_source {
        source.verify_unchanged()?;
    }
    let mut group = Vec::new();
    for (file, requested) in [
        (primary.take(), command.output.as_ref()),
        (positional.take(), command.options.position_output.as_ref()),
        (
            annotated.take(),
            command.options.annotated_variants.as_ref(),
        ),
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
    input_session.verify()?;
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
        } else if *flag == "--annotated-variants" {
            "annotated-variants"
        } else if command.panel && *flag == "-o" {
            "panel-summary"
        } else {
            "exact-evidence"
        };
        let format = if *flag == "--annotated-variants" {
            match rosalind::variant_annotation::annotation_format(
                command.options.annotated_variants.as_ref().unwrap(),
            )
            .expect("validated annotation format")
            {
                rosalind::variant_io::VariantFormat::Vcf => "vcf",
                rosalind::variant_io::VariantFormat::VcfGz => "vcf.gz",
                rosalind::variant_io::VariantFormat::Bcf => "bcf",
            }
        } else if *flag == "--position-output" || command.options.format == FeatureFormat::ArrowIpc
        {
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
