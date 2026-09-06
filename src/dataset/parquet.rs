//! Bounded, deterministic Parquet exports with transactional directory publication.

use super::persisted::{hash_budgeted, runtime_check};
use super::*;
use crate::evidence::{
    evidence_record_batch, EvidenceBatch, EvidenceRequirements, EVIDENCE_ARROW_BATCH_ROWS,
};
use crate::util::atomic_directory::AtomicDirectory;
use ::parquet::arrow::ArrowWriter;
use ::parquet::basic::Compression;
use ::parquet::file::properties::{EnabledStatistics, WriterProperties};
use arrow_array::RecordBatch;

/// Receipt filename within a successful Parquet export directory.
pub const PARQUET_EXPORT_MANIFEST_NAME: &str = "parquet-export.manifest.json";
/// Deterministic serialization recipe, independent of source partition/batch sizes.
pub const PARQUET_EXPORT_SEMANTICS: &str = "evidence-parquet-export-v1";

/// Fixed-envelope physical export controls; these never alter scientific rows.
#[derive(Debug, Clone, Copy)]
pub struct ParquetExportOptions {
    /// Maximum output rows per file, between 1 and 16,384. Default 16,384.
    pub rows_per_file: usize,
    /// Maximum derived receipt bytes. Default 32 MiB. File-inventory state and
    /// canonical serialization transients are admitted before export begins.
    pub max_receipt_bytes: usize,
}
impl Default for ParquetExportOptions {
    fn default() -> Self {
        Self {
            rows_per_file: CANONICAL_TILE_BASES as usize,
            max_receipt_bytes: 32 << 20,
        }
    }
}

/// A complete published export. File paths and receipt identities refer to the
/// final destination; the source inputs are verified dataset artifacts, never
/// falsely revalidated original BAM/reference provenance claims.
#[derive(Debug, Clone)]
pub struct ParquetExportOutcome {
    /// Atomically published directory.
    pub directory: PathBuf,
    /// Derived receipt in the published directory.
    pub manifest_path: PathBuf,
    /// Self-verified derived receipt, including output and source artifact hashes.
    pub manifest: RunManifest,
    /// Complete ordered Parquet files, including one schema-only file if empty.
    pub files: Vec<FileHash>,
    /// Whole-process source reader, encoder and finalization admission plan.
    pub plan: DatasetReadPlan,
    /// Query rows and verified source partition measurements.
    pub stats: EvidenceRunStats,
}

/// Plan source reading, projection, one Parquet row group, one bounded footer,
/// output inventory and derived receipt serialization without creating artifacts.
pub fn plan_parquet_export(
    dataset: &VerifiedEvidenceDataset,
    query: &DatasetQuery,
    options: &ParquetExportOptions,
    execution: &EvidenceExecution,
) -> Result<DatasetReadPlan, DatasetError> {
    plan_parquet_export_with_inputs(dataset, query, options, execution, None)
}

/// Include external query input identities in the same admission plan used by
/// `export_parquet_dataset_with_inputs`, including final receipt serialization.
pub fn plan_parquet_export_with_inputs(
    dataset: &VerifiedEvidenceDataset,
    query: &DatasetQuery,
    options: &ParquetExportOptions,
    execution: &EvidenceExecution,
    query_inputs: Option<&VerifiedInputSession>,
) -> Result<DatasetReadPlan, DatasetError> {
    let execution = reserve_query_inputs(dataset, query, options, execution, query_inputs)?;
    let execution = &execution;
    let bound = encoder_bound(dataset, query, options)?;
    let analyzer = ParquetConsumer {
        directory: None,
        fields: query.fields,
        rows_per_file: options.rows_per_file,
        bound,
        writer: None,
        current_rows: 0,
        files: Vec::new(),
        total_rows: 0,
        budget: dataset.execution_budget(execution),
    };
    dataset.plan(query, &analyzer, execution)
}

/// Export a query through a private sibling directory and publish only after all
/// data files, hashes, receipt bytes and input guards succeed. Existing output
/// directories, even empty ones, are refused. Unsorted input query selections are
/// normalized by the verified reader. Compression and dictionaries are disabled;
/// row groups never exceed 1,024 rows and each footer covers at most 16,384 rows.
pub fn export_parquet_dataset(
    dataset: &mut VerifiedEvidenceDataset,
    query: &DatasetQuery,
    destination: impl AsRef<Path>,
    options: &ParquetExportOptions,
    execution: &EvidenceExecution,
) -> Result<ParquetExportOutcome, DatasetError> {
    export_parquet_dataset_with_inputs(dataset, query, destination, options, execution, None)
}

/// Export with identities captured before parsing external VCF/BED selection
/// inputs. Their verified hashes participate in lineage and they are rechecked
/// after finalization, immediately before directory publication.
pub fn export_parquet_dataset_with_inputs(
    dataset: &mut VerifiedEvidenceDataset,
    query: &DatasetQuery,
    destination: impl AsRef<Path>,
    options: &ParquetExportOptions,
    execution: &EvidenceExecution,
    query_inputs: Option<&VerifiedInputSession>,
) -> Result<ParquetExportOutcome, DatasetError> {
    export_impl(
        dataset,
        query,
        destination.as_ref(),
        options,
        execution,
        query_inputs,
        || {},
    )
}

fn export_impl(
    dataset: &mut VerifiedEvidenceDataset,
    query: &DatasetQuery,
    destination: &Path,
    options: &ParquetExportOptions,
    execution: &EvidenceExecution,
    query_inputs: Option<&VerifiedInputSession>,
    before_publish: impl FnOnce(),
) -> Result<ParquetExportOutcome, DatasetError> {
    if let Some(inputs) = query_inputs {
        inputs.verify()?;
    }
    let execution = reserve_query_inputs(dataset, query, options, execution, query_inputs)?;
    let execution = &execution;
    let plan = plan_parquet_export(dataset, query, options, execution)?;
    let destination = absolute_destination(destination)?;
    let source_root = Path::new(&dataset.source_hashes()[0].path)
        .parent()
        .expect("verified manifest directory");
    if destination.starts_with(source_root) {
        return Err(DatasetError::Incompatible(
            "Parquet destination must be outside its source dataset directory".into(),
        ));
    }
    let stage = AtomicDirectory::create(&destination)?;
    let bound = encoder_bound(dataset, query, options)?;
    let mut consumer = ParquetConsumer {
        directory: Some(stage.path().to_path_buf()),
        fields: query.fields,
        rows_per_file: options.rows_per_file,
        bound,
        writer: None,
        current_rows: 0,
        files: Vec::new(),
        total_rows: 0,
        budget: dataset.execution_budget(execution),
    };
    let mut stats = dataset.visit_batches(query, &mut consumer, execution)?;
    if consumer.total_rows != plan.selected_loci {
        return Err(DatasetError::Corrupt(
            "Parquet export row count differs from the admitted query".into(),
        ));
    }
    let files: Vec<_> = consumer
        .files
        .iter()
        .map(|file| FileHash {
            path: destination.join(&file.path).display().to_string(),
            blake3: file.blake3.clone(),
        })
        .collect();
    let mut manifest = RunManifest::new("dataset export parquet");
    manifest.tool_version = env!("CARGO_PKG_VERSION").into();
    manifest.inputs = dataset.source_hashes().to_vec();
    manifest.inputs.extend(dataset.verified_partition_hashes()?);
    if let Some(inputs) = query_inputs {
        manifest
            .inputs
            .extend(inputs.identities().iter().map(|input| FileHash {
                path: input.path.clone(),
                blake3: input.blake3.clone(),
            }));
    }
    manifest.outputs = files.clone();
    manifest.params.extend([
        (
            "dataset.export_semantics".into(),
            PARQUET_EXPORT_SEMANTICS.into(),
        ),
        (
            "dataset.source_namespace".into(),
            dataset.descriptor().dataset_namespace.clone(),
        ),
        (
            "dataset.compatibility_blake3".into(),
            dataset.descriptor().compatibility_blake3.clone(),
        ),
        (
            "dataset.query_blake3".into(),
            dataset_query_digest(query, dataset.contigs())?,
        ),
        (
            "dataset.query".into(),
            canonical_dataset_query(query, dataset.contigs())?,
        ),
        ("replay.supported".into(), "false".into()),
        (
            "replay.reason".into(),
            "derived export has no executable replay recipe".into(),
        ),
        (
            "evidence.semantics".into(),
            dataset.descriptor().semantics.clone(),
        ),
        ("evidence.fields".into(), query.fields.bits().to_string()),
        (
            "evidence.fields_version".into(),
            query.fields.mask_version().to_string(),
        ),
        (
            "evidence.schema".into(),
            query.fields.schema_version().to_string(),
        ),
        (
            "evidence.source_fields".into(),
            dataset.fields().bits().to_string(),
        ),
        (
            "evidence.original_sources".into(),
            "provenance-claims-not-rehashed".into(),
        ),
        ("parquet.compression".into(), "uncompressed".into()),
        ("parquet.dictionary".into(), "false".into()),
        (
            "parquet.row_group_rows".into(),
            EVIDENCE_ARROW_BATCH_ROWS.to_string(),
        ),
        (
            "parquet.rows_per_file".into(),
            options.rows_per_file.to_string(),
        ),
        ("parquet.file_count".into(), files.len().to_string()),
        ("outcome.rows".into(), stats.emitted_loci.to_string()),
        ("run_status".into(), "completed".into()),
    ]);
    for (index, _) in files.iter().enumerate() {
        manifest.params.insert(
            format!("artifact.output.{index}.role"),
            "derived-evidence".into(),
        );
        manifest
            .params
            .insert(format!("artifact.output.{index}.format"), "parquet".into());
    }
    manifest.measurements.insert(
        "predicted_peak_rss_bytes".into(),
        plan.predicted_peak_rss_bytes.to_string(),
    );
    let query_digest = dataset_query_digest(query, dataset.contigs())?;
    let science = blake3::hash(
        format!(
            "{PARQUET_EXPORT_SEMANTICS}\n{}\n{query_digest}\n",
            dataset.descriptor().compatibility_blake3
        )
        .as_bytes(),
    )
    .to_hex()
    .to_string();
    manifest
        .params
        .insert("evidence.science_blake3".into(), science);
    manifest.measurements.insert(
        "peak_rss_bytes".into(),
        crate::util::rss::peak_rss_bytes().to_string(),
    );
    manifest.finalize();
    let receipt = manifest.to_canonical_json();
    if receipt.len() > options.max_receipt_bytes {
        return Err(DatasetError::Incompatible(
            "Parquet receipt exceeds its admitted metadata envelope".into(),
        ));
    }
    if manifest.self_hash_ok() != Some(true) || manifest.measurement_hash_ok() == Some(false) {
        return Err(DatasetError::Corrupt(
            "derived Parquet receipt failed its self-check".into(),
        ));
    }
    let receipt_path = stage.path().join(PARQUET_EXPORT_MANIFEST_NAME);
    let mut file = File::create(&receipt_path)?;
    file.write_all(receipt.as_bytes())?;
    file.sync_all()?;
    drop(file);
    drop(receipt);
    runtime_check(dataset.execution_budget(execution))?;
    // Observe receipt encoding and sync, then persist the final observation.
    stats.peak_rss_bytes = crate::util::rss::peak_rss_bytes();
    manifest
        .measurements
        .insert("peak_rss_bytes".into(), stats.peak_rss_bytes.to_string());
    manifest.finalize();
    let receipt = manifest.to_canonical_json();
    if receipt.len() > options.max_receipt_bytes {
        return Err(DatasetError::Incompatible(
            "Parquet receipt exceeds its admitted metadata envelope".into(),
        ));
    }
    let mut file = File::create(&receipt_path)?;
    file.write_all(receipt.as_bytes())?;
    file.sync_all()?;
    drop(file);
    drop(receipt);
    before_publish();
    dataset.verify_unchanged()?;
    if let Some(inputs) = query_inputs {
        inputs.verify()?;
    }
    runtime_check(dataset.execution_budget(execution))?;
    let directory = stage.commit()?;
    Ok(ParquetExportOutcome {
        manifest_path: directory.join(PARQUET_EXPORT_MANIFEST_NAME),
        directory,
        manifest,
        files,
        plan,
        stats,
    })
}

struct ParquetConsumer {
    directory: Option<PathBuf>,
    fields: EvidenceFields,
    rows_per_file: usize,
    bound: u64,
    writer: Option<ArrowWriter<File>>,
    current_rows: usize,
    files: Vec<FileHash>,
    total_rows: u64,
    budget: Option<u64>,
}
impl ParquetConsumer {
    fn open(&mut self, batch: &RecordBatch) -> Result<(), EvidenceError> {
        let properties = WriterProperties::builder()
            .set_created_by("Rosalind evidence-parquet-export-v1 parquet-55.2.0".into())
            .set_compression(Compression::UNCOMPRESSED)
            .set_dictionary_enabled(false)
            .set_statistics_enabled(EnabledStatistics::Chunk)
            .set_max_row_group_size(EVIDENCE_ARROW_BATCH_ROWS)
            .set_write_batch_size(EVIDENCE_ARROW_BATCH_ROWS)
            .build();
        let path = self
            .directory
            .as_ref()
            .expect("export directory exists")
            .join(part_name(self.files.len()));
        let file = fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(path)?;
        self.writer = Some(
            ArrowWriter::try_new(file, batch.schema(), Some(properties)).map_err(parquet_error)?,
        );
        self.current_rows = 0;
        Ok(())
    }
    fn finish_file(&mut self) -> Result<(), EvidenceError> {
        if let Some(writer) = self.writer.take() {
            // into_inner writes the footer and propagates write/flush failures.
            let file = writer.into_inner().map_err(parquet_error)?;
            runtime_check(self.budget)?;
            file.sync_all()?;
            drop(file);
            let name = part_name(self.files.len());
            let path = self
                .directory
                .as_ref()
                .expect("export directory exists")
                .join(&name);
            self.files.push(FileHash {
                path: name,
                blake3: hash_budgeted(&path, self.budget)?,
            });
            self.current_rows = 0;
        }
        Ok(())
    }
}
impl EvidenceAnalyzer for ParquetConsumer {
    fn requirements(&self) -> EvidenceRequirements {
        EvidenceRequirements {
            fields: self.fields,
            requires_reference: false,
            context_bases: 0,
            retained_bytes: Some(self.bound),
        }
    }
    fn on_batch(&mut self, batch: &EvidenceBatch) -> Result<(), EvidenceError> {
        runtime_check(self.budget)?;
        let batch = evidence_record_batch(batch)?;
        let mut offset = 0;
        while offset < batch.num_rows() {
            if self.writer.is_none() {
                self.open(&batch)?;
            }
            let count = (batch.num_rows() - offset).min(self.rows_per_file - self.current_rows);
            self.writer
                .as_mut()
                .expect("writer opened")
                .write(&batch.slice(offset, count))
                .map_err(parquet_error)?;
            self.current_rows += count;
            self.total_rows = self
                .total_rows
                .checked_add(count as u64)
                .ok_or(EvidenceError::CounterOverflow)?;
            offset += count;
            if self.current_rows == self.rows_per_file {
                self.finish_file()?;
            }
            runtime_check(self.budget)?;
        }
        Ok(())
    }
    fn finish(&mut self) -> Result<(), EvidenceError> {
        if self.writer.is_none() && self.files.is_empty() {
            // Empty queries still carry their physical schema for downstream readers.
            let batch =
                evidence_record_batch(&EvidenceBatch::new(0, "", 0, self.fields, Vec::new()))?;
            self.open(&batch)?;
        }
        self.finish_file()?;
        runtime_check(self.budget)?;
        Ok(())
    }
    fn additional_memory_bytes(&self) -> Option<u64> {
        Some(self.bound)
    }
}

fn reserve_query_inputs(
    dataset: &VerifiedEvidenceDataset,
    query: &DatasetQuery,
    options: &ParquetExportOptions,
    execution: &EvidenceExecution,
    query_inputs: Option<&VerifiedInputSession>,
) -> Result<EvidenceExecution, DatasetError> {
    if let Some(inputs) = query_inputs {
        inputs.verify()?;
    }
    let mut execution = execution.clone();
    let bound = encoder_bound(dataset, query, options)?;
    let identity_bytes = query_inputs.map_or(0, |inputs| {
        inputs.identities().iter().fold(0u64, |sum, input| {
            sum.saturating_add(input.path.len() as u64 + input.role.len() as u64 + 256)
        })
    });
    execution.analyzer_bytes = execution
        .analyzer_bytes
        .max(bound.saturating_add(identity_bytes.saturating_mul(16)));
    Ok(execution)
}

fn encoder_bound(
    dataset: &VerifiedEvidenceDataset,
    query: &DatasetQuery,
    options: &ParquetExportOptions,
) -> Result<u64, DatasetError> {
    if options.rows_per_file == 0
        || options.rows_per_file > CANONICAL_TILE_BASES as usize
        || options.max_receipt_bytes == 0
    {
        return Err(DatasetError::Incompatible(
            "Parquet rows_per_file must be in 1..=16384 and metadata envelope must be positive"
                .into(),
        ));
    }
    let (intervals, sites) = query.selection.normalize(dataset.contigs())?;
    let rows: u64 = intervals
        .iter()
        .map(|interval| u64::from(interval.end - interval.start))
        .sum();
    let files = rows.div_ceil(options.rows_per_file as u64).max(1);
    let source_root_bytes = dataset.source_hashes()[0].path.len() as u64;
    let inputs = dataset
        .descriptor()
        .partitions
        .iter()
        .fold(0u64, |sum, part| {
            sum.saturating_add(
                2 * source_root_bytes
                    + part.arrow.path.len() as u64
                    + part.receipt.path.len() as u64
                    + 512,
            )
        });
    // Final absolute output paths are capped at 4KiB; indexed artifact roles and
    // duplicate inventories/canonical receipt buffers are included explicitly.
    let receipt_bytes = (64u64 << 10)
        .saturating_add(inputs)
        .saturating_add((intervals.len() as u64).saturating_mul(128))
        .saturating_add((sites.len() as u64).saturating_mul(256))
        .saturating_add(files.saturating_mul(4608));
    if receipt_bytes > options.max_receipt_bytes as u64 {
        return Err(DatasetError::Incompatible("Parquet file inventory exceeds the declared receipt envelope; increase rows_per_file or max_receipt_bytes".into()));
    }
    let max_contig = dataset
        .contigs()
        .iter()
        .map(|c| c.name.len() as u64)
        .max()
        .unwrap_or(0);
    let row_group_bytes = query
        .fields
        .storage_bytes_per_locus()
        .saturating_add(max_contig)
        .saturating_add(128)
        .saturating_mul(EVIDENCE_ARROW_BATCH_ROWS as u64);
    // Independent Arrow arrays, Parquet encoders and serialized pages coexist;
    // fixed space covers at most 16 row-group footers, levels and schema state.
    Ok((32u64 << 20)
        .saturating_add(row_group_bytes.saturating_mul(6))
        .saturating_add(receipt_bytes.saturating_mul(16)))
}
fn absolute_destination(path: &Path) -> Result<PathBuf, DatasetError> {
    let name = path.file_name().ok_or_else(|| {
        DatasetError::Incompatible("Parquet destination requires a new directory name".into())
    })?;
    let parent = path
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let path = fs::canonicalize(parent)?.join(name);
    if path.as_os_str().len() > 4096 || path.to_str().is_none() {
        return Err(DatasetError::Incompatible(
            "Parquet destination exceeds its UTF-8 path envelope".into(),
        ));
    }
    Ok(path)
}
fn parquet_error(error: ::parquet::errors::ParquetError) -> EvidenceError {
    EvidenceError::Analyzer(format!("Parquet encoding failed: {error}"))
}
fn part_name(index: usize) -> String {
    format!("part-{index:08}.parquet")
}

#[cfg(test)]
mod tests {
    use super::super::persisted::tests::Fixture;
    use super::*;
    use ::parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
    use arrow_array::{Array, FixedSizeListArray, UInt32Array, UInt64Array};
    use arrow_schema::DataType;

    #[test]
    fn roundtrip_preserves_unsigned_scalars_fixed_lists_and_bounded_deterministic_files() {
        let fixture = Fixture::new(false);
        let mut dataset = fixture.open();
        let query = DatasetQuery {
            selection: dataset.selection(),
            fields: EvidenceFields::ALL_SUPPORTED,
        };
        let options = ParquetExportOptions {
            rows_per_file: 1500,
            ..Default::default()
        };
        let output = export_parquet_dataset(
            &mut dataset,
            &query,
            fixture.root.join("parquet-a"),
            &options,
            &EvidenceExecution::default(),
        )
        .unwrap();
        assert_eq!(output.stats.emitted_loci, 2305);
        assert_eq!(output.files.len(), 2);
        assert_eq!(output.manifest.self_hash_ok(), Some(true));
        assert_eq!(output.manifest.params["replay.supported"], "false");
        assert!(output.manifest.params["dataset.query"].contains("selection"));
        let mut rows = 0usize;
        let mut nonzero = Vec::new();
        for file in &output.files {
            assert_eq!(blake3_file(Path::new(&file.path)).unwrap(), file.blake3);
            let builder =
                ParquetRecordBatchReaderBuilder::try_new(File::open(&file.path).unwrap()).unwrap();
            assert!(builder.metadata().file_metadata().num_rows() <= 1500);
            for group in builder.metadata().row_groups() {
                assert!(group.num_rows() <= 1024);
            }
            assert_eq!(
                builder
                    .schema()
                    .field_with_name("callable_depth")
                    .unwrap()
                    .data_type(),
                &DataType::UInt64
            );
            assert!(
                matches!(builder.schema().field_with_name("base_quality_histogram").unwrap().data_type(),DataType::FixedSizeList(field,94) if field.data_type()==&DataType::UInt64)
            );
            for batch in builder.with_batch_size(1024).build().unwrap() {
                let batch = batch.unwrap();
                rows += batch.num_rows();
                let pos = batch
                    .column_by_name("pos")
                    .unwrap()
                    .as_any()
                    .downcast_ref::<UInt32Array>()
                    .unwrap();
                let depths = batch
                    .column_by_name("callable_depth")
                    .unwrap()
                    .as_any()
                    .downcast_ref::<UInt64Array>()
                    .unwrap();
                let hist = batch
                    .column_by_name("base_quality_histogram")
                    .unwrap()
                    .as_any()
                    .downcast_ref::<FixedSizeListArray>()
                    .unwrap();
                for row in 0..batch.num_rows() {
                    if depths.value(row) > 0 {
                        let values = hist.value(row);
                        let values = values.as_any().downcast_ref::<UInt64Array>().unwrap();
                        assert_eq!(values.values().iter().sum::<u64>(), 1);
                        nonzero.push(pos.value(row));
                    }
                }
            }
        }
        assert_eq!(rows, 2305);
        assert_eq!(nonzero, vec![2, 1026, CANONICAL_TILE_BASES + 2]);
        let second = export_parquet_dataset(
            &mut dataset,
            &query,
            fixture.root.join("parquet-b"),
            &options,
            &EvidenceExecution::default(),
        )
        .unwrap();
        assert_eq!(
            output.files.iter().map(|f| &f.blake3).collect::<Vec<_>>(),
            second.files.iter().map(|f| &f.blake3).collect::<Vec<_>>()
        );
    }

    #[test]
    fn empty_projection_export_is_schema_bearing_and_existing_destinations_are_preserved() {
        let fixture = Fixture::new(true);
        let mut dataset = fixture.open();
        let query = DatasetQuery {
            selection: EvidenceSelection::Intervals(Vec::new()),
            fields: EvidenceFields::DEPTHS,
        };
        let target = fixture.root.join("empty");
        let output = export_parquet_dataset(
            &mut dataset,
            &query,
            &target,
            &ParquetExportOptions::default(),
            &EvidenceExecution::default(),
        )
        .unwrap();
        assert_eq!(output.files.len(), 1);
        assert_eq!(output.stats.emitted_loci, 0);
        let builder =
            ParquetRecordBatchReaderBuilder::try_new(File::open(&output.files[0].path).unwrap())
                .unwrap();
        assert_eq!(builder.metadata().file_metadata().num_rows(), 0);
        assert!(builder.schema().field_with_name("callable_depth").is_ok());
        assert!(builder
            .schema()
            .field_with_name("base_quality_sum")
            .is_err());
        assert!(export_parquet_dataset(
            &mut dataset,
            &query,
            &target,
            &ParquetExportOptions::default(),
            &EvidenceExecution::default()
        )
        .is_err());
        assert_eq!(
            blake3_file(Path::new(&output.files[0].path)).unwrap(),
            output.files[0].blake3
        );
    }

    #[test]
    fn corruption_late_query_mutation_and_nested_destinations_never_publish_success() {
        let fixture = Fixture::new(true);
        let mut dataset = fixture.open();
        let query = DatasetQuery {
            selection: dataset.selection(),
            fields: EvidenceFields::DEPTHS,
        };
        let root = fixture.manifest.parent().unwrap();
        assert!(export_parquet_dataset(
            &mut dataset,
            &query,
            root.join("nested"),
            &ParquetExportOptions::default(),
            &EvidenceExecution::default()
        )
        .is_err());
        assert!(!root.join("nested").exists());
        let input = fixture.root.join("query.bed");
        fs::write(&input, b"chr1\t1\t2\n").unwrap();
        let session = VerifiedInputSession::open([("regions".into(), input.clone())]).unwrap();
        let output = fixture.root.join("late-mutation");
        let result = export_impl(
            &mut dataset,
            &query,
            &output,
            &ParquetExportOptions::default(),
            &EvidenceExecution::default(),
            Some(&session),
            || fs::write(&input, b"chr1\t1\t3\n").unwrap(),
        );
        assert!(result.is_err());
        assert!(!output.exists());
        let artifact = root.join(&dataset.descriptor().partitions[0].arrow.path);
        fs::write(&artifact, b"corrupt").unwrap();
        let target = fixture.root.join("corrupt");
        let mut dataset = fixture.open();
        assert!(export_parquet_dataset(
            &mut dataset,
            &query,
            &target,
            &ParquetExportOptions::default(),
            &EvidenceExecution::default()
        )
        .is_err());
        assert!(!target.exists());
        assert!(!fs::read_dir(&fixture.root).unwrap().any(|entry| entry
            .unwrap()
            .file_name()
            .to_string_lossy()
            .starts_with(".rosalind-dataset-")));
    }

    #[test]
    fn encoder_and_receipt_inventory_are_admitted_before_creating_directory() {
        let fixture = Fixture::new(false);
        let dataset = fixture.open();
        let query = DatasetQuery {
            selection: dataset.selection(),
            fields: EvidenceFields::ALL_SUPPORTED,
        };
        let plan = plan_parquet_export(
            &dataset,
            &query,
            &ParquetExportOptions::default(),
            &EvidenceExecution::default(),
        )
        .unwrap();
        assert!(plan.analyzer_bytes > 32 << 20);
        let execution = EvidenceExecution {
            memory_budget_bytes: Some(plan.predicted_peak_rss_bytes - 1),
            ..Default::default()
        };
        assert!(matches!(
            plan_parquet_export(
                &dataset,
                &query,
                &ParquetExportOptions::default(),
                &execution
            ),
            Err(DatasetError::Evidence(EvidenceError::Refused { .. }))
        ));
        assert!(plan_parquet_export(
            &dataset,
            &query,
            &ParquetExportOptions {
                max_receipt_bytes: 1,
                ..Default::default()
            },
            &EvidenceExecution::default()
        )
        .is_err());
    }

    #[test]
    fn parquet_preserves_unsigned_values_above_signed_integer_range() {
        let fixture = Fixture::new(true);
        let directory = fixture.root.join("unsigned-encoding");
        fs::create_dir(&directory).unwrap();
        let mut batch = EvidenceBatch::new(
            0,
            "chr1",
            0,
            EvidenceFields::ALLELES,
            vec![crate::evidence::EvidenceLocus {
                position: 0,
                reference: b'A',
                requested_alts: Vec::new(),
            }],
        );
        batch.row_mut(0).unwrap().alleles.unwrap().allele_counts[0] = u64::MAX;
        let mut consumer = ParquetConsumer {
            directory: Some(directory.clone()),
            fields: EvidenceFields::ALLELES,
            rows_per_file: 16384,
            bound: 64 << 20,
            writer: None,
            current_rows: 0,
            files: Vec::new(),
            total_rows: 0,
            budget: None,
        };
        consumer.on_batch(&batch).unwrap();
        consumer.finish().unwrap();
        let mut reader = ParquetRecordBatchReaderBuilder::try_new(
            File::open(directory.join(&consumer.files[0].path)).unwrap(),
        )
        .unwrap()
        .build()
        .unwrap();
        let batch = reader.next().unwrap().unwrap();
        let counts = batch
            .column_by_name("a")
            .unwrap()
            .as_any()
            .downcast_ref::<UInt64Array>()
            .unwrap();
        assert_eq!(counts.value(0), u64::MAX);
    }
}
