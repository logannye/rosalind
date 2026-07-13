//! Receipt-driven canonical merge for first-party shard artifacts.

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use arrow_array::RecordBatch;
use arrow_ipc::reader::StreamReader;
use arrow_ipc::writer::{IpcWriteOptions, StreamWriter};
use arrow_ipc::MetadataVersion;
use arrow_select::concat::concat_batches;
use thiserror::Error;

use crate::provenance::{blake3_file, FileHash, RunManifest};
use crate::util::atomic::{write_atomic, AtomicFile};

const BATCH_ROWS: usize = 65_536;

/// Successful merged artifact and receipt paths.
#[derive(Debug, Clone)]
pub struct MergeOutcome {
    /// Canonical merged artifact.
    pub output: PathBuf,
    /// Schema-5 merge receipt.
    pub manifest: PathBuf,
    /// Number of complete shards.
    pub shard_count: u32,
}

/// Canonical merge refusal or I/O failure.
#[derive(Debug, Error)]
pub enum MergeError {
    /// Missing input or destination collision.
    #[error("merge input/output error: {0}")]
    Io(#[from] io::Error),
    /// Parent receipt failed integrity verification.
    #[error("tampered shard receipt: {0}")]
    Tampered(String),
    /// Shards are incomplete or incompatible.
    #[error("incompatible shard set: {0}")]
    Incompatible(String),
    /// Artifact format has no first-party canonical codec.
    #[error("unsupported merge codec {0:?}")]
    Unsupported(String),
}

impl MergeError {
    /// Stable CLI exit code.
    pub fn exit_code(&self) -> i32 {
        match self {
            Self::Io(_) => 2,
            Self::Incompatible(_) | Self::Unsupported(_) => 3,
            Self::Tampered(_) => 5,
        }
    }
}

#[derive(Debug)]
struct Shard {
    index: u32,
    count: u32,
    receipt_path: PathBuf,
    artifact_path: PathBuf,
    manifest: RunManifest,
}

/// Validate and merge a complete first-party shard set transactionally.
pub fn merge_shards(
    manifests: &[PathBuf],
    input_roots: &[PathBuf],
    output: &Path,
    manifest_output: Option<&Path>,
    replace: bool,
) -> Result<MergeOutcome, MergeError> {
    if manifests.is_empty() {
        return Err(MergeError::Incompatible(
            "no shard manifests supplied".into(),
        ));
    }
    let receipt_path = manifest_output
        .map(Path::to_path_buf)
        .unwrap_or_else(|| PathBuf::from(format!("{}.manifest.json", output.display())));
    if output == receipt_path {
        return Err(MergeError::Io(io::Error::new(
            io::ErrorKind::InvalidInput,
            "merge artifact and receipt destinations must differ",
        )));
    }
    if !replace && receipt_path.exists() {
        return Err(MergeError::Io(io::Error::new(
            io::ErrorKind::AlreadyExists,
            format!("merge receipt already exists: {}", receipt_path.display()),
        )));
    }

    let mut shards = Vec::with_capacity(manifests.len());
    for receipt_path in manifests {
        let text = fs::read_to_string(receipt_path)?;
        let manifest = RunManifest::from_canonical_json(&text)
            .map_err(|error| MergeError::Tampered(error.to_string()))?;
        if manifest.self_hash_ok() != Some(true) || manifest.measurement_hash_ok() == Some(false) {
            return Err(MergeError::Tampered(receipt_path.display().to_string()));
        }
        if manifest.params.get("partition.kind").map(String::as_str) != Some("shard") {
            return Err(MergeError::Incompatible(format!(
                "{} is not a shard receipt",
                receipt_path.display()
            )));
        }
        if manifest
            .params
            .get("run_status")
            .is_some_and(|status| status != "completed")
        {
            return Err(MergeError::Incompatible(format!(
                "{} did not complete successfully",
                receipt_path.display()
            )));
        }
        let index = parse_u32(&manifest, "partition.shard_index")?;
        let count = parse_u32(&manifest, "partition.shard_count")?;
        let output_entry = manifest.outputs.first().ok_or_else(|| {
            MergeError::Incompatible(format!("{} has no output artifact", receipt_path.display()))
        })?;
        let artifact_path = locate_artifact(output_entry, receipt_path, input_roots)?;
        shards.push(Shard {
            index,
            count,
            receipt_path: receipt_path.clone(),
            artifact_path,
            manifest,
        });
    }
    shards.sort_by_key(|shard| shard.index);
    let count = shards[0].count;
    if count as usize != shards.len()
        || shards
            .iter()
            .enumerate()
            .any(|(index, shard)| shard.count != count || shard.index != index as u32)
    {
        return Err(MergeError::Incompatible(format!(
            "expected exactly shard indices 0..{count}"
        )));
    }
    let expected = compatibility_claim(&shards[0].manifest);
    for shard in &shards[1..] {
        if compatibility_claim(&shard.manifest) != expected {
            return Err(MergeError::Incompatible(format!(
                "shard {} analyzer/reference/parameter claims differ",
                shard.index
            )));
        }
    }
    let format = shards[0]
        .manifest
        .params
        .get("artifact.output.0.format")
        .or_else(|| shards[0].manifest.params.get("artifact.format"))
        .cloned()
        .ok_or_else(|| MergeError::Unsupported("unknown".into()))?;

    let mut atomic = AtomicFile::create(output)?;
    match format.as_str() {
        "tsv" | "coverage-tsv" => merge_text(&shards, atomic.file_mut(), TextCodec::Tsv)?,
        "vcf-sites" => merge_text(&shards, atomic.file_mut(), TextCodec::Vcf)?,
        "gvcf" => merge_text(&shards, atomic.file_mut(), TextCodec::Gvcf)?,
        "arrow-ipc" => merge_arrow(&shards, atomic.file_mut())?,
        other => return Err(MergeError::Unsupported(other.to_string())),
    }
    atomic.commit(replace)?;

    let mut receipt = RunManifest::new("merge");
    let mut parent_claims = Vec::with_capacity(shards.len());
    for shard in &shards {
        receipt.inputs.push(FileHash {
            path: shard.receipt_path.display().to_string(),
            blake3: blake3_file(&shard.receipt_path)?,
        });
        receipt.inputs.push(FileHash {
            path: shard.artifact_path.display().to_string(),
            blake3: blake3_file(&shard.artifact_path)?,
        });
        if let Some(claim) = shard.manifest.params.get("manifest_blake3") {
            parent_claims.push(claim.clone());
        }
    }
    receipt.outputs.push(FileHash {
        path: output.display().to_string(),
        blake3: blake3_file(output)?,
    });
    receipt.params.insert("merge.codec".into(), format);
    receipt
        .params
        .insert("merge.shard_count".into(), count.to_string());
    receipt.params.insert(
        "merge.parent_claims".into(),
        format!(
            "[{}]",
            parent_claims
                .iter()
                .map(|claim| format!("\"{claim}\""))
                .collect::<Vec<_>>()
                .join(",")
        ),
    );
    receipt.params.insert(
        "artifact.output.0.format".into(),
        receipt.params["merge.codec"].clone(),
    );
    receipt
        .params
        .insert("artifact.output.0.role".into(), "merged-analysis".into());
    receipt
        .params
        .insert("run_status".into(), "completed".into());
    receipt.finalize();
    write_atomic(
        &receipt_path,
        receipt.to_canonical_json().as_bytes(),
        replace,
    )?;
    Ok(MergeOutcome {
        output: output.to_path_buf(),
        manifest: receipt_path,
        shard_count: count,
    })
}

fn parse_u32(manifest: &RunManifest, key: &str) -> Result<u32, MergeError> {
    manifest
        .params
        .get(key)
        .and_then(|value| value.parse().ok())
        .ok_or_else(|| MergeError::Incompatible(format!("missing or invalid {key}")))
}

fn compatibility_claim(manifest: &RunManifest) -> (Vec<String>, BTreeMap<String, String>) {
    let inputs = manifest
        .inputs
        .iter()
        .take(2)
        .map(|input| input.blake3.clone())
        .collect();
    let ignored: BTreeSet<&str> = [
        "command",
        "command_argv",
        "manifest_blake3",
        "measurement_blake3",
        "feature_rows",
        "predicted_working_set_bytes",
        "shard_count",
        "shard_index",
        "target_triple",
        "rustc_version",
        "code_dirty",
    ]
    .into_iter()
    .collect();
    let params = manifest
        .params
        .iter()
        .filter(|(key, _)| !key.starts_with("partition.") && !ignored.contains(key.as_str()))
        .map(|(key, value)| (key.clone(), value.clone()))
        .collect();
    (inputs, params)
}

fn locate_artifact(
    output: &FileHash,
    receipt: &Path,
    roots: &[PathBuf],
) -> Result<PathBuf, MergeError> {
    let recorded = PathBuf::from(&output.path);
    let sibling = receipt
        .parent()
        .unwrap_or_else(|| Path::new("."))
        .join(&output.path);
    for candidate in [recorded, sibling] {
        if candidate.is_file() && blake3_file(&candidate)? == output.blake3 {
            return Ok(candidate);
        }
    }
    for root in roots {
        if let Some(found) = find_hash(root, &output.blake3)? {
            return Ok(found);
        }
    }
    Err(MergeError::Io(io::Error::new(
        io::ErrorKind::NotFound,
        format!("cannot content-locate shard artifact {}", output.blake3),
    )))
}

fn find_hash(root: &Path, expected: &str) -> io::Result<Option<PathBuf>> {
    if root.is_file() {
        return Ok((blake3_file(root)? == expected).then(|| root.to_path_buf()));
    }
    for entry in fs::read_dir(root)? {
        let path = entry?.path();
        if path.is_dir() {
            if let Some(found) = find_hash(&path, expected)? {
                return Ok(Some(found));
            }
        } else if path.is_file() && blake3_file(&path)? == expected {
            return Ok(Some(path));
        }
    }
    Ok(None)
}

#[derive(Clone, Copy)]
enum TextCodec {
    Tsv,
    Vcf,
    Gvcf,
}

fn merge_text(
    shards: &[Shard],
    output: &mut impl Write,
    codec: TextCodec,
) -> Result<(), MergeError> {
    let mut writer = BufWriter::new(output);
    let mut expected_header: Option<Vec<String>> = None;
    let mut pending_gvcf: Option<String> = None;
    for shard in shards {
        let reader = BufReader::new(fs::File::open(&shard.artifact_path)?);
        let mut header = Vec::new();
        let mut records = Vec::new();
        for line in reader.lines() {
            let line = line?;
            if line.starts_with('#') {
                header.push(line);
            } else if !line.is_empty() {
                records.push(line);
            }
        }
        if let Some(expected) = &expected_header {
            if expected != &header {
                return Err(MergeError::Incompatible("text headers differ".into()));
            }
        } else {
            for line in &header {
                writeln!(writer, "{line}")?;
            }
            expected_header = Some(header);
        }
        for record in records {
            if matches!(codec, TextCodec::Gvcf) {
                if let Some(previous) = pending_gvcf.take() {
                    if let Some(merged) = coalesce_gvcf(&previous, &record) {
                        pending_gvcf = Some(merged);
                    } else {
                        writeln!(writer, "{previous}")?;
                        pending_gvcf = Some(record);
                    }
                } else {
                    pending_gvcf = Some(record);
                }
            } else {
                writeln!(writer, "{record}")?;
            }
        }
    }
    if let Some(record) = pending_gvcf {
        writeln!(writer, "{record}")?;
    }
    writer.flush()?;
    Ok(())
}

fn coalesce_gvcf(left: &str, right: &str) -> Option<String> {
    let mut left_fields: Vec<String> = left.split('\t').map(str::to_string).collect();
    let right_fields: Vec<&str> = right.split('\t').collect();
    if left_fields.len() != 10
        || right_fields.len() != 10
        || left_fields[4] != "<NON_REF>"
        || right_fields[4] != "<NON_REF>"
        || left_fields[0] != right_fields[0]
    {
        return None;
    }
    let left_end = left_fields[7].strip_prefix("END=")?.parse::<u64>().ok()?;
    let right_start = right_fields[1].parse::<u64>().ok()?;
    if right_start != left_end + 1 || left_fields[8] != right_fields[8] {
        return None;
    }
    let left_sample: Vec<String> = left_fields[9].split(':').map(str::to_string).collect();
    let right_sample: Vec<&str> = right_fields[9].split(':').collect();
    if left_sample.len() != 4
        || right_sample.len() != 4
        || left_sample[0] != "0/0"
        || right_sample[0] != "0/0"
    {
        return None;
    }
    let left_gq = left_sample[2].parse::<u8>().ok()?;
    let right_gq = right_sample[2].parse::<u8>().ok()?;
    if left_gq / 10 != right_gq / 10 {
        return None;
    }
    let min_gq = left_gq.min(right_gq);
    let min_dp = left_sample[3]
        .parse::<u32>()
        .ok()?
        .min(right_sample[3].parse::<u32>().ok()?);
    left_fields[7] = right_fields[7].to_string();
    left_fields[9] = format!(
        "{}:{}:{}:{}",
        left_sample[0], left_sample[1], min_gq, min_dp
    );
    Some(left_fields.join("\t"))
}

fn merge_arrow(shards: &[Shard], output: &mut fs::File) -> Result<(), MergeError> {
    let mut writer: Option<StreamWriter<&mut fs::File>> = None;
    let mut pending: Vec<RecordBatch> = Vec::new();
    let mut pending_rows = 0usize;
    let mut expected_schema = None;
    for shard in shards {
        let file = fs::File::open(&shard.artifact_path)?;
        let reader = StreamReader::try_new(file, None)
            .map_err(|error| MergeError::Incompatible(error.to_string()))?;
        let schema = reader.schema();
        if expected_schema
            .as_ref()
            .is_some_and(|expected| expected != &schema)
        {
            return Err(MergeError::Incompatible("Arrow schemas differ".into()));
        }
        if expected_schema.is_none() {
            let options = IpcWriteOptions::try_new(8, false, MetadataVersion::V5)
                .map_err(|error| MergeError::Incompatible(error.to_string()))?;
            writer = Some(
                StreamWriter::try_new_with_options(&mut *output, &schema, options)
                    .map_err(|error| MergeError::Incompatible(error.to_string()))?,
            );
            expected_schema = Some(schema.clone());
        }
        for batch in reader {
            let batch = batch.map_err(|error| MergeError::Incompatible(error.to_string()))?;
            pending_rows += batch.num_rows();
            pending.push(batch);
            while pending_rows >= BATCH_ROWS {
                let schema = expected_schema.as_ref().expect("initialized schema");
                let combined = concat_batches(schema, &pending)
                    .map_err(|error| MergeError::Incompatible(error.to_string()))?;
                writer
                    .as_mut()
                    .expect("initialized writer")
                    .write(&combined.slice(0, BATCH_ROWS))
                    .map_err(|error| MergeError::Incompatible(error.to_string()))?;
                let remainder = combined.slice(BATCH_ROWS, combined.num_rows() - BATCH_ROWS);
                pending_rows = remainder.num_rows();
                pending = if remainder.num_rows() == 0 {
                    Vec::new()
                } else {
                    vec![remainder]
                };
            }
        }
    }
    if let Some(writer) = writer.as_mut() {
        if pending_rows != 0 {
            let combined = concat_batches(expected_schema.as_ref().unwrap(), &pending)
                .map_err(|error| MergeError::Incompatible(error.to_string()))?;
            writer
                .write(&combined)
                .map_err(|error| MergeError::Incompatible(error.to_string()))?;
        }
        writer
            .finish()
            .map_err(|error| MergeError::Incompatible(error.to_string()))?;
    }
    Ok(())
}
