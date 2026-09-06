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
        if manifest.params.get("pileup.semantics").map(String::as_str) != Some("exact-or-fail-v1") {
            return Err(MergeError::Incompatible(format!(
                "{} lacks exact-or-fail-v1 pileup semantics; historical sampled shards cannot be canonically merged; reproduce them with their historical producer or rerun all shards with the current exact producer",
                receipt_path.display()
            )));
        }
        if manifest
            .params
            .get("partition.algorithm")
            .map(String::as_str)
            != Some("reference-span-v1")
        {
            return Err(MergeError::Incompatible(format!(
                "{} has an unsupported partition algorithm",
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
        if manifest.outputs.len() != 1 {
            return Err(MergeError::Incompatible(format!(
                "{} must name exactly one first-party output artifact",
                receipt_path.display()
            )));
        }
        let output_entry = &manifest.outputs[0];
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
        .map(|input| input.blake3.clone())
        .collect();
    let ignored: BTreeSet<&str> = [
        "command",
        "command_argv",
        "manifest_blake3",
        "measurement_blake3",
        "feature_rows",
        "over_max_depth",
        "reads_skipped_total",
        "predicted_working_set_bytes",
        "max_depth",
        "max_read_len",
        "memory_budget_mb",
        "enforce",
        "require_os_limit",
        "contract.assurance",
        "os.memory_limit_bytes",
        "analyzer.memory_model",
        "analyzer.max_additional_bytes",
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
        .filter(|(key, _)| {
            (!key.starts_with("partition.") || key.as_str() == "partition.algorithm")
                && !key.starts_with("outcome.")
                && !ignored.contains(key.as_str())
        })
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

// Retain a constant-size header identity, not the header or shard contents.
// Text merge memory is bounded by the largest record plus the writer buffer.
type HeaderIdentity = (blake3::Hash, u64);

fn merge_text(
    shards: &[Shard],
    output: &mut impl Write,
    codec: TextCodec,
) -> Result<(), MergeError> {
    let mut writer = BufWriter::new(output);
    let mut expected_header = None;
    let mut pending_gvcf = None;
    for shard in shards {
        let reader = BufReader::new(fs::File::open(&shard.artifact_path)?);
        merge_text_stream(
            reader,
            &mut writer,
            codec,
            &mut expected_header,
            &mut pending_gvcf,
        )?;
    }
    if let Some(record) = pending_gvcf {
        writeln!(writer, "{record}")?;
    }
    writer.flush()?;
    Ok(())
}

fn merge_text_stream(
    mut reader: impl BufRead,
    writer: &mut impl Write,
    codec: TextCodec,
    expected_header: &mut Option<HeaderIdentity>,
    pending_gvcf: &mut Option<String>,
) -> Result<(), MergeError> {
    let first_shard = expected_header.is_none();
    let mut header_hash = blake3::Hasher::new();
    let mut header_lines = 0u64;
    let mut header_finished = false;
    let mut line = String::new();
    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        let record = line.trim_end_matches(['\r', '\n']);
        if record.is_empty() {
            continue;
        }
        if record.starts_with('#') {
            if header_finished {
                return Err(MergeError::Incompatible(
                    "text header encountered after records".into(),
                ));
            }
            header_hash.update(&(record.len() as u64).to_le_bytes());
            header_hash.update(record.as_bytes());
            header_lines += 1;
            if first_shard {
                writeln!(writer, "{record}")?;
            }
            continue;
        }
        if !header_finished {
            check_text_header(expected_header, &header_hash, header_lines)?;
            header_finished = true;
        }
        if matches!(codec, TextCodec::Gvcf) {
            if let Some(previous) = pending_gvcf.take() {
                if let Some(merged) = coalesce_gvcf(&previous, record) {
                    *pending_gvcf = Some(merged);
                } else {
                    writeln!(writer, "{previous}")?;
                    *pending_gvcf = Some(record.to_string());
                }
            } else {
                *pending_gvcf = Some(record.to_string());
            }
        } else {
            writeln!(writer, "{record}")?;
        }
    }
    if !header_finished {
        check_text_header(expected_header, &header_hash, header_lines)?;
    }
    Ok(())
}

fn check_text_header(
    expected: &mut Option<HeaderIdentity>,
    hash: &blake3::Hasher,
    lines: u64,
) -> Result<(), MergeError> {
    let actual = (hash.finalize(), lines);
    if expected
        .as_ref()
        .is_some_and(|identity| identity != &actual)
    {
        return Err(MergeError::Incompatible("text headers differ".into()));
    }
    *expected = Some(actual);
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::Cell;
    use std::io::{Cursor, Read};
    use std::rc::Rc;

    #[test]
    fn historical_sampled_shards_are_refused_before_publication() {
        let dir = std::env::temp_dir().join(format!(
            "rosalind-historical-merge-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&dir).unwrap();
        let artifact = dir.join("shard.tsv");
        fs::write(&artifact, "#contig\tpos\tdepth\nchr1\t1\t1\n").unwrap();
        for semantics in [None, Some("sampled-v0")] {
            let mut receipt = RunManifest::new("features");
            receipt.params.extend([
                ("partition.kind".into(), "shard".into()),
                ("partition.algorithm".into(), "reference-span-v1".into()),
                ("partition.shard_count".into(), "1".into()),
                ("partition.shard_index".into(), "0".into()),
                ("artifact.format".into(), "tsv".into()),
                ("run_status".into(), "completed".into()),
            ]);
            if let Some(semantics) = semantics {
                receipt
                    .params
                    .insert("pileup.semantics".into(), semantics.into());
            }
            receipt.outputs.push(FileHash {
                path: artifact.display().to_string(),
                blake3: blake3_file(&artifact).unwrap(),
            });
            receipt.finalize();
            let manifest = dir.join("shard.json");
            fs::write(&manifest, receipt.to_canonical_json()).unwrap();
            let output = dir.join("merged.tsv");
            let error = merge_shards(&[manifest], &[], &output, None, false).unwrap_err();
            assert!(error.to_string().contains("historical sampled shards"));
            assert!(!output.exists());
        }
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn outcome_counters_and_execution_capacity_do_not_change_compatibility() {
        let mut a = RunManifest::new("features");
        a.params.extend([
            ("pileup.semantics".into(), "exact-or-fail-v1".into()),
            ("partition.algorithm".into(), "reference-span-v1".into()),
            ("min_mapq".into(), "20".into()),
            ("max_depth".into(), "100".into()),
            ("reads_skipped_total".into(), "2".into()),
        ]);
        let mut b = a.clone();
        b.params.insert("reads_skipped_total".into(), "13".into());
        b.params.insert("max_depth".into(), "1000".into());
        b.params
            .insert("outcome.covered_bases".into(), "999".into());
        assert_eq!(compatibility_claim(&a), compatibility_claim(&b));
        b.params.insert("min_mapq".into(), "30".into());
        assert_ne!(compatibility_claim(&a), compatibility_claim(&b));
        b.params.insert("min_mapq".into(), "20".into());
        b.params
            .insert("pileup.semantics".into(), "sampled-v0".into());
        assert_ne!(compatibility_claim(&a), compatibility_claim(&b));
        b.params
            .insert("pileup.semantics".into(), "exact-or-fail-v1".into());
        b.params
            .insert("partition.algorithm".into(), "unknown".into());
        assert_ne!(compatibility_claim(&a), compatibility_claim(&b));
    }

    #[test]
    fn text_merge_emits_records_before_reading_the_complete_shard() {
        struct ObservedRead {
            input: Cursor<Vec<u8>>,
            written: Rc<Cell<usize>>,
        }
        impl Read for ObservedRead {
            fn read(&mut self, output: &mut [u8]) -> io::Result<usize> {
                if self.input.position() >= 32_768 && self.written.get() == 0 {
                    return Err(io::Error::other("shard was buffered before output"));
                }
                self.input.read(output)
            }
        }
        struct ObservedWrite(Rc<Cell<usize>>);
        impl Write for ObservedWrite {
            fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
                self.0.set(self.0.get() + bytes.len());
                Ok(bytes.len())
            }
            fn flush(&mut self) -> io::Result<()> {
                Ok(())
            }
        }
        let bytes = format!("#contig\tpos\tdepth\n{}", "chr1\t1\t20\n".repeat(20_000)).into_bytes();
        let expected_len = bytes.len();
        let written = Rc::new(Cell::new(0));
        let reader = BufReader::new(ObservedRead {
            input: Cursor::new(bytes),
            written: Rc::clone(&written),
        });
        let mut writer = BufWriter::new(ObservedWrite(Rc::clone(&written)));
        merge_text_stream(reader, &mut writer, TextCodec::Tsv, &mut None, &mut None).unwrap();
        writer.flush().unwrap();
        assert_eq!(written.get(), expected_len);
    }

    #[test]
    fn empty_shard_headers_are_checked_and_late_headers_are_rejected() {
        let mut expected = None;
        let mut pending = None;
        let mut output = Vec::new();
        for input in ["#a\tb\n", "#a\tb\n1\t2\n", "#a\tb\n"] {
            merge_text_stream(
                Cursor::new(input),
                &mut output,
                TextCodec::Tsv,
                &mut expected,
                &mut pending,
            )
            .unwrap();
        }
        assert_eq!(output, b"#a\tb\n1\t2\n");
        assert!(merge_text_stream(
            Cursor::new("#different\n"),
            &mut output,
            TextCodec::Tsv,
            &mut expected,
            &mut pending
        )
        .is_err());
        assert!(merge_text_stream(
            Cursor::new("#a\tb\n3\t4\n#late\n"),
            &mut output,
            TextCodec::Tsv,
            &mut expected,
            &mut pending
        )
        .is_err());
    }

    #[test]
    fn gvcf_boundary_coalescing_retains_only_one_pending_record() {
        let mut expected = None;
        let mut pending = None;
        let mut output = Vec::new();
        let inputs = [
            "#header\nchr1\t1\t.\tA\t<NON_REF>\t.\tPASS\tEND=2\tGT:PL:GQ:MIN_DP\t0/0:0,30,40:35:9\n",
            "#header\nchr1\t3\t.\tA\t<NON_REF>\t.\tPASS\tEND=4\tGT:PL:GQ:MIN_DP\t0/0:0,30,40:32:7\n",
        ];
        for input in inputs {
            merge_text_stream(
                Cursor::new(input),
                &mut output,
                TextCodec::Gvcf,
                &mut expected,
                &mut pending,
            )
            .unwrap();
        }
        assert_eq!(output, b"#header\n");
        assert_eq!(
            pending.as_deref(),
            Some("chr1\t1\t.\tA\t<NON_REF>\t.\tPASS\tEND=4\tGT:PL:GQ:MIN_DP\t0/0:0,30,40:32:7")
        );
    }
}
