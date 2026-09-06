//! Streaming, verified differences between exact per-locus evidence datasets.
//!
//! Two bounded row queues normalize Arrow and TSV into the same canonical
//! representation. Reference dictionaries come from content-verified alignment
//! headers; execution settings and file encoding do not change locus identity.

use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::sync::mpsc::{self, Receiver, SyncSender};

use rust_htslib::bam::{self, Read};
use thiserror::Error;

use crate::core::ContigSet;
use crate::evidence::{
    read_evidence_batches, EvidenceAnalyzer, EvidenceReference, EvidenceTsvWriter,
};
use crate::provenance::{command_template_tokens, verify_receipt, VerifyOpts};
use crate::util::atomic::AtomicFile;

const MAX_ROW_BYTES: usize = 64 * 1024;

/// Comparison integrity, schema, or local artifact failure.
#[derive(Debug, Error)]
pub enum DatasetDiffError {
    /// A receipt or its recorded files failed verification.
    #[error("evidence diff integrity failure: {0}")]
    Integrity(String),
    /// The datasets do not expose compatible exact per-locus schemas.
    #[error("incompatible evidence datasets: {0}")]
    Incompatible(String),
    /// A recorded artifact could not be read, or the new output could not be written.
    #[error("evidence diff I/O failed: {0}")]
    Io(#[from] io::Error),
    /// A stream decoder failed unexpectedly.
    #[error("evidence diff decoder terminated unexpectedly")]
    WorkerPanic,
}
impl DatasetDiffError {
    /// Integrity failures exit 5, incompatible requests 3, local failures 2.
    pub fn exit_code(&self) -> i32 {
        match self {
            Self::Integrity(_) => 5,
            Self::Incompatible(_) => 3,
            Self::Io(_) | Self::WorkerPanic => 2,
        }
    }
}

/// Counts in a successfully published locus-delta table.
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct DatasetDiffSummary {
    /// Loci whose shared metric values differ.
    pub changed_loci: u64,
    /// Loci present only in the second dataset.
    pub added_loci: u64,
    /// Loci present only in the first dataset.
    pub removed_loci: u64,
    /// Data rows written, one for each differing metric (including added/removed).
    pub metric_changes: u64,
}

struct VerifiedDataset {
    path: PathBuf,
    arrow: bool,
    contigs: ContigSet,
}

fn operand<'a>(tokens: &'a [String], flags: &[&str]) -> Option<&'a str> {
    tokens
        .windows(2)
        .find(|pair| flags.contains(&pair[0].as_str()))
        .map(|pair| pair[1].as_str())
}
fn verified_dataset(receipt_path: &Path) -> Result<VerifiedDataset, DatasetDiffError> {
    let text = fs::read_to_string(receipt_path)?;
    let verified = verify_receipt(
        &text,
        &VerifyOpts {
            rehash_files: true,
            ..VerifyOpts::default()
        },
    );
    if !verified.ok {
        return Err(DatasetDiffError::Integrity(format!(
            "{}: {}",
            receipt_path.display(),
            verified.problems.join("; ")
        )));
    }
    let manifest = verified.manifest.unwrap();
    let compatible = [
        ("evidence.schema", "1"),
        ("evidence.sampling", "none"),
        ("run_status", "completed"),
    ];
    if compatible
        .iter()
        .any(|(key, value)| manifest.params.get(*key).map(String::as_str) != Some(*value))
    {
        return Err(DatasetDiffError::Incompatible(
            "requires completed exact evidence schema version 1 receipts".into(),
        ));
    }
    let tokens = command_template_tokens(&manifest).map_err(DatasetDiffError::Integrity)?;
    let panel = manifest.params.get("analyzer.id").map(String::as_str) == Some("panel-qc");
    let output_flags: &[&str] = if panel {
        &["--position-output"]
    } else {
        &["-o", "--output"]
    };
    let output_hash = operand(&tokens, output_flags)
        .and_then(|token| token.strip_prefix("@out:"))
        .ok_or_else(|| {
            DatasetDiffError::Incompatible(
                if panel {
                    "panel comparison requires a recorded --position-output Arrow artifact"
                } else {
                    "receipt has no recorded per-locus output"
                }
                .into(),
            )
        })?;
    let output = manifest
        .outputs
        .iter()
        .find(|file| file.blake3 == output_hash)
        .ok_or_else(|| {
            DatasetDiffError::Integrity("output operand is absent from recorded outputs".into())
        })?;
    let alignment_hash = operand(&tokens, &["--alignments"])
        .and_then(|token| token.strip_prefix("@in:"))
        .ok_or_else(|| {
            DatasetDiffError::Incompatible("receipt has no verified alignment dictionary".into())
        })?;
    let alignment = manifest
        .inputs
        .iter()
        .find(|file| file.blake3 == alignment_hash)
        .ok_or_else(|| {
            DatasetDiffError::Integrity("alignment operand is absent from recorded inputs".into())
        })?;
    let reader = bam::Reader::from_path(&alignment.path).map_err(|error| {
        DatasetDiffError::Integrity(format!(
            "cannot read verified alignment dictionary: {error}"
        ))
    })?;
    let header = reader.header();
    let mut contigs = ContigSet::new();
    for tid in 0..header.target_count() {
        let name = std::str::from_utf8(header.tid2name(tid)).map_err(|_| {
            DatasetDiffError::Incompatible("alignment contig name is not UTF-8".into())
        })?;
        if name.len() > 4096 || name.chars().any(char::is_control) {
            return Err(DatasetDiffError::Incompatible(
                "unsupported contig name in dictionary".into(),
            ));
        }
        let length = header
            .target_len(tid)
            .and_then(|value| u32::try_from(value).ok())
            .ok_or_else(|| {
                DatasetDiffError::Incompatible(
                    "alignment contig length is absent or exceeds u32".into(),
                )
            })?;
        if contigs.iter().any(|contig| contig.name.as_ref() == name) {
            return Err(DatasetDiffError::Incompatible(
                "alignment dictionary contains duplicate contigs".into(),
            ));
        }
        contigs.push(name.to_owned(), length);
    }
    // The extractor orders by the analysis reference, which may differ from
    // BAM target order. Reference-free coverage instead uses the BAM dictionary.
    if let Some(reference_operand) = operand(&tokens, &["--reference", "--cram-reference"]) {
        let reference_hash = reference_operand.strip_prefix("@in:").ok_or_else(|| {
            DatasetDiffError::Integrity("reference is not a content operand".into())
        })?;
        let reference_file = manifest
            .inputs
            .iter()
            .find(|file| file.blake3 == reference_hash)
            .ok_or_else(|| {
                DatasetDiffError::Integrity("reference operand is absent from inputs".into())
            })?;
        let fai_flag = if operand(&tokens, &["--reference"]).is_some() {
            "--reference-fai"
        } else {
            "--cram-reference-fai"
        };
        let fai = operand(&tokens, &[fai_flag])
            .map(|token| {
                let hash = token.strip_prefix("@in:").ok_or_else(|| {
                    DatasetDiffError::Integrity("FAI is not a content operand".into())
                })?;
                manifest
                    .inputs
                    .iter()
                    .find(|file| file.blake3 == hash)
                    .map(|file| Path::new(&file.path))
                    .ok_or_else(|| {
                        DatasetDiffError::Integrity("FAI operand is absent from inputs".into())
                    })
            })
            .transpose()?;
        let reference =
            EvidenceReference::open_with_fai(&reference_file.path, fai).map_err(|error| {
                DatasetDiffError::Integrity(format!(
                    "invalid recorded reference dictionary: {error}"
                ))
            })?;
        if reference.contigs().len() != contigs.len()
            || reference.contigs().iter().any(|contig| {
                contigs
                    .by_name(&contig.name)
                    .is_none_or(|other| other.length != contig.length)
            })
        {
            return Err(DatasetDiffError::Incompatible(
                "reference and alignment dictionaries disagree".into(),
            ));
        }
        contigs = reference.contigs().clone();
    }
    let arrow = if panel {
        true
    } else {
        match operand(&tokens, &["--format"]).unwrap_or("tsv") {
            "tsv" => false,
            "arrow-ipc" => true,
            other => {
                return Err(DatasetDiffError::Incompatible(format!(
                    "unsupported per-locus encoding {other}"
                )))
            }
        }
    };
    Ok(VerifiedDataset {
        path: PathBuf::from(&output.path),
        arrow,
        contigs,
    })
}

/// Verify both source receipts and all recorded artifacts, then atomically create
/// a new TSV delta table. Existing output files are never replaced. Added/removed
/// loci emit one delta for every metric; unchanged datasets contain only a header.
pub fn diff_evidence_datasets(
    a: &Path,
    b: &Path,
    output: &Path,
) -> Result<DatasetDiffSummary, DatasetDiffError> {
    if output.exists() {
        return Err(io::Error::new(
            io::ErrorKind::AlreadyExists,
            format!("{} already exists", output.display()),
        )
        .into());
    }
    let first = verified_dataset(a)?;
    let second = verified_dataset(b)?;
    if !first
        .contigs
        .iter()
        .map(|c| (&c.name, c.length))
        .eq(second.contigs.iter().map(|c| (&c.name, c.length)))
    {
        return Err(DatasetDiffError::Incompatible(
            "reference dictionaries differ in contig names, lengths, or order".into(),
        ));
    }
    let mut header_writer = EvidenceTsvWriter::new(Vec::new());
    header_writer
        .finish()
        .map_err(|error| DatasetDiffError::Integrity(error.to_string()))?;
    let header = String::from_utf8(header_writer.into_inner()).unwrap();
    let metrics: Vec<&str> = header.trim_end().split('\t').skip(2).collect();
    let dictionary: BTreeMap<String, (u32, u32)> = first
        .contigs
        .iter()
        .map(|contig| (contig.name.to_string(), (contig.id, contig.length)))
        .collect();
    let mut staged = AtomicFile::create(output)?;
    let mut writer = io::BufWriter::new(staged.file_mut());
    writeln!(writer, "contig\tpos\tchange\tmetric\tbefore\tafter")?;
    let summary = std::thread::scope(|scope| {
        let (send_a, recv_a) = mpsc::sync_channel(2);
        let (send_b, recv_b) = mpsc::sync_channel(2);
        let first_handle = scope.spawn(|| send_rows(first, send_a));
        let second_handle = scope.spawn(|| send_rows(second, send_b));
        let compared = (|| {
            let mut first = Rows::new(&recv_a, &header, &dictionary)?;
            let mut second = Rows::new(&recv_b, &header, &dictionary)?;
            compare_rows(&mut first, &mut second, &metrics, &mut writer)
        })();
        // Dropping receivers unblocks both bounded producers after any failure.
        drop(recv_a);
        drop(recv_b);
        let joined_a = first_handle.join();
        let joined_b = second_handle.join();
        if joined_a.is_err() || joined_b.is_err() {
            return Err(DatasetDiffError::WorkerPanic);
        }
        compared
    })?;
    writer.flush()?;
    drop(writer);
    staged.commit(false)?;
    Ok(summary)
}

/// A writer that forwards complete canonical lines through a bounded row queue.
struct LineSender {
    sender: SyncSender<Result<String, DatasetDiffError>>,
    pending: Vec<u8>,
}
impl Write for LineSender {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        for piece in bytes.split_inclusive(|byte| *byte == b'\n') {
            if self.pending.len() + piece.len() > MAX_ROW_BYTES {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "canonical evidence row exceeds bound",
                ));
            }
            self.pending.extend_from_slice(piece);
            if piece.last() == Some(&b'\n') {
                let bytes = std::mem::take(&mut self.pending);
                let line = String::from_utf8(bytes)
                    .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))?;
                self.sender.send(Ok(line)).map_err(|_| {
                    io::Error::new(io::ErrorKind::BrokenPipe, "diff consumer stopped")
                })?;
            }
        }
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
fn send_rows(dataset: VerifiedDataset, sender: SyncSender<Result<String, DatasetDiffError>>) {
    let result = (|| {
        let file = File::open(&dataset.path)?;
        if dataset.arrow {
            let mut writer = EvidenceTsvWriter::new(LineSender {
                sender: sender.clone(),
                pending: Vec::new(),
            });
            read_evidence_batches(BufReader::new(file), &dataset.contigs, |batch| {
                writer.on_batch(batch)
            })
            .map_err(|error| {
                DatasetDiffError::Integrity(format!("invalid Arrow evidence: {error}"))
            })?;
            writer
                .finish()
                .map_err(|error| DatasetDiffError::Integrity(error.to_string()))?;
        } else {
            let mut reader = BufReader::new(file);
            loop {
                let mut bytes = Vec::new();
                let mut limited = std::io::Read::take(&mut reader, (MAX_ROW_BYTES + 1) as u64);
                if limited.read_until(b'\n', &mut bytes)? == 0 {
                    break;
                }
                if bytes.len() > MAX_ROW_BYTES || bytes.last() != Some(&b'\n') {
                    return Err(DatasetDiffError::Integrity(
                        "TSV evidence row exceeds bound or is unterminated".into(),
                    ));
                }
                let line = String::from_utf8(bytes)
                    .map_err(|error| DatasetDiffError::Integrity(error.to_string()))?;
                if sender.send(Ok(line)).is_err() {
                    return Ok(());
                }
            }
        }
        Ok::<_, DatasetDiffError>(())
    })();
    if let Err(error) = result {
        let _ = sender.send(Err(error));
    }
}
struct Row {
    key: (u32, u32),
    contig: String,
    values: Vec<String>,
}
struct Rows<'a> {
    receiver: &'a Receiver<Result<String, DatasetDiffError>>,
    dictionary: &'a BTreeMap<String, (u32, u32)>,
    previous: Option<(u32, u32)>,
}
impl<'a> Rows<'a> {
    fn new(
        receiver: &'a Receiver<Result<String, DatasetDiffError>>,
        header: &str,
        dictionary: &'a BTreeMap<String, (u32, u32)>,
    ) -> Result<Self, DatasetDiffError> {
        let actual = receiver
            .recv()
            .map_err(|_| DatasetDiffError::Integrity("missing evidence TSV header".into()))??;
        if actual != header {
            return Err(DatasetDiffError::Incompatible(
                "per-locus TSV schema differs from exact evidence version 1".into(),
            ));
        }
        Ok(Self {
            receiver,
            dictionary,
            previous: None,
        })
    }
    fn next(&mut self) -> Result<Option<Row>, DatasetDiffError> {
        let line = match self.receiver.recv() {
            Ok(line) => line?,
            Err(_) => return Ok(None),
        };
        let fields: Vec<&str> = line.trim_end_matches('\n').split('\t').collect();
        if fields.len() != 34 {
            return Err(DatasetDiffError::Integrity(
                "invalid number of TSV evidence fields".into(),
            ));
        }
        let &(contig, length) = self.dictionary.get(fields[0]).ok_or_else(|| {
            DatasetDiffError::Integrity("evidence contig absent from dictionary".into())
        })?;
        let position = fields[1]
            .parse::<u32>()
            .ok()
            .filter(|pos| *pos > 0 && *pos <= length)
            .ok_or_else(|| DatasetDiffError::Integrity("invalid evidence coordinate".into()))?;
        if position.to_string() != fields[1] {
            return Err(DatasetDiffError::Integrity(
                "noncanonical evidence coordinate".into(),
            ));
        }
        let key = (contig, position);
        if self.previous.is_some_and(|previous| previous >= key) {
            return Err(DatasetDiffError::Integrity(
                "evidence loci must be unique and sorted by dictionary".into(),
            ));
        }
        self.previous = Some(key);
        if fields[2].len() != 1
            || !b"ACGTN".contains(&fields[2].as_bytes()[0])
            || (fields[3] != "."
                && (fields[3]
                    .split(',')
                    .any(|alt| alt.len() != 1 || !b"ACGT".contains(&alt.as_bytes()[0]))
                    || fields[3].len() > 7))
        {
            return Err(DatasetDiffError::Integrity(
                "invalid evidence reference or ALT".into(),
            ));
        }
        if fields[4..32].iter().any(|value| {
            value
                .parse::<u64>()
                .map_or(true, |number| number.to_string() != *value)
        }) {
            return Err(DatasetDiffError::Integrity(
                "invalid evidence scalar counter".into(),
            ));
        }
        for (field, bins) in [(fields[32], 94usize), (fields[33], 255)] {
            if field == "." {
                continue;
            }
            let mut previous = None;
            for item in field.split(',') {
                let (bin, count) = item.split_once(':').ok_or_else(|| {
                    DatasetDiffError::Integrity("invalid evidence histogram".into())
                })?;
                let bin = bin
                    .parse::<usize>()
                    .ok()
                    .filter(|bin| *bin < bins)
                    .ok_or_else(|| {
                        DatasetDiffError::Integrity("histogram bin out of range".into())
                    })?;
                if count
                    .parse::<u64>()
                    .ok()
                    .filter(|count| *count > 0)
                    .is_none()
                    || previous.is_some_and(|before| before >= bin)
                {
                    return Err(DatasetDiffError::Integrity(
                        "noncanonical evidence histogram".into(),
                    ));
                }
                previous = Some(bin);
            }
        }
        Ok(Some(Row {
            key,
            contig: fields[0].to_owned(),
            values: fields[2..]
                .iter()
                .map(|value| (*value).to_owned())
                .collect(),
        }))
    }
}
fn emit_delta(
    writer: &mut dyn Write,
    row: &Row,
    change: &str,
    metrics: &[&str],
    other: Option<&Row>,
) -> Result<u64, DatasetDiffError> {
    let mut count = 0;
    for (index, metric) in metrics.iter().enumerate() {
        let (before, after) = match (change, other) {
            ("added", _) => (".", row.values[index].as_str()),
            ("removed", _) => (row.values[index].as_str(), "."),
            (_, Some(other)) => (row.values[index].as_str(), other.values[index].as_str()),
            _ => unreachable!(),
        };
        if change == "changed" && before == after {
            continue;
        }
        writeln!(
            writer,
            "{}\t{}\t{change}\t{metric}\t{before}\t{after}",
            row.contig, row.key.1
        )?;
        count += 1;
    }
    Ok(count)
}
fn compare_rows(
    a: &mut Rows<'_>,
    b: &mut Rows<'_>,
    metrics: &[&str],
    writer: &mut dyn Write,
) -> Result<DatasetDiffSummary, DatasetDiffError> {
    let mut summary = DatasetDiffSummary::default();
    let mut first = a.next()?;
    let mut second = b.next()?;
    while first.is_some() || second.is_some() {
        let order = match (&first, &second) {
            (Some(a), Some(b)) => a.key.cmp(&b.key),
            (Some(_), None) => Ordering::Less,
            (None, Some(_)) => Ordering::Greater,
            _ => break,
        };
        let count = match order {
            Ordering::Less => {
                let count = emit_delta(writer, first.as_ref().unwrap(), "removed", metrics, None)?;
                summary.removed_loci += 1;
                first = a.next()?;
                count
            }
            Ordering::Greater => {
                let count = emit_delta(writer, second.as_ref().unwrap(), "added", metrics, None)?;
                summary.added_loci += 1;
                second = b.next()?;
                count
            }
            Ordering::Equal => {
                let count = emit_delta(
                    writer,
                    first.as_ref().unwrap(),
                    "changed",
                    metrics,
                    second.as_ref(),
                )?;
                if count > 0 {
                    summary.changed_loci += 1;
                }
                first = a.next()?;
                second = b.next()?;
                count
            }
        };
        summary.metric_changes += count;
    }
    Ok(summary)
}
