//! Checked CRAM 3.0 allocation envelope, before htslib opens a decoder.
//!
//! The native record cap remains cooperative: malformed/oversize decoded records
//! can fail after native allocation. For successful records inside that cap, the
//! plan charges complete container storage, not just the requested genomic tile.
use super::{EvidenceError, EvidenceExecution, EvidenceRequest};
use crate::dataset::InputSnapshot;
use flate2::{read::MultiGzDecoder, Crc};
use rust_htslib::bam::{self, Read as BamRead};
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Cursor, Read};
use std::path::{Path, PathBuf};

const METADATA_CAP: u64 = 8 << 20;
const BLOCK_CAP: u64 = 64 << 20;
const CONTAINER_CAP: u64 = 256 << 20;
const RECORDS_CAP: u64 = 1_000_000;
const BLOCKS_CAP: u64 = 4096;
const SCAN_SCRATCH: u64 = 256 << 10;
const FAI_LINE_CAP: u64 = 8192;

/// Observed maxima from checked CRAM container metadata, before native decode.
#[derive(Debug, Clone, Default)]
pub struct CramEnvelope {
    /// Versioned supported-codec and allocation-estimator contract.
    pub model: &'static str,
    /// Number of data containers in the inspected file.
    pub containers: u64,
    /// Largest complete container's record count.
    pub max_container_records: u64,
    /// Largest complete container's decoded sequence-base count.
    pub max_container_bases: u64,
    /// Largest number of complete slices in one container.
    pub max_container_slices: u64,
    /// Largest block inventory in one container, including headers.
    pub max_container_blocks: u64,
    /// Largest sum of compressed block bytes within one container.
    pub max_compressed_bytes: u64,
    /// Largest sum of declared and independently checked uncompressed block sizes.
    pub max_uncompressed_bytes: u64,
    /// Largest compression-header metadata size.
    pub max_compression_header_bytes: u64,
    /// Largest SAM header block, including its length prefix.
    pub header_bytes: u64,
    /// Decoder reference coexistence allowance, including FASTA line endings.
    pub reference_bytes: u64,
    /// Compressed/text CRAI storage, parsed nodes and temporary expansion.
    pub index_bytes: u64,
    /// Records checked by the complete sequential validation pass.
    pub validated_records: u64,
    /// Largest payload among all records, including off-target records.
    pub validated_max_record_bytes: u64,
    /// Largest read among all records, including off-target records.
    pub validated_max_read_length: u64,
    /// Wall time of whole-file native validation in microseconds.
    pub validation_wall_micros: u64,
    /// Process high-water RSS after whole-file validation.
    pub validation_peak_rss_bytes: u64,
    /// Total records declared by checked container metadata.
    pub declared_records: u64,
    /// Total sequence bases from the complete validation pass.
    pub validated_bases: u64,
    /// Total sequence bases declared by checked containers.
    pub declared_bases: u64,
}

impl CramEnvelope {
    /// Additional per-decoder reservation beyond ordinary record/allocator slack.
    pub fn additional_bytes(&self, _execution: &EvidenceExecution) -> Result<u64, EvidenceError> {
        // A complete native validation pass has checked every record, including
        // records outside future indexed queries. Its maximum is source-bound.
        // htslib retains whole-slice crecs plus generated names/aux/CIGAR/seq/qual.
        // CIGAR capacity doubles; other blocks grow by 25%+800. Three times the
        // complete successful BAM payload envelope covers old/new reallocations.
        // Sequence is expanded (not BAM's packed representation), so charge four
        // bytes per declared base separately. crecs are <256 bytes on the pinned
        // supported 64-bit ABI; two arrays cover container transition lifetime.
        // Compression codecs are <1KiB state per encoded metadata byte/symbol;
        // this also covers nested byte-array/Huffman tables in CRAM3.0 headers.
        self.structural_bytes()?
            .checked_add(
                self.max_container_records
                    .checked_mul(self.validated_max_record_bytes)
                    .and_then(|n| n.checked_mul(3))
                    .ok_or(EvidenceError::CounterOverflow)?,
            )
            .ok_or(EvidenceError::CounterOverflow)
    }

    fn structural_bytes(&self) -> Result<u64, EvidenceError> {
        // Cross-container totals do not prove each container's declared bases.
        // Once validated, every container is also bounded by N * global max RL.
        let base_bound = self
            .max_container_records
            .checked_mul(self.validated_max_read_length)
            .ok_or(EvidenceError::CounterOverflow)?
            .max(self.max_container_bases);
        let parts = [
            self.max_container_records.checked_mul(512),
            base_bound.checked_mul(4),
            self.max_container_slices.checked_mul(16 << 10),
            self.max_container_blocks.checked_mul(512),
            self.max_compressed_bytes.checked_mul(2),
            self.max_uncompressed_bytes.checked_mul(2),
            self.max_compression_header_bytes.checked_mul(1024),
            self.header_bytes.checked_mul(16),
            Some(self.reference_bytes),
            Some(self.index_bytes),
            Some(SCAN_SCRATCH),
        ];
        parts.into_iter().try_fold(0u64, |total, value| {
            total
                .checked_add(value.ok_or(EvidenceError::CounterOverflow)?)
                .ok_or(EvidenceError::CounterOverflow)
        })
    }
}

#[derive(Debug)]
pub(crate) struct CramPreflight {
    pub envelope: CramEnvelope,
    snapshot: InputSnapshot,
    pub index: PathBuf,
}
impl CramPreflight {
    pub fn verify(&self) -> Result<(), EvidenceError> {
        self.snapshot.verify()
    }

    pub fn inspect(request: &EvidenceRequest) -> Result<Self, EvidenceError> {
        let fasta = request
            .cram_reference
            .as_ref()
            .or(request.reference.as_ref())
            .ok_or_else(|| unsupported("an explicit local FASTA and adjacent FAI are required"))?;
        let adjacent = PathBuf::from(format!("{}.fai", fasta.display()));
        if let Some(explicit) = request
            .cram_reference_fai
            .as_ref()
            .or(request.reference_fai.as_ref())
        {
            if std::fs::canonicalize(explicit)? != std::fs::canonicalize(&adjacent)? {
                return Err(unsupported(
                    "decoder FAI must be adjacent to the verified FASTA",
                ));
            }
        }
        let mut paths = vec![request.alignments.clone(), fasta.clone(), adjacent.clone()];
        let index = if let Some(index) = &request.alignment_index {
            index.clone()
        } else {
            [
                PathBuf::from(format!("{}.crai", request.alignments.display())),
                request.alignments.with_extension("crai"),
            ]
            .into_iter()
            .find(|path| path.is_file())
            .ok_or_else(|| invalid("local CRAI index is required"))?
        };
        paths.push(index.clone());
        let snapshot = InputSnapshot::capture(paths)?;
        admit(request.execution.memory_budget_bytes, SCAN_SCRATCH)?;
        let fasta_length = std::fs::metadata(fasta)?.len();
        let layouts = read_fai(
            &adjacent,
            fasta_length,
            request.execution.memory_budget_bytes,
        )?;
        let index_bytes = read_crai(
            &index,
            layouts.len(),
            std::fs::metadata(&request.alignments)?.len(),
            request.execution.memory_budget_bytes,
        )?;
        let mut envelope = scan(
            &request.alignments,
            &layouts,
            &request.execution,
            Some(&index),
        )?;
        envelope.index_bytes = index_bytes;
        snapshot.verify()?;
        // This phase has checked structural bounds but has not yet proven a
        // per-record payload bound. Its full-file native decoding is cooperative:
        // htslib can allocate before the next record checkpoint. Never describe
        // it as a hard prospective allocation guarantee.
        admit(
            request.execution.memory_budget_bytes,
            envelope
                .structural_bytes()?
                .checked_add(
                    (request.execution.max_record_bytes as u64)
                        .checked_mul(3)
                        .ok_or(EvidenceError::CounterOverflow)?,
                )
                .ok_or(EvidenceError::CounterOverflow)?
                .checked_add(8 << 20)
                .ok_or(EvidenceError::CounterOverflow)?,
        )?;
        validate_records(request, fasta, &mut envelope)?;
        snapshot.verify()?;
        let validated = Self {
            envelope,
            snapshot,
            index,
        };
        validated.admit_decoder(&request.execution)?;
        Ok(validated)
    }

    pub fn admit_decoder(&self, execution: &EvidenceExecution) -> Result<(), EvidenceError> {
        if self.envelope.validated_max_record_bytes > execution.max_record_bytes as u64
            || self.envelope.validated_max_read_length > execution.max_read_len as u64
        {
            return Err(EvidenceError::RecordLimit(
                "worker limits exclude records in the validated whole-file CRAM session".into(),
            ));
        }
        admit(
            execution.memory_budget_bytes,
            self.envelope
                .additional_bytes(execution)?
                .checked_add(
                    (execution.max_record_bytes as u64)
                        .checked_mul(3)
                        .ok_or(EvidenceError::CounterOverflow)?,
                )
                .ok_or(EvidenceError::CounterOverflow)?
                .checked_add(8 << 20)
                .ok_or(EvidenceError::CounterOverflow)?,
        )
    }
}

fn validate_records(
    request: &EvidenceRequest,
    fasta: &Path,
    envelope: &mut CramEnvelope,
) -> Result<(), EvidenceError> {
    let started = std::time::Instant::now();
    let mut reader = bam::Reader::from_path(&request.alignments)
        .map_err(|error| invalid(&format!("open whole-file validation reader: {error}")))?;
    reader
        .set_reference(fasta)
        .map_err(|error| invalid(&format!("set explicit validation reference: {error}")))?;
    check_peak(request.execution.memory_budget_bytes)?;
    let mut record = bam::Record::new();
    while let Some(result) = reader.read(&mut record) {
        result.map_err(|error| invalid(&format!("whole-file validation decode: {error}")))?;
        crate::core::governor::checkpoint()?;
        let bytes = record.inner().l_data;
        if bytes < 0 || bytes as usize > request.execution.max_record_bytes {
            return Err(EvidenceError::RecordLimit(format!("whole-file CRAM record payload {bytes} exceeds max_record_bytes {} (including off-target records)", request.execution.max_record_bytes)));
        }
        let length = record.seq_len();
        if length > request.execution.max_read_len {
            return Err(EvidenceError::RecordLimit(format!("whole-file CRAM read length {length} exceeds max_read_len {} (including off-target records)", request.execution.max_read_len)));
        }
        envelope.validated_records = envelope
            .validated_records
            .checked_add(1)
            .ok_or(EvidenceError::CounterOverflow)?;
        envelope.validated_max_record_bytes = envelope.validated_max_record_bytes.max(bytes as u64);
        envelope.validated_max_read_length = envelope.validated_max_read_length.max(length as u64);
        envelope.validated_bases = envelope
            .validated_bases
            .checked_add(length as u64)
            .ok_or(EvidenceError::CounterOverflow)?;
        check_peak(request.execution.memory_budget_bytes)?;
    }
    if envelope.validated_records != envelope.declared_records {
        return Err(invalid(
            "whole-file decoded record count differs from containers",
        ));
    }
    if envelope.validated_bases != envelope.declared_bases {
        return Err(invalid(
            "whole-file decoded base count differs from containers",
        ));
    }
    drop(reader);
    check_peak(request.execution.memory_budget_bytes)?;
    envelope.validation_wall_micros = started.elapsed().as_micros().try_into().unwrap_or(u64::MAX);
    envelope.validation_peak_rss_bytes = crate::util::rss::peak_rss_bytes();
    Ok(())
}
fn check_peak(budget: Option<u64>) -> Result<(), EvidenceError> {
    crate::core::governor::checkpoint()?;
    if let Some(budget) = budget {
        let needed = crate::util::rss::peak_rss_bytes();
        if needed > budget {
            return Err(crate::core::CoreError::BudgetExceeded { needed, budget }.into());
        }
    }
    Ok(())
}

#[derive(Debug)]
struct FaiLayout {
    length: u64,
    line_bases: u64,
    line_bytes: u64,
}
impl FaiLayout {
    fn bytes(&self, span: u64) -> Result<u64, EvidenceError> {
        // A requested portion can start/end on arbitrary line boundaries.
        span.checked_add(
            (span / self.line_bases + 2)
                .checked_mul(self.line_bytes - self.line_bases)
                .ok_or(EvidenceError::CounterOverflow)?,
        )
        .ok_or(EvidenceError::CounterOverflow)
    }
}
fn read_fai(
    path: &Path,
    fasta_length: u64,
    budget: Option<u64>,
) -> Result<BTreeMap<String, FaiLayout>, EvidenceError> {
    let mut reader = BufReader::new(File::open(path)?);
    let mut result = BTreeMap::new();
    let mut total = 0u64;
    loop {
        // Admit both the next bounded line and conservative retained map/string
        // overhead before allocating. Even malformed tab-heavy input stays small.
        admit(
            budget,
            total
                .checked_mul(16)
                .and_then(|n| n.checked_add(SCAN_SCRATCH + FAI_LINE_CAP * 2))
                .ok_or(EvidenceError::CounterOverflow)?,
        )?;
        let mut line = String::new();
        let count = reader
            .by_ref()
            .take(FAI_LINE_CAP + 1)
            .read_line(&mut line)?;
        if count == 0 {
            break;
        }
        if count as u64 > FAI_LINE_CAP {
            return Err(unsupported("FAI line exceeds 8192 bytes"));
        }
        total += count as u64;
        if total > METADATA_CAP {
            return Err(unsupported("FAI metadata exceeds 8MiB"));
        }
        let fields: Vec<_> = line
            .trim_end_matches(['\r', '\n'])
            .splitn(6, '\t')
            .collect();
        if fields.len() != 5 {
            return Err(invalid("invalid FASTA index"));
        }
        let number = |i: usize| {
            fields[i]
                .parse::<u64>()
                .map_err(|_| invalid("invalid FASTA index integer"))
        };
        let length = number(1)?;
        let offset = number(2)?;
        let layout = FaiLayout {
            length,
            line_bases: number(3)?,
            line_bytes: number(4)?,
        };
        let last_base = length.checked_sub(1).and_then(|last| {
            last.checked_div(layout.line_bases)
                .and_then(|line| line.checked_mul(layout.line_bytes))
                .and_then(|base| base.checked_add(last.checked_rem(layout.line_bases)?))
                .and_then(|base| base.checked_add(1))
        });
        if fields[0].is_empty()
            || length == 0
            || layout.line_bases == 0
            || layout.line_bytes < layout.line_bases
            || offset >= fasta_length
            || last_base.is_none_or(|extent| extent > fasta_length - offset)
            || result.insert(fields[0].to_owned(), layout).is_some()
        {
            return Err(invalid("invalid or duplicate FASTA index entry"));
        }
    }
    if result.is_empty() {
        return Err(invalid("empty FASTA index"));
    }
    Ok(result)
}

fn read_crai(
    path: &Path,
    contigs: usize,
    alignment_bytes: u64,
    budget: Option<u64>,
) -> Result<u64, EvidenceError> {
    let file = File::open(path)?;
    let compressed = file.metadata()?.len();
    if compressed > BLOCK_CAP {
        return Err(unsupported("CRAI file exceeds 64MiB"));
    }
    let mut file = BufReader::new(file);
    let gzip = file.fill_buf()?.starts_with(&[31, 139]);
    let input: Box<dyn Read> = if gzip {
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };
    let mut input = BufReader::new(input);
    let mut decoded = 0u64;
    let mut entries = 0u64;
    let mut containment = Vec::<(i64, i64, i64)>::with_capacity(65);
    loop {
        let reservation = compressed
            .checked_mul(2)
            .and_then(|n| n.checked_add(decoded.checked_mul(4)?))
            .and_then(|n| n.checked_add(entries.checked_mul(4096)?))
            .and_then(|n| n.checked_add((contigs as u64).checked_mul(256)?))
            .and_then(|n| n.checked_add(SCAN_SCRATCH))
            .ok_or(EvidenceError::CounterOverflow)?;
        admit(budget, reservation)?;
        let mut line = String::new();
        let count = input.by_ref().take(257).read_line(&mut line)?;
        if count == 0 {
            return Ok(reservation);
        }
        decoded += count as u64;
        entries += 1;
        if count > 256 || decoded > BLOCK_CAP || entries > RECORDS_CAP {
            return Err(unsupported(
                "CRAI text exceeds bounded line/byte/entry envelope",
            ));
        }
        let mut values = [0i64; 6];
        let mut fields = line.trim_end_matches(['\r', '\n']).splitn(7, '\t');
        for value in &mut values {
            *value = fields
                .next()
                .ok_or_else(|| invalid("truncated CRAI entry"))?
                .parse()
                .map_err(|_| invalid("invalid CRAI integer"))?;
        }
        if fields.next().is_some()
            || values[0] < -1
            || values[0] >= contigs as i64
            || values[1] < 0
            || values[1] > i32::MAX as i64
            || values[2] < 0
            || values[2] > i32::MAX as i64 - values[1]
            || values[3] < 0
            || values[3] as u64 >= alignment_bytes
            || values[4] < 0
            || values[4] as u64 > CONTAINER_CAP
            || values[5] <= 0
            || values[5] as u64 > CONTAINER_CAP
        {
            return Err(invalid("CRAI entry exceeds dictionary/file envelope"));
        }
        let end = values[1] + values[2] - 1;
        while containment.last().is_some_and(|&(id, start, last)| {
            id != values[0] || values[1] < start || end > last || (start == 0 && id == -1)
        }) {
            containment.pop();
        }
        if containment.len() >= 64 {
            return Err(unsupported("CRAI containment depth exceeds 64"));
        }
        containment.push((values[0], values[1], end));
    }
}

fn open_crai(path: &Path) -> Result<BufReader<Box<dyn Read>>, EvidenceError> {
    let mut file = BufReader::new(File::open(path)?);
    let input: Box<dyn Read> = if file.fill_buf()?.starts_with(&[31, 139]) {
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };
    Ok(BufReader::new(input))
}
fn next_crai(reader: &mut BufReader<Box<dyn Read>>) -> Result<Option<[i64; 6]>, EvidenceError> {
    crate::core::governor::checkpoint()?;
    let mut line = String::new();
    let n = reader.by_ref().take(257).read_line(&mut line)?;
    if n == 0 {
        return Ok(None);
    }
    if n > 256 {
        return Err(invalid("CRAI changed or contains oversized line"));
    }
    let mut fields = line.trim_end_matches(['\r', '\n']).splitn(7, '\t');
    let mut values = [0; 6];
    for value in &mut values {
        *value = fields
            .next()
            .ok_or_else(|| invalid("truncated CRAI"))?
            .parse()
            .map_err(|_| invalid("invalid CRAI integer"))?;
    }
    if fields.next().is_some() {
        return Err(invalid("extra CRAI fields"));
    }
    Ok(Some(values))
}

struct Checked<R> {
    inner: R,
    crc: Crc,
    offset: u64,
}
impl<R: Read> Read for Checked<R> {
    fn read(&mut self, bytes: &mut [u8]) -> std::io::Result<usize> {
        let count = self.inner.read(bytes)?;
        self.crc.update(&bytes[..count]);
        self.offset += count as u64;
        Ok(count)
    }
}
impl<R: Read> Checked<R> {
    fn byte(&mut self) -> Result<u8, EvidenceError> {
        let mut b = [0];
        self.read_exact(&mut b)?;
        Ok(b[0])
    }
    fn itf8(&mut self) -> Result<i32, EvidenceError> {
        let b = self.byte()?;
        let n = if b < 128 {
            0
        } else if b < 192 {
            1
        } else if b < 224 {
            2
        } else if b < 240 {
            3
        } else {
            4
        };
        let mut value = u32::from(if n == 4 {
            b & 15
        } else {
            b & ((1 << (7 - n)) - 1)
        });
        for i in 0..n {
            let next = self.byte()?;
            value = if n == 4 && i == 3 {
                (value << 4) | u32::from(next & 15)
            } else {
                (value << 8) | u32::from(next)
            };
        }
        Ok(value as i32)
    }
    fn ltf8(&mut self) -> Result<u64, EvidenceError> {
        let b = self.byte()?;
        let n = b.leading_ones();
        let mut value = if n == 8 {
            0
        } else {
            u64::from(b & ((1 << (7 - n)) - 1))
        };
        for _ in 0..n {
            value = (value << 8) | u64::from(self.byte()?);
        }
        Ok(value)
    }
    fn size(&mut self, cap: u64, label: &str) -> Result<u64, EvidenceError> {
        let value = self.itf8()?;
        if value < 0 || value as u64 > cap {
            return Err(unsupported(label));
        }
        Ok(value as u64)
    }
    fn check_crc(&mut self) -> Result<(), EvidenceError> {
        let expected = self.crc.sum();
        let mut raw = [0; 4];
        self.read_exact(&mut raw)?;
        if u32::from_le_bytes(raw) != expected {
            return Err(invalid("CRC32 mismatch"));
        }
        Ok(())
    }
}

#[derive(Debug)]
struct Block {
    kind: u8,
    compressed: u64,
    uncompressed: u64,
    metadata: Vec<u8>,
}
fn block<R: Read>(reader: &mut Checked<R>, budget: Option<u64>) -> Result<Block, EvidenceError> {
    reader.crc.reset();
    let method = reader.byte()?;
    let kind = reader.byte()?;
    let _id = reader.itf8()?;
    let compressed = reader.size(BLOCK_CAP, "compressed block exceeds 64MiB")?;
    let uncompressed = reader.size(BLOCK_CAP, "uncompressed block exceeds 64MiB")?;
    if !matches!(kind, 0..=5) {
        return Err(invalid("unknown block content type"));
    }
    if !matches!(method, 0 | 1 | 4) {
        return Err(unsupported("only RAW, GZIP and rANS4 codecs are admitted"));
    }
    let is_metadata = kind <= 3;
    if is_metadata && (uncompressed > METADATA_CAP || method == 4) {
        return Err(unsupported(
            "metadata blocks require RAW/GZIP and <=8MiB decoded bytes",
        ));
    }
    let retained = if is_metadata { uncompressed } else { 0 };
    let factor = if kind == 1 { 1024 } else { 16 };
    admit(
        budget,
        retained
            .checked_mul(factor)
            .and_then(|n| n.checked_add(SCAN_SCRATCH))
            .ok_or(EvidenceError::CounterOverflow)?,
    )?;
    let mut metadata = Vec::with_capacity(retained as usize);
    let mut scratch = [0u8; 64 * 1024];
    let mut input = reader.take(compressed);
    match method {
        0 => {
            if compressed != uncompressed {
                return Err(invalid("RAW block sizes disagree"));
            }
            loop {
                crate::core::governor::checkpoint()?;
                let n = input.read(&mut scratch)?;
                if n == 0 {
                    break;
                }
                if is_metadata {
                    metadata.extend_from_slice(&scratch[..n]);
                }
            }
        }
        1 => {
            let mut decoder = MultiGzDecoder::new(&mut input);
            let mut size = 0u64;
            loop {
                crate::core::governor::checkpoint()?;
                let n = decoder.read(&mut scratch)?;
                if n == 0 {
                    break;
                }
                size += n as u64;
                if size > uncompressed {
                    return Err(invalid("GZIP expansion exceeds declared block size"));
                }
                if is_metadata {
                    metadata.extend_from_slice(&scratch[..n]);
                }
            }
            if size != uncompressed {
                return Err(invalid("GZIP decoded length disagrees with block"));
            }
        }
        4 => {
            if compressed < 9 {
                return Err(invalid("truncated rANS prefix"));
            }
            let mut prefix = [0; 9];
            input.read_exact(&mut prefix)?;
            if prefix[0] > 1
                || u64::from(u32::from_le_bytes(prefix[1..5].try_into().unwrap())) != compressed - 9
                || u64::from(u32::from_le_bytes(prefix[5..9].try_into().unwrap())) != uncompressed
            {
                return Err(invalid("rANS order or in-band sizes disagree with block"));
            }
            loop {
                crate::core::governor::checkpoint()?;
                let n = input.read(&mut scratch)?;
                if n == 0 {
                    break;
                }
            }
        }
        _ => unreachable!(),
    }
    if input.limit() != 0 {
        return Err(invalid("truncated compressed block"));
    }
    input.into_inner().check_crc()?;
    Ok(Block {
        kind,
        compressed,
        uncompressed,
        metadata,
    })
}

fn scan(
    path: &Path,
    layouts: &BTreeMap<String, FaiLayout>,
    execution: &EvidenceExecution,
    index: Option<&Path>,
) -> Result<CramEnvelope, EvidenceError> {
    let mut index = index.map(open_crai).transpose()?;
    let file = File::open(path)?;
    let file_length = file.metadata()?.len();
    let mut reader = Checked {
        inner: BufReader::new(file),
        crc: Crc::new(),
        offset: 0,
    };
    let mut magic = [0; 26];
    reader.read_exact(&mut magic)?;
    if &magic[..6] != b"CRAM\x03\x00" {
        return Err(unsupported("only CRAM3.0 has a validated decoder envelope"));
    }
    let mut envelope = CramEnvelope {
        model: "cram-container-envelope-v1",
        ..Default::default()
    };
    let mut dictionary = Vec::<(String, u64)>::new();
    let mut max_partial_ref = 0;
    let mut max_cached_ref = 0;
    let mut seen_header = false;
    let mut seen_eof = false;
    while reader.offset < file_length {
        crate::core::governor::checkpoint()?;
        if seen_eof {
            return Err(invalid("data follows CRAM EOF container"));
        }
        let container_offset = reader.offset;
        reader.crc.reset();
        let mut raw = [0; 4];
        reader.read_exact(&mut raw)?;
        let length = u32::from_le_bytes(raw) as u64;
        if length > CONTAINER_CAP {
            return Err(unsupported("container bytes exceed 256MiB"));
        }
        let ref_id = reader.itf8()?;
        let start = reader.size(u32::MAX as u64, "invalid container start")?;
        let span = reader.size(u32::MAX as u64, "invalid container span")?;
        let records = reader.size(RECORDS_CAP, "container exceeds one million records")?;
        let _counter = reader.ltf8()?;
        let bases = reader.ltf8()?;
        let blocks = reader.size(BLOCKS_CAP, "container exceeds 4096 blocks")?;
        let landmarks = reader.size(1024, "container exceeds 1024 slices")?;
        let mut previous = None;
        let mut slice_offsets = Vec::with_capacity(landmarks as usize);
        for _ in 0..landmarks {
            let point = reader.size(length, "invalid landmark offset")?;
            if previous.is_some_and(|p| p >= point) {
                return Err(invalid("unordered container landmarks"));
            }
            previous = Some(point);
            slice_offsets.push(point);
        }
        reader.check_crc()?;
        let end = reader
            .offset
            .checked_add(length)
            .ok_or(EvidenceError::CounterOverflow)?;
        let container_start = reader.offset;
        if end > file_length {
            return Err(invalid("container extends past file"));
        }
        if ref_id < -1 {
            return Err(unsupported(
                "multi-reference containers are not admitted by envelope v1",
            ));
        }
        if records > 0 && bases > records.saturating_mul(execution.max_read_len as u64) {
            return Err(EvidenceError::RecordLimit(
                "CRAM container mean read length exceeds declared max_read_len".into(),
            ));
        }
        let mut compressed = 0u64;
        let mut uncompressed = 0u64;
        let mut slice_records = 0u64;
        let mut slices = 0;
        let mut slice_remaining = 0;
        let mut seen_compression = false;
        for _ in 0..blocks {
            let block_start = reader.offset - container_start;
            let decoded = block(&mut reader, execution.memory_budget_bytes)?;
            compressed = compressed
                .checked_add(decoded.compressed)
                .ok_or(EvidenceError::CounterOverflow)?;
            uncompressed = uncompressed
                .checked_add(decoded.uncompressed)
                .ok_or(EvidenceError::CounterOverflow)?;
            if uncompressed > CONTAINER_CAP {
                return Err(unsupported("container uncompressed blocks exceed 256MiB"));
            }
            match decoded.kind {
                0 => {
                    if records != 0 || envelope.containers != 0 {
                        return Err(invalid("misplaced SAM header block"));
                    }
                    envelope.header_bytes = envelope
                        .header_bytes
                        .checked_add(decoded.uncompressed)
                        .ok_or(EvidenceError::CounterOverflow)?;
                    if !seen_header {
                        seen_header = true;
                        dictionary = parse_sam_header(&decoded.metadata, layouts)?;
                    }
                }
                1 => {
                    if seen_compression || slices > 0 || slice_remaining > 0 {
                        return Err(invalid("misplaced or duplicate compression header"));
                    }
                    seen_compression = true;
                    compression_header(&decoded.metadata)?;
                    envelope.max_compression_header_bytes = envelope
                        .max_compression_header_bytes
                        .max(decoded.uncompressed);
                }
                2 | 3 => {
                    let (id, position, width, count, inventory) =
                        slice_header(&decoded.metadata, decoded.kind)?;
                    if !seen_compression
                        || slice_remaining != 0
                        || slice_offsets.get(slices as usize) != Some(&block_start)
                    {
                        return Err(invalid("slice block inventory or landmark disagrees"));
                    }
                    slice_remaining = inventory;
                    if let Some(index) = &mut index {
                        let slice_end = slice_offsets
                            .get(slices as usize + 1)
                            .copied()
                            .unwrap_or(length);
                        let expected = [
                            i64::from(id),
                            position as i64,
                            width as i64,
                            container_offset as i64,
                            block_start as i64,
                            (slice_end - block_start) as i64,
                        ];
                        if next_crai(index)? != Some(expected) {
                            return Err(invalid(
                                "CRAI entry disagrees with checked slice coordinates or offsets",
                            ));
                        }
                    }
                    if id != ref_id
                        || count == 0
                        || count > records
                        || (id >= 0 && (position < start || position + width > start + span))
                    {
                        return Err(invalid("slice metadata exceeds its container"));
                    }
                    slice_records += count;
                    slices += 1;
                    if id >= 0 {
                        let (name, contig_length) =
                            dictionary.get(id as usize).ok_or_else(|| {
                                invalid("slice reference id is absent from SAM header")
                            })?;
                        if position == 0 || width == 0 || position - 1 + width > *contig_length {
                            return Err(invalid("slice reference range exceeds dictionary"));
                        }
                        let layout = &layouts[name];
                        // htslib promotes spans over half a contig to a cached
                        // full reference. That cache can coexist with old/new
                        // small reference buffers on later single-ref queries.
                        if width.saturating_sub(1) >= contig_length / 2 {
                            max_cached_ref = max_cached_ref.max(layout.bytes(*contig_length)?);
                        } else {
                            max_partial_ref = max_partial_ref.max(layout.bytes(width)?);
                        }
                    }
                }
                4 | 5 if records > 0 => {
                    if slice_remaining == 0 {
                        return Err(invalid("data block is outside a slice"));
                    }
                    slice_remaining -= 1;
                }
                _ => {}
            }
        }
        if records == 0 && !seen_eof && ref_id >= 0 && reader.offset < end {
            let padding = end - reader.offset;
            if padding > METADATA_CAP {
                return Err(unsupported("SAM header padding exceeds 8MiB"));
            }
            envelope.header_bytes = envelope
                .header_bytes
                .checked_add(padding)
                .ok_or(EvidenceError::CounterOverflow)?;
            let mut remaining = reader.by_ref().take(padding);
            std::io::copy(&mut remaining, &mut std::io::sink())?;
        }
        if reader.offset != end {
            return Err(invalid("container length or block count disagrees"));
        }
        // Header and EOF containers are decoded too. Their compressed padding
        // may dominate a tiny data container and must be admitted before open.
        envelope.max_compressed_bytes = envelope.max_compressed_bytes.max(compressed);
        envelope.max_uncompressed_bytes = envelope.max_uncompressed_bytes.max(uncompressed);
        envelope.max_container_blocks = envelope.max_container_blocks.max(blocks);
        envelope.max_container_slices = envelope.max_container_slices.max(landmarks);
        if records > 0 {
            if !seen_header
                || slice_records != records
                || slices != landmarks
                || slice_remaining != 0
            {
                return Err(invalid("container/slice record counts disagree"));
            }
            envelope.containers += 1;
            envelope.declared_records = envelope
                .declared_records
                .checked_add(records)
                .ok_or(EvidenceError::CounterOverflow)?;
            envelope.declared_bases = envelope
                .declared_bases
                .checked_add(bases)
                .ok_or(EvidenceError::CounterOverflow)?;
            envelope.max_container_records = envelope.max_container_records.max(records);
            envelope.max_container_bases = envelope.max_container_bases.max(bases);
        } else if ref_id == -1 && start == 0x454f46 {
            seen_eof = true;
        }
    }
    if !seen_header || !seen_eof {
        return Err(invalid("missing SAM header or CRAM EOF container"));
    }
    if index
        .as_mut()
        .map(next_crai)
        .transpose()?
        .flatten()
        .is_some()
    {
        return Err(invalid("CRAI has entries absent from CRAM"));
    }
    envelope.reference_bytes = max_cached_ref
        .checked_add(
            max_partial_ref
                .checked_mul(2)
                .ok_or(EvidenceError::CounterOverflow)?,
        )
        .ok_or(EvidenceError::CounterOverflow)?;
    Ok(envelope)
}

// These parsers do not instantiate native codecs. In particular htslib's Huffman
// constructor allocates its symbol table before validating the remaining bytes;
// the byte-proportional codec allowance is valid only after this checked pass.
struct Metadata<'a> {
    bytes: &'a [u8],
}
impl<'a> Metadata<'a> {
    fn take(&mut self, count: usize) -> Result<&'a [u8], EvidenceError> {
        let part = self
            .bytes
            .get(..count)
            .ok_or_else(|| invalid("truncated compression metadata"))?;
        self.bytes = &self.bytes[count..];
        Ok(part)
    }
    fn integer(&mut self) -> Result<i32, EvidenceError> {
        let mut reader = Checked {
            inner: Cursor::new(self.bytes),
            crc: Crc::new(),
            offset: 0,
        };
        let value = reader.itf8()?;
        self.take(reader.offset as usize)?;
        Ok(value)
    }
    fn count(&mut self, cap: usize) -> Result<usize, EvidenceError> {
        let value = self.integer()?;
        if value < 0 || value as usize > cap {
            return Err(invalid(
                "compression metadata count exceeds available bytes",
            ));
        }
        Ok(value as usize)
    }
    fn section(&mut self) -> Result<Metadata<'a>, EvidenceError> {
        let size = self.count(self.bytes.len())?;
        Ok(Metadata {
            bytes: self.take(size)?,
        })
    }
    fn finish(self) -> Result<(), EvidenceError> {
        if self.bytes.is_empty() {
            Ok(())
        } else {
            Err(invalid("trailing compression metadata"))
        }
    }
}

fn codec(encoding: i32, bytes: &[u8], depth: usize) -> Result<(), EvidenceError> {
    if depth > 16 {
        return Err(unsupported("codec nesting exceeds 16 levels"));
    }
    let mut metadata = Metadata { bytes };
    match encoding {
        0 => {} // NULL is admitted only with no parameters.
        1 => {
            metadata.integer()?;
        } // EXTERNAL
        3 => {
            let count = metadata.count(metadata.bytes.len().saturating_sub(1) / 2)?;
            for _ in 0..count {
                metadata.integer()?;
            }
            if metadata.count(count)? != count {
                return Err(invalid("Huffman symbol and code-length counts differ"));
            }
            let mut largest = 0;
            for _ in 0..count {
                let length = metadata.count(31)?;
                largest = largest.max(length);
            }
            if count > 0 && largest >= count {
                return Err(invalid("invalid Huffman code lengths"));
            }
        }
        4 => {
            for _ in 0..2 {
                let child = metadata.integer()?;
                if child == 0 {
                    return Err(invalid("NULL nested byte-array codec"));
                }
                let params = metadata.section()?;
                codec(child, params.bytes, depth + 1)?;
            }
        }
        5 => {
            metadata.take(1)?;
            metadata.integer()?;
        }
        6 => {
            metadata.integer()?;
            metadata.count(32)?;
        }
        7 => {
            metadata.integer()?;
            metadata.count(31)?;
        }
        9 => {
            metadata.integer()?;
        }
        _ => return Err(unsupported("unvalidated CRAM encoding codec")),
    }
    metadata.finish()
}

fn compression_header(bytes: &[u8]) -> Result<(), EvidenceError> {
    let mut all = Metadata { bytes };
    let mut preservation = all.section()?;
    let entries = preservation.count(preservation.bytes.len() / 3)?;
    let mut keys = std::collections::BTreeSet::new();
    for _ in 0..entries {
        let key: [u8; 2] = preservation.take(2)?.try_into().unwrap();
        if !keys.insert(key) {
            return Err(invalid("duplicate preservation key"));
        }
        match &key {
            b"MI" | b"UI" | b"PI" | b"RN" | b"AP" | b"RR" | b"QO" => {
                if preservation.take(1)?[0] > 1 {
                    return Err(invalid("invalid preservation boolean"));
                }
            }
            b"SM" => {
                preservation.take(5)?;
            }
            b"TD" => {
                let tags = preservation.section()?;
                if !tags.bytes.is_empty() && tags.bytes.last() != Some(&0) {
                    return Err(invalid("unterminated tag dictionary"));
                }
                for row in tags.bytes.split(|byte| *byte == 0) {
                    if row.len() % 3 != 0 {
                        return Err(invalid("invalid tag dictionary tuple"));
                    }
                }
            }
            _ => return Err(unsupported("unknown preservation map key")),
        }
    }
    preservation.finish()?;
    let mut records = all.section()?;
    let count = records.count(records.bytes.len() / 4)?;
    keys.clear();
    for _ in 0..count {
        let key: [u8; 2] = records.take(2)?.try_into().unwrap();
        if !keys.insert(key) {
            return Err(invalid("duplicate record codec key"));
        }
        if !matches!(
            &key,
            b"BF"
                | b"CF"
                | b"RI"
                | b"RL"
                | b"AP"
                | b"RG"
                | b"MF"
                | b"NS"
                | b"NP"
                | b"TS"
                | b"NF"
                | b"TC"
                | b"TN"
                | b"FN"
                | b"FC"
                | b"FP"
                | b"BS"
                | b"IN"
                | b"SC"
                | b"DL"
                | b"BA"
                | b"BB"
                | b"RS"
                | b"PD"
                | b"HC"
                | b"MQ"
                | b"RN"
                | b"QS"
                | b"QQ"
                | b"TL"
        ) {
            return Err(unsupported("unknown record encoding key"));
        }
        let encoding = records.integer()?;
        let params = records.section()?;
        codec(encoding, params.bytes, 0)?;
    }
    records.finish()?;
    let mut tags = all.section()?;
    let count = tags.count(tags.bytes.len() / 3)?;
    let mut tag_keys = std::collections::BTreeSet::new();
    for _ in 0..count {
        let key = tags.count(0x00ff_ffff)?;
        if !tag_keys.insert(key) {
            return Err(invalid("duplicate tag encoding key"));
        }
        let encoding = tags.integer()?;
        if encoding == 0 {
            return Err(invalid("NULL tag encoding"));
        }
        let params = tags.section()?;
        codec(encoding, params.bytes, 0)?;
    }
    tags.finish()?;
    all.finish()
}

fn parse_sam_header(
    bytes: &[u8],
    layouts: &BTreeMap<String, FaiLayout>,
) -> Result<Vec<(String, u64)>, EvidenceError> {
    let size = bytes
        .get(..4)
        .ok_or_else(|| invalid("truncated SAM header"))?;
    let length = u32::from_le_bytes(size.try_into().unwrap()) as usize;
    if length > bytes.len() - 4 {
        return Err(invalid("SAM header length differs from block"));
    }
    let text = std::str::from_utf8(&bytes[4..4 + length])
        .map_err(|_| invalid("SAM header is not UTF-8"))?;
    let mut dictionary = Vec::new();
    let mut coordinate = false;
    let mut names = std::collections::BTreeSet::new();
    for line in text.lines() {
        if line.starts_with("@HD\t") {
            coordinate = line.split('\t').any(|f| f == "SO:coordinate");
        }
        if !line.starts_with("@SQ\t") {
            continue;
        }
        let mut name = None;
        let mut length = None;
        for field in line.split('\t') {
            if let Some(value) = field.strip_prefix("SN:") {
                if name.replace(value).is_some() {
                    return Err(invalid("duplicate SQ name field"));
                }
            }
            if let Some(value) = field.strip_prefix("LN:") {
                if length
                    .replace(
                        value
                            .parse::<u64>()
                            .map_err(|_| invalid("invalid SQ length"))?,
                    )
                    .is_some()
                {
                    return Err(invalid("duplicate SQ length field"));
                }
            }
        }
        let name = name.ok_or_else(|| invalid("SQ name absent"))?;
        let length = length.ok_or_else(|| invalid("SQ length absent"))?;
        if !names.insert(name)
            || layouts
                .get(name)
                .is_none_or(|layout| layout.length != length)
        {
            return Err(invalid("SAM/FAI dictionaries disagree"));
        }
        dictionary.push((name.to_owned(), length));
    }
    if !coordinate {
        return Err(unsupported("SAM header must declare SO:coordinate"));
    }
    if dictionary.len() != layouts.len() {
        return Err(invalid("SAM/FAI dictionaries differ"));
    }
    Ok(dictionary)
}
fn slice_header(bytes: &[u8], kind: u8) -> Result<(i32, u64, u64, u64, u64), EvidenceError> {
    let mut reader = Checked {
        inner: Cursor::new(bytes),
        crc: Crc::new(),
        offset: 0,
    };
    let (id, start, span) = if kind == 2 {
        (
            reader.itf8()?,
            reader.size(u32::MAX as u64, "invalid slice start")?,
            reader.size(u32::MAX as u64, "invalid slice span")?,
        )
    } else {
        (-1, 0, 0)
    };
    if id < -1 {
        return Err(unsupported(
            "multi-reference slices are not admitted by envelope v1",
        ));
    }
    let count = reader.size(RECORDS_CAP, "slice record count exceeds envelope")?;
    let _counter = reader.ltf8()?;
    let blocks = reader.size(BLOCKS_CAP, "slice block count exceeds envelope")?;
    let ids = reader.size(BLOCKS_CAP, "slice content-id count exceeds envelope")?;
    if ids == 0 || blocks < ids {
        return Err(invalid("invalid slice block inventory"));
    }
    for _ in 0..ids {
        reader.itf8()?;
    }
    if kind == 2 {
        reader.itf8()?;
    }
    let mut md5 = [0; 16];
    reader.read_exact(&mut md5)?;
    Ok((id, start, span, count, blocks))
}
fn admit(budget: Option<u64>, additional: u64) -> Result<(), EvidenceError> {
    crate::core::governor::checkpoint()?;
    if let Some(budget) = budget {
        let needed = crate::util::rss::peak_rss_bytes().saturating_add(additional);
        if needed > budget {
            return Err(EvidenceError::Refused { needed, budget });
        }
    }
    Ok(())
}
fn invalid(message: &str) -> EvidenceError {
    EvidenceError::InvalidInput(format!("CRAM preflight: {message}"))
}
fn unsupported(message: &str) -> EvidenceError {
    EvidenceError::InvalidRequest(format!(
        "CRAM decoder envelope v1 refuses this input: {message}"
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use std::sync::atomic::{AtomicU64, Ordering};
    static NEXT: AtomicU64 = AtomicU64::new(0);
    struct Temp(PathBuf);
    impl Temp {
        fn new(bytes: &[u8]) -> Self {
            let path = std::env::temp_dir().join(format!(
                "rosalind-cram-parser-{}-{}",
                std::process::id(),
                NEXT.fetch_add(1, Ordering::Relaxed)
            ));
            std::fs::write(&path, bytes).unwrap();
            Self(path)
        }
    }
    impl Drop for Temp {
        fn drop(&mut self) {
            let _ = std::fs::remove_file(&self.0);
        }
    }
    fn integer(value: i32, out: &mut Vec<u8>) {
        let n = value as u32;
        if n < 1 << 7 {
            out.push(n as u8);
        } else if n < 1 << 14 {
            out.extend([0x80 | (n >> 8) as u8, n as u8]);
        } else if n < 1 << 21 {
            out.extend([0xc0 | (n >> 16) as u8, (n >> 8) as u8, n as u8]);
        } else if n < 1 << 28 {
            out.extend([
                0xe0 | (n >> 24) as u8,
                (n >> 16) as u8,
                (n >> 8) as u8,
                n as u8,
            ]);
        } else {
            out.extend([
                0xf0 | (n >> 28) as u8,
                (n >> 20) as u8,
                (n >> 12) as u8,
                (n >> 4) as u8,
                (n & 15) as u8,
            ]);
        }
    }
    fn section(bytes: &[u8], out: &mut Vec<u8>) {
        integer(bytes.len() as i32, out);
        out.extend(bytes);
    }
    fn crc(bytes: &mut Vec<u8>) {
        let mut crc = Crc::new();
        crc.update(bytes);
        bytes.extend(crc.sum().to_le_bytes());
    }
    fn raw_block(method: u8, kind: u8, payload: &[u8], size: i32) -> Vec<u8> {
        let mut bytes = vec![method, kind, 0];
        integer(payload.len() as i32, &mut bytes);
        integer(size, &mut bytes);
        bytes.extend(payload);
        crc(&mut bytes);
        bytes
    }
    fn container(id: i32, start: i32, blocks: &[Vec<u8>]) -> Vec<u8> {
        let body: Vec<u8> = blocks.iter().flatten().copied().collect();
        let mut out = (body.len() as u32).to_le_bytes().to_vec();
        for n in [id, start, 0, 0, 0, 0, blocks.len() as i32, 0] {
            integer(n, &mut out);
        }
        crc(&mut out);
        out.extend(body);
        out
    }
    fn parse_block(bytes: &[u8]) -> Result<Block, EvidenceError> {
        block(
            &mut Checked {
                inner: Cursor::new(bytes),
                crc: Crc::new(),
                offset: 0,
            },
            None,
        )
    }
    #[test]
    fn forged_huffman_counts_and_recursive_codec_metadata_fail_before_native_allocation() {
        let mut forged = Vec::new();
        integer(i32::MAX, &mut forged);
        assert!(codec(3, &forged, 0).is_err());
        // One constant symbol with zero-bit coding is valid.
        assert!(codec(3, &[1, 65, 1, 0], 0).is_ok());
        let mut records = vec![1, b'R', b'L', 3];
        section(&forged, &mut records);
        let mut header = vec![1, 0];
        section(&records, &mut header);
        header.extend([1, 0]);
        assert!(compression_header(&header).is_err());
        let mut nested = vec![1, 1, 0, 1, 1, 0];
        for _ in 0..18 {
            let mut outer = vec![4];
            section(&nested, &mut outer);
            outer.extend([1, 1, 0]);
            nested = outer;
        }
        assert!(codec(4, &nested, 0).is_err());
        assert!(compression_header(&[1, 0, 1, 0, 1, 0]).is_ok());
    }
    #[test]
    fn gzip_expansion_rans_internal_lengths_and_crc_are_checked() {
        let mut gzip = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        gzip.write_all(&[0; 8192]).unwrap();
        let payload = gzip.finish().unwrap();
        assert!(parse_block(&raw_block(1, 5, &payload, 1))
            .unwrap_err()
            .to_string()
            .contains("expansion"));
        let mut rans = vec![0];
        rans.extend(0u32.to_le_bytes());
        rans.extend(1000u32.to_le_bytes());
        assert!(parse_block(&raw_block(4, 5, &rans, 1)).is_err());
        let mut bad = raw_block(0, 5, b"abc", 3);
        *bad.last_mut().unwrap() ^= 1;
        assert!(parse_block(&bad).unwrap_err().to_string().contains("CRC32"));
        assert!(parse_block(&raw_block(2, 5, b"abc", 3)).is_err());
    }
    #[test]
    fn fai_and_crai_refuse_unbounded_lines_and_bad_declared_coordinates() {
        let huge = Temp::new(&vec![b'\t'; FAI_LINE_CAP as usize + 1]);
        assert!(matches!(
            read_fai(&huge.0, u64::MAX, Some(1)),
            Err(EvidenceError::Refused { .. })
        ));
        assert!(read_fai(&huge.0, u64::MAX, None).is_err());
        let extent = Temp::new(b"chr1\t100\t5\t10\t11\n");
        assert!(read_fai(&extent.0, 105, None).is_err());
        let forged = Temp::new(b"2147483647\t1\t10\t1\t1\t1\n");
        assert!(read_crai(&forged.0, 1, 100, None).is_err());
    }
    #[test]
    fn header_padding_and_all_container_blocks_are_reserved() {
        let sam = b"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:10\n";
        let mut sam_bytes = (sam.len() as u32).to_le_bytes().to_vec();
        sam_bytes.extend(sam);
        let padding = vec![0; 8192];
        let mut file = b"CRAM\x03\x00".to_vec();
        file.extend([0; 20]);
        file.extend(container(
            0,
            0,
            &[
                raw_block(0, 0, &sam_bytes, sam_bytes.len() as i32),
                raw_block(0, 5, &padding, padding.len() as i32),
            ],
        ));
        file.extend(container(
            -1,
            0x454f46,
            &[raw_block(0, 1, &[1, 0, 1, 0, 1, 0], 6)],
        ));
        let path = Temp::new(&file);
        let layouts = BTreeMap::from([(
            "chr1".into(),
            FaiLayout {
                length: 10,
                line_bases: 10,
                line_bytes: 11,
            },
        )]);
        let envelope = scan(&path.0, &layouts, &EvidenceExecution::default(), None).unwrap();
        assert!(envelope.max_compressed_bytes >= 8192);
        assert!(envelope.max_uncompressed_bytes >= 8192);
        assert_eq!(envelope.max_container_blocks, 2);
        assert!(
            envelope
                .additional_bytes(&EvidenceExecution::default())
                .unwrap()
                >= 4 * 8192
        );
    }
    #[test]
    fn repeated_decoder_bound_uses_global_verified_maxima_and_checked_arithmetic() {
        let mut envelope = CramEnvelope {
            max_container_records: 100,
            validated_max_record_bytes: 300,
            validated_max_read_length: 150,
            ..Default::default()
        };
        let first = envelope
            .additional_bytes(&EvidenceExecution::default())
            .unwrap();
        envelope.max_container_records = 200;
        assert_eq!(
            envelope
                .additional_bytes(&EvidenceExecution::default())
                .unwrap()
                - first,
            100 * (3 * 300 + 512 + 4 * 150)
        );
        envelope.validated_max_record_bytes = u64::MAX;
        assert!(matches!(
            envelope.additional_bytes(&EvidenceExecution::default()),
            Err(EvidenceError::CounterOverflow)
        ));
    }
}
