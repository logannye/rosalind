use super::*;
use arrow_array::{
    ArrayRef, FixedSizeListArray, RecordBatch, StringArray, UInt32Array, UInt64Array,
};
use arrow_ipc::{
    writer::{IpcWriteOptions, StreamWriter},
    MetadataVersion,
};
use arrow_schema::{DataType, Field, Schema};
use std::io::Write;
use std::sync::Arc;

/// Canonical Arrow batch size, independent of computation tile size or budget.
pub const EVIDENCE_ARROW_BATCH_ROWS: usize = 1024;

const SCALAR_NAMES: [&str; 28] = [
    "prefilter_depth",
    "aligned_depth",
    "callable_depth",
    "a",
    "c",
    "g",
    "t",
    "a_fwd",
    "a_rev",
    "c_fwd",
    "c_rev",
    "g_fwd",
    "g_rev",
    "t_fwd",
    "t_rev",
    "base_quality_sum",
    "mapping_quality_sum",
    "read_position_sum",
    "read_length_sum",
    "filtered_secondary",
    "filtered_supplementary",
    "filtered_qc_fail",
    "filtered_duplicate",
    "filtered_unavailable_mapq",
    "filtered_low_mapq",
    "filtered_unavailable_base_quality",
    "filtered_low_base_quality",
    "filtered_ambiguous_base",
];
const ALLELE_SUM_NAMES: [&str; 4] = [
    "allele_base_quality_sum",
    "allele_mapping_quality_sum",
    "allele_read_position_sum",
    "allele_read_length_sum",
];
fn allele_sum<'a>(row: EvidenceRowRef<'a>, index: usize) -> &'a [u64; 4] {
    let values = row.allele_quality.expect("field capability checked");
    match index {
        0 => &values.base_quality_sum,
        1 => &values.mapping_quality_sum,
        2 => &values.read_position_sum,
        3 => &values.read_length_sum,
        _ => unreachable!("known allele summary"),
    }
}
fn scalar_group(index: usize) -> EvidenceFields {
    match index {
        0..=2 | 19..=27 => EvidenceFields::DEPTHS,
        3..=6 => EvidenceFields::ALLELES,
        7..=14 => EvidenceFields::STRANDS,
        15..=16 => EvidenceFields::QUALITY_SUMS,
        17..=18 => EvidenceFields::READ_POSITION,
        _ => unreachable!("known scalar index"),
    }
}
fn scalar(row: EvidenceRowRef<'_>, index: usize) -> u64 {
    match index {
        0 => row.depths.unwrap().prefilter_depth,
        1 => row.depths.unwrap().aligned_depth,
        2 => row.depths.unwrap().callable_depth,
        3..=6 => row.alleles.unwrap().allele_counts[index - 3],
        7..=14 => row.strands.unwrap().strand_counts[(index - 7) / 2][(index - 7) % 2],
        15 => row.quality_sums.unwrap().base_quality_sum,
        16 => row.quality_sums.unwrap().mapping_quality_sum,
        17 => row.read_position.unwrap().read_position_sum,
        18 => row.read_position.unwrap().read_length_sum,
        19..=27 => {
            let f = row.depths.unwrap().filters;
            [
                f.secondary,
                f.supplementary,
                f.qc_fail,
                f.duplicate,
                f.unavailable_mapq,
                f.low_mapq,
                f.unavailable_base_quality,
                f.low_base_quality,
                f.ambiguous_base,
            ][index - 19]
        }
        _ => unreachable!("known scalar index"),
    }
}
fn requirements(fields: EvidenceFields, retained_bytes: Option<u64>) -> EvidenceRequirements {
    EvidenceRequirements {
        fields,
        requires_reference: false,
        context_bases: 0,
        retained_bytes,
    }
}
fn require_fields(actual: EvidenceFields, expected: EvidenceFields) -> Result<(), EvidenceError> {
    if actual != expected {
        return Err(EvidenceError::InvalidInput(
            "evidence batch field mask differs from configured encoder".into(),
        ));
    }
    Ok(())
}
fn write_histogram(out: &mut dyn Write, bins: &[u64]) -> std::io::Result<()> {
    let mut first = true;
    for (bin, count) in bins.iter().enumerate().filter(|(_, count)| **count > 0) {
        if !first {
            write!(out, ",")?;
        }
        write!(out, "{bin}:{count}")?;
        first = false;
    }
    if first {
        write!(out, ".")?;
    }
    Ok(())
}

/// Deterministic, streaming TSV writer. Histograms use sparse `quality:count`
/// pairs in ascending quality order; a dot means an all-zero histogram.
#[derive(Debug)]
pub struct EvidenceTsvWriter<W: Write> {
    out: W,
    started: bool,
    fields: EvidenceFields,
}
impl<W: Write> EvidenceTsvWriter<W> {
    /// Construct the configured value without starting evidence extraction.
    pub fn new(out: W) -> Self {
        Self::with_fields(out, EvidenceFields::ALL)
    }
    /// Configure the physical field groups before any rows, including empty output.
    pub fn with_fields(out: W, fields: EvidenceFields) -> Self {
        Self {
            out,
            started: false,
            fields,
        }
    }
    /// Finish encoding when necessary and return the underlying output writer.
    pub fn into_inner(self) -> W {
        self.out
    }
    fn header(&mut self) -> Result<(), EvidenceError> {
        if !self.started {
            write!(self.out, "#contig\tpos\tref\trequested_alts")?;
            for (index, name) in SCALAR_NAMES.iter().enumerate() {
                if self.fields.contains(scalar_group(index)) {
                    write!(self.out, "\t{name}")?;
                }
            }
            if self.fields.contains(EvidenceFields::QUALITY_HISTOGRAMS) {
                write!(
                    self.out,
                    "\tbase_quality_histogram\tmapping_quality_histogram"
                )?;
            }
            if self.fields.contains(EvidenceFields::ALLELE_QUALITY) {
                for name in ALLELE_SUM_NAMES {
                    write!(self.out, "\t{name}")?;
                }
            }
            writeln!(self.out)?;
            self.started = true;
        }
        Ok(())
    }
}
impl<W: Write> EvidenceAnalyzer for EvidenceTsvWriter<W> {
    fn requirements(&self) -> EvidenceRequirements {
        requirements(self.fields, self.additional_memory_bytes())
    }
    fn on_batch(&mut self, batch: &EvidenceBatch) -> Result<(), EvidenceError> {
        require_fields(batch.fields(), self.fields)?;
        self.header()?;
        for row in batch.rows() {
            write!(
                self.out,
                "{}\t{}\t{}\t",
                batch.contig,
                row.position as u64 + 1,
                char::from(row.reference)
            )?;
            if row.requested_alts.is_empty() {
                write!(self.out, ".")?;
            } else {
                for (i, alternate) in row.requested_alts.iter().enumerate() {
                    if i > 0 {
                        write!(self.out, ",")?;
                    }
                    write!(self.out, "{}", char::from(*alternate))?;
                }
            }
            for index in 0..SCALAR_NAMES.len() {
                if self.fields.contains(scalar_group(index)) {
                    write!(self.out, "\t{}", scalar(row, index))?;
                }
            }
            if let Some(histograms) = row.quality_histograms {
                write!(self.out, "\t")?;
                write_histogram(&mut self.out, &histograms.base_quality_histogram)?;
                write!(self.out, "\t")?;
                write_histogram(&mut self.out, &histograms.mapping_quality_histogram)?;
            }
            if self.fields.contains(EvidenceFields::ALLELE_QUALITY) {
                for index in 0..4 {
                    let values = allele_sum(row, index);
                    write!(
                        self.out,
                        "\t{},{},{},{}",
                        values[0], values[1], values[2], values[3]
                    )?;
                }
            }
            writeln!(self.out)?;
        }
        Ok(())
    }
    fn finish(&mut self) -> Result<(), EvidenceError> {
        self.header()?;
        self.out.flush()?;
        Ok(())
    }
    fn additional_memory_bytes(&self) -> Option<u64> {
        Some(16 << 10)
    }
}

fn schema(selected: EvidenceFields) -> Schema {
    let mut fields = vec![
        Field::new("contig", DataType::Utf8, false),
        Field::new("pos", DataType::UInt32, false),
        Field::new("ref", DataType::Utf8, false),
        Field::new("requested_alts", DataType::Utf8, false),
    ];
    fields.extend(
        SCALAR_NAMES
            .iter()
            .enumerate()
            .filter(|(index, _)| selected.contains(scalar_group(*index)))
            .map(|(_, name)| Field::new(*name, DataType::UInt64, false)),
    );
    if selected.contains(EvidenceFields::QUALITY_HISTOGRAMS) {
        fields.push(Field::new(
            "base_quality_histogram",
            DataType::FixedSizeList(
                Arc::new(Field::new("item", DataType::UInt64, false)),
                BASE_QUALITY_BINS as i32,
            ),
            false,
        ));
        fields.push(Field::new(
            "mapping_quality_histogram",
            DataType::FixedSizeList(
                Arc::new(Field::new("item", DataType::UInt64, false)),
                MAPPING_QUALITY_BINS as i32,
            ),
            false,
        ));
    }
    if selected.contains(EvidenceFields::ALLELE_QUALITY) {
        for name in ALLELE_SUM_NAMES {
            fields.push(Field::new(
                name,
                DataType::FixedSizeList(Arc::new(Field::new("item", DataType::UInt64, false)), 4),
                false,
            ));
        }
    }
    Schema::new(fields).with_metadata(std::collections::HashMap::from([(
        "rosalind.evidence.schema".into(),
        if selected == EvidenceFields::ALL {
            "1".into()
        } else {
            format!("2;fields-v{}={}", selected.mask_version(), selected.bits())
        },
    )]))
}

/// Convert a bounded native batch to Arrow arrays without serializing IPC.
/// The returned arrays allocate independently; callers must reserve their memory.
/// At most 1,024 rows are accepted, matching canonical output batching.
pub fn evidence_record_batch(batch: &EvidenceBatch) -> Result<RecordBatch, EvidenceError> {
    if batch.len() > EVIDENCE_ARROW_BATCH_ROWS {
        return Err(EvidenceError::InvalidRequest(
            "Arrow conversion requires at most 1,024 evidence rows".into(),
        ));
    }
    encode_rows(
        batch,
        std::iter::repeat_n(batch.contig.as_str(), batch.len()),
    )
}

fn encode_rows<'a>(
    rows: &EvidenceBatch,
    contigs: impl Iterator<Item = &'a str>,
) -> Result<RecordBatch, EvidenceError> {
    let mut arrays: Vec<ArrayRef> = Vec::with_capacity(34);
    arrays.push(Arc::new(StringArray::from_iter_values(contigs)));
    arrays.push(Arc::new(UInt32Array::from_iter_values(
        rows.rows().map(|row| row.position + 1),
    )));
    arrays.push(Arc::new(StringArray::from_iter_values(
        rows.rows().map(|row| char::from(row.reference).to_string()),
    )));
    arrays.push(Arc::new(StringArray::from_iter_values(rows.rows().map(
        |row| String::from_utf8(row.requested_alts.to_vec()).expect("validated nucleotide ALT"),
    ))));
    for field in 0..SCALAR_NAMES.len() {
        if !rows.fields().contains(scalar_group(field)) {
            continue;
        }
        arrays.push(Arc::new(UInt64Array::from_iter_values(
            rows.rows().map(|row| scalar(row, field)),
        )));
    }
    if rows.fields().contains(EvidenceFields::QUALITY_HISTOGRAMS) {
        for (bins, base_quality) in [(BASE_QUALITY_BINS, true), (MAPPING_QUALITY_BINS, false)] {
            let values: ArrayRef =
                Arc::new(UInt64Array::from_iter_values(rows.rows().flat_map(|row| {
                    let histogram = row.quality_histograms.unwrap();
                    let values: &[u64] = if base_quality {
                        &histogram.base_quality_histogram
                    } else {
                        &histogram.mapping_quality_histogram
                    };
                    values.iter().copied()
                })));
            arrays.push(Arc::new(
                FixedSizeListArray::try_new(
                    Arc::new(Field::new("item", DataType::UInt64, false)),
                    bins as i32,
                    values,
                    None,
                )
                .map_err(arrow_error)?,
            ));
        }
    }
    if rows.fields().contains(EvidenceFields::ALLELE_QUALITY) {
        for index in 0..4 {
            let values: ArrayRef = Arc::new(UInt64Array::from_iter_values(
                rows.rows()
                    .flat_map(|row| allele_sum(row, index).iter().copied()),
            ));
            arrays.push(Arc::new(
                FixedSizeListArray::try_new(
                    Arc::new(Field::new("item", DataType::UInt64, false)),
                    4,
                    values,
                    None,
                )
                .map_err(arrow_error)?,
            ));
        }
    }
    RecordBatch::try_new(Arc::new(schema(rows.fields())), arrays).map_err(arrow_error)
}

/// Canonical Arrow IPC writer. Buffering is limited to 1024 evidence rows plus
/// encoding scratch. Construction performs no writes, so admission can run first.
pub struct EvidenceArrowWriter<W: Write> {
    output: Option<W>,
    writer: Option<StreamWriter<W>>,
    rows: EvidenceBatch,
    contigs: Vec<Arc<str>>,
    fields: EvidenceFields,
    last_contig: Option<Arc<str>>,
    finished: bool,
}
impl<W: Write> std::fmt::Debug for EvidenceArrowWriter<W> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("EvidenceArrowWriter")
            .field("buffered_rows", &self.rows.len())
            .field("finished", &self.finished)
            .finish()
    }
}
impl<W: Write> EvidenceArrowWriter<W> {
    /// Construct the configured value without starting evidence extraction.
    pub fn new(output: W) -> Self {
        Self::with_fields(output, EvidenceFields::ALL)
    }
    /// Configure the physical field groups before the stream header is written.
    pub fn with_fields(output: W, fields: EvidenceFields) -> Self {
        Self {
            output: Some(output),
            writer: None,
            rows: EvidenceBatch::new(0, "", 0, fields, Vec::new()),
            contigs: Vec::new(),
            fields,
            last_contig: None,
            finished: false,
        }
    }
    fn ensure_writer(&mut self) -> Result<(), EvidenceError> {
        if self.writer.is_none() {
            let output = self
                .output
                .take()
                .ok_or_else(|| EvidenceError::Analyzer("Arrow output unavailable".into()))?;
            let options =
                IpcWriteOptions::try_new(8, false, MetadataVersion::V5).map_err(arrow_error)?;
            self.writer = Some(
                StreamWriter::try_new_with_options(output, &schema(self.fields), options)
                    .map_err(arrow_error)?,
            );
        }
        Ok(())
    }
    fn flush_batch(&mut self) -> Result<(), EvidenceError> {
        self.ensure_writer()?;
        if self.rows.is_empty() {
            return Ok(());
        }
        let batch = encode_rows(
            &self.rows,
            self.contigs.iter().map(|contig| contig.as_ref()),
        )?;
        self.writer
            .as_mut()
            .unwrap()
            .write(&batch)
            .map_err(arrow_error)?;
        self.rows.clear();
        self.contigs.clear();
        Ok(())
    }
    /// Finish encoding when necessary and return the underlying output writer.
    pub fn into_inner(mut self) -> Result<W, EvidenceError> {
        self.finish()?;
        self.writer
            .take()
            .unwrap()
            .into_inner()
            .map_err(arrow_error)
    }
}
impl<W: Write> EvidenceAnalyzer for EvidenceArrowWriter<W> {
    fn requirements(&self) -> EvidenceRequirements {
        requirements(self.fields, self.additional_memory_bytes())
    }
    fn on_batch(&mut self, batch: &EvidenceBatch) -> Result<(), EvidenceError> {
        require_fields(batch.fields(), self.fields)?;
        if self.finished {
            return Err(EvidenceError::Analyzer(
                "cannot append to a finished Arrow stream".into(),
            ));
        }
        if batch.contig.len() > 4096 {
            return Err(EvidenceError::InvalidInput(
                "Arrow evidence contig name exceeds 4096 bytes".into(),
            ));
        }
        let contig = match &self.last_contig {
            Some(contig) if contig.as_ref() == batch.contig => Arc::clone(contig),
            _ => {
                let contig: Arc<str> = Arc::from(batch.contig.as_str());
                self.last_contig = Some(Arc::clone(&contig));
                contig
            }
        };
        for row in batch.rows() {
            if self.contigs.capacity() == 0 {
                self.rows.reserve_exact(EVIDENCE_ARROW_BATCH_ROWS);
                self.contigs.reserve_exact(EVIDENCE_ARROW_BATCH_ROWS);
            }
            self.rows.push_row(row)?;
            self.contigs.push(Arc::clone(&contig));
            if self.rows.len() == EVIDENCE_ARROW_BATCH_ROWS {
                self.flush_batch()?;
            }
        }
        Ok(())
    }
    fn finish(&mut self) -> Result<(), EvidenceError> {
        if !self.finished {
            self.flush_batch()?;
            self.writer
                .as_mut()
                .unwrap()
                .finish()
                .map_err(arrow_error)?;
            self.finished = true;
        }
        Ok(())
    }
    fn additional_memory_bytes(&self) -> Option<u64> {
        // Row storage, Arrow arrays, IPC scratch, plus worst-case names for a
        // batch spanning many short contigs. Old feature encoding is unchanged.
        Some(
            self.fields.storage_bytes_per_locus() * EVIDENCE_ARROW_BATCH_ROWS as u64 * 4
                + (16 << 20),
        )
    }
}
fn arrow_error(error: impl std::fmt::Display) -> EvidenceError {
    EvidenceError::Analyzer(error.to_string())
}

/// Physical schema descriptor, available even for an empty evidence stream.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EvidenceArtifactMetadata {
    /// Actual stored field groups.
    pub fields: EvidenceFields,
    /// Evidence schema version, independent of receipt-envelope version.
    pub schema_version: u32,
}

fn body_envelope(fields: EvidenceFields) -> u64 {
    // At most 1024 contig names of 4096 bytes, small REF/ALT/offset arrays,
    // selected counters, alignment padding, and bounded Arrow bookkeeping.
    (EVIDENCE_ARROW_BATCH_ROWS as u64)
        .saturating_mul(4160 + fields.storage_bytes_per_locus())
        .saturating_add(65_536)
        .min(MAX_IPC_BODY_BYTES)
}
/// Conservative retained IPC buffers, decoded arrays, projected callback storage,
/// and metadata. Does not include memory retained by the downstream analyzer.
/// Use ALL_SUPPORTED when the input field set is unknown, or enforce the expected set with
/// `read_evidence_batches_expected_fields` before admitting a smaller reservation.
pub fn evidence_reader_memory_bytes(fields: EvidenceFields) -> u64 {
    2 * body_envelope(fields)
        + fields.storage_bytes_per_locus() * EVIDENCE_ARROW_BATCH_ROWS as u64
        + (1 << 20)
}

/// Stream verified canonical evidence into bounded callbacks.
pub fn read_evidence_batches<R: std::io::Read>(
    reader: R,
    contigs: &crate::core::ContigSet,
    on_batch: impl FnMut(&EvidenceBatch) -> Result<(), EvidenceError>,
) -> Result<(), EvidenceError> {
    read_evidence_batches_with_metadata(reader, contigs, on_batch).map(|_| ())
}

/// Read physical fields without synthesizing omitted metrics; return the schema
/// descriptor even when the stream contains no loci.
pub fn read_evidence_batches_with_metadata<R: std::io::Read>(
    reader: R,
    contigs: &crate::core::ContigSet,
    on_batch: impl FnMut(&EvidenceBatch) -> Result<(), EvidenceError>,
) -> Result<EvidenceArtifactMetadata, EvidenceError> {
    read_evidence_batches_impl(reader, contigs, None, on_batch)
}

/// Check a known field mask before decoding any body buffers. Budgeted readers
/// can reserve `evidence_reader_memory_bytes(expected)` without accepting a
/// larger, differently projected artifact first. Empty streams are also checked.
pub fn read_evidence_batches_expected_fields<R: std::io::Read>(
    reader: R,
    contigs: &crate::core::ContigSet,
    expected: EvidenceFields,
    on_batch: impl FnMut(&EvidenceBatch) -> Result<(), EvidenceError>,
) -> Result<EvidenceArtifactMetadata, EvidenceError> {
    read_evidence_batches_impl(reader, contigs, Some(expected), on_batch)
}

fn read_evidence_batches_impl<R: std::io::Read>(
    reader: R,
    contigs: &crate::core::ContigSet,
    expected: Option<EvidenceFields>,
    mut on_batch: impl FnMut(&EvidenceBatch) -> Result<(), EvidenceError>,
) -> Result<EvidenceArtifactMetadata, EvidenceError> {
    use arrow_array::Array;
    let invalid = |message: &str| EvidenceError::InvalidInput(message.into());
    let mut reader = arrow_ipc::reader::StreamReader::try_new(BoundedIpcReader::new(reader), None)
        .map_err(arrow_error)?;
    let metadata = reader
        .schema()
        .metadata()
        .get("rosalind.evidence.schema")
        .cloned()
        .ok_or_else(|| invalid("missing evidence schema metadata"))?;
    let fields = if metadata == "1" {
        EvidenceFields::ALL
    } else {
        let bits = metadata
            .strip_prefix("2;fields-v")
            .and_then(|value| value.split_once('='))
            .and_then(|(_, bits)| bits.parse::<u32>().ok())
            .ok_or_else(|| invalid("unsupported evidence schema or field mask version"))?;
        EvidenceFields::from_bits(bits)?
    };
    if reader.schema().as_ref() != &schema(fields) {
        return Err(invalid("cached Arrow evidence schema mismatch"));
    }
    if expected.is_some_and(|expected| expected != fields) {
        return Err(invalid(
            "evidence schema field mask differs from expected fields",
        ));
    }
    reader.get_mut().max_body_bytes = body_envelope(fields);
    reader.get_mut().fields = Some(fields);
    let descriptor = EvidenceArtifactMetadata {
        fields,
        schema_version: fields.schema_version(),
    };
    let scalar_indices: Vec<usize> = (0..SCALAR_NAMES.len())
        .filter(|index| fields.contains(scalar_group(*index)))
        .collect();
    let mut previous = None;
    for encoded in &mut reader {
        let encoded = encoded.map_err(arrow_error)?;
        if encoded.num_rows() > EVIDENCE_ARROW_BATCH_ROWS {
            return Err(invalid(
                "cached Arrow evidence batch exceeds canonical 1024 rows",
            ));
        }
        if encoded
            .columns()
            .iter()
            .any(|array| array.null_count() != 0)
        {
            return Err(invalid("cached evidence cannot contain nulls"));
        }
        let strings = |column: usize| {
            encoded
                .column(column)
                .as_any()
                .downcast_ref::<StringArray>()
                .ok_or_else(|| invalid("invalid evidence string array"))
        };
        let names = strings(0)?;
        let positions = encoded
            .column(1)
            .as_any()
            .downcast_ref::<UInt32Array>()
            .ok_or_else(|| invalid("invalid position array"))?;
        let references = strings(2)?;
        let alts = strings(3)?;
        let mut numeric: [Option<&UInt64Array>; 28] = [None; 28];
        for (column, index) in scalar_indices.iter().enumerate() {
            numeric[*index] = Some(
                encoded
                    .column(column + 4)
                    .as_any()
                    .downcast_ref::<UInt64Array>()
                    .ok_or_else(|| invalid("invalid scalar evidence array"))?,
            );
        }
        let histogram = |offset: usize| -> Result<Option<&FixedSizeListArray>, EvidenceError> {
            if !fields.contains(EvidenceFields::QUALITY_HISTOGRAMS) {
                return Ok(None);
            }
            let array = encoded
                .column(4 + scalar_indices.len() + offset)
                .as_any()
                .downcast_ref::<FixedSizeListArray>()
                .ok_or_else(|| invalid("invalid histogram array"))?;
            if array.values().null_count() != 0 {
                return Err(invalid("histograms cannot contain nulls"));
            }
            Ok(Some(array))
        };
        let bq = histogram(0)?;
        let mq = histogram(1)?;
        let mut allele_arrays = Vec::new();
        if fields.contains(EvidenceFields::ALLELE_QUALITY) {
            let start = 4
                + scalar_indices.len()
                + if fields.contains(EvidenceFields::QUALITY_HISTOGRAMS) {
                    2
                } else {
                    0
                };
            for column in start..start + 4 {
                let array = encoded
                    .column(column)
                    .as_any()
                    .downcast_ref::<FixedSizeListArray>()
                    .ok_or_else(|| invalid("invalid per-allele summary array"))?;
                if array.values().null_count() != 0 {
                    return Err(invalid("per-allele summaries cannot contain nulls"));
                }
                allele_arrays.push(array);
            }
        }
        let mut index = 0;
        while index < encoded.num_rows() {
            let start = index;
            let contig = contigs
                .by_name(names.value(index))
                .ok_or_else(|| invalid("cached evidence contig not in dictionary"))?;
            let first_position = positions
                .value(index)
                .checked_sub(1)
                .ok_or_else(|| invalid("evidence POS must be 1-based"))?;
            let canonical_start = first_position / CANONICAL_TILE_BASES * CANONICAL_TILE_BASES;
            let mut loci = Vec::new();
            while index < encoded.num_rows() {
                let position = positions
                    .value(index)
                    .checked_sub(1)
                    .ok_or_else(|| invalid("evidence POS must be 1-based"))?;
                if names.value(index) != contig.name.as_ref()
                    || position / CANONICAL_TILE_BASES * CANONICAL_TILE_BASES != canonical_start
                {
                    break;
                }
                if names.value(index).len() > 4096 || alts.value(index).len() > 3 {
                    return Err(invalid("cached evidence string exceeds canonical bounds"));
                }
                if position >= contig.length {
                    return Err(invalid("cached evidence position out of bounds"));
                }
                let locus = (contig.id, position);
                if previous.is_some_and(|before| locus <= before) {
                    return Err(invalid(
                        "cached evidence rows are not uniquely coordinate sorted",
                    ));
                }
                previous = Some(locus);
                let reference = references.value(index).as_bytes();
                if reference.len() != 1 || !b"ACGTN".contains(&reference[0]) {
                    return Err(invalid("invalid cached reference base"));
                }
                let requested_alts = alts.value(index).as_bytes().to_vec();
                if requested_alts.iter().any(|base| !b"ACGT".contains(base)) {
                    return Err(invalid("invalid cached SNV ALT"));
                }
                loci.push(EvidenceLocus {
                    position,
                    reference: reference[0],
                    requested_alts,
                });
                index += 1;
            }
            let mut batch = EvidenceBatch::new(
                contig.id,
                contig.name.to_string(),
                canonical_start,
                fields,
                loci,
            );
            for offset in 0..batch.len() {
                let source = start + offset;
                let v = |column: usize| {
                    numeric[column]
                        .expect("schema field invariant")
                        .value(source)
                };
                let row = batch.row_mut(offset).unwrap();
                if let Some(depths) = row.depths {
                    *depths = EvidenceDepths {
                        prefilter_depth: v(0),
                        aligned_depth: v(1),
                        callable_depth: v(2),
                        filters: EvidenceFilterCounts {
                            secondary: v(19),
                            supplementary: v(20),
                            qc_fail: v(21),
                            duplicate: v(22),
                            unavailable_mapq: v(23),
                            low_mapq: v(24),
                            unavailable_base_quality: v(25),
                            low_base_quality: v(26),
                            ambiguous_base: v(27),
                        },
                    };
                }
                if let Some(alleles) = row.alleles {
                    alleles.allele_counts = [v(3), v(4), v(5), v(6)];
                }
                if let Some(strands) = row.strands {
                    strands.strand_counts =
                        [[v(7), v(8)], [v(9), v(10)], [v(11), v(12)], [v(13), v(14)]];
                }
                if let Some(quality) = row.quality_sums {
                    quality.base_quality_sum = v(15);
                    quality.mapping_quality_sum = v(16);
                }
                if let Some(position) = row.read_position {
                    position.read_position_sum = v(17);
                    position.read_length_sum = v(18);
                }
                if let Some(histograms) = row.quality_histograms {
                    for (array, destination) in [
                        (
                            bq.unwrap(),
                            histograms.base_quality_histogram.as_mut_slice(),
                        ),
                        (
                            mq.unwrap(),
                            histograms.mapping_quality_histogram.as_mut_slice(),
                        ),
                    ] {
                        let value = array.value(source);
                        let values = value
                            .as_any()
                            .downcast_ref::<UInt64Array>()
                            .ok_or_else(|| invalid("invalid histogram values"))?;
                        destination.copy_from_slice(values.values());
                    }
                }
                if let Some(sums) = row.allele_quality {
                    for (column, target) in [
                        &mut sums.base_quality_sum,
                        &mut sums.mapping_quality_sum,
                        &mut sums.read_position_sum,
                        &mut sums.read_length_sum,
                    ]
                    .into_iter()
                    .enumerate()
                    {
                        let array = allele_arrays[column].value(source);
                        let values = array
                            .as_any()
                            .downcast_ref::<UInt64Array>()
                            .ok_or_else(|| invalid("invalid per-allele summary values"))?;
                        target.copy_from_slice(values.values());
                    }
                }
                validate_row(batch.row(offset).unwrap())?;
            }
            on_batch(&batch)?;
        }
    }
    reader.get_mut().validate_complete()?;
    Ok(descriptor)
}

fn validate_row(row: EvidenceRowRef<'_>) -> Result<(), EvidenceError> {
    let invalid = || {
        EvidenceError::InvalidInput(
            "cached evidence counters or ALT annotations are inconsistent".into(),
        )
    };
    let sum = |values: &[u64]| values.iter().map(|value| u128::from(*value)).sum::<u128>();
    let weighted = |values: &[u64]| {
        values
            .iter()
            .enumerate()
            .map(|(quality, count)| quality as u128 * u128::from(*count))
            .sum::<u128>()
    };
    if row
        .requested_alts
        .windows(2)
        .any(|bases| bases[0] >= bases[1])
        || row.requested_alts.contains(&row.reference)
        || (!row.requested_alts.is_empty() && row.reference == b'N')
    {
        return Err(invalid());
    }
    let mut depth = None;
    let mut check_depth = |value: u128| -> Result<(), EvidenceError> {
        if value > u64::MAX as u128 || depth.is_some_and(|before| before != value) {
            return Err(invalid());
        }
        depth = Some(value);
        Ok(())
    };
    if let Some(d) = row.depths {
        let f = d.filters;
        let read_filtered = sum(&[
            f.secondary,
            f.supplementary,
            f.qc_fail,
            f.duplicate,
            f.unavailable_mapq,
            f.low_mapq,
        ]);
        let base_filtered = sum(&[
            f.unavailable_base_quality,
            f.low_base_quality,
            f.ambiguous_base,
        ]);
        if u128::from(d.prefilter_depth) != u128::from(d.aligned_depth) + read_filtered
            || u128::from(d.aligned_depth) != u128::from(d.callable_depth) + base_filtered
        {
            return Err(invalid());
        }
        check_depth(u128::from(d.callable_depth))?;
    }
    if let Some(a) = row.alleles {
        check_depth(sum(&a.allele_counts))?;
    }
    if let Some(s) = row.strands {
        check_depth(s.strand_counts.iter().map(|pair| sum(pair)).sum())?;
        if let Some(a) = row.alleles {
            if s.strand_counts
                .iter()
                .zip(a.allele_counts)
                .any(|(pair, allele)| sum(pair) != u128::from(allele))
            {
                return Err(invalid());
            }
        }
    }
    if let Some(h) = row.quality_histograms {
        check_depth(sum(&h.base_quality_histogram))?;
        check_depth(sum(&h.mapping_quality_histogram))?;
        if let Some(q) = row.quality_sums {
            if weighted(&h.base_quality_histogram) != u128::from(q.base_quality_sum)
                || weighted(&h.mapping_quality_histogram) != u128::from(q.mapping_quality_sum)
            {
                return Err(invalid());
            }
        }
    }
    if let Some(q) = row.quality_sums {
        if depth.is_some_and(|d| {
            u128::from(q.base_quality_sum) > d * 93 || u128::from(q.mapping_quality_sum) > d * 254
        }) {
            return Err(invalid());
        }
    }
    if let Some(a) = row.allele_quality {
        let total_bq = sum(&a.base_quality_sum);
        let total_mq = sum(&a.mapping_quality_sum);
        let total_position = sum(&a.read_position_sum);
        let total_length = sum(&a.read_length_sum);
        if depth.is_some_and(|count| {
            total_bq > count * 93
                || total_mq > count * 254
                || total_position + count > total_length
                || (count == 0 && total_length != 0)
        }) {
            return Err(invalid());
        }
        if let Some(h) = row.quality_histograms {
            if total_bq != weighted(&h.base_quality_histogram)
                || total_mq != weighted(&h.mapping_quality_histogram)
            {
                return Err(invalid());
            }
        }
        if let Some(q) = row.quality_sums {
            if sum(&a.base_quality_sum) != u128::from(q.base_quality_sum)
                || sum(&a.mapping_quality_sum) != u128::from(q.mapping_quality_sum)
            {
                return Err(invalid());
            }
        }
        if let Some(p) = row.read_position {
            if sum(&a.read_position_sum) != u128::from(p.read_position_sum)
                || sum(&a.read_length_sum) != u128::from(p.read_length_sum)
            {
                return Err(invalid());
            }
        }
        for allele in 0..4 {
            let count = row
                .alleles
                .map(|counts| u128::from(counts.allele_counts[allele]))
                .or_else(|| {
                    row.strands
                        .map(|strands| sum(&strands.strand_counts[allele]))
                });
            let length = u128::from(a.read_length_sum[allele]);
            // Even without counts, each contributing read has positive length
            // and its zero-based position is strictly less than that length.
            if u128::from(a.base_quality_sum[allele]) > length * 93
                || u128::from(a.mapping_quality_sum[allele]) > length * 254
                || (length > 0 && u128::from(a.read_position_sum[allele]) >= length)
            {
                return Err(invalid());
            }
            if count.is_some_and(|count| {
                u128::from(a.base_quality_sum[allele]) > count * 93
                    || u128::from(a.mapping_quality_sum[allele]) > count * 254
                    || (count == 0 && a.read_length_sum[allele] != 0)
            }) || u128::from(a.read_position_sum[allele]) + count.unwrap_or(0)
                > u128::from(a.read_length_sum[allele])
            {
                return Err(invalid());
            }
        }
    }
    if let Some(p) = row.read_position {
        if u128::from(p.read_position_sum) + depth.unwrap_or(0) > u128::from(p.read_length_sum) {
            return Err(invalid());
        }
    }
    Ok(())
}

// Validate each frame before Arrow sees any attacker-controlled allocation
// length. Canonical evidence uses uncompressed record batches, no dictionaries,
// <=1024 rows, and contig strings <=4096 bytes. An 8 MiB body envelope covers
// the resulting histogram, scalar and string buffers without buffering a body.
const MAX_IPC_METADATA_BYTES: usize = 65_536;
const MAX_IPC_BODY_BYTES: u64 = 8 << 20;
struct BoundedIpcReader<R> {
    source: R,
    prefix: std::io::Cursor<Vec<u8>>,
    body_remaining: u64,
    max_body_bytes: u64,
    fields: Option<EvidenceFields>,
    ended: bool,
}
impl<R: std::io::Read> BoundedIpcReader<R> {
    fn new(source: R) -> Self {
        Self {
            source,
            prefix: std::io::Cursor::new(Vec::new()),
            body_remaining: 0,
            max_body_bytes: MAX_IPC_BODY_BYTES,
            fields: None,
            ended: false,
        }
    }
    fn validate_complete(&mut self) -> std::io::Result<()> {
        if !self.ended {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "canonical evidence IPC is missing its end marker",
            ));
        }
        let mut trailing = [0u8; 1];
        if self.source.read(&mut trailing)? != 0 {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "canonical evidence IPC has trailing bytes",
            ));
        }
        Ok(())
    }
    fn next_frame(&mut self) -> std::io::Result<()> {
        let invalid = |message: &str| std::io::Error::new(std::io::ErrorKind::InvalidData, message);
        let mut first = [0u8; 4];
        // Arrow requires a complete end marker; premature EOF remains an error.
        self.source.read_exact(&mut first)?;
        let mut prefix = first.to_vec();
        let length = if first == [255; 4] {
            let mut size = [0u8; 4];
            self.source.read_exact(&mut size)?;
            prefix.extend_from_slice(&size);
            i32::from_le_bytes(size)
        } else {
            i32::from_le_bytes(first)
        };
        if length < 0 || length as usize > MAX_IPC_METADATA_BYTES {
            return Err(invalid(
                "evidence IPC metadata exceeds the bounded frame envelope",
            ));
        }
        if length == 0 {
            self.ended = true;
            self.prefix = std::io::Cursor::new(prefix);
            return Ok(());
        }
        let offset = prefix.len();
        prefix.resize(offset + length as usize, 0);
        self.source.read_exact(&mut prefix[offset..])?;
        let message = arrow_ipc::root_as_message(&prefix[offset..])
            .map_err(|_| invalid("invalid evidence IPC metadata"))?;
        let body = message.bodyLength();
        if body < 0 || body as u64 > self.max_body_bytes {
            return Err(invalid(
                "evidence IPC body exceeds the bounded frame envelope",
            ));
        }
        match message.header_type() {
            arrow_ipc::MessageHeader::Schema if body == 0 => {}
            arrow_ipc::MessageHeader::RecordBatch => {
                let batch = message
                    .header_as_record_batch()
                    .ok_or_else(|| invalid("invalid evidence IPC record batch"))?;
                if batch.length() < 0 || batch.length() as usize > EVIDENCE_ARROW_BATCH_ROWS {
                    return Err(invalid(
                        "evidence IPC record batch exceeds canonical 1024 rows",
                    ));
                }
                if batch.compression().is_some() {
                    return Err(invalid(
                        "canonical evidence IPC does not permit compressed buffers",
                    ));
                }
                let fields = self
                    .fields
                    .ok_or_else(|| invalid("record batch precedes verified schema"))?;
                let scalar_count = (0..SCALAR_NAMES.len())
                    .filter(|index| fields.contains(scalar_group(*index)))
                    .count();
                let histograms = fields.contains(EvidenceFields::QUALITY_HISTOGRAMS);
                let nodes = batch
                    .nodes()
                    .ok_or_else(|| invalid("evidence IPC field nodes missing"))?;
                let buffers = batch
                    .buffers()
                    .ok_or_else(|| invalid("evidence IPC buffers missing"))?;
                let allele_sums = fields.contains(EvidenceFields::ALLELE_QUALITY);
                let histogram_nodes = if histograms { 4 } else { 0 };
                let expected_nodes =
                    4 + scalar_count + histogram_nodes + if allele_sums { 8 } else { 0 };
                let expected_buffers = 11
                    + scalar_count * 2
                    + if histograms { 6 } else { 0 }
                    + if allele_sums { 12 } else { 0 };
                if nodes.len() != expected_nodes || buffers.len() != expected_buffers {
                    return Err(invalid(
                        "evidence IPC node/buffer count differs from projected schema",
                    ));
                }
                if batch
                    .variadicBufferCounts()
                    .is_some_and(|counts| !counts.is_empty())
                {
                    return Err(invalid("canonical evidence IPC has no variadic buffers"));
                }
                let rows = batch.length();
                for (index, node) in nodes.iter().enumerate() {
                    let expected_length = match index.checked_sub(4 + scalar_count) {
                        Some(1) if histograms => rows * BASE_QUALITY_BINS as i64,
                        Some(3) if histograms => rows * MAPPING_QUALITY_BINS as i64,
                        Some(index)
                            if allele_sums
                                && index >= histogram_nodes
                                && (index - histogram_nodes) % 2 == 1 =>
                        {
                            rows * 4
                        }
                        _ => rows,
                    };
                    if node.length() != expected_length || node.null_count() != 0 {
                        return Err(invalid(
                            "evidence IPC array length/null count differs from projected schema",
                        ));
                    }
                }
                // Arrow may copy every unaligned buffer before validating the
                // completed RecordBatch. Canonical disjoint, aligned ranges
                // prevent overlapping slices from multiplying those allocations,
                // and checked bounds prevent Arrow's buffer slicing from panicking.
                let mut previous_end = 0;
                for buffer in buffers {
                    let offset = buffer.offset();
                    let length = buffer.length();
                    let end = offset
                        .checked_add(length)
                        .ok_or_else(|| invalid("evidence IPC buffer range overflow"))?;
                    if offset < 0 || length < 0 || end > body {
                        return Err(invalid("evidence IPC buffer outside declared body"));
                    }
                    if offset % 8 != 0 {
                        return Err(invalid(
                            "canonical evidence IPC buffers require 8-byte alignment",
                        ));
                    }
                    if length != 0 {
                        if offset < previous_end {
                            return Err(invalid(
                                "canonical evidence IPC buffers overlap or are out of order",
                            ));
                        }
                        previous_end = end;
                    }
                }
            }
            _ => {
                return Err(invalid(
                    "canonical evidence IPC only permits schema and record batch messages",
                ))
            }
        }
        self.body_remaining = body as u64;
        self.prefix = std::io::Cursor::new(prefix);
        Ok(())
    }
}
impl<R: std::io::Read> std::io::Read for BoundedIpcReader<R> {
    fn read(&mut self, buffer: &mut [u8]) -> std::io::Result<usize> {
        if buffer.is_empty() {
            return Ok(0);
        }
        loop {
            if self.prefix.position() < self.prefix.get_ref().len() as u64 {
                return self.prefix.read(buffer);
            }
            if self.body_remaining != 0 {
                let allowed = buffer.len().min(self.body_remaining as usize);
                let count = self.source.read(&mut buffer[..allowed])?;
                if count == 0 {
                    return Err(std::io::Error::new(
                        std::io::ErrorKind::UnexpectedEof,
                        "truncated evidence IPC body",
                    ));
                }
                self.body_remaining -= count as u64;
                return Ok(count);
            }
            if self.ended {
                return Ok(0);
            }
            self.next_frame()?;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture(fields: EvidenceFields, rows: u32) -> (Vec<u8>, crate::core::ContigSet) {
        let mut contigs = crate::core::ContigSet::new();
        contigs.push("chr1", 100);
        let batch = EvidenceBatch::new(
            0,
            "chr1",
            0,
            fields,
            (0..rows)
                .map(|position| EvidenceLocus {
                    position,
                    reference: b'A',
                    requested_alts: Vec::new(),
                })
                .collect(),
        );
        let mut writer = EvidenceArrowWriter::with_fields(Vec::new(), fields);
        writer.on_batch(&batch).unwrap();
        (writer.into_inner().unwrap(), contigs)
    }

    fn indirect(bytes: &[u8], offset: usize) -> usize {
        offset + u32::from_le_bytes(bytes[offset..offset + 4].try_into().unwrap()) as usize
    }
    fn table_field(bytes: &[u8], table: usize, slot: u16) -> usize {
        let vtable = (table as isize
            - i32::from_le_bytes(bytes[table..table + 4].try_into().unwrap()) as isize)
            as usize;
        let relative = u16::from_le_bytes(
            bytes[vtable + slot as usize..vtable + slot as usize + 2]
                .try_into()
                .unwrap(),
        ) as usize;
        assert_ne!(relative, 0);
        table + relative
    }

    #[test]
    fn malformed_buffer_ranges_are_refused_from_metadata_before_body_reads() {
        let (bytes, contigs) = fixture(EvidenceFields::ALL, 1);
        let schema_end = 8 + u32::from_le_bytes(bytes[4..8].try_into().unwrap()) as usize;
        let metadata_start = schema_end + 8;
        let metadata_end = metadata_start
            + u32::from_le_bytes(bytes[schema_end + 4..schema_end + 8].try_into().unwrap())
                as usize;
        let message = indirect(&bytes, metadata_start);
        let batch = indirect(
            &bytes,
            table_field(&bytes, message, arrow_ipc::Message::VT_HEADER),
        );
        let buffers = indirect(
            &bytes,
            table_field(&bytes, batch, arrow_ipc::RecordBatch::VT_BUFFERS),
        );
        let nodes = indirect(
            &bytes,
            table_field(&bytes, batch, arrow_ipc::RecordBatch::VT_NODES),
        );
        let body_length = i64::from_le_bytes(
            bytes[table_field(&bytes, message, arrow_ipc::Message::VT_BODYLENGTH)..][..8]
                .try_into()
                .unwrap(),
        );
        let first_nonempty = buffers + 4 + 16; // contig offsets, following its empty validity bitmap
        let next_nonempty = buffers + 4 + 32; // contig UTF-8 bytes
        let first_offset = i64::from_le_bytes(
            bytes[first_nonempty..first_nonempty + 8]
                .try_into()
                .unwrap(),
        );
        let cases = [
            (first_nonempty, -8, "outside declared body"),
            (first_nonempty + 8, -1, "outside declared body"),
            (first_nonempty + 8, body_length + 1, "outside declared body"),
            (first_nonempty, 1, "8-byte alignment"),
            (next_nonempty, first_offset, "overlap or are out of order"),
            (nodes + 4, 2, "array length/null count"),
            (nodes + 12, 1, "array length/null count"),
        ];
        for (offset, value, expected) in cases {
            // No body is available at all: these errors must arise while
            // validating metadata, before Arrow allocates or reads body data.
            let mut malformed = bytes[..metadata_end].to_vec();
            malformed[offset..offset + 8].copy_from_slice(&value.to_le_bytes());
            let error =
                read_evidence_batches(malformed.as_slice(), &contigs, |_| Ok(())).unwrap_err();
            assert!(
                error.to_string().contains(expected),
                "expected {expected}, got {error}"
            );
        }
        for vector in [buffers, nodes] {
            let mut malformed = bytes[..metadata_end].to_vec();
            let count = u32::from_le_bytes(malformed[vector..vector + 4].try_into().unwrap());
            malformed[vector..vector + 4].copy_from_slice(&(count - 1).to_le_bytes());
            let error =
                read_evidence_batches(malformed.as_slice(), &contigs, |_| Ok(())).unwrap_err();
            assert!(error.to_string().contains("node/buffer count"), "{error}");
        }
    }

    #[test]
    fn expected_projection_is_checked_before_any_body_and_on_empty_streams() {
        let (full, contigs) = fixture(EvidenceFields::ALL, 1);
        let header_end = 8 + u32::from_le_bytes(full[4..8].try_into().unwrap()) as usize;
        let error = read_evidence_batches_expected_fields(
            &full[..header_end],
            &contigs,
            EvidenceFields::DEPTHS,
            |_| panic!("unexpected callback"),
        )
        .unwrap_err();
        assert!(
            error.to_string().contains("differs from expected fields"),
            "{error}"
        );
        for bits in 0..=EvidenceFields::ALL_SUPPORTED.bits() {
            let fields = EvidenceFields::from_bits(bits).unwrap();
            for rows in [0, 2] {
                let (bytes, contigs) = fixture(fields, rows);
                let mut observed = 0;
                let metadata = read_evidence_batches_expected_fields(
                    bytes.as_slice(),
                    &contigs,
                    fields,
                    |batch| {
                        assert_eq!(batch.fields(), fields);
                        observed += batch.len();
                        Ok(())
                    },
                )
                .unwrap();
                assert_eq!(observed, rows as usize);
                assert_eq!(
                    metadata,
                    EvidenceArtifactMetadata {
                        fields,
                        schema_version: fields.schema_version()
                    }
                );
            }
        }
        let (empty, contigs) = fixture(EvidenceFields::ALL, 0);
        assert!(read_evidence_batches_expected_fields(
            empty.as_slice(),
            &contigs,
            EvidenceFields::DEPTHS,
            |_| Ok(())
        )
        .is_err());
    }

    #[test]
    fn allele_sums_validate_against_each_independently_available_count_source() {
        for (bits, case) in [(64, 0), (64, 1), (65, 0), (66, 0), (68, 0), (80, 0)] {
            let fields = EvidenceFields::from_bits(bits).unwrap();
            let mut batch = EvidenceBatch::new(
                0,
                "chr1",
                0,
                fields,
                vec![EvidenceLocus {
                    position: 0,
                    reference: b'A',
                    requested_alts: vec![b'C'],
                }],
            );
            let row = batch.row_mut(0).unwrap();
            let sums = row.allele_quality.unwrap();
            match bits {
                64 if case == 0 => sums.base_quality_sum[0] = 1,
                64 => {
                    sums.read_position_sum[0] = 1;
                    sums.read_length_sum[0] = 1;
                }
                // No callable observations can contribute quality or length.
                65 => sums.base_quality_sum[0] = 1,
                66 => sums.read_length_sum[0] = 1,
                // A strand count independently proves no C observations exist.
                68 => {
                    row.strands.unwrap().strand_counts[0][0] = 1;
                    sums.read_length_sum[0] = 1;
                    sums.base_quality_sum[1] = 1;
                }
                // Pooled histograms constrain sums even without QUALITY_SUMS.
                80 => {
                    let hist = row.quality_histograms.unwrap();
                    hist.base_quality_histogram[30] = 1;
                    hist.mapping_quality_histogram[60] = 1;
                    sums.base_quality_sum[0] = 31;
                    sums.mapping_quality_sum[0] = 60;
                    sums.read_length_sum[0] = 1;
                }
                _ => unreachable!(),
            }
            let mut writer = EvidenceArrowWriter::with_fields(Vec::new(), fields);
            writer.on_batch(&batch).unwrap();
            let bytes = writer.into_inner().unwrap();
            let mut contigs = crate::core::ContigSet::new();
            contigs.push("chr1", 1);
            assert!(
                read_evidence_batches_expected_fields(bytes.as_slice(), &contigs, fields, |_| {
                    panic!("invalid values must fail before consumer access")
                })
                .is_err(),
                "mask={bits}"
            );
        }
    }

    #[test]
    fn final_buffered_flush_error_is_returned_after_arrow_end_of_stream() {
        #[derive(Default)]
        struct FailingFinalFlush(Vec<u8>);
        impl Write for FailingFinalFlush {
            fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
                self.0.extend_from_slice(bytes);
                Ok(bytes.len())
            }
            fn flush(&mut self) -> std::io::Result<()> {
                if self.0.ends_with(&[255, 255, 255, 255, 0, 0, 0, 0]) {
                    Err(std::io::Error::other("final EOS flush failed"))
                } else {
                    Ok(())
                }
            }
        }
        let mut writer =
            EvidenceArrowWriter::new(std::io::BufWriter::new(FailingFinalFlush::default()));
        // Arrow's continuation helper flushes EOS through the BufWriter. A
        // final I/O failure must escape before callers publish a complete file.
        let error = writer.finish().unwrap_err();
        assert!(error.to_string().contains("final EOS flush failed"));
    }
}
