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
fn scalars(row: &EvidenceRow) -> [u64; 28] {
    [
        row.prefilter_depth,
        row.aligned_depth,
        row.callable_depth,
        row.allele_counts[0],
        row.allele_counts[1],
        row.allele_counts[2],
        row.allele_counts[3],
        row.strand_counts[0][0],
        row.strand_counts[0][1],
        row.strand_counts[1][0],
        row.strand_counts[1][1],
        row.strand_counts[2][0],
        row.strand_counts[2][1],
        row.strand_counts[3][0],
        row.strand_counts[3][1],
        row.base_quality_sum,
        row.mapping_quality_sum,
        row.read_position_sum,
        row.read_length_sum,
        row.filters.secondary,
        row.filters.supplementary,
        row.filters.qc_fail,
        row.filters.duplicate,
        row.filters.unavailable_mapq,
        row.filters.low_mapq,
        row.filters.unavailable_base_quality,
        row.filters.low_base_quality,
        row.filters.ambiguous_base,
    ]
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
}
impl<W: Write> EvidenceTsvWriter<W> {
    /// Construct the configured value without starting evidence extraction.
    pub fn new(out: W) -> Self {
        Self {
            out,
            started: false,
        }
    }
    /// Finish encoding when necessary and return the underlying output writer.
    pub fn into_inner(self) -> W {
        self.out
    }
    fn header(&mut self) -> Result<(), EvidenceError> {
        if !self.started {
            write!(self.out, "#contig\tpos\tref\trequested_alts")?;
            for name in SCALAR_NAMES {
                write!(self.out, "\t{name}")?;
            }
            writeln!(
                self.out,
                "\tbase_quality_histogram\tmapping_quality_histogram"
            )?;
            self.started = true;
        }
        Ok(())
    }
}
impl<W: Write> EvidenceAnalyzer for EvidenceTsvWriter<W> {
    fn on_batch(&mut self, batch: &EvidenceBatch) -> Result<(), EvidenceError> {
        self.header()?;
        for row in &batch.rows {
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
            for value in scalars(row) {
                write!(self.out, "\t{value}")?;
            }
            write!(self.out, "\t")?;
            write_histogram(&mut self.out, &row.base_quality_histogram)?;
            write!(self.out, "\t")?;
            write_histogram(&mut self.out, &row.mapping_quality_histogram)?;
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

fn schema() -> Schema {
    let mut fields = vec![
        Field::new("contig", DataType::Utf8, false),
        Field::new("pos", DataType::UInt32, false),
        Field::new("ref", DataType::Utf8, false),
        Field::new("requested_alts", DataType::Utf8, false),
    ];
    fields.extend(
        SCALAR_NAMES
            .iter()
            .map(|name| Field::new(*name, DataType::UInt64, false)),
    );
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
    Schema::new(fields).with_metadata(std::collections::HashMap::from([(
        "rosalind.evidence.schema".into(),
        EVIDENCE_SCHEMA_VERSION.to_string(),
    )]))
}

/// Canonical Arrow IPC writer. Buffering is limited to 1024 evidence rows plus
/// encoding scratch. Construction performs no writes, so admission can run first.
pub struct EvidenceArrowWriter<W: Write> {
    output: Option<W>,
    writer: Option<StreamWriter<W>>,
    rows: Vec<(Arc<str>, EvidenceRow)>,
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
        Self {
            output: Some(output),
            writer: None,
            rows: Vec::new(),
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
                StreamWriter::try_new_with_options(output, &schema(), options)
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
        let mut arrays: Vec<ArrayRef> = Vec::with_capacity(34);
        arrays.push(Arc::new(StringArray::from_iter_values(
            self.rows.iter().map(|(contig, _)| contig.as_ref()),
        )));
        arrays.push(Arc::new(UInt32Array::from_iter_values(
            self.rows.iter().map(|(_, row)| row.position + 1),
        )));
        arrays.push(Arc::new(StringArray::from_iter_values(
            self.rows
                .iter()
                .map(|(_, row)| char::from(row.reference).to_string()),
        )));
        arrays.push(Arc::new(StringArray::from_iter_values(
            self.rows.iter().map(|(_, row)| {
                String::from_utf8(row.requested_alts.clone()).expect("validated nucleotide ALT")
            }),
        )));
        for field in 0..SCALAR_NAMES.len() {
            arrays.push(Arc::new(UInt64Array::from_iter_values(
                self.rows.iter().map(|(_, row)| scalars(row)[field]),
            )));
        }
        let bq: ArrayRef = Arc::new(UInt64Array::from_iter_values(
            self.rows
                .iter()
                .flat_map(|(_, row)| row.base_quality_histogram.iter().copied()),
        ));
        let mq: ArrayRef =
            Arc::new(UInt64Array::from_iter_values(self.rows.iter().flat_map(
                |(_, row)| row.mapping_quality_histogram.iter().copied(),
            )));
        arrays.push(Arc::new(
            FixedSizeListArray::try_new(
                Arc::new(Field::new("item", DataType::UInt64, false)),
                BASE_QUALITY_BINS as i32,
                bq,
                None,
            )
            .map_err(arrow_error)?,
        ));
        arrays.push(Arc::new(
            FixedSizeListArray::try_new(
                Arc::new(Field::new("item", DataType::UInt64, false)),
                MAPPING_QUALITY_BINS as i32,
                mq,
                None,
            )
            .map_err(arrow_error)?,
        ));
        let batch = RecordBatch::try_new(Arc::new(schema()), arrays).map_err(arrow_error)?;
        self.writer
            .as_mut()
            .unwrap()
            .write(&batch)
            .map_err(arrow_error)?;
        self.rows.clear();
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
    fn on_batch(&mut self, batch: &EvidenceBatch) -> Result<(), EvidenceError> {
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
        for row in &batch.rows {
            if self.rows.is_empty() && self.rows.capacity() == 0 {
                self.rows.reserve_exact(EVIDENCE_ARROW_BATCH_ROWS);
            }
            self.rows.push((Arc::clone(&contig), row.clone()));
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
            (std::mem::size_of::<EvidenceRow>() as u64) * EVIDENCE_ARROW_BATCH_ROWS as u64 * 4
                + (8 << 20),
        )
    }
}
fn arrow_error(error: impl std::fmt::Display) -> EvidenceError {
    EvidenceError::Analyzer(error.to_string())
}

/// Stream a canonical evidence Arrow artifact into bounded callbacks. Cached
/// input is schema-checked, rejects nulls/oversized batches, and is required to
/// be uniquely sorted. `contigs` supplies stable IDs from the verified reference.
pub fn read_evidence_batches<R: std::io::Read>(
    reader: R,
    contigs: &crate::core::ContigSet,
    mut on_batch: impl FnMut(&EvidenceBatch) -> Result<(), EvidenceError>,
) -> Result<(), EvidenceError> {
    use arrow_array::Array;
    let mut reader = arrow_ipc::reader::StreamReader::try_new(BoundedIpcReader::new(reader), None)
        .map_err(arrow_error)?;
    if reader.schema().as_ref() != &schema() {
        return Err(EvidenceError::InvalidInput(
            "cached Arrow evidence schema mismatch".into(),
        ));
    }
    let mut previous = None;
    for encoded in &mut reader {
        let encoded = encoded.map_err(arrow_error)?;
        if encoded.num_rows() > EVIDENCE_ARROW_BATCH_ROWS {
            return Err(EvidenceError::InvalidInput(
                "cached Arrow evidence batch exceeds canonical 1024 rows".into(),
            ));
        }
        if encoded
            .columns()
            .iter()
            .any(|array| array.null_count() != 0)
        {
            return Err(EvidenceError::InvalidInput(
                "cached evidence cannot contain nulls".into(),
            ));
        }
        let names = encoded
            .column(0)
            .as_any()
            .downcast_ref::<StringArray>()
            .ok_or_else(|| EvidenceError::InvalidInput("invalid contig array".into()))?;
        let positions = encoded
            .column(1)
            .as_any()
            .downcast_ref::<UInt32Array>()
            .ok_or_else(|| EvidenceError::InvalidInput("invalid position array".into()))?;
        let references = encoded
            .column(2)
            .as_any()
            .downcast_ref::<StringArray>()
            .ok_or_else(|| EvidenceError::InvalidInput("invalid reference array".into()))?;
        let alts = encoded
            .column(3)
            .as_any()
            .downcast_ref::<StringArray>()
            .ok_or_else(|| EvidenceError::InvalidInput("invalid ALT array".into()))?;
        let mut numeric = Vec::new();
        for column in 4..32 {
            numeric.push(
                encoded
                    .column(column)
                    .as_any()
                    .downcast_ref::<UInt64Array>()
                    .ok_or_else(|| {
                        EvidenceError::InvalidInput("invalid scalar evidence array".into())
                    })?,
            );
        }
        let bq = encoded
            .column(32)
            .as_any()
            .downcast_ref::<FixedSizeListArray>()
            .ok_or_else(|| EvidenceError::InvalidInput("invalid BQ histogram".into()))?;
        let mq = encoded
            .column(33)
            .as_any()
            .downcast_ref::<FixedSizeListArray>()
            .ok_or_else(|| EvidenceError::InvalidInput("invalid MAPQ histogram".into()))?;
        let mut batch: Option<EvidenceBatch> = None;
        for index in 0..encoded.num_rows() {
            if names.value(index).len() > 4096 || alts.value(index).len() > 3 {
                return Err(EvidenceError::InvalidInput(
                    "cached evidence string exceeds canonical bounds".into(),
                ));
            }
            let contig = contigs.by_name(names.value(index)).ok_or_else(|| {
                EvidenceError::InvalidInput("cached evidence contig not in dictionary".into())
            })?;
            let position = positions.value(index).checked_sub(1).ok_or_else(|| {
                EvidenceError::InvalidInput("evidence POS must be 1-based".into())
            })?;
            if position >= contig.length {
                return Err(EvidenceError::InvalidInput(
                    "cached evidence position out of bounds".into(),
                ));
            }
            let locus = (contig.id, position);
            if previous.is_some_and(|previous| locus <= previous) {
                return Err(EvidenceError::InvalidInput(
                    "cached evidence rows are not uniquely coordinate sorted".into(),
                ));
            }
            previous = Some(locus);
            let ref_text = references.value(index);
            if ref_text.len() != 1 || !b"ACGTN".contains(&ref_text.as_bytes()[0]) {
                return Err(EvidenceError::InvalidInput(
                    "invalid cached reference base".into(),
                ));
            }
            let v = |column: usize| numeric[column].value(index);
            let mut row = EvidenceRow {
                position,
                reference: ref_text.as_bytes()[0],
                prefilter_depth: v(0),
                aligned_depth: v(1),
                callable_depth: v(2),
                allele_counts: [v(3), v(4), v(5), v(6)],
                strand_counts: [[v(7), v(8)], [v(9), v(10)], [v(11), v(12)], [v(13), v(14)]],
                base_quality_sum: v(15),
                mapping_quality_sum: v(16),
                read_position_sum: v(17),
                read_length_sum: v(18),
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
                requested_alts: alts.value(index).as_bytes().to_vec(),
                ..EvidenceRow::default()
            };
            if row
                .requested_alts
                .iter()
                .any(|base| !b"ACGT".contains(base))
            {
                return Err(EvidenceError::InvalidInput("invalid cached SNV ALT".into()));
            }
            let bq_value = bq.value(index);
            let mq_value = mq.value(index);
            let bq_values = bq_value
                .as_any()
                .downcast_ref::<UInt64Array>()
                .ok_or_else(|| EvidenceError::InvalidInput("invalid BQ histogram values".into()))?;
            let mq_values = mq_value
                .as_any()
                .downcast_ref::<UInt64Array>()
                .ok_or_else(|| {
                    EvidenceError::InvalidInput("invalid MAPQ histogram values".into())
                })?;
            if bq_values.null_count() != 0 || mq_values.null_count() != 0 {
                return Err(EvidenceError::InvalidInput(
                    "histograms cannot contain nulls".into(),
                ));
            }
            row.base_quality_histogram
                .copy_from_slice(bq_values.values());
            row.mapping_quality_histogram
                .copy_from_slice(mq_values.values());
            validate_row(&row)?;
            let canonical_start = position / CANONICAL_TILE_BASES * CANONICAL_TILE_BASES;
            if batch.as_ref().is_some_and(|batch| {
                batch.contig_id != contig.id || batch.canonical_tile_start != canonical_start
            }) {
                on_batch(batch.as_ref().unwrap())?;
                batch = None;
            }
            batch
                .get_or_insert_with(|| EvidenceBatch {
                    contig_id: contig.id,
                    contig: contig.name.to_string(),
                    canonical_tile_start: canonical_start,
                    rows: Vec::new(),
                })
                .rows
                .push(row);
        }
        if let Some(batch) = batch {
            on_batch(&batch)?;
        }
    }
    reader.get_mut().validate_complete()?;
    Ok(())
}

fn validate_row(row: &EvidenceRow) -> Result<(), EvidenceError> {
    let sum = |values: &[u64]| values.iter().map(|value| u128::from(*value)).sum::<u128>();
    let weighted = |values: &[u64]| {
        values
            .iter()
            .enumerate()
            .map(|(quality, count)| quality as u128 * u128::from(*count))
            .sum::<u128>()
    };
    let f = &row.filters;
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
    if row
        .requested_alts
        .windows(2)
        .any(|bases| bases[0] >= bases[1])
        || row.requested_alts.contains(&row.reference)
        || (!row.requested_alts.is_empty() && row.reference == b'N')
        || u128::from(row.prefilter_depth) != u128::from(row.aligned_depth) + read_filtered
        || u128::from(row.aligned_depth) != u128::from(row.callable_depth) + base_filtered
        || sum(&row.allele_counts) != u128::from(row.callable_depth)
        || row
            .strand_counts
            .iter()
            .zip(row.allele_counts)
            .any(|(strands, allele)| sum(strands) != u128::from(allele))
        || sum(&row.base_quality_histogram) != u128::from(row.callable_depth)
        || sum(&row.mapping_quality_histogram) != u128::from(row.callable_depth)
        || weighted(&row.base_quality_histogram) != u128::from(row.base_quality_sum)
        || weighted(&row.mapping_quality_histogram) != u128::from(row.mapping_quality_sum)
        || u128::from(row.read_position_sum) + u128::from(row.callable_depth)
            > u128::from(row.read_length_sum)
    {
        return Err(EvidenceError::InvalidInput(
            "cached evidence counters or ALT annotations are inconsistent".into(),
        ));
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
    ended: bool,
}
impl<R: std::io::Read> BoundedIpcReader<R> {
    fn new(source: R) -> Self {
        Self {
            source,
            prefix: std::io::Cursor::new(Vec::new()),
            body_remaining: 0,
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
        if body < 0 || body as u64 > MAX_IPC_BODY_BYTES {
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
                if batch.nodes().is_some_and(|nodes| {
                    nodes.iter().any(|node| {
                        node.length() < 0
                            || node.length() as usize
                                > EVIDENCE_ARROW_BATCH_ROWS * MAPPING_QUALITY_BINS
                    })
                }) {
                    return Err(invalid(
                        "evidence IPC array exceeds the bounded row envelope",
                    ));
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
