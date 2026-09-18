//! Canonical bounded cohort output. Unmeasured cells use actual Arrow nulls;
//! TSV uses a dot for nulls and preserves integer counts without float conversion.

use super::runtime::{CohortConsumer, CohortRow};
use super::summary::{required_fields, CandidateObservation, CandidateReducer, CandidateSummary};
use super::{CohortError, Result};
use crate::core::ContigSet;
use crate::evidence::{
    evidence_record_batch, EvidenceBatch, EvidenceFields, SnvSite, CANONICAL_TILE_BASES,
};
use arrow_array::{
    Array, ArrayRef, BooleanArray, FixedSizeListArray, RecordBatch, StringArray, UInt32Array,
    UInt64Array,
};
use arrow_ipc::{writer::IpcWriteOptions, writer::StreamWriter, MetadataVersion};
use arrow_schema::{DataType, Field, Schema};
use std::io::Write;
use std::sync::Arc;

pub(crate) const OUTPUT_BATCH_ROWS: usize = 1024;
const MAX_TEXT_BYTES: usize = 4096;
const MAX_MEMBER_ID_BYTES: usize = 256;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum CohortOutputFormat {
    Arrow,
    Tsv,
}

/// Includes retained native rows, nullable Arrow projection, IPC staging, and keys.
/// Callers must also reserve their writer (e.g. file buffer) and managed receipt.
pub(crate) fn memory_bytes(fields: EvidenceFields, max_text_length: usize) -> u64 {
    (OUTPUT_BATCH_ROWS as u64)
        .saturating_mul(
            fields
                .storage_bytes_per_locus()
                .saturating_mul(6)
                .saturating_add((max_text_length as u64).saturating_mul(8))
                .saturating_add(4096),
        )
        .saturating_add(1 << 20)
}

/// Reducers are bounded by ALT-expanded candidates in one genomic window.
pub(crate) fn summary_memory_bytes(max_window_candidates: usize, max_text_length: usize) -> u64 {
    (max_window_candidates as u64)
        .saturating_mul((std::mem::size_of::<WindowCandidate>() as u64).saturating_mul(2))
        .saturating_add(
            (OUTPUT_BATCH_ROWS as u64).saturating_mul(
                (max_text_length as u64)
                    .saturating_mul(6)
                    .saturating_add(8192),
            ),
        )
        .saturating_add(1 << 20)
}

pub(crate) fn max_text_length(contigs: &ContigSet) -> usize {
    contigs
        .iter()
        .map(|contig| contig.name.len())
        .max()
        .unwrap_or(0)
        .max(MAX_MEMBER_ID_BYTES)
}

fn validate_text(value: &str, format: CohortOutputFormat, maximum: usize) -> Result<()> {
    if value.len() > maximum {
        return Err(CohortError::Limit(format!(
            "cohort output text exceeds {maximum} bytes"
        )));
    }
    if format == CohortOutputFormat::Tsv
        && value
            .bytes()
            .any(|byte| matches!(byte, b'\t' | b'\n' | b'\r'))
    {
        return Err(CohortError::Incompatible(
            "TSV identifiers cannot contain tab, newline or carriage return".into(),
        ));
    }
    Ok(())
}

fn validate_contigs(contigs: &ContigSet, format: CohortOutputFormat) -> Result<()> {
    for contig in contigs.iter() {
        validate_text(&contig.name, format, MAX_TEXT_BYTES)?;
    }
    Ok(())
}

fn arrow_error(error: arrow_schema::ArrowError) -> CohortError {
    CohortError::Corrupt(format!("cohort Arrow encoding: {error}"))
}

struct BatchSink<W: Write> {
    output: Option<W>,
    arrow: Option<StreamWriter<W>>,
    format: CohortOutputFormat,
    started: bool,
    finished: bool,
}
impl<W: Write> BatchSink<W> {
    fn new(output: W, format: CohortOutputFormat) -> Self {
        Self {
            output: Some(output),
            arrow: None,
            format,
            started: false,
            finished: false,
        }
    }
    fn start(&mut self, schema: &Schema) -> Result<()> {
        if self.started {
            return Ok(());
        }
        match self.format {
            CohortOutputFormat::Arrow => {
                let output = self.output.take().ok_or_else(|| {
                    CohortError::Corrupt("cohort encoder output unavailable".into())
                })?;
                let options =
                    IpcWriteOptions::try_new(8, false, MetadataVersion::V5).map_err(arrow_error)?;
                self.arrow = Some(
                    StreamWriter::try_new_with_options(output, schema, options)
                        .map_err(arrow_error)?,
                );
            }
            CohortOutputFormat::Tsv => {
                let output = self.output.as_mut().expect("unstarted output");
                for (index, field) in schema.fields().iter().enumerate() {
                    if index != 0 {
                        write!(output, "\t")?;
                    }
                    write!(output, "{}", field.name())?;
                }
                writeln!(output)?;
            }
        }
        self.started = true;
        Ok(())
    }
    fn write(&mut self, batch: &RecordBatch) -> Result<()> {
        if self.finished {
            return Err(CohortError::Corrupt(
                "cannot append to a finished cohort output".into(),
            ));
        }
        self.start(batch.schema().as_ref())?;
        if batch.num_rows() == 0 {
            return Ok(());
        }
        match self.format {
            CohortOutputFormat::Arrow => self
                .arrow
                .as_mut()
                .expect("started Arrow output")
                .write(batch)
                .map_err(arrow_error)?,
            CohortOutputFormat::Tsv => {
                let out = self.output.as_mut().expect("TSV output");
                for row in 0..batch.num_rows() {
                    for (column, array) in batch.columns().iter().enumerate() {
                        if column != 0 {
                            write!(out, "\t")?;
                        }
                        write_tsv_value(out, array.as_ref(), row)?;
                    }
                    writeln!(out)?;
                }
            }
        }
        Ok(())
    }
    fn finish(&mut self, schema: &Schema) -> Result<()> {
        if !self.finished {
            self.start(schema)?;
            match self.format {
                CohortOutputFormat::Arrow => self
                    .arrow
                    .as_mut()
                    .expect("started Arrow output")
                    .finish()
                    .map_err(arrow_error)?,
                CohortOutputFormat::Tsv => self.output.as_mut().unwrap().flush()?,
            }
            self.finished = true;
        }
        Ok(())
    }
    fn into_inner(mut self) -> Result<W> {
        if !self.finished {
            return Err(CohortError::Corrupt("cohort output is unfinished".into()));
        }
        match self.format {
            CohortOutputFormat::Arrow => {
                self.arrow.take().unwrap().into_inner().map_err(arrow_error)
            }
            CohortOutputFormat::Tsv => Ok(self.output.take().unwrap()),
        }
    }
}

fn write_tsv_value(out: &mut dyn Write, array: &dyn Array, row: usize) -> Result<()> {
    if array.is_null(row) {
        write!(out, ".")?;
        return Ok(());
    }
    match array.data_type() {
        DataType::Utf8 => {
            let value = array
                .as_any()
                .downcast_ref::<StringArray>()
                .unwrap()
                .value(row);
            validate_text(value, CohortOutputFormat::Tsv, MAX_TEXT_BYTES)?;
            write!(out, "{value}")?;
        }
        DataType::UInt32 => write!(
            out,
            "{}",
            array
                .as_any()
                .downcast_ref::<UInt32Array>()
                .unwrap()
                .value(row)
        )?,
        DataType::UInt64 => write!(
            out,
            "{}",
            array
                .as_any()
                .downcast_ref::<UInt64Array>()
                .unwrap()
                .value(row)
        )?,
        DataType::Boolean => write!(
            out,
            "{}",
            array
                .as_any()
                .downcast_ref::<BooleanArray>()
                .unwrap()
                .value(row)
        )?,
        DataType::FixedSizeList(_, _) => {
            let values = array
                .as_any()
                .downcast_ref::<FixedSizeListArray>()
                .unwrap()
                .value(row);
            for index in 0..values.len() {
                if index != 0 {
                    write!(out, ",")?;
                }
                write_tsv_value(out, values.as_ref(), index)?;
            }
        }
        other => {
            return Err(CohortError::Corrupt(format!(
                "unsupported cohort TSV type {other:?}"
            )))
        }
    }
    Ok(())
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct CandidateKey {
    contig: u32,
    position: u32,
    reference: u8,
    alternate: u8,
}
impl CandidateKey {
    fn from_row(row: CohortRow<'_>) -> Self {
        Self {
            contig: row.contig,
            position: row.position,
            reference: row.reference,
            alternate: row.alternate,
        }
    }
    fn validate(self, contigs: &ContigSet) -> Result<()> {
        let contig = contigs
            .by_id(self.contig)
            .ok_or_else(|| CohortError::Corrupt("cohort output refers to unknown contig".into()))?;
        if self.position as u64 >= contig.length as u64
            || !b"ACGT".contains(&self.reference)
            || !b"ACGT".contains(&self.alternate)
            || self.reference == self.alternate
        {
            return Err(CohortError::Corrupt("invalid output candidate SNV".into()));
        }
        Ok(())
    }
}

fn candidate_fields() -> Vec<Field> {
    vec![
        Field::new("contig", DataType::Utf8, false),
        Field::new("pos", DataType::UInt32, false),
        Field::new("ref", DataType::Utf8, false),
        Field::new("alt", DataType::Utf8, false),
    ]
}
fn candidate_arrays(keys: &[CandidateKey], contigs: &ContigSet) -> Vec<ArrayRef> {
    vec![
        Arc::new(StringArray::from_iter_values(keys.iter().map(|key| {
            contigs
                .by_id(key.contig)
                .expect("validated candidate contig")
                .name
                .as_ref()
        }))),
        Arc::new(UInt32Array::from_iter_values(
            keys.iter().map(|key| key.position + 1),
        )),
        Arc::new(StringArray::from_iter_values(
            keys.iter().map(|key| char::from(key.reference).to_string()),
        )),
        Arc::new(StringArray::from_iter_values(
            keys.iter().map(|key| char::from(key.alternate).to_string()),
        )),
    ]
}
fn output_schema(fields: Vec<Field>, kind: &str) -> Schema {
    Schema::new(fields).with_metadata(std::collections::HashMap::from([(
        "rosalind.cohort.schema".into(),
        format!("1;{kind};counting-unit=read-observations"),
    )]))
}

struct ExtractCell {
    key: CandidateKey,
    member_id: String,
    observation: CandidateObservation,
    evidence_index: Option<u32>,
}

pub(crate) struct ExtractEncoder<'a, W: Write> {
    sink: BatchSink<W>,
    contigs: &'a ContigSet,
    fields: EvidenceFields,
    min_callable_depth: u64,
    cells: Vec<ExtractCell>,
    observed: EvidenceBatch,
}
impl<'a, W: Write> ExtractEncoder<'a, W> {
    pub fn new(
        output: W,
        format: CohortOutputFormat,
        contigs: &'a ContigSet,
        fields: EvidenceFields,
        min_callable_depth: u64,
    ) -> Result<Self> {
        validate_contigs(contigs, format)?;
        CandidateObservation::unmeasured(min_callable_depth)?;
        if !fields.contains(required_fields()) {
            return Err(CohortError::Incompatible(
                "candidate extraction requires depths and allele counts".into(),
            ));
        }
        Ok(Self {
            sink: BatchSink::new(output, format),
            contigs,
            fields,
            min_callable_depth,
            cells: Vec::new(),
            observed: EvidenceBatch::new(0, "", 0, fields, Vec::new()),
        })
    }
    fn record_batch(&self) -> Result<RecordBatch> {
        let evidence = evidence_record_batch(&self.observed)?;
        let indices = UInt32Array::from_iter(self.cells.iter().map(|cell| cell.evidence_index));
        let mut fields = vec![Field::new("sample_id", DataType::Utf8, false)];
        fields.extend(candidate_fields());
        fields.push(Field::new("status", DataType::Utf8, false));
        let keys: Vec<_> = self.cells.iter().map(|cell| cell.key).collect();
        let mut arrays: Vec<ArrayRef> = vec![Arc::new(StringArray::from_iter_values(
            self.cells.iter().map(|cell| cell.member_id.as_str()),
        ))];
        arrays.extend(candidate_arrays(&keys, self.contigs));
        arrays.push(Arc::new(StringArray::from_iter_values(
            self.cells.iter().map(|cell| {
                if cell.evidence_index.is_some() {
                    "observed"
                } else {
                    "unmeasured"
                }
            }),
        )));
        for (field, array) in evidence
            .schema()
            .fields()
            .iter()
            .zip(evidence.columns())
            .skip(4)
        {
            fields.push(field.as_ref().clone().with_nullable(true));
            arrays.push(
                arrow_select::take::take(array.as_ref(), &indices, None).map_err(arrow_error)?,
            );
        }
        for name in [
            "alt_count",
            "observed_alt_fraction_numerator",
            "observed_alt_fraction_denominator",
        ] {
            fields.push(Field::new(name, DataType::UInt64, true));
        }
        arrays.push(Arc::new(UInt64Array::from_iter(
            self.cells.iter().map(|cell| cell.observation.alt_count),
        )));
        arrays.push(Arc::new(UInt64Array::from_iter(self.cells.iter().map(
            |cell| {
                cell.observation
                    .observed_alt_fraction
                    .map(|fraction| fraction.numerator)
            },
        ))));
        arrays.push(Arc::new(UInt64Array::from_iter(self.cells.iter().map(
            |cell| {
                cell.observation
                    .observed_alt_fraction
                    .map(|fraction| fraction.denominator)
            },
        ))));
        for name in ["depth_eligible", "alt_supported"] {
            fields.push(Field::new(name, DataType::Boolean, true));
        }
        arrays.push(Arc::new(BooleanArray::from_iter(
            self.cells
                .iter()
                .map(|cell| cell.observation.depth_eligible),
        )));
        arrays.push(Arc::new(BooleanArray::from_iter(
            self.cells.iter().map(|cell| cell.observation.alt_supported),
        )));
        fields.push(Field::new("min_callable_depth", DataType::UInt64, false));
        arrays.push(Arc::new(UInt64Array::from_iter_values(
            self.cells.iter().map(|_| self.min_callable_depth),
        )));
        RecordBatch::try_new(Arc::new(output_schema(fields, "candidate-extract")), arrays)
            .map_err(arrow_error)
    }
    fn flush_batch(&mut self) -> Result<()> {
        let batch = self.record_batch()?;
        self.sink.write(&batch)?;
        self.cells.clear();
        self.observed.clear();
        Ok(())
    }
    pub fn into_inner(mut self) -> Result<W> {
        self.finish()?;
        self.sink.into_inner()
    }
}
impl<W: Write> CohortConsumer for ExtractEncoder<'_, W> {
    fn retained_bytes(&self) -> Option<u64> {
        Some(memory_bytes(self.fields, max_text_length(self.contigs)))
    }
    fn on_row(&mut self, row: CohortRow<'_>) -> Result<()> {
        if self.sink.finished {
            return Err(CohortError::Corrupt(
                "cannot append to finished cohort output".into(),
            ));
        }
        let key = CandidateKey::from_row(row);
        key.validate(self.contigs)?;
        validate_text(row.member_id, self.sink.format, MAX_MEMBER_ID_BYTES)?;
        let (observation, evidence_index) = match row.evidence {
            Some(evidence) => {
                if evidence.position != row.position || evidence.reference != row.reference {
                    return Err(CohortError::Corrupt(
                        "candidate key differs from evidence row".into(),
                    ));
                }
                let observation = CandidateObservation::from_row(
                    evidence,
                    row.alternate,
                    self.min_callable_depth,
                )?;
                let index = self.observed.len() as u32;
                self.observed.push_row(evidence)?;
                (observation, Some(index))
            }
            None => (
                CandidateObservation::unmeasured(self.min_callable_depth)?,
                None,
            ),
        };
        self.cells.push(ExtractCell {
            key,
            member_id: row.member_id.to_owned(),
            observation,
            evidence_index,
        });
        if self.cells.len() == OUTPUT_BATCH_ROWS {
            self.flush_batch()?;
        }
        Ok(())
    }
    fn finish(&mut self) -> Result<()> {
        if !self.sink.finished {
            self.flush_batch()?;
            self.sink.finish(self.record_batch()?.schema().as_ref())?;
        }
        Ok(())
    }
}

struct WindowCandidate {
    key: CandidateKey,
    reducer: CandidateReducer,
    last_member: Option<usize>,
}
struct SummaryCell {
    key: CandidateKey,
    summary: CandidateSummary,
}

pub(crate) struct SummaryEncoder<'a, W: Write> {
    sink: BatchSink<W>,
    contigs: &'a ContigSet,
    n_requested: u64,
    min_callable_depth: u64,
    max_window_candidates: usize,
    window: Vec<WindowCandidate>,
    window_open: bool,
    cells: Vec<SummaryCell>,
}
impl<'a, W: Write> SummaryEncoder<'a, W> {
    pub fn new(
        output: W,
        format: CohortOutputFormat,
        contigs: &'a ContigSet,
        n_requested: u64,
        min_callable_depth: u64,
        max_window_candidates: usize,
    ) -> Result<Self> {
        validate_contigs(contigs, format)?;
        CandidateReducer::new(n_requested, min_callable_depth)?;
        if max_window_candidates > 3 * CANONICAL_TILE_BASES as usize {
            return Err(CohortError::Limit(
                "summary window exceeds canonical tile candidates".into(),
            ));
        }
        Ok(Self {
            sink: BatchSink::new(output, format),
            contigs,
            n_requested,
            min_callable_depth,
            max_window_candidates,
            window: Vec::new(),
            window_open: false,
            cells: Vec::new(),
        })
    }
    fn record_batch(&self) -> Result<RecordBatch> {
        let mut fields = candidate_fields();
        let keys: Vec<_> = self.cells.iter().map(|cell| cell.key).collect();
        let mut arrays = candidate_arrays(&keys, self.contigs);
        macro_rules! count {
            ($name:ident) => {
                fields.push(Field::new(stringify!($name), DataType::UInt64, false));
                arrays.push(Arc::new(UInt64Array::from_iter_values(
                    self.cells.iter().map(|cell| cell.summary.$name),
                )) as ArrayRef);
            };
        }
        count!(min_callable_depth);
        count!(n_requested);
        count!(n_observed);
        count!(n_depth_eligible);
        count!(n_alt_supported);
        count!(callable_total);
        count!(alt_total);
        count!(eligible_callable_total);
        count!(eligible_alt_total);
        macro_rules! fraction {
            ($name:ident) => {
                for name in [
                    concat!(stringify!($name), "_numerator"),
                    concat!(stringify!($name), "_denominator"),
                ] {
                    fields.push(Field::new(name, DataType::UInt64, true));
                }
                arrays.push(Arc::new(UInt64Array::from_iter(
                    self.cells
                        .iter()
                        .map(|cell| cell.summary.$name.map(|fraction| fraction.numerator)),
                )) as ArrayRef);
                arrays.push(Arc::new(UInt64Array::from_iter(
                    self.cells
                        .iter()
                        .map(|cell| cell.summary.$name.map(|fraction| fraction.denominator)),
                )) as ArrayRef);
            };
        }
        fraction!(depth_eligible_support_fraction);
        fraction!(observed_alt_fraction);
        fraction!(eligible_alt_fraction);
        RecordBatch::try_new(Arc::new(output_schema(fields, "candidate-summary")), arrays)
            .map_err(arrow_error)
    }
    fn flush_batch(&mut self) -> Result<()> {
        self.sink.write(&self.record_batch()?)?;
        self.cells.clear();
        Ok(())
    }
    pub fn into_inner(mut self) -> Result<W> {
        self.finish()?;
        self.sink.into_inner()
    }
}
impl<W: Write> CohortConsumer for SummaryEncoder<'_, W> {
    fn retained_bytes(&self) -> Option<u64> {
        Some(summary_memory_bytes(
            self.max_window_candidates,
            max_text_length(self.contigs),
        ))
    }
    fn begin_window(&mut self, sites: &[SnvSite]) -> Result<()> {
        if self.window_open || self.sink.finished {
            return Err(CohortError::Corrupt(
                "summary window lifecycle is invalid".into(),
            ));
        }
        let count = sites
            .iter()
            .try_fold(0usize, |total, site| {
                total.checked_add(site.alternates.len())
            })
            .ok_or_else(|| CohortError::Limit("summary candidate count overflow".into()))?;
        if count > self.max_window_candidates {
            return Err(CohortError::Limit(
                "summary window exceeds declared reducer reservation".into(),
            ));
        }
        self.window.clear();
        self.window.reserve_exact(count);
        for site in sites {
            for &alternate in &site.alternates {
                let key = CandidateKey {
                    contig: site.contig,
                    position: site.position,
                    reference: site.reference,
                    alternate,
                };
                key.validate(self.contigs)?;
                if self
                    .window
                    .last()
                    .is_some_and(|previous| previous.key >= key)
                {
                    return Err(CohortError::Corrupt(
                        "summary candidates must be unique and canonically ordered".into(),
                    ));
                }
                self.window.push(WindowCandidate {
                    key,
                    reducer: CandidateReducer::new(self.n_requested, self.min_callable_depth)?,
                    last_member: None,
                });
            }
        }
        self.window_open = true;
        Ok(())
    }
    fn on_row(&mut self, row: CohortRow<'_>) -> Result<()> {
        if !self.window_open {
            return Err(CohortError::Corrupt(
                "summary row outside active window".into(),
            ));
        }
        let key = CandidateKey::from_row(row);
        let index = self
            .window
            .binary_search_by_key(&key, |candidate| candidate.key)
            .map_err(|_| {
                CohortError::Corrupt("summary row is not in the active candidate window".into())
            })?;
        let candidate = &mut self.window[index];
        if candidate
            .last_member
            .is_some_and(|previous| previous >= row.member_index)
        {
            return Err(CohortError::Corrupt(
                "summary member cells are repeated or out of order".into(),
            ));
        }
        let observation = match row.evidence {
            Some(evidence) => {
                if evidence.position != row.position || evidence.reference != row.reference {
                    return Err(CohortError::Corrupt(
                        "candidate key differs from evidence row".into(),
                    ));
                }
                CandidateObservation::from_row(evidence, row.alternate, self.min_callable_depth)?
            }
            None => CandidateObservation::unmeasured(self.min_callable_depth)?,
        };
        candidate.reducer.push(&observation)?;
        candidate.last_member = Some(row.member_index);
        Ok(())
    }
    fn end_window(&mut self) -> Result<()> {
        if !self.window_open {
            return Err(CohortError::Corrupt("no active summary window".into()));
        }
        // Check all denominators before emitting any summaries from this window.
        for candidate in &self.window {
            candidate.reducer.finish()?;
        }
        for index in 0..self.window.len() {
            let candidate = &self.window[index];
            self.cells.push(SummaryCell {
                key: candidate.key,
                summary: candidate.reducer.finish()?,
            });
            if self.cells.len() == OUTPUT_BATCH_ROWS {
                self.flush_batch()?;
            }
        }
        self.window.clear();
        self.window_open = false;
        Ok(())
    }
    fn finish(&mut self) -> Result<()> {
        if self.window_open {
            return Err(CohortError::Corrupt(
                "cannot finish an active summary window".into(),
            ));
        }
        if !self.sink.finished {
            self.flush_batch()?;
            self.sink.finish(self.record_batch()?.schema().as_ref())?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::evidence::{EvidenceAlleles, EvidenceDepths, EvidenceLocus, EvidenceRowRef};
    use arrow_ipc::reader::StreamReader;
    use std::io::Cursor;

    fn contigs() -> ContigSet {
        let mut contigs = ContigSet::new();
        contigs.push("chr1", 100_000);
        contigs
    }
    fn batches(bytes: Vec<u8>) -> Vec<RecordBatch> {
        StreamReader::try_new(Cursor::new(bytes), None)
            .unwrap()
            .map(|batch| batch.unwrap())
            .collect()
    }
    fn u64_column<'a>(batch: &'a RecordBatch, name: &str) -> &'a UInt64Array {
        batch
            .column_by_name(name)
            .unwrap()
            .as_any()
            .downcast_ref()
            .unwrap()
    }
    fn row<'a>(
        position: u32,
        depth: &'a EvidenceDepths,
        alleles: &'a EvidenceAlleles,
    ) -> EvidenceRowRef<'a> {
        EvidenceRowRef {
            position,
            reference: b'A',
            requested_alts: b"CG",
            depths: Some(depth),
            alleles: Some(alleles),
            strands: None,
            quality_sums: None,
            quality_histograms: None,
            read_position: None,
            allele_quality: None,
        }
    }
    fn cell<'a>(
        member_index: usize,
        position: u32,
        evidence: Option<EvidenceRowRef<'a>>,
    ) -> CohortRow<'a> {
        CohortRow {
            member_index,
            member_id: if member_index == 0 {
                "sample-A"
            } else {
                "sample-B"
            },
            contig: 0,
            position,
            reference: b'A',
            alternate: b'C',
            evidence,
        }
    }
    fn sites(start: u32, end: u32) -> Vec<SnvSite> {
        (start..end)
            .map(|position| SnvSite {
                contig: 0,
                position,
                reference: b'A',
                alternates: vec![b'C'],
            })
            .collect()
    }

    #[test]
    fn extraction_preserves_null_zero_and_uint64_in_arrow_and_tsv() {
        let contigs = contigs();
        let zero = EvidenceDepths::default();
        let no_alleles = EvidenceAlleles::default();
        let depth = EvidenceDepths {
            callable_depth: u64::MAX,
            ..EvidenceDepths::default()
        };
        let alleles = EvidenceAlleles {
            allele_counts: [0, u64::MAX, 0, 0],
        };
        let mut arrow = ExtractEncoder::new(
            Vec::new(),
            CohortOutputFormat::Arrow,
            &contigs,
            required_fields(),
            10,
        )
        .unwrap();
        let mut tsv = ExtractEncoder::new(
            Vec::new(),
            CohortOutputFormat::Tsv,
            &contigs,
            required_fields(),
            10,
        )
        .unwrap();
        for encoder in [&mut arrow, &mut tsv] {
            encoder.on_row(cell(0, 0, None)).unwrap();
            encoder
                .on_row(cell(0, 1, Some(row(1, &zero, &no_alleles))))
                .unwrap();
            encoder
                .on_row(cell(0, 2, Some(row(2, &depth, &alleles))))
                .unwrap();
        }
        let records = batches(arrow.into_inner().unwrap());
        assert_eq!(records.len(), 1);
        let record = &records[0];
        let callable = u64_column(record, "callable_depth");
        assert!(callable.is_null(0));
        assert_eq!(callable.value(1), 0);
        assert_eq!(callable.value(2), u64::MAX);
        let denominator = u64_column(record, "observed_alt_fraction_denominator");
        assert!(denominator.is_null(0));
        assert!(denominator.is_null(1));
        assert_eq!(denominator.value(2), u64::MAX);
        let eligible = record
            .column_by_name("depth_eligible")
            .unwrap()
            .as_any()
            .downcast_ref::<BooleanArray>()
            .unwrap();
        assert!(eligible.is_null(0));
        assert!(!eligible.value(1));
        assert!(eligible.value(2));
        let text = String::from_utf8(tsv.into_inner().unwrap()).unwrap();
        let lines: Vec<_> = text
            .lines()
            .map(|line| line.split('\t').collect::<Vec<_>>())
            .collect();
        let index = |name: &str| lines[0].iter().position(|field| *field == name).unwrap();
        assert_eq!(lines[1][index("status")], "unmeasured");
        assert_eq!(lines[1][index("callable_depth")], ".");
        assert_eq!(lines[2][index("callable_depth")], "0");
        assert_eq!(lines[2][index("observed_alt_fraction_denominator")], ".");
        assert_eq!(lines[3][index("alt_count")], "18446744073709551615");
        assert_eq!(
            lines[3][index("observed_alt_fraction_denominator")],
            "18446744073709551615"
        );
    }

    #[test]
    fn projected_optional_groups_are_nullable_including_all_unmeasured_lists() {
        let contigs = contigs();
        for observed in [false, true] {
            let fields = EvidenceFields::ALL_SUPPORTED;
            let batch = EvidenceBatch::new(
                0,
                "chr1",
                0,
                fields,
                vec![EvidenceLocus {
                    position: 1,
                    reference: b'A',
                    requested_alts: b"C".to_vec(),
                }],
            );
            let mut encoder =
                ExtractEncoder::new(Vec::new(), CohortOutputFormat::Arrow, &contigs, fields, 10)
                    .unwrap();
            encoder.on_row(cell(0, 0, None)).unwrap();
            if observed {
                encoder.on_row(cell(0, 1, batch.row(0))).unwrap();
            }
            let records = batches(encoder.into_inner().unwrap());
            let record = &records[0];
            for name in [
                "a_fwd",
                "mapping_quality_sum",
                "read_length_sum",
                "base_quality_histogram",
                "allele_base_quality_sum",
            ] {
                let array = record.column_by_name(name).unwrap();
                assert!(array.is_null(0), "{name}");
                assert!(record.schema().field_with_name(name).unwrap().is_nullable());
                if observed {
                    assert!(!array.is_null(1), "{name}");
                }
            }
        }
        let mut minimal = ExtractEncoder::new(
            Vec::new(),
            CohortOutputFormat::Arrow,
            &contigs,
            required_fields(),
            10,
        )
        .unwrap();
        minimal.on_row(cell(0, 0, None)).unwrap();
        assert!(batches(minimal.into_inner().unwrap())[0]
            .column_by_name("a_fwd")
            .is_none());
    }

    #[test]
    fn canonical_extract_batches_ignore_runtime_window_boundaries() {
        let contigs = contigs();
        let zero = EvidenceDepths::default();
        let alleles = EvidenceAlleles::default();
        for format in [CohortOutputFormat::Arrow, CohortOutputFormat::Tsv] {
            let mut outputs = Vec::new();
            for width in [1, 31, 2048] {
                let mut encoder =
                    ExtractEncoder::new(Vec::new(), format, &contigs, required_fields(), 10)
                        .unwrap();
                for start in (0..2053).step_by(width) {
                    let window = sites(start, (start + width as u32).min(2053));
                    encoder.begin_window(&window).unwrap();
                    for site in &window {
                        let evidence =
                            (site.position % 3 != 0).then(|| row(site.position, &zero, &alleles));
                        encoder.on_row(cell(0, site.position, evidence)).unwrap();
                    }
                    encoder.end_window().unwrap();
                }
                outputs.push(encoder.into_inner().unwrap());
            }
            assert_eq!(outputs[0], outputs[1]);
            assert_eq!(outputs[0], outputs[2]);
            if format == CohortOutputFormat::Arrow {
                assert_eq!(
                    batches(outputs.remove(0))
                        .iter()
                        .map(RecordBatch::num_rows)
                        .collect::<Vec<_>>(),
                    vec![1024, 1024, 5]
                );
            }
        }
    }

    #[test]
    fn summary_oracle_and_canonical_batches_ignore_window_width() {
        let contigs = contigs();
        let depth = EvidenceDepths {
            callable_depth: 10,
            ..EvidenceDepths::default()
        };
        let alleles = EvidenceAlleles {
            allele_counts: [9, 1, 0, 0],
        };
        for format in [CohortOutputFormat::Arrow, CohortOutputFormat::Tsv] {
            let mut outputs = Vec::new();
            for width in [1, 17, 2048] {
                let mut encoder =
                    SummaryEncoder::new(Vec::new(), format, &contigs, 2, 10, 2048).unwrap();
                for start in (0..2053).step_by(width) {
                    let window = sites(start, (start + width as u32).min(2053));
                    encoder.begin_window(&window).unwrap();
                    for member in 0..2 {
                        for site in &window {
                            encoder
                                .on_row(cell(
                                    member,
                                    site.position,
                                    (member == 0).then(|| row(site.position, &depth, &alleles)),
                                ))
                                .unwrap();
                        }
                    }
                    encoder.end_window().unwrap();
                }
                outputs.push(encoder.into_inner().unwrap());
            }
            assert_eq!(outputs[0], outputs[1]);
            assert_eq!(outputs[0], outputs[2]);
            if format == CohortOutputFormat::Arrow {
                let records = batches(outputs.remove(0));
                assert_eq!(
                    records
                        .iter()
                        .map(RecordBatch::num_rows)
                        .collect::<Vec<_>>(),
                    vec![1024, 1024, 5]
                );
                for (name, expected) in [
                    ("n_requested", 2),
                    ("n_observed", 1),
                    ("n_depth_eligible", 1),
                    ("n_alt_supported", 1),
                    ("callable_total", 10),
                    ("alt_total", 1),
                    ("observed_alt_fraction_numerator", 1),
                    ("observed_alt_fraction_denominator", 10),
                    ("depth_eligible_support_fraction_denominator", 1),
                ] {
                    assert_eq!(u64_column(&records[0], name).value(0), expected, "{name}");
                }
            } else {
                let text = String::from_utf8(outputs.remove(0)).unwrap();
                let mut lines = text.lines();
                let names: Vec<_> = lines.next().unwrap().split('\t').collect();
                let values: Vec<_> = lines.next().unwrap().split('\t').collect();
                for (name, value) in [
                    ("n_requested", "2"),
                    ("n_observed", "1"),
                    ("callable_total", "10"),
                    ("alt_total", "1"),
                ] {
                    assert_eq!(
                        values[names.iter().position(|field| *field == name).unwrap()],
                        value
                    );
                }
            }
        }
    }

    #[test]
    fn empty_selection_and_zero_members_have_defined_schemas_and_null_fractions() {
        let contigs = contigs();
        let extract = ExtractEncoder::new(
            Vec::new(),
            CohortOutputFormat::Arrow,
            &contigs,
            required_fields(),
            10,
        )
        .unwrap()
        .into_inner()
        .unwrap();
        let reader = StreamReader::try_new(Cursor::new(extract), None).unwrap();
        assert!(reader
            .schema()
            .field_with_name("callable_depth")
            .unwrap()
            .is_nullable());
        assert_eq!(reader.count(), 0);
        let empty = SummaryEncoder::new(Vec::new(), CohortOutputFormat::Arrow, &contigs, 0, 10, 0)
            .unwrap()
            .into_inner()
            .unwrap();
        assert!(batches(empty).is_empty());
        let mut summary =
            SummaryEncoder::new(Vec::new(), CohortOutputFormat::Arrow, &contigs, 0, 10, 1).unwrap();
        summary.begin_window(&sites(0, 1)).unwrap();
        summary.end_window().unwrap();
        let output = batches(summary.into_inner().unwrap());
        assert_eq!(u64_column(&output[0], "n_requested").value(0), 0);
        assert_eq!(u64_column(&output[0], "n_observed").value(0), 0);
        assert!(u64_column(&output[0], "observed_alt_fraction_denominator").is_null(0));
    }

    #[test]
    fn encoders_refuse_unsafe_text_unreserved_windows_and_repeated_cells() {
        let mut unsafe_contigs = ContigSet::new();
        unsafe_contigs.push("chr1\tother", 10);
        assert!(ExtractEncoder::new(
            Vec::new(),
            CohortOutputFormat::Tsv,
            &unsafe_contigs,
            required_fields(),
            10
        )
        .is_err());
        let contigs = contigs();
        let mut writer = Vec::new();
        {
            let mut encoder = ExtractEncoder::new(
                &mut writer,
                CohortOutputFormat::Tsv,
                &contigs,
                required_fields(),
                10,
            )
            .unwrap();
            let mut unsafe_row = cell(0, 0, None);
            unsafe_row.member_id = "line\nbreak";
            assert!(encoder.on_row(unsafe_row).is_err());
        }
        assert!(
            writer.is_empty(),
            "constructors and rejected rows write nothing"
        );
        let mut summary =
            SummaryEncoder::new(Vec::new(), CohortOutputFormat::Arrow, &contigs, 2, 10, 1).unwrap();
        assert!(summary.begin_window(&sites(0, 2)).is_err());
        summary.begin_window(&sites(0, 1)).unwrap();
        summary.on_row(cell(0, 0, None)).unwrap();
        assert!(summary.on_row(cell(0, 0, None)).is_err());
        assert!(
            summary.end_window().is_err(),
            "incomplete denominator cannot publish a summary"
        );
        summary.on_row(cell(1, 0, None)).unwrap();
        summary.end_window().unwrap();
        assert_eq!(batches(summary.into_inner().unwrap())[0].num_rows(), 1);
    }
}
