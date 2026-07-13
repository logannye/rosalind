//! Canonical Arrow IPC feature output.

use std::io::{self, Write};
use std::sync::Arc;

use arrow_array::{ArrayRef, Float64Array, RecordBatch, StringArray, UInt32Array};
use arrow_ipc::writer::{IpcWriteOptions, StreamWriter};
use arrow_ipc::MetadataVersion;
use arrow_schema::{DataType, Field, Schema};

use crate::call::ColumnAnalyzer;
use crate::pileup::PileupColumn;

/// Versioned Arrow feature schema.
pub const FEATURE_ARROW_SCHEMA_VERSION: u32 = 1;
/// Canonical record-batch row count.
pub const FEATURE_ARROW_BATCH_ROWS: usize = 65_536;

fn schema() -> Schema {
    let mut fields = vec![
        Field::new("contig", DataType::Utf8, false),
        Field::new("pos", DataType::UInt32, false),
        Field::new("ref", DataType::Utf8, false),
    ];
    for name in [
        "depth",
        "raw_depth",
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
    ] {
        fields.push(Field::new(name, DataType::UInt32, false));
    }
    fields.push(Field::new("mean_bq", DataType::Float64, false));
    fields.push(Field::new("mean_mapq", DataType::Float64, false));
    Schema::new(fields)
}

#[derive(Debug, Default)]
struct ChunkWriter {
    bytes: Vec<u8>,
}

impl Write for ChunkWriter {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.bytes.extend_from_slice(bytes);
        Ok(bytes.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

#[derive(Debug)]
struct FeatureRow {
    contig: String,
    pos: u32,
    reference: String,
    counts: [u32; 14],
    mean_bq: f64,
    mean_mapq: f64,
}

/// First-party canonical Arrow IPC feature analyzer.
pub struct FeatureArrowAnalyzer {
    writer: StreamWriter<ChunkWriter>,
    rows: Vec<FeatureRow>,
    emitted_rows: u64,
    finished: bool,
}

impl std::fmt::Debug for FeatureArrowAnalyzer {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("FeatureArrowAnalyzer")
            .field("buffered_rows", &self.rows.len())
            .field("emitted_rows", &self.emitted_rows)
            .field("finished", &self.finished)
            .finish()
    }
}

impl FeatureArrowAnalyzer {
    /// Construct an encoder with a schema header already buffered.
    pub fn new() -> io::Result<Self> {
        let options =
            IpcWriteOptions::try_new(8, false, MetadataVersion::V5).map_err(arrow_error)?;
        let writer = StreamWriter::try_new_with_options(ChunkWriter::default(), &schema(), options)
            .map_err(arrow_error)?;
        Ok(Self {
            writer,
            rows: Vec::with_capacity(FEATURE_ARROW_BATCH_ROWS),
            emitted_rows: 0,
            finished: false,
        })
    }

    /// Number of feature rows accepted.
    pub fn rows(&self) -> u64 {
        self.emitted_rows + self.rows.len() as u64
    }

    fn flush_bytes(&mut self, out: &mut dyn Write) -> io::Result<()> {
        let bytes = std::mem::take(&mut self.writer.get_mut().bytes);
        out.write_all(&bytes)
    }

    fn flush_batch(&mut self, out: &mut dyn Write) -> io::Result<()> {
        self.flush_bytes(out)?;
        if self.rows.is_empty() {
            return Ok(());
        }
        let rows = std::mem::take(&mut self.rows);
        let mut arrays: Vec<ArrayRef> = Vec::with_capacity(19);
        arrays.push(Arc::new(StringArray::from_iter_values(
            rows.iter().map(|row| row.contig.as_str()),
        )));
        arrays.push(Arc::new(UInt32Array::from_iter_values(
            rows.iter().map(|row| row.pos),
        )));
        arrays.push(Arc::new(StringArray::from_iter_values(
            rows.iter().map(|row| row.reference.as_str()),
        )));
        for index in 0..14 {
            arrays.push(Arc::new(UInt32Array::from_iter_values(
                rows.iter().map(|row| row.counts[index]),
            )));
        }
        arrays.push(Arc::new(Float64Array::from_iter_values(
            rows.iter().map(|row| row.mean_bq),
        )));
        arrays.push(Arc::new(Float64Array::from_iter_values(
            rows.iter().map(|row| row.mean_mapq),
        )));
        let batch = RecordBatch::try_new(Arc::new(schema()), arrays).map_err(arrow_error)?;
        self.writer.write(&batch).map_err(arrow_error)?;
        self.emitted_rows += rows.len() as u64;
        self.rows = Vec::with_capacity(FEATURE_ARROW_BATCH_ROWS);
        self.flush_bytes(out)
    }
}

impl Default for FeatureArrowAnalyzer {
    fn default() -> Self {
        Self::new().expect("canonical Arrow schema and options are valid")
    }
}

impl ColumnAnalyzer for FeatureArrowAnalyzer {
    fn params(&self) -> std::collections::BTreeMap<String, String> {
        std::collections::BTreeMap::from([
            ("feature_rows".into(), self.rows().to_string()),
            (
                "feature.schema".into(),
                FEATURE_ARROW_SCHEMA_VERSION.to_string(),
            ),
            (
                "feature.batch_rows".into(),
                FEATURE_ARROW_BATCH_ROWS.to_string(),
            ),
            ("artifact.format".into(), "arrow-ipc".into()),
        ])
    }

    fn on_column(
        &mut self,
        column: &PileupColumn,
        contig: &str,
        out: &mut dyn Write,
    ) -> io::Result<()> {
        let depth = column.depth();
        let allele = column.allele_counts();
        let strand = column.strand_counts();
        let (sum_bq, sum_mapq) = column.obs.iter().fold((0u64, 0u64), |(bq, mapq), obs| {
            (bq + obs.base_qual as u64, mapq + obs.mapq as u64)
        });
        let mean = |sum: u64| {
            if depth == 0 {
                0.0
            } else {
                let scaled = (sum * 100 + depth as u64 / 2) / depth as u64;
                scaled as f64 / 100.0
            }
        };
        self.rows.push(FeatureRow {
            contig: contig.to_string(),
            pos: column.locus.pos.0 + 1,
            reference: char::from(column.ref_base).to_string(),
            counts: [
                depth,
                column.raw_depth,
                allele[0],
                allele[1],
                allele[2],
                allele[3],
                strand[0][0],
                strand[0][1],
                strand[1][0],
                strand[1][1],
                strand[2][0],
                strand[2][1],
                strand[3][0],
                strand[3][1],
            ],
            mean_bq: mean(sum_bq),
            mean_mapq: mean(sum_mapq),
        });
        if self.rows.len() == FEATURE_ARROW_BATCH_ROWS {
            self.flush_batch(out)?;
        }
        Ok(())
    }

    fn finish(&mut self, out: &mut dyn Write) -> io::Result<()> {
        if self.finished {
            return Ok(());
        }
        self.flush_batch(out)?;
        self.writer.finish().map_err(arrow_error)?;
        self.flush_bytes(out)?;
        self.finished = true;
        Ok(())
    }
}

fn arrow_error(error: impl std::fmt::Display) -> io::Error {
    io::Error::other(error.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::{Locus, Position};
    use crate::pileup::Obs;
    use arrow_ipc::reader::StreamReader;

    fn encode_rows(rows: usize) -> Vec<u8> {
        let mut analyzer = FeatureArrowAnalyzer::new().unwrap();
        let mut bytes = Vec::new();
        for position in 0..rows {
            let column = PileupColumn {
                locus: Locus {
                    contig: 0,
                    pos: Position(position as u32),
                },
                ref_base: b'A',
                raw_depth: 1,
                obs: vec![Obs {
                    allele: 0,
                    base_qual: 30,
                    mapq: 60,
                    reverse: false,
                }],
            };
            analyzer.on_column(&column, "chr1", &mut bytes).unwrap();
        }
        analyzer.finish(&mut bytes).unwrap();
        bytes
    }

    #[test]
    fn empty_stream_has_schema_and_no_batches() {
        let bytes = encode_rows(0);
        let reader = StreamReader::try_new(std::io::Cursor::new(bytes), None).unwrap();
        assert_eq!(reader.schema().fields().len(), 19);
        assert_eq!(reader.count(), 0);
    }

    #[test]
    fn schema_field_order_types_and_nullability_are_frozen() {
        let fields = schema().fields().to_vec();
        let expected = [
            ("contig", DataType::Utf8),
            ("pos", DataType::UInt32),
            ("ref", DataType::Utf8),
            ("depth", DataType::UInt32),
            ("raw_depth", DataType::UInt32),
            ("a", DataType::UInt32),
            ("c", DataType::UInt32),
            ("g", DataType::UInt32),
            ("t", DataType::UInt32),
            ("a_fwd", DataType::UInt32),
            ("a_rev", DataType::UInt32),
            ("c_fwd", DataType::UInt32),
            ("c_rev", DataType::UInt32),
            ("g_fwd", DataType::UInt32),
            ("g_rev", DataType::UInt32),
            ("t_fwd", DataType::UInt32),
            ("t_rev", DataType::UInt32),
            ("mean_bq", DataType::Float64),
            ("mean_mapq", DataType::Float64),
        ];
        assert_eq!(fields.len(), expected.len());
        for (field, (name, data_type)) in fields.iter().zip(expected) {
            assert_eq!(field.name(), name);
            assert_eq!(field.data_type(), &data_type);
            assert!(!field.is_nullable());
        }
        assert!(schema().metadata().is_empty());
    }

    #[test]
    fn fixed_batch_boundaries_and_repeat_bytes() {
        for rows in [0, 1, 65_535, 65_536, 65_537, 131_072] {
            let first = encode_rows(rows);
            assert_eq!(
                first,
                encode_rows(rows),
                "{rows} rows must repeat byte-for-byte"
            );
            let batches = StreamReader::try_new(std::io::Cursor::new(first), None)
                .unwrap()
                .map(|batch| batch.unwrap().num_rows())
                .collect::<Vec<_>>();
            let full = rows / FEATURE_ARROW_BATCH_ROWS;
            let tail = rows % FEATURE_ARROW_BATCH_ROWS;
            let mut expected = vec![FEATURE_ARROW_BATCH_ROWS; full];
            if tail != 0 {
                expected.push(tail);
            }
            assert_eq!(batches, expected);
        }
    }
}
