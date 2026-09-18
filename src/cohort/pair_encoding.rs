//! Bounded paired reports. Each side keeps its measurement and technical screen;
//! differences use exact sign/magnitude arithmetic, never floating point.

use super::encoding::{
    arrow_error, candidate_arrays, candidate_fields, max_text_length, output_schema,
    validate_contigs, validate_text, BatchSink, CandidateKey, CohortOutputFormat,
    OUTPUT_BATCH_ROWS,
};
use super::pairs::{compare_observations, PairedObservation};
use super::runtime::{CohortConsumer, CohortRow, PairConsumer, PairWindow};
use super::summary::{CandidateObservation, ObservationStatus};
use super::{CohortError, Result};
use crate::core::ContigSet;
use crate::evidence::{SnvSite, CANONICAL_TILE_BASES};
use arrow_array::{ArrayRef, BooleanArray, RecordBatch, StringArray, UInt64Array};
use arrow_schema::{DataType, Field};
use std::io::Write;
use std::sync::Arc;

struct WindowCandidate {
    key: CandidateKey,
    left: Option<CandidateObservation>,
    right: Option<CandidateObservation>,
}
struct PairIdentity {
    id: String,
    left: String,
    right: String,
    left_index: usize,
    right_index: usize,
}
struct PairCell {
    key: CandidateKey,
    pair_id: String,
    left_id: String,
    right_id: String,
    value: PairedObservation,
}

/// Retained one-window side observations plus the fixed output batch, nullable
/// Arrow arrays and IPC staging. The caller separately reserves file/receipt and
/// pair-table metadata. No reservation scales as pairs times candidate loci.
pub(crate) fn pair_memory_bytes(max_window_candidates: usize, max_text: usize) -> u64 {
    (max_window_candidates as u64)
        .saturating_mul((std::mem::size_of::<WindowCandidate>() as u64).saturating_mul(2))
        .saturating_add(
            (OUTPUT_BATCH_ROWS as u64)
                .saturating_mul((max_text as u64).saturating_mul(16).saturating_add(16_384)),
        )
        .saturating_add(1 << 20)
}

pub(crate) struct PairEncoder<'a, W: Write> {
    sink: BatchSink<W>,
    format: CohortOutputFormat,
    contigs: &'a ContigSet,
    min_callable_depth: u64,
    max_window_candidates: usize,
    active: Option<PairIdentity>,
    next_left: usize,
    next_right: usize,
    window: Vec<WindowCandidate>,
    cells: Vec<PairCell>,
    finished: bool,
}
impl<'a, W: Write> PairEncoder<'a, W> {
    pub(crate) fn new(
        output: W,
        format: CohortOutputFormat,
        contigs: &'a ContigSet,
        min_callable_depth: u64,
        max_window_candidates: usize,
    ) -> Result<Self> {
        validate_contigs(contigs, format)?;
        CandidateObservation::unmeasured(min_callable_depth)?;
        if max_window_candidates > 3 * CANONICAL_TILE_BASES as usize {
            return Err(CohortError::Limit(
                "paired window exceeds canonical candidate envelope".into(),
            ));
        }
        Ok(Self {
            sink: BatchSink::new(output, format),
            format,
            contigs,
            min_callable_depth,
            max_window_candidates,
            active: None,
            next_left: 0,
            next_right: 0,
            window: Vec::new(),
            cells: Vec::new(),
            finished: false,
        })
    }

    fn record_batch(&self) -> Result<RecordBatch> {
        let mut fields = vec![
            Field::new("pair_id", DataType::Utf8, false),
            Field::new("left_member_id", DataType::Utf8, false),
            Field::new("right_member_id", DataType::Utf8, false),
        ];
        let mut arrays: Vec<ArrayRef> = vec![
            Arc::new(StringArray::from_iter_values(
                self.cells.iter().map(|cell| cell.pair_id.as_str()),
            )),
            Arc::new(StringArray::from_iter_values(
                self.cells.iter().map(|cell| cell.left_id.as_str()),
            )),
            Arc::new(StringArray::from_iter_values(
                self.cells.iter().map(|cell| cell.right_id.as_str()),
            )),
        ];
        fields.extend(candidate_fields());
        let keys: Vec<_> = self.cells.iter().map(|cell| cell.key).collect();
        arrays.extend(candidate_arrays(&keys, self.contigs));
        for (prefix, left) in [("left", true), ("right", false)] {
            let observations: Vec<_> = self
                .cells
                .iter()
                .map(|cell| {
                    if left {
                        cell.value.left
                    } else {
                        cell.value.right
                    }
                })
                .collect();
            fields.push(Field::new(
                format!("{prefix}_status"),
                DataType::Utf8,
                false,
            ));
            arrays.push(Arc::new(StringArray::from_iter_values(
                observations.iter().map(|value| match value.status {
                    ObservationStatus::Observed => "observed",
                    ObservationStatus::Unmeasured => "unmeasured",
                }),
            )));
            for name in ["callable_depth", "alt_count"] {
                fields.push(Field::new(
                    format!("{prefix}_{name}"),
                    DataType::UInt64,
                    true,
                ));
            }
            arrays.push(Arc::new(UInt64Array::from_iter(
                observations.iter().map(|value| value.callable_depth),
            )));
            arrays.push(Arc::new(UInt64Array::from_iter(
                observations.iter().map(|value| value.alt_count),
            )));
            for name in ["depth_eligible", "alt_supported"] {
                fields.push(Field::new(
                    format!("{prefix}_{name}"),
                    DataType::Boolean,
                    true,
                ));
            }
            arrays.push(Arc::new(BooleanArray::from_iter(
                observations.iter().map(|value| value.depth_eligible),
            )));
            arrays.push(Arc::new(BooleanArray::from_iter(
                observations.iter().map(|value| value.alt_supported),
            )));
            for name in ["alt_fraction_numerator", "alt_fraction_denominator"] {
                fields.push(Field::new(
                    format!("{prefix}_{name}"),
                    DataType::UInt64,
                    true,
                ));
            }
            arrays.push(Arc::new(UInt64Array::from_iter(observations.iter().map(
                |value| {
                    value
                        .observed_alt_fraction
                        .map(|fraction| fraction.numerator)
                },
            ))));
            arrays.push(Arc::new(UInt64Array::from_iter(observations.iter().map(
                |value| {
                    value
                        .observed_alt_fraction
                        .map(|fraction| fraction.denominator)
                },
            ))));
        }
        fields.push(Field::new("min_callable_depth", DataType::UInt64, false));
        arrays.push(Arc::new(UInt64Array::from_iter_values(
            self.cells
                .iter()
                .map(|cell| cell.value.left.min_callable_depth),
        )));
        fields.push(Field::new("both_depth_eligible", DataType::Boolean, true));
        arrays.push(Arc::new(BooleanArray::from_iter(
            self.cells.iter().map(|cell| cell.value.both_depth_eligible),
        )));
        fields.push(Field::new("difference_negative", DataType::Boolean, true));
        arrays.push(Arc::new(BooleanArray::from_iter(self.cells.iter().map(
            |cell| cell.value.right_minus_left.map(|value| value.negative),
        ))));
        // u64 depth products need up to 39 decimal digits and exceed both i64
        // and i128. UTF8 decimal strings preserve their exact unsigned values.
        for name in ["difference_numerator", "difference_denominator"] {
            fields.push(Field::new(name, DataType::Utf8, true));
        }
        arrays.push(Arc::new(StringArray::from_iter(self.cells.iter().map(
            |cell| {
                cell.value
                    .right_minus_left
                    .map(|value| value.magnitude.to_string())
            },
        ))));
        arrays.push(Arc::new(StringArray::from_iter(self.cells.iter().map(
            |cell| {
                cell.value
                    .right_minus_left
                    .map(|value| value.denominator.to_string())
            },
        ))));
        RecordBatch::try_new(Arc::new(output_schema(fields, "paired-candidates")), arrays)
            .map_err(arrow_error)
    }

    fn flush_batch(&mut self) -> Result<()> {
        if !self.cells.is_empty() {
            self.sink.write(&self.record_batch()?)?;
            self.cells.clear();
        }
        Ok(())
    }

    pub(crate) fn into_inner(mut self) -> Result<W> {
        self.finish()?;
        self.sink.into_inner()
    }
}

impl<W: Write> PairConsumer for PairEncoder<'_, W> {
    fn begin_pair_window(&mut self, pair: PairWindow<'_>, sites: &[SnvSite]) -> Result<()> {
        if self.active.is_some() || self.finished {
            return Err(CohortError::Corrupt(
                "paired window lifecycle is invalid".into(),
            ));
        }
        for id in [pair.pair_id, pair.left_member_id, pair.right_member_id] {
            validate_text(id, self.format, 256)?;
            if id.is_empty() || id.chars().any(char::is_control) {
                return Err(CohortError::Incompatible(
                    "paired identifiers must be nonempty and have no controls".into(),
                ));
            }
        }
        if pair.left_member_index == pair.right_member_index
            || pair.left_member_id == pair.right_member_id
        {
            return Err(CohortError::Corrupt(
                "paired window requires distinct members".into(),
            ));
        }
        let count = sites
            .iter()
            .try_fold(0usize, |count, site| {
                count.checked_add(site.alternates.len())
            })
            .ok_or_else(|| CohortError::Limit("paired candidate count overflow".into()))?;
        if count > self.max_window_candidates {
            return Err(CohortError::Limit(
                "paired window exceeds declared reservation".into(),
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
                        "paired candidates must be unique and canonically ordered".into(),
                    ));
                }
                self.window.push(WindowCandidate {
                    key,
                    left: None,
                    right: None,
                });
            }
        }
        self.active = Some(PairIdentity {
            id: pair.pair_id.into(),
            left: pair.left_member_id.into(),
            right: pair.right_member_id.into(),
            left_index: pair.left_member_index,
            right_index: pair.right_member_index,
        });
        self.next_left = 0;
        self.next_right = 0;
        Ok(())
    }
}

impl<W: Write> CohortConsumer for PairEncoder<'_, W> {
    fn retained_bytes(&self) -> Option<u64> {
        Some(pair_memory_bytes(
            self.max_window_candidates,
            max_text_length(self.contigs),
        ))
    }

    fn on_row(&mut self, row: CohortRow<'_>) -> Result<()> {
        let pair = self
            .active
            .as_ref()
            .ok_or_else(|| CohortError::Corrupt("paired row outside active window".into()))?;
        let left = if row.member_index == pair.left_index && row.member_id == pair.left {
            true
        } else if row.member_index == pair.right_index && row.member_id == pair.right {
            false
        } else {
            return Err(CohortError::Corrupt(
                "paired row belongs to another member".into(),
            ));
        };
        if (!left && self.next_left != self.window.len()) || (left && self.next_right != 0) {
            return Err(CohortError::Corrupt(
                "paired sides must arrive serially left then right".into(),
            ));
        }
        let index = if left {
            self.next_left
        } else {
            self.next_right
        };
        let candidate = self
            .window
            .get_mut(index)
            .filter(|candidate| candidate.key == CandidateKey::from_row(row))
            .ok_or_else(|| {
                CohortError::Corrupt(
                    "paired rows are repeated, extra, or out of candidate order".into(),
                )
            })?;
        let observation = match row.evidence {
            Some(evidence) => {
                if evidence.position != row.position || evidence.reference != row.reference {
                    return Err(CohortError::Corrupt(
                        "paired candidate differs from evidence row".into(),
                    ));
                }
                CandidateObservation::from_row(evidence, row.alternate, self.min_callable_depth)?
            }
            None => CandidateObservation::unmeasured(self.min_callable_depth)?,
        };
        if left {
            candidate.left = Some(observation);
            self.next_left += 1;
        } else {
            candidate.right = Some(observation);
            self.next_right += 1;
        }
        Ok(())
    }

    fn end_window(&mut self) -> Result<()> {
        if self.active.is_none()
            || self.next_left != self.window.len()
            || self.next_right != self.window.len()
        {
            return Err(CohortError::Corrupt(
                "paired window is absent or lacks one side".into(),
            ));
        }
        // Check the complete window before any of its cells reach the sink.
        for candidate in &self.window {
            compare_observations(candidate.left.unwrap(), candidate.right.unwrap())?;
        }
        for index in 0..self.window.len() {
            let pair = self.active.as_ref().unwrap();
            let candidate = &self.window[index];
            self.cells.push(PairCell {
                key: candidate.key,
                pair_id: pair.id.clone(),
                left_id: pair.left.clone(),
                right_id: pair.right.clone(),
                value: compare_observations(candidate.left.unwrap(), candidate.right.unwrap())?,
            });
            if self.cells.len() == OUTPUT_BATCH_ROWS {
                self.flush_batch()?;
            }
        }
        self.window.clear();
        self.active = None;
        Ok(())
    }

    fn finish(&mut self) -> Result<()> {
        if self.active.is_some() {
            return Err(CohortError::Corrupt(
                "cannot finish active paired window".into(),
            ));
        }
        if !self.finished {
            self.flush_batch()?;
            self.sink.finish(self.record_batch()?.schema().as_ref())?;
            self.finished = true;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::evidence::{EvidenceAlleles, EvidenceDepths, EvidenceRowRef};
    use arrow_array::{Array, BooleanArray};
    use arrow_ipc::reader::StreamReader;
    use std::io::Cursor;

    fn contigs() -> ContigSet {
        let mut value = ContigSet::new();
        value.push("chr1", 100_000);
        value
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
    fn pair(reverse: bool) -> PairWindow<'static> {
        PairWindow {
            pair_id: if reverse { "reversed" } else { "forward" },
            left_member_id: if reverse { "B" } else { "A" },
            right_member_id: if reverse { "A" } else { "B" },
            left_member_index: if reverse { 1 } else { 0 },
            right_member_index: if reverse { 0 } else { 1 },
        }
    }
    fn push(
        encoder: &mut PairEncoder<'_, Vec<u8>>,
        index: usize,
        position: u32,
        value: Option<(u64, u64)>,
    ) -> Result<()> {
        let depth = EvidenceDepths {
            callable_depth: value.map_or(0, |v| v.0),
            ..EvidenceDepths::default()
        };
        let alleles = EvidenceAlleles {
            allele_counts: [
                value.map_or(0, |v| v.0 - v.1),
                value.map_or(0, |v| v.1),
                0,
                0,
            ],
        };
        encoder.on_row(CohortRow {
            member_index: index,
            member_id: if index == 0 { "A" } else { "B" },
            contig: 0,
            position,
            reference: b'A',
            alternate: b'C',
            evidence: value.map(|_| EvidenceRowRef {
                position,
                reference: b'A',
                requested_alts: b"C",
                depths: Some(&depth),
                alleles: Some(&alleles),
                strands: None,
                quality_sums: None,
                quality_histograms: None,
                read_position: None,
                allele_quality: None,
            }),
        })
    }
    fn batches(bytes: Vec<u8>) -> Vec<RecordBatch> {
        StreamReader::try_new(Cursor::new(bytes), None)
            .unwrap()
            .map(|batch| batch.unwrap())
            .collect()
    }
    fn string<'a>(batch: &'a RecordBatch, name: &str) -> &'a StringArray {
        batch
            .column_by_name(name)
            .unwrap()
            .as_any()
            .downcast_ref()
            .unwrap()
    }
    fn boolean<'a>(batch: &'a RecordBatch, name: &str) -> &'a BooleanArray {
        batch
            .column_by_name(name)
            .unwrap()
            .as_any()
            .downcast_ref()
            .unwrap()
    }

    #[test]
    fn paired_outputs_preserve_exact_large_differences_missingness_and_direction() {
        let contigs = contigs();
        for format in [CohortOutputFormat::Arrow, CohortOutputFormat::Tsv] {
            let mut encoder = PairEncoder::new(Vec::new(), format, &contigs, 10, 4).unwrap();
            for reverse in [false, true] {
                let context = pair(reverse);
                encoder.begin_pair_window(context, &sites(0, 4)).unwrap();
                for side in [context.left_member_index, context.right_member_index] {
                    for (position, value) in if side == 0 {
                        [
                            Some((u64::MAX, u64::MAX - 1)),
                            None,
                            Some((0, 0)),
                            Some((2, 1)),
                        ]
                    } else {
                        [
                            Some((u64::MAX - 1, u64::MAX - 2)),
                            Some((10, 3)),
                            Some((10, 3)),
                            Some((4, 3)),
                        ]
                    }
                    .into_iter()
                    .enumerate()
                    {
                        push(&mut encoder, side, position as u32, value).unwrap();
                    }
                }
                encoder.end_window().unwrap();
            }
            let bytes = encoder.into_inner().unwrap();
            if format == CohortOutputFormat::Arrow {
                let output = batches(bytes);
                assert_eq!(output.len(), 1);
                let output = &output[0];
                assert_eq!(output.num_rows(), 8);
                let numerator = string(output, "difference_numerator");
                let denominator = string(output, "difference_denominator");
                let negative = boolean(output, "difference_negative");
                assert_eq!(numerator.value(0), "1");
                assert_eq!(
                    denominator.value(0),
                    "340282366920938463408034375210639556610"
                );
                assert!(negative.value(0));
                assert!(!negative.value(4));
                for index in [1, 2, 5, 6] {
                    assert!(numerator.is_null(index));
                    assert!(negative.is_null(index));
                }
                assert_eq!(numerator.value(3), "2");
                assert_eq!(denominator.value(3), "8");
                assert!(!negative.value(3));
                assert!(negative.value(7));
                assert!(boolean(output, "both_depth_eligible").is_null(1));
                assert!(!boolean(output, "both_depth_eligible").value(2));
                assert!(!boolean(output, "both_depth_eligible").value(3));
                assert_eq!(string(output, "left_status").value(1), "unmeasured");
                assert_eq!(string(output, "left_status").value(2), "observed");
            } else {
                let text = String::from_utf8(bytes).unwrap();
                let rows: Vec<_> = text
                    .lines()
                    .map(|line| line.split('\t').collect::<Vec<_>>())
                    .collect();
                let column = |name| rows[0].iter().position(|value| *value == name).unwrap();
                assert_eq!(
                    rows[1][column("difference_denominator")],
                    "340282366920938463408034375210639556610"
                );
                assert_eq!(rows[1][column("difference_negative")], "true");
                assert_eq!(rows[5][column("difference_negative")], "false");
                assert_eq!(rows[2][column("left_callable_depth")], ".");
                assert_eq!(rows[3][column("left_callable_depth")], "0");
                assert_eq!(rows[3][column("difference_numerator")], ".");
            }
        }
    }

    #[test]
    fn canonical_paired_batches_cross_window_and_pair_boundaries() {
        let contigs = contigs();
        for format in [CohortOutputFormat::Arrow, CohortOutputFormat::Tsv] {
            let mut outputs = Vec::new();
            for width in [1, 31, 2048] {
                let mut encoder = PairEncoder::new(Vec::new(), format, &contigs, 10, 2048).unwrap();
                for reverse in [false, true] {
                    let context = pair(reverse);
                    for start in (0..1101).step_by(width) {
                        let window = sites(start, (start + width as u32).min(1101));
                        encoder.begin_pair_window(context, &window).unwrap();
                        for side in [context.left_member_index, context.right_member_index] {
                            for site in &window {
                                push(
                                    &mut encoder,
                                    side,
                                    site.position,
                                    (site.position % 3 != 0).then_some((12, side as u64)),
                                )
                                .unwrap();
                            }
                        }
                        encoder.end_window().unwrap();
                    }
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
                    vec![1024, 1024, 154]
                );
            }
        }
    }

    #[test]
    fn paired_lifecycle_rejects_unreserved_incomplete_reordered_and_wrong_member_rows() {
        let contigs = contigs();
        let mut encoder =
            PairEncoder::new(Vec::new(), CohortOutputFormat::Arrow, &contigs, 10, 1).unwrap();
        assert!(encoder
            .begin_pair_window(pair(false), &sites(0, 2))
            .is_err());
        encoder
            .begin_pair_window(pair(false), &sites(0, 1))
            .unwrap();
        assert!(push(&mut encoder, 1, 0, None).is_err());
        assert!(push(&mut encoder, 2, 0, None).is_err());
        assert!(push(&mut encoder, 0, 1, None).is_err());
        push(&mut encoder, 0, 0, None).unwrap();
        assert!(push(&mut encoder, 0, 0, None).is_err());
        assert!(encoder.end_window().is_err());
        assert!(encoder.finish().is_err());
        push(&mut encoder, 1, 0, Some((0, 0))).unwrap();
        encoder.end_window().unwrap();
        assert_eq!(batches(encoder.into_inner().unwrap())[0].num_rows(), 1);
        let empty = PairEncoder::new(Vec::new(), CohortOutputFormat::Arrow, &contigs, 10, 0)
            .unwrap()
            .into_inner()
            .unwrap();
        let reader = StreamReader::try_new(Cursor::new(empty), None).unwrap();
        assert!(reader
            .schema()
            .field_with_name("difference_denominator")
            .unwrap()
            .is_nullable());
        assert_eq!(reader.count(), 0);
    }
}
