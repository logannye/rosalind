use super::*;

/// Identity and requested alleles for one selected locus, independent of metrics.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EvidenceLocus {
    /// Zero-based coordinate within the batch contig.
    pub position: u32,
    /// Reference base, or N when no analysis sequence is available.
    pub reference: u8,
    /// Requested SNV alternatives in canonical A/C/G/T order.
    pub requested_alts: Vec<u8>,
}

/// Depths and mutually exclusive first-failure counters.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct EvidenceDepths {
    /// Matched observations before profile filters.
    pub prefilter_depth: u64,
    /// Observations passing read-level filters.
    pub aligned_depth: u64,
    /// A/C/G/T observations passing every filter.
    pub callable_depth: u64,
    /// Exclusive rejection reasons under the declared profile.
    pub filters: EvidenceFilterCounts,
}

/// Callable observations by nucleotide.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct EvidenceAlleles {
    /// Counts in A/C/G/T order.
    pub allele_counts: [u64; 4],
}

/// Callable observations by nucleotide and alignment strand.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct EvidenceStrands {
    /// A/C/G/T counts, each ordered forward then reverse.
    pub strand_counts: [[u64; 2]; 4],
}

/// Exact quality sufficient statistics over callable observations.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct EvidenceQualitySums {
    /// Sum of callable base qualities.
    pub base_quality_sum: u64,
    /// Sum of callable mapping qualities.
    pub mapping_quality_sum: u64,
}

/// Exact quality distributions over callable observations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EvidenceQualityHistograms {
    /// Counts for base qualities 0 through 93.
    pub base_quality_histogram: [u64; BASE_QUALITY_BINS],
    /// Counts for mapping qualities 0 through 254.
    pub mapping_quality_histogram: [u64; MAPPING_QUALITY_BINS],
}
impl Default for EvidenceQualityHistograms {
    fn default() -> Self {
        Self {
            base_quality_histogram: [0; BASE_QUALITY_BINS],
            mapping_quality_histogram: [0; MAPPING_QUALITY_BINS],
        }
    }
}

/// Exact stored-read-position sufficient statistics over callable observations.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct EvidenceReadPosition {
    /// Zero-based stored-SEQ offsets, reversed back on reverse alignments.
    pub read_position_sum: u64,
    /// Full stored sequence lengths, including soft-clipped bases.
    pub read_length_sum: u64,
}

/// Selected loci with physically present, contiguous summary groups. An absent
/// group has no allocation and does not imply zero measurements. Group lengths
/// and availability are maintained by the batch constructors and append methods.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EvidenceBatch {
    /// Stable identifier from the reference dictionary.
    pub contig_id: u32,
    /// Canonical contig name.
    pub contig: String,
    /// Start of the immutable 16,384-base ownership tile.
    pub canonical_tile_start: u32,
    fields: EvidenceFields,
    loci: Vec<EvidenceLocus>,
    depths: Option<Vec<EvidenceDepths>>,
    alleles: Option<Vec<EvidenceAlleles>>,
    strands: Option<Vec<EvidenceStrands>>,
    quality_sums: Option<Vec<EvidenceQualitySums>>,
    quality_histograms: Option<Vec<EvidenceQualityHistograms>>,
    read_position: Option<Vec<EvidenceReadPosition>>,
}

/// A borrowed locus. Check a group's availability before reading its metrics.
/// Materializing an owned full row is explicit and fails for a projection.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EvidenceRowRef<'a> {
    /// Zero-based coordinate within the contig.
    pub position: u32,
    /// Reference base or N when unavailable.
    pub reference: u8,
    /// Requested alternatives; empty for interval selection.
    pub requested_alts: &'a [u8],
    /// Present only when DEPTHS was requested.
    pub depths: Option<&'a EvidenceDepths>,
    /// Present only when ALLELES was requested.
    pub alleles: Option<&'a EvidenceAlleles>,
    /// Present only when STRANDS was requested.
    pub strands: Option<&'a EvidenceStrands>,
    /// Present only when QUALITY_SUMS was requested.
    pub quality_sums: Option<&'a EvidenceQualitySums>,
    /// Present only when QUALITY_HISTOGRAMS was requested.
    pub quality_histograms: Option<&'a EvidenceQualityHistograms>,
    /// Present only when READ_POSITION was requested.
    pub read_position: Option<&'a EvidenceReadPosition>,
}
impl EvidenceRowRef<'_> {
    /// Physical groups available on this row.
    pub fn fields(&self) -> EvidenceFields {
        let mut fields = EvidenceFields::from_bits(0).expect("empty field mask");
        for (present, group) in [
            (self.depths.is_some(), EvidenceFields::DEPTHS),
            (self.alleles.is_some(), EvidenceFields::ALLELES),
            (self.strands.is_some(), EvidenceFields::STRANDS),
            (self.quality_sums.is_some(), EvidenceFields::QUALITY_SUMS),
            (
                self.quality_histograms.is_some(),
                EvidenceFields::QUALITY_HISTOGRAMS,
            ),
            (self.read_position.is_some(), EvidenceFields::READ_POSITION),
        ] {
            if present {
                fields = fields.union(group);
            }
        }
        fields
    }
    /// Copy a complete row for legacy full-evidence consumers. Missing metrics
    /// are never synthesized; the caller must explicitly request ALL first.
    pub fn try_to_full_row(&self) -> Result<EvidenceRow, EvidenceError> {
        if self.fields() != EvidenceFields::ALL {
            return Err(EvidenceError::InvalidRequest(
                "cannot materialize a full row from projected evidence".into(),
            ));
        }
        let depths = self.depths.expect("batch field invariant");
        let quality = self.quality_sums.expect("batch field invariant");
        let histograms = self.quality_histograms.expect("batch field invariant");
        let position = self.read_position.expect("batch field invariant");
        Ok(EvidenceRow {
            position: self.position,
            reference: self.reference,
            requested_alts: self.requested_alts.to_vec(),
            prefilter_depth: depths.prefilter_depth,
            aligned_depth: depths.aligned_depth,
            callable_depth: depths.callable_depth,
            filters: depths.filters,
            allele_counts: self.alleles.expect("batch field invariant").allele_counts,
            strand_counts: self.strands.expect("batch field invariant").strand_counts,
            base_quality_sum: quality.base_quality_sum,
            mapping_quality_sum: quality.mapping_quality_sum,
            base_quality_histogram: histograms.base_quality_histogram,
            mapping_quality_histogram: histograms.mapping_quality_histogram,
            read_position_sum: position.read_position_sum,
            read_length_sum: position.read_length_sum,
        })
    }
}

/// Internal mutation borrows preserve the set and length of allocated groups.
pub(crate) struct EvidenceRowMut<'a> {
    pub position: u32,
    pub reference: u8,
    pub depths: Option<&'a mut EvidenceDepths>,
    pub alleles: Option<&'a mut EvidenceAlleles>,
    pub strands: Option<&'a mut EvidenceStrands>,
    pub quality_sums: Option<&'a mut EvidenceQualitySums>,
    pub quality_histograms: Option<&'a mut EvidenceQualityHistograms>,
    pub read_position: Option<&'a mut EvidenceReadPosition>,
}

impl EvidenceBatch {
    /// Allocate zero-initialized counters only for requested groups. Every locus
    /// remains present, including selected positions with no observations.
    pub fn new(
        contig_id: u32,
        contig: impl Into<String>,
        canonical_tile_start: u32,
        fields: EvidenceFields,
        loci: Vec<EvidenceLocus>,
    ) -> Self {
        let len = loci.len();
        Self {
            contig_id,
            contig: contig.into(),
            canonical_tile_start,
            fields,
            loci,
            depths: fields
                .contains(EvidenceFields::DEPTHS)
                .then(|| vec![EvidenceDepths::default(); len]),
            alleles: fields
                .contains(EvidenceFields::ALLELES)
                .then(|| vec![EvidenceAlleles::default(); len]),
            strands: fields
                .contains(EvidenceFields::STRANDS)
                .then(|| vec![EvidenceStrands::default(); len]),
            quality_sums: fields
                .contains(EvidenceFields::QUALITY_SUMS)
                .then(|| vec![EvidenceQualitySums::default(); len]),
            quality_histograms: fields
                .contains(EvidenceFields::QUALITY_HISTOGRAMS)
                .then(|| vec![EvidenceQualityHistograms::default(); len]),
            read_position: fields
                .contains(EvidenceFields::READ_POSITION)
                .then(|| vec![EvidenceReadPosition::default(); len]),
        }
    }

    /// Explicit full-schema adapter for existing fixtures and owned consumers.
    /// The engine and encoders use projected storage directly.
    pub fn from_full_rows(
        contig_id: u32,
        contig: impl Into<String>,
        canonical_tile_start: u32,
        rows: Vec<EvidenceRow>,
    ) -> Self {
        let mut batch = Self::new(
            contig_id,
            contig,
            canonical_tile_start,
            EvidenceFields::ALL,
            Vec::new(),
        );
        batch.reserve_exact(rows.len());
        for row in rows {
            batch.loci.push(EvidenceLocus {
                position: row.position,
                reference: row.reference,
                requested_alts: row.requested_alts,
            });
            batch.depths.as_mut().unwrap().push(EvidenceDepths {
                prefilter_depth: row.prefilter_depth,
                aligned_depth: row.aligned_depth,
                callable_depth: row.callable_depth,
                filters: row.filters,
            });
            batch.alleles.as_mut().unwrap().push(EvidenceAlleles {
                allele_counts: row.allele_counts,
            });
            batch.strands.as_mut().unwrap().push(EvidenceStrands {
                strand_counts: row.strand_counts,
            });
            batch
                .quality_sums
                .as_mut()
                .unwrap()
                .push(EvidenceQualitySums {
                    base_quality_sum: row.base_quality_sum,
                    mapping_quality_sum: row.mapping_quality_sum,
                });
            batch
                .quality_histograms
                .as_mut()
                .unwrap()
                .push(EvidenceQualityHistograms {
                    base_quality_histogram: row.base_quality_histogram,
                    mapping_quality_histogram: row.mapping_quality_histogram,
                });
            batch
                .read_position
                .as_mut()
                .unwrap()
                .push(EvidenceReadPosition {
                    read_position_sum: row.read_position_sum,
                    read_length_sum: row.read_length_sum,
                });
        }
        batch
    }

    /// Number of selected loci in this batch.
    pub fn len(&self) -> usize {
        self.loci.len()
    }
    /// Whether no loci are present. The declared field set remains meaningful.
    pub fn is_empty(&self) -> bool {
        self.loci.is_empty()
    }
    /// Physical groups present in every row.
    pub fn fields(&self) -> EvidenceFields {
        self.fields
    }
    /// Borrow the locus identities without materializing metrics.
    pub fn loci(&self) -> &[EvidenceLocus] {
        &self.loci
    }
    /// Borrow one selected row by its index.
    pub fn row(&self, index: usize) -> Option<EvidenceRowRef<'_>> {
        let locus = self.loci.get(index)?;
        Some(EvidenceRowRef {
            position: locus.position,
            reference: locus.reference,
            requested_alts: &locus.requested_alts,
            depths: self.depths.as_ref().map(|values| &values[index]),
            alleles: self.alleles.as_ref().map(|values| &values[index]),
            strands: self.strands.as_ref().map(|values| &values[index]),
            quality_sums: self.quality_sums.as_ref().map(|values| &values[index]),
            quality_histograms: self
                .quality_histograms
                .as_ref()
                .map(|values| &values[index]),
            read_position: self.read_position.as_ref().map(|values| &values[index]),
        })
    }
    /// Iterate borrowed rows in their existing coordinate order.
    pub fn rows(&self) -> impl ExactSizeIterator<Item = EvidenceRowRef<'_>> + DoubleEndedIterator {
        (0..self.len()).map(|index| self.row(index).expect("batch index invariant"))
    }
    pub(crate) fn row_mut(&mut self, index: usize) -> Option<EvidenceRowMut<'_>> {
        let locus = self.loci.get(index)?;
        Some(EvidenceRowMut {
            position: locus.position,
            reference: locus.reference,
            depths: self.depths.as_mut().map(|values| &mut values[index]),
            alleles: self.alleles.as_mut().map(|values| &mut values[index]),
            strands: self.strands.as_mut().map(|values| &mut values[index]),
            quality_sums: self.quality_sums.as_mut().map(|values| &mut values[index]),
            quality_histograms: self
                .quality_histograms
                .as_mut()
                .map(|values| &mut values[index]),
            read_position: self.read_position.as_mut().map(|values| &mut values[index]),
        })
    }
    pub(crate) fn reserve_exact(&mut self, additional: usize) {
        self.loci.reserve_exact(additional);
        if let Some(values) = &mut self.depths {
            values.reserve_exact(additional);
        }
        if let Some(values) = &mut self.alleles {
            values.reserve_exact(additional);
        }
        if let Some(values) = &mut self.strands {
            values.reserve_exact(additional);
        }
        if let Some(values) = &mut self.quality_sums {
            values.reserve_exact(additional);
        }
        if let Some(values) = &mut self.quality_histograms {
            values.reserve_exact(additional);
        }
        if let Some(values) = &mut self.read_position {
            values.reserve_exact(additional);
        }
    }
    pub(crate) fn clear(&mut self) {
        self.loci.clear();
        if let Some(values) = &mut self.depths {
            values.clear();
        }
        if let Some(values) = &mut self.alleles {
            values.clear();
        }
        if let Some(values) = &mut self.strands {
            values.clear();
        }
        if let Some(values) = &mut self.quality_sums {
            values.clear();
        }
        if let Some(values) = &mut self.quality_histograms {
            values.clear();
        }
        if let Some(values) = &mut self.read_position {
            values.clear();
        }
    }
    /// Replace locus identities while retaining each summary allocation. The
    /// caller reuses the returned identity vector as its next-window scratch.
    pub(crate) fn reset_loci(&mut self, loci: Vec<EvidenceLocus>) -> Vec<EvidenceLocus> {
        let len = loci.len();
        macro_rules! reset {
            ($values:expr, $ty:ty) => {
                if let Some(values) = $values {
                    values.resize(len, <$ty>::default());
                    values.fill(<$ty>::default());
                }
            };
        }
        reset!(&mut self.depths, EvidenceDepths);
        reset!(&mut self.alleles, EvidenceAlleles);
        reset!(&mut self.strands, EvidenceStrands);
        reset!(&mut self.quality_sums, EvidenceQualitySums);
        reset!(&mut self.quality_histograms, EvidenceQualityHistograms);
        reset!(&mut self.read_position, EvidenceReadPosition);
        std::mem::replace(&mut self.loci, loci)
    }
    pub(crate) fn push_row(&mut self, row: EvidenceRowRef<'_>) -> Result<(), EvidenceError> {
        if !row.fields().contains(self.fields) {
            return Err(EvidenceError::InvalidRequest(
                "source row lacks requested evidence groups".into(),
            ));
        }
        self.loci.push(EvidenceLocus {
            position: row.position,
            reference: row.reference,
            requested_alts: row.requested_alts.to_vec(),
        });
        if let Some(values) = &mut self.depths {
            values.push(*row.depths.expect("source field invariant"));
        }
        if let Some(values) = &mut self.alleles {
            values.push(*row.alleles.expect("source field invariant"));
        }
        if let Some(values) = &mut self.strands {
            values.push(*row.strands.expect("source field invariant"));
        }
        if let Some(values) = &mut self.quality_sums {
            values.push(*row.quality_sums.expect("source field invariant"));
        }
        if let Some(values) = &mut self.quality_histograms {
            values.push(*row.quality_histograms.expect("source field invariant"));
        }
        if let Some(values) = &mut self.read_position {
            values.push(*row.read_position.expect("source field invariant"));
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn loci(count: u32) -> Vec<EvidenceLocus> {
        (0..count)
            .map(|position| EvidenceLocus {
                position,
                reference: b'A',
                requested_alts: Vec::new(),
            })
            .collect()
    }

    #[test]
    fn every_mask_allocates_exactly_its_requested_groups() {
        for bits in 0..=EvidenceFields::ALL.bits() {
            let fields = EvidenceFields::from_bits(bits).unwrap();
            let batch = EvidenceBatch::new(0, "chr1", 0, fields, loci(3));
            assert_eq!(batch.len(), 3);
            assert_eq!(batch.fields(), fields);
            let mut bytes = batch.loci.capacity() * std::mem::size_of::<EvidenceLocus>();
            macro_rules! check {
                ($member:ident, $group:ident, $ty:ty) => {
                    assert_eq!(
                        batch.$member.is_some(),
                        fields.contains(EvidenceFields::$group)
                    );
                    if let Some(values) = &batch.$member {
                        assert_eq!(values.len(), 3);
                        assert_eq!(values.capacity(), 3);
                        bytes += values.capacity() * std::mem::size_of::<$ty>();
                    }
                };
            }
            check!(depths, DEPTHS, EvidenceDepths);
            check!(alleles, ALLELES, EvidenceAlleles);
            check!(strands, STRANDS, EvidenceStrands);
            check!(quality_sums, QUALITY_SUMS, EvidenceQualitySums);
            check!(
                quality_histograms,
                QUALITY_HISTOGRAMS,
                EvidenceQualityHistograms
            );
            check!(read_position, READ_POSITION, EvidenceReadPosition);
            assert_eq!(bytes as u64, fields.storage_bytes_per_locus() * 3);
            for row in batch.rows() {
                assert_eq!(row.fields(), fields);
                assert_eq!(row.try_to_full_row().is_ok(), fields == EvidenceFields::ALL);
            }
            assert!(batch.row(3).is_none());
        }
        let panel = EvidenceFields::DEPTHS.union(EvidenceFields::QUALITY_SUMS);
        assert!(
            panel.storage_bytes_per_locus() * 10 < EvidenceFields::ALL.storage_bytes_per_locus()
        );
    }

    #[test]
    fn owned_full_adapter_and_projected_append_preserve_nonzero_values() {
        let mut row = EvidenceRow {
            position: 5,
            reference: b'A',
            requested_alts: vec![b'T'],
            prefilter_depth: 3,
            aligned_depth: 2,
            callable_depth: 1,
            allele_counts: [0, 0, 0, 1],
            strand_counts: [[0, 0], [0, 0], [0, 0], [0, 1]],
            base_quality_sum: 30,
            mapping_quality_sum: 60,
            read_position_sum: 4,
            read_length_sum: 10,
            ..EvidenceRow::default()
        };
        row.filters.duplicate = 1;
        row.filters.low_base_quality = 1;
        row.base_quality_histogram[30] = 1;
        row.mapping_quality_histogram[60] = 1;
        let full = EvidenceBatch::from_full_rows(0, "chr1", 0, vec![row.clone()]);
        assert_eq!(full.row(0).unwrap().try_to_full_row().unwrap(), row);
        let mask = EvidenceFields::DEPTHS.union(EvidenceFields::QUALITY_SUMS);
        let mut projected = EvidenceBatch::new(0, "chr1", 0, mask, Vec::new());
        projected.reserve_exact(4);
        projected.push_row(full.row(0).unwrap()).unwrap();
        assert_eq!(projected.row(0).unwrap().depths.unwrap().callable_depth, 1);
        assert_eq!(
            projected
                .row(0)
                .unwrap()
                .quality_sums
                .unwrap()
                .mapping_quality_sum,
            60
        );
        assert!(projected.row(0).unwrap().quality_histograms.is_none());
        let mut destination = EvidenceBatch::new(0, "chr1", 0, EvidenceFields::ALL, Vec::new());
        assert!(destination.push_row(projected.row(0).unwrap()).is_err());
        assert!(destination.is_empty());
        projected.clear();
        assert!(projected.is_empty());
        assert_eq!(projected.fields(), mask);
        assert!(projected.quality_histograms.is_none());
        let mut borrowed = full.row(0).unwrap();
        borrowed.depths = None;
        assert!(borrowed.try_to_full_row().is_err());
    }
}
