use super::EvidenceError;
use crate::core::ContigSet;
use crate::selection::{GenomicInterval, IntervalSet};
use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

type ResolvedSelection = (Vec<GenomicInterval>, BTreeMap<(u32, u32), SnvSite>);

/// One requested SNV locus; multiple ALT alleles are represented once.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SnvSite {
    /// Reference contig identifier or canonical contig name for a batch.
    pub contig: u32,
    /// Zero-based reference coordinate within the contig.
    pub position: u32,
    /// Normalized reference base at this locus, or N when unavailable.
    pub reference: u8,
    /// Requested one-base ALT alleles, normalized and deduplicated on resolution.
    pub alternates: Vec<u8>,
}
#[derive(Debug, Clone, Default)]
/// Requested genomic loci; normalization preserves annotations and removes duplicate output positions.
pub enum EvidenceSelection {
    #[default]
    /// Emit every position of the ordered reference dictionary.
    WholeGenome,
    /// Emit the normalized union of half-open intervals, including zero-depth loci.
    Intervals(Vec<GenomicInterval>),
    /// Emit unique requested SNV loci after validating their REF bases.
    Sites(Vec<SnvSite>),
}
impl EvidenceSelection {
    /// Parse local BED coordinates against the supplied reference dictionary.
    pub fn from_bed(path: impl AsRef<Path>, contigs: &ContigSet) -> Result<Self, EvidenceError> {
        let set = IntervalSet::from_bed(path, contigs)
            .map_err(|error| EvidenceError::InvalidRequest(error.to_string()))?;
        Ok(Self::Intervals(set.intervals().to_vec()))
    }
    /// Parse VCF, compressed VCF, or BCF using the default record envelopes.
    /// Duplicate loci union ALT alleles and must agree on REF. No index is needed.
    pub fn from_vcf(path: impl AsRef<Path>, contigs: &ContigSet) -> Result<Self, EvidenceError> {
        Self::from_variants(path, contigs, crate::variant_io::VariantLimits::default())
    }

    /// Parse typed SNV records with explicit cooperative header/record envelopes.
    /// The normalized selection is independent of record order and compression;
    /// use [`crate::variant_io::VariantReader`] when original records are needed.
    pub fn from_variants(
        path: impl AsRef<Path>,
        contigs: &ContigSet,
        limits: crate::variant_io::VariantLimits,
    ) -> Result<Self, EvidenceError> {
        let mut reader = crate::variant_io::VariantReader::open(path, limits)?;
        let mut sites: BTreeMap<(u32, u32), (u8, BTreeSet<u8>)> = BTreeMap::new();
        while let Some(record) = reader.read()? {
            let site = crate::variant_io::parse_snv_record(record, contigs)?;
            let entry = sites
                .entry((site.contig, site.position))
                .or_insert((site.reference, BTreeSet::new()));
            if entry.0 != site.reference {
                return Err(EvidenceError::InvalidInput(
                    "duplicate variant locus has conflicting REF".into(),
                ));
            }
            entry.1.extend(site.alternates);
        }
        Ok(Self::Sites(
            sites
                .into_iter()
                .map(|((contig, position), (reference, alternates))| SnvSite {
                    contig,
                    position,
                    reference,
                    alternates: alternates.into_iter().collect(),
                })
                .collect(),
        ))
    }
    pub(crate) fn normalize(
        &self,
        contigs: &ContigSet,
    ) -> Result<ResolvedSelection, EvidenceError> {
        let mut sites = BTreeMap::new();
        let intervals = match self {
            Self::WholeGenome => contigs
                .iter()
                .map(|c| GenomicInterval {
                    contig: c.id,
                    start: 0,
                    end: c.length,
                })
                .collect(),
            Self::Intervals(intervals) => intervals.clone(),
            Self::Sites(input) => {
                for site in input {
                    if !b"ACGT".contains(&site.reference)
                        || site.alternates.is_empty()
                        || site
                            .alternates
                            .iter()
                            .any(|b| !b"ACGT".contains(b) || *b == site.reference)
                    {
                        return Err(EvidenceError::InvalidRequest(
                            "site requires A/C/G/T SNV REF/ALT".into(),
                        ));
                    }
                    let existing = sites
                        .entry((site.contig, site.position))
                        .or_insert_with(|| site.clone());
                    if existing.reference != site.reference {
                        return Err(EvidenceError::InvalidRequest(
                            "conflicting REF at duplicated site".into(),
                        ));
                    }
                    existing.alternates.extend_from_slice(&site.alternates);
                    existing.alternates.sort_unstable();
                    existing.alternates.dedup();
                }
                sites
                    .values()
                    .map(|site| {
                        Ok(GenomicInterval {
                            contig: site.contig,
                            start: site.position,
                            end: site
                                .position
                                .checked_add(1)
                                .ok_or(EvidenceError::CounterOverflow)?,
                        })
                    })
                    .collect::<Result<Vec<_>, EvidenceError>>()?
            }
        };
        let normalized = IntervalSet::new(intervals, contigs)
            .map_err(|error| EvidenceError::InvalidRequest(error.to_string()))?;
        Ok((normalized.intervals().to_vec(), sites))
    }
}
