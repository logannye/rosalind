//! Multi-contig genome index: a [`BlockedFMIndex`] over the concatenated
//! reference paired with the [`ContigSet`] that maps global offsets back to
//! `(contig, position)`. Exact-match queries return [`Locus`]es and reject any
//! hit that would straddle a contig boundary (the bwa "bns" rule). The FM-index
//! kernel is unchanged — multi-contig awareness lives entirely in the coordinate
//! mapping. Persistence + a zero-copy mmap view are Phase B3; wiring the
//! seed/chain/extend aligner onto this is Phase B4.

use std::sync::Arc;

use thiserror::Error;

use crate::core::{ContigSet, Locus};
use crate::genomics::{BlockedFMIndex, FMIndexError};

/// Largest concatenated-genome length supported. The suffix array stores
/// positions as `u32` over a text of `len + 1` symbols (one sentinel), so the
/// concatenation must satisfy `len + 1 <= u32::MAX`, i.e. `len <= u32::MAX - 1`
/// (~4.29 Gbp — covers the human genome and the vast majority of edge targets).
pub const MAX_GENOME_LEN: u64 = u32::MAX as u64 - 1;

/// Errors from building or querying a [`GenomeIndex`].
#[derive(Debug, Error)]
pub enum GenomeIndexError {
    /// The concatenated genome exceeds the `u32` suffix-array ceiling.
    #[error("concatenated genome length {len} exceeds the {max}-byte limit (u32 suffix array)")]
    GenomeTooLarge {
        /// Concatenated length in bytes.
        len: u64,
        /// Maximum supported length.
        max: u64,
    },
    /// The reference was empty (no contigs / no bases).
    #[error("genome index requires a non-empty reference")]
    EmptyGenome,
    /// Failure constructing the underlying FM-index (e.g. an unsupported base).
    #[error("fm-index error: {0}")]
    FmIndex(#[from] FMIndexError),
}

/// Validate that a concatenated-genome length fits the `u32` suffix array.
fn validate_concat_len(len: u64) -> Result<(), GenomeIndexError> {
    if len == 0 {
        return Err(GenomeIndexError::EmptyGenome);
    }
    if len > MAX_GENOME_LEN {
        return Err(GenomeIndexError::GenomeTooLarge {
            len,
            max: MAX_GENOME_LEN,
        });
    }
    Ok(())
}

/// An FM-index over a concatenated multi-contig reference plus the contig map.
#[derive(Debug)]
pub struct GenomeIndex {
    fm: BlockedFMIndex,
    contigs: ContigSet,
    reference: Arc<[u8]>,
}

impl GenomeIndex {
    /// Build an index over `reference` (the concatenation of every contig in
    /// `contigs`, in id order). The reference is normalised to uppercase on
    /// store so `reference()` agrees with the FM-index's internal (uppercased)
    /// keys; non-A/C/G/T/N bytes surface as [`FMIndexError`]. Errors if the
    /// concatenation is empty or exceeds [`MAX_GENOME_LEN`].
    pub fn build(contigs: ContigSet, reference: Arc<[u8]>) -> Result<Self, GenomeIndexError> {
        validate_concat_len(reference.len() as u64)?;
        // Normalise to uppercase on store (matching BWTAligner) so bytes from
        // `reference()` agree with the FM-index's uppercased keys — this avoids
        // a silent mismatch when the aligner extracts windows (Phase B4).
        let reference: Arc<[u8]> = Arc::from(
            reference
                .iter()
                .map(|b| b.to_ascii_uppercase())
                .collect::<Vec<u8>>()
                .into_boxed_slice(),
        );
        // Mirror BWTAligner's block-size heuristic.
        let block_size = ((reference.len() as f64).sqrt().ceil() as usize).max(64);
        let fm = BlockedFMIndex::build(&reference, block_size)?;
        Ok(Self {
            fm,
            contigs,
            reference,
        })
    }

    /// Build from per-contig `(name, sequence)` pairs: assembles the
    /// [`ContigSet`] and the concatenated reference (in the given order), then
    /// delegates to [`GenomeIndex::build`]. Sequences should be uppercase
    /// A/C/G/T/N.
    pub fn from_named_sequences(named: &[(String, Vec<u8>)]) -> Result<Self, GenomeIndexError> {
        let mut contigs = ContigSet::new();
        let total: usize = named.iter().map(|(_, seq)| seq.len()).sum();
        let mut concat = Vec::with_capacity(total);
        for (name, seq) in named {
            contigs.push(name.clone(), seq.len() as u32);
            concat.extend_from_slice(seq);
        }
        Self::build(contigs, Arc::from(concat.into_boxed_slice()))
    }

    /// The contig map for resolving `Locus`es.
    pub fn contigs(&self) -> &ContigSet {
        &self.contigs
    }

    /// The concatenated forward reference (uppercase A/C/G/T/N).
    pub fn reference(&self) -> &[u8] {
        &self.reference
    }

    /// The underlying FM-index (consumed by the aligner in Phase B4).
    pub fn fm(&self) -> &BlockedFMIndex {
        &self.fm
    }

    /// Locate exact occurrences of `pattern`, returning up to `max_hits`
    /// [`Locus`]es with boundary-straddling hits removed, sorted by
    /// `(contig, pos)` for a deterministic, kernel-independent order.
    ///
    /// Note: `max_hits` bounds the *candidate* suffixes located before boundary
    /// filtering, so the returned count may be smaller.
    pub fn locate_exact(&self, pattern: &[u8], max_hits: usize) -> Vec<Locus> {
        if pattern.is_empty() {
            return Vec::new();
        }
        let pattern: Vec<u8> = pattern.iter().map(|b| b.to_ascii_uppercase()).collect();
        let interval = self.fm.backward_search(&pattern);
        let len = pattern.len() as u32;
        let mut loci: Vec<Locus> = self
            .fm
            .locate_interval(interval, max_hits)
            .into_iter()
            .filter_map(|g| self.contigs.resolve_span(g as u64, len))
            .collect();
        loci.sort_unstable();
        loci
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::Position;

    // chr1 = ACGTACGT (global 0..8); chr2 = TTTTGGGG (global 8..16).
    // Concatenation: ACGTACGTTTTTGGGG
    fn two_contig_index() -> GenomeIndex {
        let mut contigs = ContigSet::new();
        contigs.push("chr1", 8);
        contigs.push("chr2", 8);
        let reference: Arc<[u8]> = Arc::from(b"ACGTACGTTTTTGGGG".to_vec().into_boxed_slice());
        GenomeIndex::build(contigs, reference).expect("build should succeed")
    }

    // Naive reference scan resolving each within-a-contig match to a Locus.
    fn naive_loci(reference: &[u8], pattern: &[u8], contigs: &ContigSet) -> Vec<Locus> {
        let mut out = Vec::new();
        if pattern.is_empty() || pattern.len() > reference.len() {
            return out;
        }
        for start in 0..=reference.len() - pattern.len() {
            if &reference[start..start + pattern.len()] == pattern {
                if let Some(locus) = contigs.resolve_span(start as u64, pattern.len() as u32) {
                    out.push(locus);
                }
            }
        }
        out.sort_unstable();
        out
    }

    #[test]
    fn stores_contigs_and_reference() {
        let idx = two_contig_index();
        assert_eq!(idx.contigs().len(), 2);
        assert_eq!(idx.contigs().by_name("chr2").unwrap().global_offset, 8);
        assert_eq!(idx.reference(), b"ACGTACGTTTTTGGGG");
    }

    #[test]
    fn locates_a_pattern_unique_to_chr2() {
        let idx = two_contig_index();
        assert_eq!(
            idx.locate_exact(b"GGGG", 16),
            vec![Locus {
                contig: 1,
                pos: Position(4)
            }]
        );
    }

    #[test]
    fn locates_a_pattern_unique_to_chr1() {
        let idx = two_contig_index();
        assert_eq!(
            idx.locate_exact(b"ACGTAC", 16),
            vec![Locus {
                contig: 0,
                pos: Position(0)
            }]
        );
    }

    #[test]
    fn rejects_a_boundary_straddling_match() {
        let idx = two_contig_index();
        // "GTTTTT" matches only at global 6 (chr1[6..8]="GT" + chr2[0..4]="TTTT"),
        // which crosses the chr1/chr2 boundary, so no Locus is returned.
        let loci = idx.locate_exact(b"GTTTTT", 16);
        assert!(
            loci.is_empty(),
            "boundary-straddling hit must be rejected, got {loci:?}"
        );
    }

    #[test]
    fn exact_match_matches_a_naive_scan() {
        let idx = two_contig_index();
        let reference = idx.reference().to_vec();
        for pattern in [
            b"ACGT".as_slice(),
            b"GT".as_slice(),
            b"TTTT".as_slice(),
            b"T".as_slice(),
            b"CG".as_slice(),
            b"GTTTTT".as_slice(),
        ] {
            let expected = naive_loci(&reference, pattern, idx.contigs());
            let got = idx.locate_exact(pattern, 1024);
            assert_eq!(
                got,
                expected,
                "mismatch for pattern {:?}",
                std::str::from_utf8(pattern).unwrap()
            );
        }
    }

    #[test]
    fn locate_exact_is_deterministic() {
        let idx = two_contig_index();
        assert_eq!(idx.locate_exact(b"T", 1024), idx.locate_exact(b"T", 1024));
    }

    #[test]
    fn empty_pattern_returns_nothing() {
        let idx = two_contig_index();
        assert!(idx.locate_exact(b"", 16).is_empty());
    }

    #[test]
    fn validate_concat_len_enforces_the_u32_ceiling() {
        assert!(matches!(
            validate_concat_len(0),
            Err(GenomeIndexError::EmptyGenome)
        ));
        assert!(validate_concat_len(1).is_ok());
        assert!(validate_concat_len(MAX_GENOME_LEN).is_ok());
        assert!(matches!(
            validate_concat_len(MAX_GENOME_LEN + 1),
            Err(GenomeIndexError::GenomeTooLarge { .. })
        ));
        assert!(matches!(
            validate_concat_len(u32::MAX as u64),
            Err(GenomeIndexError::GenomeTooLarge { .. })
        ));
    }

    #[test]
    fn max_hits_caps_candidates() {
        let idx = two_contig_index();
        // "ACGT" occurs at chr1:0 and chr1:4 (both within chr1).
        assert_eq!(idx.locate_exact(b"ACGT", 16).len(), 2);
        // Capping the candidate suffixes at 1 yields at most one locus.
        assert_eq!(idx.locate_exact(b"ACGT", 1).len(), 1);
    }

    #[test]
    fn from_named_sequences_builds_a_multi_contig_index() {
        let idx = GenomeIndex::from_named_sequences(&[
            ("chr1".to_string(), b"ACGTACGT".to_vec()),
            ("chr2".to_string(), b"TTTTGGGG".to_vec()),
            ("chr3".to_string(), b"CCCCAAAA".to_vec()),
        ])
        .expect("build should succeed");

        assert_eq!(idx.contigs().len(), 3);
        assert_eq!(idx.reference(), b"ACGTACGTTTTTGGGGCCCCAAAA");
        // "CCCC" is unique to chr3 at pos 0 (global 16).
        assert_eq!(
            idx.locate_exact(b"CCCC", 16),
            vec![Locus {
                contig: 2,
                pos: Position(0)
            }]
        );
    }

    #[test]
    fn from_named_sequences_rejects_an_empty_genome() {
        assert!(matches!(
            GenomeIndex::from_named_sequences(&[]),
            Err(GenomeIndexError::EmptyGenome)
        ));
    }
}
