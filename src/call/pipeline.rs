//! Orchestration: drive a `ReadSource` through the `PileupEngine` and call a
//! germline genotype at every covered position, collecting the variant sites
//! (hom-ref columns are abstained on by `call_germline` and never appear).

use std::ops::Range;
use std::sync::Arc;

use crate::call::{call_germline, GermlineCall, GermlineParams};
use crate::core::{CoreError, Locus};
use crate::pileup::{PileupEngine, PileupParams, ReadSource};

/// Call germline variants across `region` of `contig`. Returns one entry per
/// emitted site as `(locus, ref_base, call)`; the caller pairs these into VCF
/// rows. Hom-ref / no-evidence positions are abstained on (absent from output).
pub fn call_germline_region<S: ReadSource>(
    source: S,
    reference: Arc<[u8]>,
    contig: u32,
    region: Range<u32>,
    pileup_params: PileupParams,
    germline_params: &GermlineParams,
) -> Result<Vec<(Locus, u8, GermlineCall)>, CoreError> {
    let engine = PileupEngine::new(source, reference, contig, region, pileup_params);
    let mut out = Vec::new();
    for column in engine {
        let column = column?;
        if let Some(call) = call_germline(&column, germline_params) {
            out.push((column.locus, column.ref_base, call));
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::call::Genotype;
    use crate::core::{AlignedRead, CigarOp, CigarOpKind, Position, SamFlags};
    use crate::pileup::SliceSource;

    fn read(pos: u32, seq: &[u8], reverse: bool) -> AlignedRead {
        let flags = if reverse {
            SamFlags(SamFlags::REVERSE)
        } else {
            SamFlags::default()
        };
        AlignedRead {
            contig: 0,
            pos: Position(pos),
            mapq: 60,
            flags,
            cigar: vec![CigarOp::new(CigarOpKind::Match, seq.len() as u32)],
            seq: Arc::from(seq.to_vec().into_boxed_slice()),
            qual: Arc::from(vec![35u8; seq.len()].into_boxed_slice()),
        }
    }

    #[test]
    fn calls_a_het_snv_including_a_reverse_strand_read() {
        // Reference AAAA over [0,4). At position 1, reads split A (ref) / C (alt).
        // One alt read is reverse-strand — proving the end-to-end reverse-strand
        // fix (it must contribute a C, not a complemented G).
        let reference: Arc<[u8]> = Arc::from(b"AAAA".to_vec().into_boxed_slice());
        let reads = vec![
            read(0, b"AAAA", false), // all ref
            read(0, b"ACAA", false), // alt C at pos 1 (forward)
            read(0, b"ACAA", true),  // alt C at pos 1 (reverse — SEQ already forward)
            read(0, b"ACAA", false),
            read(0, b"AAAA", false),
            read(0, b"ACAA", true),
        ];
        let source = SliceSource::new(reads);
        let calls = call_germline_region(
            source,
            reference,
            0,
            0..4,
            PileupParams::default(),
            &GermlineParams::default(),
        )
        .unwrap();

        // Exactly one variant site, at position 1, alt C, heterozygous.
        assert_eq!(calls.len(), 1);
        let (locus, ref_base, call) = &calls[0];
        assert_eq!(locus.pos, Position(1));
        assert_eq!(*ref_base, b'A');
        assert_eq!(call.alt_base, b'C');
        assert_eq!(call.genotype, Genotype::Het);
        assert!(call.dp >= 6);
    }

    #[test]
    fn pure_reference_yields_no_calls() {
        let reference: Arc<[u8]> = Arc::from(b"ACGT".to_vec().into_boxed_slice());
        let reads = vec![read(0, b"ACGT", false), read(0, b"ACGT", false)];
        let source = SliceSource::new(reads);
        let calls = call_germline_region(
            source,
            reference,
            0,
            0..4,
            PileupParams::default(),
            &GermlineParams::default(),
        )
        .unwrap();
        assert!(calls.is_empty());
    }
}
