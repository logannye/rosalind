//! The canonical aligned-read record and CIGAR projection.
//!
//! SEQ is stored forward-reference-oriented (as in SAM/BAM); strand is
//! metadata and is never applied to the bytes. Records are read-length-
//! agnostic — short Illumina reads and long Nanopore/PacBio reads are handled
//! identically.

use std::sync::Arc;

use crate::core::locus::Position;

/// CIGAR operation kind (SAM). Consuming semantics follow the SAM spec.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CigarOpKind {
    /// Alignment match or mismatch (`M`/`=`/`X`): consumes ref and read.
    Match,
    /// Insertion to the reference (`I`): consumes read only.
    Insertion,
    /// Deletion from the reference (`D`): consumes ref only.
    Deletion,
    /// Skipped reference region (`N`, e.g. an intron): consumes ref only.
    RefSkip,
    /// Soft clip (`S`): read bases present but unaligned; consumes read only.
    SoftClip,
    /// Hard clip (`H`): trimmed bases absent from the read; consumes neither.
    HardClip,
}

impl CigarOpKind {
    /// Whether this op consumes reference bases.
    pub fn consumes_ref(self) -> bool {
        matches!(self, CigarOpKind::Match | CigarOpKind::Deletion | CigarOpKind::RefSkip)
    }

    /// Whether this op consumes read/query bases.
    pub fn consumes_read(self) -> bool {
        matches!(self, CigarOpKind::Match | CigarOpKind::Insertion | CigarOpKind::SoftClip)
    }
}

/// A single CIGAR operation with its run length.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CigarOp {
    /// The operation kind.
    pub kind: CigarOpKind,
    /// Number of bases the operation spans.
    pub len: u32,
}

impl CigarOp {
    /// Construct a CIGAR operation.
    pub fn new(kind: CigarOpKind, len: u32) -> Self {
        Self { kind, len }
    }
}

/// SAM flag bitset (the subset Rosalind uses).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct SamFlags(pub u16);

impl SamFlags {
    /// Template has multiple segments (paired).
    pub const PAIRED: u16 = 0x1;
    /// Each segment properly aligned (proper pair).
    pub const PROPER_PAIR: u16 = 0x2;
    /// Segment unmapped.
    pub const UNMAPPED: u16 = 0x4;
    /// Segment reverse-strand.
    pub const REVERSE: u16 = 0x10;
    /// Secondary alignment.
    pub const SECONDARY: u16 = 0x100;
    /// PCR or optical duplicate.
    pub const DUPLICATE: u16 = 0x400;
    /// Supplementary alignment.
    pub const SUPPLEMENTARY: u16 = 0x800;

    /// Whether the given flag bit is set.
    pub fn contains(self, bit: u16) -> bool {
        self.0 & bit != 0
    }

    /// Segment is unmapped.
    pub fn is_unmapped(self) -> bool {
        self.contains(Self::UNMAPPED)
    }

    /// Segment maps to the reverse strand (metadata only — SEQ stays forward).
    pub fn is_reverse(self) -> bool {
        self.contains(Self::REVERSE)
    }

    /// Alignment is secondary.
    pub fn is_secondary(self) -> bool {
        self.contains(Self::SECONDARY)
    }

    /// Alignment is supplementary.
    pub fn is_supplementary(self) -> bool {
        self.contains(Self::SUPPLEMENTARY)
    }

    /// Read is a PCR/optical duplicate.
    pub fn is_duplicate(self) -> bool {
        self.contains(Self::DUPLICATE)
    }
}

/// One projected base: a read base aligned to a specific reference position.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RefBase {
    /// 0-based reference position.
    pub ref_pos: u32,
    /// Offset of the base within the read's `seq`/`qual`.
    pub read_offset: usize,
}

/// A read aligned to a single contig.
///
/// `seq` is uppercase ASCII in forward-reference orientation; `pos` is the
/// 0-based leftmost reference coordinate. Read-length-agnostic.
#[derive(Debug, Clone)]
pub struct AlignedRead {
    /// Contig id (index into the originating `ContigSet`).
    pub contig: u32,
    /// 0-based leftmost reference coordinate.
    pub pos: Position,
    /// Mapping quality (Phred-scaled).
    pub mapq: u8,
    /// SAM flags.
    pub flags: SamFlags,
    /// CIGAR describing the alignment.
    pub cigar: Vec<CigarOp>,
    /// Read sequence, uppercase ASCII, forward orientation.
    pub seq: Arc<[u8]>,
    /// Per-base Phred qualities (parallel to `seq`).
    pub qual: Arc<[u8]>,
}

impl AlignedRead {
    /// Reference bases spanned by the alignment (sum of ref-consuming ops).
    pub fn ref_span(&self) -> u32 {
        self.cigar.iter().filter(|o| o.kind.consumes_ref()).map(|o| o.len).sum()
    }

    /// Half-open reference end coordinate. CIGAR-derived — never `pos + seq_len`.
    pub fn end(&self) -> u32 {
        self.pos.0 + self.ref_span()
    }

    /// Project each `Match` base to its reference position by walking the CIGAR.
    /// Insertions/soft-clips consume read only; deletions/ref-skips consume ref
    /// only; hard-clips consume neither. This is the single CIGAR projection the
    /// pileup kernel consumes.
    pub fn projected_bases(&self) -> Vec<RefBase> {
        let mut out = Vec::new();
        let mut ref_pos = self.pos.0;
        let mut read_off: usize = 0;
        for op in &self.cigar {
            match op.kind {
                CigarOpKind::Match => {
                    for _ in 0..op.len {
                        out.push(RefBase { ref_pos, read_offset: read_off });
                        ref_pos += 1;
                        read_off += 1;
                    }
                }
                CigarOpKind::Insertion | CigarOpKind::SoftClip => {
                    read_off += op.len as usize;
                }
                CigarOpKind::Deletion | CigarOpKind::RefSkip => {
                    ref_pos += op.len;
                }
                CigarOpKind::HardClip => {}
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn op(kind: CigarOpKind, len: u32) -> CigarOp {
        CigarOp::new(kind, len)
    }

    fn read(pos: u32, cigar: Vec<CigarOp>, seq: &[u8]) -> AlignedRead {
        AlignedRead {
            contig: 0,
            pos: Position(pos),
            mapq: 60,
            flags: SamFlags::default(),
            cigar,
            seq: Arc::from(seq.to_vec().into_boxed_slice()),
            qual: Arc::from(vec![30u8; seq.len()].into_boxed_slice()),
        }
    }

    #[test]
    fn ref_span_and_end_are_cigar_derived_not_seq_len() {
        // 3M1D2M consumes 6 reference bases from 5 read bases.
        let r = read(100, vec![op(CigarOpKind::Match, 3), op(CigarOpKind::Deletion, 1), op(CigarOpKind::Match, 2)], b"ACGTT");
        assert_eq!(r.ref_span(), 6);
        assert_eq!(r.end(), 106);
    }

    #[test]
    fn projection_handles_softclip_insertion_deletion_refskip() {
        // 2S 3M 1I 2M 1D 2M; ref starts at 50.
        // read offsets [0,1]=softclip, [2,3,4]=M, [5]=I, [6,7]=M, (1D no read), [8,9]=M
        let r = read(
            50,
            vec![
                op(CigarOpKind::SoftClip, 2),
                op(CigarOpKind::Match, 3),
                op(CigarOpKind::Insertion, 1),
                op(CigarOpKind::Match, 2),
                op(CigarOpKind::Deletion, 1),
                op(CigarOpKind::Match, 2),
            ],
            b"NNACGTTGAA", // 10 bases: 2 clipped + 3 + 1 ins + 2 + 2
        );
        let proj = r.projected_bases();
        let pairs: Vec<(u32, usize)> = proj.iter().map(|p| (p.ref_pos, p.read_offset)).collect();
        assert_eq!(
            pairs,
            vec![
                (50, 2), (51, 3), (52, 4), // first 3M
                (53, 6), (54, 7),          // 2M after the 1I (read offset skips 5)
                (56, 8), (57, 9),          // 2M after the 1D (ref skips 55)
            ]
        );
        // Soft-clipped and inserted read bases never appear in the projection.
        assert!(!pairs.iter().any(|(_, off)| *off == 0 || *off == 1 || *off == 5));
    }

    #[test]
    fn refskip_consumes_reference_only_like_a_long_intron() {
        // 2M 100N 2M (spliced long read): ref advances over the skip, no bases there.
        let r = read(0, vec![op(CigarOpKind::Match, 2), op(CigarOpKind::RefSkip, 100), op(CigarOpKind::Match, 2)], b"ACGT");
        let proj = r.projected_bases();
        assert_eq!(proj.first().unwrap().ref_pos, 0);
        assert_eq!(proj.last().unwrap().ref_pos, 103);
        assert_eq!(proj.len(), 4);
        assert_eq!(r.ref_span(), 104);
    }

    #[test]
    fn long_read_projects_every_match_base() {
        // Read-length-agnostic: a 5000-base full-match "long read".
        let seq = vec![b'A'; 5000];
        let r = read(1000, vec![op(CigarOpKind::Match, 5000)], &seq);
        let proj = r.projected_bases();
        assert_eq!(proj.len(), 5000);
        assert_eq!(proj[0], RefBase { ref_pos: 1000, read_offset: 0 });
        assert_eq!(proj[4999], RefBase { ref_pos: 5999, read_offset: 4999 });
    }

    #[test]
    fn sam_flags_decode_common_bits() {
        let f = SamFlags(SamFlags::REVERSE | SamFlags::DUPLICATE);
        assert!(f.is_reverse());
        assert!(f.is_duplicate());
        assert!(!f.is_secondary());
        assert!(!f.is_unmapped());
    }
}
