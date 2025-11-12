use std::sync::Arc;

/// Simple CIGAR operation kinds describing how a read aligns to the reference.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CigarOpKind {
    /// Consuming match/mismatch.
    Match,
    /// Insertion relative to the reference.
    Insertion,
    /// Deletion relative to the reference.
    Deletion,
    /// Soft clipping (sequence present in read only).
    SoftClip,
    /// Hard clipping (trimmed sequence not present in read).
    HardClip,
}

/// CIGAR operation with length.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CigarOp {
    /// Operation kind.
    pub kind: CigarOpKind,
    /// Number of bases affected by the operation.
    pub len: u32,
}

impl CigarOp {
    /// Construct a new CIGAR operation.
    pub fn new(kind: CigarOpKind, len: u32) -> Self {
        Self { kind, len }
    }
}

/// Aligned read with sequence and quality information.
#[derive(Debug, Clone)]
pub struct AlignedRead {
    /// Reference contig/chromosome name.
    pub chrom: Arc<str>,
    /// 0-based leftmost reference coordinate.
    pub pos: u32,
    /// Mapping quality (Phred-scaled).
    pub mapq: u8,
    /// CIGAR describing the alignment.
    pub cigar: Vec<CigarOp>,
    /// Read sequence stored as uppercase ASCII.
    pub sequence: Arc<[u8]>,
    /// Per-base quality scores in Phred space.
    pub qualities: Arc<[u8]>,
    /// Whether the read maps to the reverse complement strand.
    pub is_reverse: bool,
}

impl AlignedRead {
    /// Construct a new aligned read wrapper.
    pub fn new(
        chrom: impl Into<Arc<str>>,
        pos: u32,
        mapq: u8,
        cigar: Vec<CigarOp>,
        sequence: impl Into<Arc<[u8]>>,
        qualities: impl Into<Arc<[u8]>>,
        is_reverse: bool,
    ) -> Self {
        Self {
            chrom: chrom.into(),
            pos,
            mapq,
            cigar,
            sequence: sequence.into(),
            qualities: qualities.into(),
            is_reverse,
        }
    }

    /// Read length inferred from the sequence.
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    /// End position (half-open) on the reference assuming contiguous match.
    pub fn end(&self) -> u32 {
        self.pos + self.len() as u32
    }

    /// Base at the provided read offset.
    pub fn base_at(&self, offset: usize) -> Option<u8> {
        self.sequence.get(offset).copied()
    }

    /// Quality score at the provided read offset.
    pub fn quality_at(&self, offset: usize) -> Option<u8> {
        self.qualities.get(offset).copied()
    }

    /// Mapping quality associated with the alignment.
    pub fn mapq(&self) -> u8 {
        self.mapq
    }
}
