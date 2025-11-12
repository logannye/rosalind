use crate::genomics::CompressedDNA;

/// Number of canonical DNA symbols tracked in rank/select (A, C, G, T, N).
pub const ALPHABET_SIZE: usize = 5;
/// Default number of bases between checkpoints.
pub const CHECKPOINT_STRIDE: usize = 256;

/// Enumeration representing base codes used for rank/select queries.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BaseCode {
    /// Adenine.
    A = 0,
    /// Cytosine.
    C = 1,
    /// Guanine.
    G = 2,
    /// Thymine/Uracil.
    T = 3,
    /// Ambiguous base (`N`).
    N = 4,
}

impl BaseCode {
    /// Attempt to parse an ASCII base into a [`BaseCode`].
    pub fn from_ascii(base: u8) -> Option<Self> {
        match base {
            b'A' | b'a' => Some(BaseCode::A),
            b'C' | b'c' => Some(BaseCode::C),
            b'G' | b'g' => Some(BaseCode::G),
            b'T' | b't' | b'U' | b'u' => Some(BaseCode::T),
            b'N' | b'n' => Some(BaseCode::N),
            _ => None,
        }
    }

    /// Convert the base code to an index into rank/select tables.
    #[inline]
    pub fn index(self) -> usize {
        self as usize
    }
}

/// Prefix-sum checkpoint for rank queries.
#[derive(Debug, Clone)]
pub struct RankSelectCheckpoint {
    /// Starting position (0-indexed) covered by this checkpoint.
    pub position: usize,
    /// Prefix counts for [A, C, G, T, N] at the start of the checkpoint.
    pub counts: [u32; ALPHABET_SIZE],
}

impl RankSelectCheckpoint {
    fn new(position: usize, counts: [u32; ALPHABET_SIZE]) -> Self {
        Self { position, counts }
    }
}

/// Rank/select index built over a [`CompressedDNA`] sequence.
#[derive(Debug, Clone)]
pub struct RankSelectIndex {
    stride: usize,
    checkpoints: Vec<RankSelectCheckpoint>,
    totals: [u32; ALPHABET_SIZE],
}

impl RankSelectIndex {
    /// Construct an index with the default stride.
    pub fn build(sequence: &CompressedDNA) -> Self {
        Self::build_with_stride(sequence, CHECKPOINT_STRIDE)
    }

    /// Construct an index with the provided stride.
    pub fn build_with_stride(sequence: &CompressedDNA, stride: usize) -> Self {
        assert!(stride > 0, "stride must be greater than zero");

        let mut checkpoints = Vec::new();
        let mut counts = [0u32; ALPHABET_SIZE];
        checkpoints.push(RankSelectCheckpoint::new(0, counts));

        for (idx, base) in sequence.iter().enumerate() {
            if idx % stride == 0 && idx != 0 {
                checkpoints.push(RankSelectCheckpoint::new(idx, counts));
            }
            let symbol = BaseCode::from_ascii(base).unwrap_or(BaseCode::N);
            counts[symbol.index()] += 1;
        }

        // Add terminal checkpoint for completeness.
        checkpoints.push(RankSelectCheckpoint::new(sequence.len(), counts));

        Self {
            stride,
            checkpoints,
            totals: counts,
        }
    }

    /// Number of bases between checkpoints.
    pub fn stride(&self) -> usize {
        self.stride
    }

    /// Total cumulative counts for each symbol.
    pub fn totals(&self) -> [u32; ALPHABET_SIZE] {
        self.totals
    }

    /// Rank query: count of `base` in `sequence[..position)`.
    pub fn rank(&self, sequence: &CompressedDNA, base: BaseCode, position: usize) -> u32 {
        let bounded = position.min(sequence.len());
        let checkpoint_idx = bounded / self.stride;
        let remainder_start = checkpoint_idx * self.stride;

        let checkpoint = &self.checkpoints[checkpoint_idx];
        let mut count = checkpoint.counts[base.index()];

        for idx in remainder_start..bounded {
            let symbol =
                BaseCode::from_ascii(sequence.base_at(idx).unwrap_or(b'N')).unwrap_or(BaseCode::N);
            if symbol == base {
                count += 1;
            }
        }

        count
    }

    /// Rank query returning counts for all bases in one pass.
    pub fn rank_all(&self, sequence: &CompressedDNA, position: usize) -> [u32; ALPHABET_SIZE] {
        let bounded = position.min(sequence.len());
        let checkpoint_idx = bounded / self.stride;
        let remainder_start = checkpoint_idx * self.stride;

        let mut counts = self.checkpoints[checkpoint_idx].counts;
        for idx in remainder_start..bounded {
            let symbol =
                BaseCode::from_ascii(sequence.base_at(idx).unwrap_or(b'N')).unwrap_or(BaseCode::N);
            counts[symbol.index()] += 1;
        }
        counts
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rank_queries_match_naive_counts() {
        let seq = b"AAACCCGGGTTTNNNAAGT";
        let compressed = CompressedDNA::compress(seq).unwrap();
        let index = RankSelectIndex::build_with_stride(&compressed, 4);

        // Compare with naive counting.
        for pos in 0..=seq.len() {
            for &base in &[
                BaseCode::A,
                BaseCode::C,
                BaseCode::G,
                BaseCode::T,
                BaseCode::N,
            ] {
                let naive = seq[..pos]
                    .iter()
                    .filter(|&&b| BaseCode::from_ascii(b).unwrap_or(BaseCode::N) == base)
                    .count() as u32;
                assert_eq!(index.rank(&compressed, base, pos), naive);
            }
        }
    }

    #[test]
    fn rank_all_returns_expected_counts() {
        let seq = b"ATCGNNATCG";
        let compressed = CompressedDNA::compress(seq).unwrap();
        let index = RankSelectIndex::build(&compressed);

        for pos in 0..=seq.len() {
            let counts = index.rank_all(&compressed, pos);
            let naive = [
                seq[..pos]
                    .iter()
                    .filter(|&&b| b == b'A' || b == b'a')
                    .count() as u32,
                seq[..pos]
                    .iter()
                    .filter(|&&b| b == b'C' || b == b'c')
                    .count() as u32,
                seq[..pos]
                    .iter()
                    .filter(|&&b| b == b'G' || b == b'g')
                    .count() as u32,
                seq[..pos]
                    .iter()
                    .filter(|&&b| b == b'T' || b == b't' || b == b'U' || b == b'u')
                    .count() as u32,
                seq[..pos]
                    .iter()
                    .filter(|&&b| b == b'N' || b == b'n')
                    .count() as u32,
            ];
            assert_eq!(counts, naive);
        }
    }

    #[test]
    fn totals_match_full_sequence() {
        let seq = b"AACCGGTTNN";
        let compressed = CompressedDNA::compress(seq).unwrap();
        let index = RankSelectIndex::build(&compressed);
        assert_eq!(
            index.totals(),
            [
                2, // A
                2, // C
                2, // G
                2, // T
                2, // N
            ]
        );
    }
}
