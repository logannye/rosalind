use crate::genomics::CompressedDNA;

/// Number of canonical DNA symbols tracked in rank/select (A, C, G, T, N).
pub const ALPHABET_SIZE: usize = 5;
/// Default number of bases per superblock in the bitvector rank structure.
///
/// Larger values reduce prefix table size; smaller values reduce the per-query scan.
pub const CHECKPOINT_STRIDE: usize = 1024;

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
///
/// This is a legacy-compatible type that describes the prefix counts at a given
/// boundary. The current implementation uses superblocks and popcount-based rank
/// queries instead of scanning raw bases.
#[derive(Debug, Clone)]
pub struct RankSelectCheckpoint {
    /// Starting position (0-indexed) covered by this checkpoint.
    pub position: usize,
    /// Prefix counts for [A, C, G, T, N] at the start of the checkpoint.
    pub counts: [u32; ALPHABET_SIZE],
}

impl RankSelectCheckpoint {
    #[allow(dead_code)]
    fn new(position: usize, counts: [u32; ALPHABET_SIZE]) -> Self {
        Self { position, counts }
    }
}

/// Rank/select index built over a [`CompressedDNA`] sequence.
///
/// This implementation stores a bitvector per base (A/C/G/T/N) and answers rank
/// queries using a two-level structure:
/// - **superblock prefix sums** (every `stride` bases)
/// - **popcount** within the superblock using masked `u64` words
#[derive(Debug, Clone)]
pub struct RankSelectIndex {
    stride: usize, // bases per superblock
    bitvectors: [Vec<u64>; ALPHABET_SIZE],
    superblocks: [Vec<u32>; ALPHABET_SIZE], // prefix counts per superblock boundary
    totals: [u32; ALPHABET_SIZE],
    len: usize,
}

impl RankSelectIndex {
    /// Construct an index with the default stride.
    pub fn build(sequence: &CompressedDNA) -> Self {
        Self::build_with_stride(sequence, CHECKPOINT_STRIDE)
    }

    /// Construct an index with the provided stride.
    pub fn build_with_stride(sequence: &CompressedDNA, stride: usize) -> Self {
        let stride = stride.max(1);
        let len = sequence.len();
        let word_len = words_for_bits(len);

        let mut bitvectors: [Vec<u64>; ALPHABET_SIZE] =
            std::array::from_fn(|_| vec![0u64; word_len]);
        let mut totals = [0u32; ALPHABET_SIZE];

        for idx in 0..len {
            let base = sequence.base_at(idx).unwrap_or(b'N');
            let symbol = BaseCode::from_ascii(base).unwrap_or(BaseCode::N);
            totals[symbol.index()] += 1;

            let (word_idx, bit_mask) = bit_position(idx);
            bitvectors[symbol.index()][word_idx] |= bit_mask;
        }

        // Build superblock prefix sums for each base.
        let num_superblocks = (len + stride - 1) / stride;
        let mut superblocks: [Vec<u32>; ALPHABET_SIZE] =
            std::array::from_fn(|_| Vec::with_capacity(num_superblocks + 1));
        for b in 0..ALPHABET_SIZE {
            superblocks[b].push(0);
        }

        for sb in 0..num_superblocks {
            let start = sb * stride;
            let end = (start + stride).min(len);
            for b in 0..ALPHABET_SIZE {
                let prev = *superblocks[b].last().unwrap();
                let add = popcount_range(&bitvectors[b], start, end);
                superblocks[b].push(prev + add);
            }
        }

        Self {
            stride,
            bitvectors,
            superblocks,
            totals,
            len,
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
        debug_assert_eq!(
            sequence.len(),
            self.len,
            "RankSelectIndex must be queried with the same sequence it was built for"
        );

        let bounded = position.min(self.len);
        let sb = bounded / self.stride;
        let within_start = sb * self.stride;

        let base_idx = base.index();
        let prefix = self.superblocks[base_idx][sb];
        let within = popcount_range(&self.bitvectors[base_idx], within_start, bounded);
        prefix + within
    }

    /// Rank query returning counts for all bases in one pass.
    pub fn rank_all(&self, sequence: &CompressedDNA, position: usize) -> [u32; ALPHABET_SIZE] {
        debug_assert_eq!(
            sequence.len(),
            self.len,
            "RankSelectIndex must be queried with the same sequence it was built for"
        );

        let bounded = position.min(self.len);
        let sb = bounded / self.stride;
        let within_start = sb * self.stride;

        let mut out = [0u32; ALPHABET_SIZE];
        for b in 0..ALPHABET_SIZE {
            let prefix = self.superblocks[b][sb];
            let within = popcount_range(&self.bitvectors[b], within_start, bounded);
            out[b] = prefix + within;
        }
        out
    }
}

#[inline]
fn words_for_bits(bits: usize) -> usize {
    if bits == 0 {
        0
    } else {
        (bits + 63) / 64
    }
}

#[inline]
fn bit_position(idx: usize) -> (usize, u64) {
    let word_idx = idx / 64;
    let bit_idx = idx % 64;
    (word_idx, 1u64 << bit_idx)
}

#[inline]
fn popcount_range(words: &[u64], start: usize, end: usize) -> u32 {
    if end <= start {
        return 0;
    }
    let start_word = start / 64;
    let end_word = (end - 1) / 64;
    let start_bit = start % 64;
    let end_bit = end % 64;

    let mut count = 0u32;

    if start_word == end_word {
        let mut w = words[start_word];
        // Mask off bits below start.
        w &= !((1u64 << start_bit) - 1);
        // Mask off bits at/above end.
        if end_bit != 0 {
            w &= (1u64 << end_bit) - 1;
        }
        count += w.count_ones();
        return count;
    }

    // First partial word.
    let first = words[start_word] & !((1u64 << start_bit) - 1);
    count += first.count_ones();

    // Full words.
    for w in &words[start_word + 1..end_word] {
        count += w.count_ones();
    }

    // Last partial word.
    let mut last = words[end_word];
    if end_bit != 0 {
        last &= (1u64 << end_bit) - 1;
    }
    count += last.count_ones();

    count
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
