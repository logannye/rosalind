//! Compact sampled suffix array for the FM-index.
//!
//! Replaces a dense `Vec<u32>` (one slot per BWT position, ~12 GB for a human
//! genome) with a sparse representation: a sampled-position **bitvector**
//! (`marks`), superblock prefix-popcounts over it for O(1) rank, and a **packed**
//! array of the sampled SA values in ascending BWT-position order. Storage is
//! `O(bwt_len bits + (bwt_len / rate) values)`. This is the structure the on-disk
//! index serialises (Phase B3b).

/// BWT positions per superblock in the rank structure over `marks`.
const RANK_STRIDE: usize = 1024;

/// A compact sampled suffix array.
#[derive(Debug, Clone)]
pub struct SampledSuffixArray {
    rate: usize,
    bwt_len: usize,
    /// Bit `i` set iff BWT position `i` carries a sampled SA value.
    marks: Vec<u64>,
    /// Prefix popcount of `marks` at each `RANK_STRIDE` boundary (`len/stride + 1`).
    superblocks: Vec<u32>,
    /// Sampled SA values, in ascending BWT-position order.
    values: Vec<u32>,
}

impl SampledSuffixArray {
    /// Build from `(bwt_position, sa_value)` pairs of the sampled positions,
    /// which MUST be yielded in **strictly ascending** `bwt_position` order.
    pub fn from_sorted_samples(
        bwt_len: usize,
        rate: usize,
        samples: impl Iterator<Item = (usize, u32)>,
    ) -> Self {
        // `(n + 63) / 64`, not `n.div_ceil(64)` (which is Rust 1.73+; MSRV is 1.72).
        let words = (bwt_len + 63) / 64;
        let mut marks = vec![0u64; words];
        let mut values = Vec::new();
        // The strictly-ascending contract is load-bearing: `values` is pushed in
        // iteration order and indexed by rank, so out-of-order input would
        // silently return the wrong SA value. Guard it in debug builds.
        #[cfg(debug_assertions)]
        let mut prev: Option<usize> = None;
        for (pos, val) in samples {
            debug_assert!(pos < bwt_len, "sample position out of range");
            #[cfg(debug_assertions)]
            {
                debug_assert!(
                    prev.map_or(true, |p| pos > p),
                    "samples must be yielded in strictly ascending BWT-position order"
                );
                prev = Some(pos);
            }
            marks[pos / 64] |= 1u64 << (pos % 64);
            values.push(val);
        }
        let superblocks = build_superblocks(&marks, bwt_len, RANK_STRIDE);
        Self {
            rate: rate.max(1),
            bwt_len,
            marks,
            superblocks,
            values,
        }
    }

    /// The sampling rate (`>= 1`).
    pub fn rate(&self) -> usize {
        self.rate
    }

    /// Number of BWT positions covered.
    pub fn len(&self) -> usize {
        self.bwt_len
    }

    /// Whether the structure covers zero positions.
    pub fn is_empty(&self) -> bool {
        self.bwt_len == 0
    }

    /// Number of sampled values stored (≈ `bwt_len / rate`).
    pub fn num_samples(&self) -> usize {
        self.values.len()
    }

    /// The sampled SA value at BWT position `bwt_idx`, or `None` if that position
    /// is not sampled. `O(1)` via superblock + word popcount.
    pub fn sample_at(&self, bwt_idx: usize) -> Option<u32> {
        debug_assert!(bwt_idx < self.bwt_len, "bwt_idx out of range");
        let word = self.marks[bwt_idx / 64];
        if word & (1u64 << (bwt_idx % 64)) == 0 {
            return None;
        }
        // rank = number of set bits strictly before bwt_idx = index into `values`.
        let sb = bwt_idx / RANK_STRIDE;
        let rank = self.superblocks[sb] + popcount_prefix(&self.marks, sb * RANK_STRIDE, bwt_idx);
        Some(self.values[rank as usize])
    }
}

/// Popcount of set bits in `[start, end)` of the `u64`-word bitvector `words`.
fn popcount_prefix(words: &[u64], start: usize, end: usize) -> u32 {
    if end <= start {
        return 0;
    }
    let start_word = start / 64;
    let end_word = (end - 1) / 64;
    let start_bit = start % 64;
    let end_bit = end % 64;

    if start_word == end_word {
        let mut w = words[start_word];
        w &= !((1u64 << start_bit) - 1);
        if end_bit != 0 {
            w &= (1u64 << end_bit) - 1;
        }
        return w.count_ones();
    }

    let mut count = (words[start_word] & !((1u64 << start_bit) - 1)).count_ones();
    for w in &words[start_word + 1..end_word] {
        count += w.count_ones();
    }
    let mut last = words[end_word];
    if end_bit != 0 {
        last &= (1u64 << end_bit) - 1;
    }
    count + last.count_ones()
}

/// Prefix popcounts of `marks` at each `stride` boundary (entry `k` = set bits in
/// `[0, k*stride)`); length `ceil(bwt_len/stride) + 1`.
fn build_superblocks(marks: &[u64], bwt_len: usize, stride: usize) -> Vec<u32> {
    let num = (bwt_len + stride - 1) / stride;
    let mut sb = Vec::with_capacity(num + 1);
    sb.push(0u32);
    let mut acc = 0u32;
    for i in 0..num {
        let start = i * stride;
        let end = ((i + 1) * stride).min(bwt_len);
        acc += popcount_prefix(marks, start, end);
        sb.push(acc);
    }
    sb
}

#[cfg(test)]
mod tests {
    use super::*;

    // bwt_len 10, rate 3; sample positions 0,3,6,9 with SA values 9,6,3,0.
    fn fixture() -> SampledSuffixArray {
        SampledSuffixArray::from_sorted_samples(
            10,
            3,
            [(0usize, 9u32), (3, 6), (6, 3), (9, 0)].into_iter(),
        )
    }

    #[test]
    fn sample_at_returns_value_at_sampled_positions() {
        let s = fixture();
        assert_eq!(s.sample_at(0), Some(9));
        assert_eq!(s.sample_at(3), Some(6));
        assert_eq!(s.sample_at(6), Some(3));
        assert_eq!(s.sample_at(9), Some(0));
    }

    #[test]
    fn sample_at_returns_none_at_unsampled_positions() {
        let s = fixture();
        for i in [1usize, 2, 4, 5, 7, 8] {
            assert_eq!(s.sample_at(i), None, "position {i} should be unsampled");
        }
    }

    #[test]
    fn metadata_is_correct() {
        let s = fixture();
        assert_eq!(s.rate(), 3);
        assert_eq!(s.len(), 10);
        assert_eq!(s.num_samples(), 4);
        assert!(!s.is_empty());
    }

    #[test]
    fn handles_a_superblock_boundary() {
        // > RANK_STRIDE positions, sampled sparsely across the boundary.
        let n = RANK_STRIDE * 2 + 5;
        let samples: Vec<(usize, u32)> = (0..n).step_by(64).map(|p| (p, p as u32)).collect();
        let s = SampledSuffixArray::from_sorted_samples(n, 64, samples.iter().copied());
        for &(p, v) in &samples {
            assert_eq!(s.sample_at(p), Some(v));
        }
        assert_eq!(s.sample_at(1), None);
        assert_eq!(s.sample_at(RANK_STRIDE + 1), None);
    }

    #[test]
    fn empty_has_no_samples() {
        let s = SampledSuffixArray::from_sorted_samples(0, 1, std::iter::empty());
        assert!(s.is_empty());
        assert_eq!(s.num_samples(), 0);
    }

    #[test]
    fn rank_is_correct_for_unaligned_positions_in_a_higher_superblock() {
        // Non-word-aligned sample positions, including one inside the 2nd
        // superblock, to exercise popcount masking beyond sb0.
        let n = RANK_STRIDE * 2;
        let samples = [(3usize, 100u32), (67, 101), (1027, 102), (2000, 103)];
        let s = SampledSuffixArray::from_sorted_samples(n, 64, samples.into_iter());
        assert_eq!(s.sample_at(3), Some(100));
        assert_eq!(s.sample_at(67), Some(101));
        assert_eq!(s.sample_at(1027), Some(102)); // bit 3 of word 16, inside sb1
        assert_eq!(s.sample_at(2000), Some(103));
        assert_eq!(s.sample_at(1026), None);
        assert_eq!(s.sample_at(1028), None);
    }
}
