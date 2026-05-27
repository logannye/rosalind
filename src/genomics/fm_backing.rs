//! One FM query algorithm over two backings.
//!
//! The blocked FM-index query algorithm (`rank`/`symbol_at`/`lf_index`/
//! `backward_search`/`sa_at`/`locate_interval`/`total`) reads a small, fixed
//! *data surface*: block-boundary cumulative counts, per-block rank/symbol, the
//! C-table, and the sampled suffix array. [`BwtBacking`] is exactly that surface.
//!
//! Making the algorithm generic over the backing is what lets a persisted,
//! memory-mapped index ([`crate::genomics::FmIndexView`]) answer queries with
//! the *same* code as the in-RAM [`crate::genomics::BlockedFMIndex`] — the view
//! can only differ in how it fetches bytes, never in the query logic, so the
//! equivalence gate (Phase B3b.2) is a complete check.
//!
//! This is an internal seam (`pub(crate)`), not part of the public API.

use crate::genomics::{BaseCode, FMInterval, FmSymbol};

/// The data surface the blocked FM query algorithm reads.
///
/// Boundary accessors are valid for `block_idx` in `0..=num_blocks()` (the
/// terminal boundary at `num_blocks()` carries the whole-BWT cumulative counts);
/// block accessors are valid for `block_idx` in `0..num_blocks()`.
pub(crate) trait BwtBacking {
    /// Length of the BWT string (reference length + 1 sentinel).
    fn bwt_len(&self) -> usize;
    /// Block size used to partition the BWT.
    fn block_size(&self) -> usize;
    /// Number of blocks (`= ceil(bwt_len / block_size)`).
    fn num_blocks(&self) -> usize;
    /// Position of the sentinel (`$`) in the BWT string.
    fn sentinel_pos(&self) -> usize;
    /// Suffix-array sample rate (`>= 1`); used only for the `sa_at` debug guard.
    fn sample_rate(&self) -> usize;
    /// The C table (`$`,A,C,G,T,N cumulative starts) used in LF-mapping.
    fn c_table(&self) -> [u32; 6];
    /// Cumulative count of base `base_index` in `BWT[.. block_idx * block_size]`.
    fn boundary_base(&self, block_idx: usize, base_index: usize) -> u32;
    /// Cumulative sentinel count before block `block_idx`.
    fn boundary_sentinel(&self, block_idx: usize) -> u32;
    /// Rank of `symbol` within block `block_idx` over its first `within` positions.
    fn block_rank(&self, block_idx: usize, symbol: FmSymbol, within: usize) -> u32;
    /// The BWT symbol at `within` inside block `block_idx`. The caller guarantees
    /// `block_idx * block_size + within != sentinel_pos` (the sentinel is handled
    /// in [`symbol_at`]), so this never needs to report `FmSymbol::Sentinel`.
    fn block_symbol(&self, block_idx: usize, within: usize) -> FmSymbol;
    /// Sampled suffix-array value at BWT position `index`, if sampled.
    fn sampled_at(&self, index: usize) -> Option<u32>;
}

/// Rank of `symbol` in `BWT[..position)`.
pub(crate) fn rank<B: BwtBacking>(b: &B, symbol: FmSymbol, position: usize) -> u32 {
    let bounded = position.min(b.bwt_len());
    let block_idx = bounded / b.block_size();

    let mut count = match symbol {
        FmSymbol::Sentinel => b.boundary_sentinel(block_idx),
        FmSymbol::Base(code) => b.boundary_base(block_idx, code.index()),
    };

    // The within-block term is skipped only at the terminal position when
    // `bwt_len` is an exact multiple of `block_size` (then `block_idx ==
    // num_blocks`, which has no block — only the terminal boundary).
    if block_idx < b.num_blocks() {
        let within = bounded - block_idx * b.block_size();
        count += b.block_rank(block_idx, symbol, within);
    }

    count
}

/// Total occurrences of `symbol` across the entire BWT string.
pub(crate) fn total<B: BwtBacking>(b: &B, symbol: FmSymbol) -> u32 {
    match symbol {
        FmSymbol::Sentinel => 1,
        FmSymbol::Base(code) => b.boundary_base(b.num_blocks(), code.index()),
    }
}

/// The symbol stored at `index` in the BWT string.
pub(crate) fn symbol_at<B: BwtBacking>(b: &B, index: usize) -> FmSymbol {
    assert!(index < b.bwt_len(), "BWT index out of range");
    if index == b.sentinel_pos() {
        return FmSymbol::Sentinel;
    }
    let block_idx = index / b.block_size();
    let within = index - block_idx * b.block_size();
    b.block_symbol(block_idx, within)
}

/// LF-mapping: the BWT index whose suffix is one position earlier.
pub(crate) fn lf_index<B: BwtBacking>(b: &B, index: usize) -> usize {
    let symbol = symbol_at(b, index);
    let occ_inclusive = rank(b, symbol, index + 1);
    let c_row = b.c_table()[symbol.order()] as usize;
    c_row + occ_inclusive as usize - 1
}

/// Exact FM-index backward search over `pattern` (ASCII A/C/G/T/N).
pub(crate) fn backward_search<B: BwtBacking>(b: &B, pattern: &[u8]) -> FMInterval {
    let mut interval = FMInterval::full(b.bwt_len());

    for &ch in pattern.iter().rev() {
        let base_code = match BaseCode::from_ascii(ch) {
            Some(code) => code,
            None => return FMInterval { lower: 0, upper: 0 },
        };
        let symbol = FmSymbol::Base(base_code);
        let c_row = b.c_table()[symbol.order()];
        let new_lower = c_row + rank(b, symbol, interval.lower as usize);
        let new_upper = c_row + rank(b, symbol, interval.upper as usize);
        interval = FMInterval {
            lower: new_lower,
            upper: new_upper,
        };
        if interval.is_empty() {
            break;
        }
    }

    interval
}

/// The suffix array value corresponding to BWT index `index`.
pub(crate) fn sa_at<B: BwtBacking>(b: &B, index: usize) -> usize {
    assert!(index < b.bwt_len(), "BWT index out of range");

    let mut current = index;
    let mut lf_steps = 0usize;

    // Following LF decreases the SA value by 1 (mod n), so within `rate` steps we
    // must reach a sampled position; the recovered SA value is the sample plus
    // the number of LF steps taken.
    loop {
        if let Some(sampled) = b.sampled_at(current) {
            return sampled as usize + lf_steps;
        }
        current = lf_index(b, current);
        lf_steps += 1;
        debug_assert!(
            lf_steps <= b.sample_rate() + 1,
            "LF steps exceeded sample rate; sampling invariant violated"
        );
    }
}

/// Locate up to `max_hits` 0-based reference positions for `interval`.
pub(crate) fn locate_interval<B: BwtBacking>(
    b: &B,
    interval: FMInterval,
    max_hits: usize,
) -> Vec<u32> {
    let max_hits = max_hits.max(1);
    let mut out = Vec::new();

    let reference_len = b.bwt_len().saturating_sub(1);
    let lower = interval.lower as usize;
    let upper = interval.upper as usize;
    for bwt_idx in lower..upper {
        if out.len() >= max_hits {
            break;
        }
        let sa = sa_at(b, bwt_idx);
        if sa < reference_len {
            out.push(sa as u32);
        }
    }

    out
}
