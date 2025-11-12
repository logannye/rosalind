use std::cmp::Ordering;

use crate::genomics::{
    BaseCode, CompressedDNA, CompressedDNAError, RankSelectIndex, ALPHABET_SIZE,
};
use thiserror::Error;

const SENTINEL_BYTE: u8 = b'$';

/// Error type returned by FM-index construction and queries.
#[derive(Debug, Error)]
pub enum FMIndexError {
    /// Reference sequence was empty.
    #[error("reference sequence must be non-empty")]
    EmptyReference,

    /// Encountered an unsupported character in the input.
    #[error("unsupported character '{ch}' at position {position}")]
    UnsupportedCharacter {
        /// Character that could not be encoded.
        ch: char,
        /// Position within the reference where the character was observed.
        position: usize,
    },

    /// Block size was zero.
    #[error("block size must be greater than zero")]
    InvalidBlockSize,

    /// Compression failure bubbling up from `CompressedDNA`.
    #[error("compression error: {0}")]
    Compression(#[from] CompressedDNAError),
}

/// Compact representation of block boundaries and cumulative counts.
#[derive(Debug, Clone)]
pub struct CompressedBoundaries {
    entries: Vec<BlockBoundary>,
}

impl CompressedBoundaries {
    fn new() -> Self {
        Self {
            entries: Vec::new(),
        }
    }

    fn push(&mut self, boundary: BlockBoundary) {
        self.entries.push(boundary);
    }

    /// Return the boundary information for `block_idx`.
    pub fn boundary(&self, block_idx: usize) -> &BlockBoundary {
        &self.entries[block_idx]
    }

    /// Iterator over block boundaries.
    pub fn iter(&self) -> impl Iterator<Item = &BlockBoundary> {
        self.entries.iter()
    }
}

/// Delimits the start of a block and carries cumulative counts at that point.
#[derive(Debug, Clone)]
pub struct BlockBoundary {
    /// Starting offset (inclusive) for the block.
    pub start: usize,
    /// Cumulative counts (A, C, G, T, N) before this block.
    pub cumulative_counts: [u32; ALPHABET_SIZE],
    /// Number of sentinel characters (`$`) encountered before this block.
    pub sentinel_count: u32,
}

/// A block of the BWT string with precomputed rank/select structure.
#[derive(Debug, Clone)]
pub struct BWTBlock {
    start: usize,
    end: usize,
    bwt: CompressedDNA,
    occ: RankSelectIndex,
    sentinel_offset: Option<usize>,
}

impl BWTBlock {
    fn len(&self) -> usize {
        self.end - self.start
    }

    fn rank_symbol(&self, symbol: FmSymbol, position: usize) -> u32 {
        let bounded = position.min(self.len());
        match symbol {
            FmSymbol::Sentinel => {
                if let Some(offset) = self.sentinel_offset {
                    if offset < bounded {
                        1
                    } else {
                        0
                    }
                } else {
                    0
                }
            }
            FmSymbol::Base(code) => {
                let mut count = self.occ.rank(&self.bwt, code, bounded);
                if code == BaseCode::N {
                    if let Some(offset) = self.sentinel_offset {
                        if offset < bounded {
                            count = count.saturating_sub(1);
                        }
                    }
                }
                count
            }
        }
    }
}

/// Symbol used in FM-index queries (includes sentinel).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FmSymbol {
    /// The unique sentinel symbol `$`.
    Sentinel,
    /// One of the canonical DNA bases.
    Base(BaseCode),
}

impl FmSymbol {
    /// Lexicographic order for cumulative C table: `$` < A < C < G < T < N.
    pub fn order(&self) -> usize {
        match self {
            FmSymbol::Sentinel => 0,
            FmSymbol::Base(BaseCode::A) => 1,
            FmSymbol::Base(BaseCode::C) => 2,
            FmSymbol::Base(BaseCode::G) => 3,
            FmSymbol::Base(BaseCode::T) => 4,
            FmSymbol::Base(BaseCode::N) => 5,
        }
    }
}

/// Blocked FM-index structure with per-block rank/select summaries.
#[derive(Debug, Clone)]
pub struct BlockedFMIndex {
    blocks: Vec<BWTBlock>,
    boundaries: CompressedBoundaries,
    c_table: [u32; 6],
    block_size: usize,
    bwt_len: usize,
    sentinel_pos: usize,
}

impl BlockedFMIndex {
    /// Build the FM-index from a reference string using a naive suffix-array
    /// construction. Suitable for demonstration and moderate input sizes.
    pub fn build(reference: &[u8], block_size: usize) -> Result<Self, FMIndexError> {
        if reference.is_empty() {
            return Err(FMIndexError::EmptyReference);
        }
        if block_size == 0 {
            return Err(FMIndexError::InvalidBlockSize);
        }

        let clean = sanitize_reference(reference)?;
        let (bwt, sentinel_pos) = build_bwt(&clean);
        let bwt_len = bwt.len();

        let mut blocks = Vec::new();
        let mut boundaries = CompressedBoundaries::new();
        let mut cumulative_counts = [0u32; ALPHABET_SIZE];
        let mut sentinel_cumulative = 0u32;

        for (block_idx, chunk) in bwt.chunks(block_size).enumerate() {
            let start = block_idx * block_size;
            let end = start + chunk.len();

            boundaries.push(BlockBoundary {
                start,
                cumulative_counts,
                sentinel_count: sentinel_cumulative,
            });

            let mut sanitized = Vec::with_capacity(chunk.len());
            let mut sentinel_offset = None;

            for (offset, &ch) in chunk.iter().enumerate() {
                if ch == SENTINEL_BYTE {
                    sentinel_offset = Some(offset);
                    sentinel_cumulative += 1;
                    sanitized.push(b'N');
                    continue;
                }
                let _code =
                    BaseCode::from_ascii(ch).ok_or_else(|| FMIndexError::UnsupportedCharacter {
                        ch: ch as char,
                        position: start + offset,
                    })?;
                sanitized.push(ch);
            }

            let bwt_compressed = CompressedDNA::compress(&sanitized)?;
            let occ = RankSelectIndex::build(&bwt_compressed);
            let mut block_counts = occ.rank_all(&bwt_compressed, chunk.len());
            if sentinel_offset.is_some() && block_counts[BaseCode::N.index()] > 0 {
                block_counts[BaseCode::N.index()] -= 1;
            }

            blocks.push(BWTBlock {
                start,
                end,
                bwt: bwt_compressed,
                occ,
                sentinel_offset,
            });

            cumulative_counts = add_counts(cumulative_counts, block_counts);
        }

        // Terminal boundary at end of BWT.
        boundaries.push(BlockBoundary {
            start: bwt_len,
            cumulative_counts,
            sentinel_count: sentinel_cumulative,
        });

        let global_totals = cumulative_counts;
        let c_table = build_c_table(global_totals);

        Ok(Self {
            blocks,
            boundaries,
            c_table,
            block_size,
            bwt_len,
            sentinel_pos,
        })
    }

    /// Length of the BWT string (reference length + 1 sentinel).
    pub fn len(&self) -> usize {
        self.bwt_len
    }

    /// Block size used for the FM-index partitioning.
    pub fn block_size(&self) -> usize {
        self.block_size
    }

    /// Number of blocks maintained by the index.
    pub fn num_blocks(&self) -> usize {
        self.blocks.len()
    }

    /// Access the C table used in LF-mapping.
    pub fn c_table(&self) -> &[u32; 6] {
        &self.c_table
    }

    /// Position of the sentinel (`$`) in the BWT string.
    pub fn sentinel_position(&self) -> usize {
        self.sentinel_pos
    }

    /// Retrieve rank of `symbol` in `BWT[..position)`.
    pub fn rank(&self, symbol: FmSymbol, position: usize) -> u32 {
        let bounded = position.min(self.bwt_len);
        let block_idx = bounded / self.block_size;
        let boundary = self.boundaries.boundary(block_idx);

        let mut count = match symbol {
            FmSymbol::Sentinel => boundary.sentinel_count,
            FmSymbol::Base(code) => boundary.cumulative_counts[code.index()],
        };

        if let Some(block) = self.blocks.get(block_idx) {
            let within = bounded - block.start;
            count += block.rank_symbol(symbol, within);
        }

        count
    }

    /// Total occurrences of `symbol` across the entire BWT string.
    pub fn total(&self, symbol: FmSymbol) -> u32 {
        match symbol {
            FmSymbol::Sentinel => 1,
            FmSymbol::Base(code) => {
                let boundary = self.boundaries.boundary(self.blocks.len());
                boundary.cumulative_counts[code.index()]
            }
        }
    }

    /// Access to the raw blocks (useful for specialized processing).
    pub fn blocks(&self) -> &[BWTBlock] {
        &self.blocks
    }

    /// Access block boundaries and cumulative counts.
    pub fn boundaries(&self) -> &CompressedBoundaries {
        &self.boundaries
    }
}

fn sanitize_reference(reference: &[u8]) -> Result<Vec<u8>, FMIndexError> {
    let mut clean = Vec::with_capacity(reference.len());
    for (idx, &ch) in reference.iter().enumerate() {
        match BaseCode::from_ascii(ch) {
            Some(code) => {
                let uppercase = match code {
                    BaseCode::A => b'A',
                    BaseCode::C => b'C',
                    BaseCode::G => b'G',
                    BaseCode::T => b'T',
                    BaseCode::N => b'N',
                };
                clean.push(uppercase);
            }
            None => {
                return Err(FMIndexError::UnsupportedCharacter {
                    ch: ch as char,
                    position: idx,
                });
            }
        }
    }
    Ok(clean)
}

fn build_bwt(reference: &[u8]) -> (Vec<u8>, usize) {
    let mut text = reference.to_vec();
    text.push(SENTINEL_BYTE);
    let sa = build_suffix_array(&text);

    let mut bwt = Vec::with_capacity(text.len());
    let mut sentinel_pos = 0usize;

    for (idx, &sa_idx) in sa.iter().enumerate() {
        let prev = if sa_idx == 0 {
            text.len() - 1
        } else {
            sa_idx - 1
        };
        let ch = text[prev];
        if sa_idx == 0 {
            sentinel_pos = idx;
        }
        bwt.push(ch);
    }

    (bwt, sentinel_pos)
}

fn build_suffix_array(text: &[u8]) -> Vec<usize> {
    let mut sa: Vec<usize> = (0..text.len()).collect();
    sa.sort_by(|&a, &b| compare_suffixes(text, a, b));
    sa
}

fn compare_suffixes(text: &[u8], lhs: usize, rhs: usize) -> Ordering {
    text[lhs..].cmp(&text[rhs..])
}

fn add_counts(lhs: [u32; ALPHABET_SIZE], rhs: [u32; ALPHABET_SIZE]) -> [u32; ALPHABET_SIZE] {
    [
        lhs[0] + rhs[0],
        lhs[1] + rhs[1],
        lhs[2] + rhs[2],
        lhs[3] + rhs[3],
        lhs[4] + rhs[4],
    ]
}

fn build_c_table(totals: [u32; ALPHABET_SIZE]) -> [u32; 6] {
    let sentinel = 1;
    let a = totals[BaseCode::A.index()];
    let c = totals[BaseCode::C.index()];
    let g = totals[BaseCode::G.index()];
    let t = totals[BaseCode::T.index()];

    [
        0,
        sentinel,
        sentinel + a,
        sentinel + a + c,
        sentinel + a + c + g,
        sentinel + a + c + g + t,
    ]
    .map(|value| value)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fm_index_builds_and_ranks() {
        let reference = b"ACGTCGTA";
        let index = BlockedFMIndex::build(reference, 4).expect("index build should succeed");

        assert_eq!(index.len(), reference.len() + 1);
        assert_eq!(index.total(FmSymbol::Base(BaseCode::A)), 2);
        assert_eq!(index.total(FmSymbol::Base(BaseCode::C)), 2);
        assert_eq!(index.total(FmSymbol::Base(BaseCode::G)), 2);
        assert_eq!(index.total(FmSymbol::Base(BaseCode::T)), 2);

        for pos in 0..=index.len() {
            let rank_a = index.rank(FmSymbol::Base(BaseCode::A), pos);
            let naive = naive_rank(reference, b'A', pos);
            assert_eq!(rank_a, naive);
        }

        // Sentinel count should always be 1.
        assert_eq!(index.total(FmSymbol::Sentinel), 1);
        assert_eq!(index.rank(FmSymbol::Sentinel, index.len()), 1);
    }

    fn naive_rank(reference: &[u8], base: u8, position: usize) -> u32 {
        // Build BWT naively for validation.
        let (bwt, _) = build_bwt(&sanitize_reference(reference).unwrap());
        let bounded = position.min(bwt.len());
        bwt[..bounded].iter().filter(|&&ch| ch == base).count() as u32
    }
}
