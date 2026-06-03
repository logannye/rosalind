use crate::genomics::fm_backing::{self, BwtBacking};
use crate::genomics::suffix_array::{sais_u32, SuffixArrayError};
use crate::genomics::FMInterval;
use crate::genomics::{
    BaseCode, CompressedDNA, CompressedDNAError, RankSelectIndex, SampledSuffixArray, ALPHABET_SIZE,
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

    /// Suffix array construction failed.
    #[error("suffix array construction failed: {0}")]
    SuffixArray(#[from] SuffixArrayError),
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

    /// Number of boundary entries (`num_blocks + 1`). (Serialization.)
    pub(crate) fn len(&self) -> usize {
        self.entries.len()
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

    /// Block start offset (inclusive) in the BWT. (Serialization.)
    pub(crate) fn start(&self) -> usize {
        self.start
    }

    /// Block end offset (exclusive) in the BWT. (Serialization.)
    pub(crate) fn end(&self) -> usize {
        self.end
    }

    /// The block's 2-bit BWT payload. (Serialization.)
    pub(crate) fn bwt(&self) -> &CompressedDNA {
        &self.bwt
    }

    /// The block's rank/select (occ) structure. (Serialization.)
    pub(crate) fn occ(&self) -> &RankSelectIndex {
        &self.occ
    }

    /// The sentinel offset within this block, if the sentinel falls here.
    /// (Serialization.)
    pub(crate) fn sentinel_offset(&self) -> Option<usize> {
        self.sentinel_offset
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
    sampled: SampledSuffixArray,
}

impl BlockedFMIndex {
    /// Build the FM-index from a reference string.
    ///
    /// Note: The FM-index is stored in blocked form for rank queries. Locating
    /// suffix array positions uses a sampled SA array, making `sa_at` O(sample_rate).
    pub fn build(reference: &[u8], block_size: usize) -> Result<Self, FMIndexError> {
        if reference.is_empty() {
            return Err(FMIndexError::EmptyReference);
        }
        if block_size == 0 {
            return Err(FMIndexError::InvalidBlockSize);
        }

        let clean = sanitize_reference(reference)?;
        // For now, pick a conservative SA sampling rate. This will become
        // configurable and/or derived from the on-disk index format.
        let sa_sample_rate = 32usize;
        let (bwt, sentinel_pos, sampled) = build_bwt_and_sa_samples(&clean, sa_sample_rate)?;
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
            sampled,
        })
    }

    /// Length of the BWT string (reference length + 1 sentinel).
    pub fn len(&self) -> usize {
        self.bwt_len
    }

    /// Whether the BWT is empty (no bases indexed).
    pub fn is_empty(&self) -> bool {
        self.bwt_len == 0
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

    /// Perform exact FM-index backward search over `pattern` (ASCII A/C/G/T/N).
    ///
    /// Returns an [`FMInterval`] over BWT indices for suffixes matching the pattern.
    pub fn backward_search(&self, pattern: &[u8]) -> FMInterval {
        fm_backing::backward_search(self, pattern)
    }

    /// Locate up to `max_hits` suffix array positions for the given interval.
    ///
    /// Returned positions are 0-based reference coordinates (excluding the sentinel).
    pub fn locate_interval(&self, interval: FMInterval, max_hits: usize) -> Vec<u32> {
        fm_backing::locate_interval(self, interval, max_hits)
    }

    /// Suffix array sample rate used for locating.
    pub fn sa_sample_rate(&self) -> usize {
        self.sampled.rate()
    }

    /// The compact sampled suffix array backing `sa_at`.
    pub fn sampled(&self) -> &SampledSuffixArray {
        &self.sampled
    }

    /// Retrieve rank of `symbol` in `BWT[..position)`.
    pub fn rank(&self, symbol: FmSymbol, position: usize) -> u32 {
        fm_backing::rank(self, symbol, position)
    }

    /// Total occurrences of `symbol` across the entire BWT string.
    pub fn total(&self, symbol: FmSymbol) -> u32 {
        fm_backing::total(self, symbol)
    }

    /// Access to the raw blocks (useful for specialized processing).
    pub fn blocks(&self) -> &[BWTBlock] {
        &self.blocks
    }

    /// Access block boundaries and cumulative counts.
    pub fn boundaries(&self) -> &CompressedBoundaries {
        &self.boundaries
    }

    /// Retrieve the symbol stored at `index` in the BWT string.
    pub fn symbol_at(&self, index: usize) -> FmSymbol {
        fm_backing::symbol_at(self, index)
    }

    /// Compute the suffix array value corresponding to the provided BWT index.
    pub fn sa_at(&self, index: usize) -> usize {
        fm_backing::sa_at(self, index)
    }
}

impl BwtBacking for BlockedFMIndex {
    fn bwt_len(&self) -> usize {
        self.bwt_len
    }

    fn block_size(&self) -> usize {
        self.block_size
    }

    fn num_blocks(&self) -> usize {
        self.blocks.len()
    }

    fn sentinel_pos(&self) -> usize {
        self.sentinel_pos
    }

    fn sample_rate(&self) -> usize {
        self.sampled.rate()
    }

    fn c_table(&self) -> [u32; 6] {
        self.c_table
    }

    fn boundary_base(&self, block_idx: usize, base_index: usize) -> u32 {
        self.boundaries.boundary(block_idx).cumulative_counts[base_index]
    }

    fn boundary_sentinel(&self, block_idx: usize) -> u32 {
        self.boundaries.boundary(block_idx).sentinel_count
    }

    fn block_rank(&self, block_idx: usize, symbol: FmSymbol, within: usize) -> u32 {
        self.blocks[block_idx].rank_symbol(symbol, within)
    }

    fn block_symbol(&self, block_idx: usize, within: usize) -> FmSymbol {
        let block = &self.blocks[block_idx];
        let base = block
            .bwt
            .base_at(within)
            .expect("BWT block should contain sequence data");
        let code = BaseCode::from_ascii(base)
            .expect("BWT symbol must be a valid DNA base except sentinel");
        FmSymbol::Base(code)
    }

    fn sampled_at(&self, index: usize) -> Option<u32> {
        self.sampled.sample_at(index)
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
                // IUPAC ambiguity codes (R,Y,S,W,K,M,B,D,H,V) appear in real
                // assemblies (GRCh38's primary assembly, many bacterial/viral
                // references). Map them to N — matching bwa/bowtie — so the index
                // step ingests a stock reference instead of aborting; the N-mask
                // carries the ambiguity, and these positions never match an
                // A/C/G/T query. Genuinely non-sequence bytes still error loudly.
                if matches!(
                    ch.to_ascii_uppercase(),
                    b'R' | b'Y' | b'S' | b'W' | b'K' | b'M' | b'B' | b'D' | b'H' | b'V'
                ) {
                    clean.push(b'N');
                } else {
                    return Err(FMIndexError::UnsupportedCharacter {
                        ch: ch as char,
                        position: idx,
                    });
                }
            }
        }
    }
    Ok(clean)
}

fn build_bwt_and_sa_samples(
    reference: &[u8],
    sa_sample_rate: usize,
) -> Result<(Vec<u8>, usize, SampledSuffixArray), FMIndexError> {
    // Map A/C/G/T/N -> 1..5, sentinel -> 0.
    let mut text: Vec<u32> = Vec::with_capacity(reference.len() + 1);
    for &b in reference {
        let code = BaseCode::from_ascii(b).expect("reference already sanitized");
        let sym = match code {
            BaseCode::A => 1u32,
            BaseCode::C => 2u32,
            BaseCode::G => 3u32,
            BaseCode::T => 4u32,
            BaseCode::N => 5u32,
        };
        text.push(sym);
    }
    text.push(0); // sentinel

    let sa = sais_u32(&text, 5)?;

    let mut bwt = Vec::with_capacity(text.len());
    let mut sentinel_pos = 0usize;

    for (bwt_idx, &sa_idx_u32) in sa.iter().enumerate() {
        let sa_idx = sa_idx_u32 as usize;
        let prev = if sa_idx == 0 {
            text.len() - 1
        } else {
            sa_idx - 1
        };
        let ch = if prev == text.len() - 1 {
            SENTINEL_BYTE
        } else {
            // Map back to ASCII reference base.
            match text[prev] {
                1 => b'A',
                2 => b'C',
                3 => b'G',
                4 => b'T',
                5 => b'N',
                _ => SENTINEL_BYTE,
            }
        };
        if sa_idx == 0 {
            sentinel_pos = bwt_idx;
        }
        bwt.push(ch);
    }

    let rate = sa_sample_rate.max(1);
    // `% == 0` (not `usize::is_multiple_of`, stable only in 1.87) to hold MSRV 1.83.
    #[allow(clippy::manual_is_multiple_of)]
    let sampled = SampledSuffixArray::from_sorted_samples(
        text.len(),
        rate,
        sa.iter().enumerate().filter_map(|(bwt_idx, &sa_idx)| {
            ((sa_idx as usize) % rate == 0).then_some((bwt_idx, sa_idx))
        }),
    );

    Ok((bwt, sentinel_pos, sampled))
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
    use crate::genomics::BWTAligner;

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
        let clean = sanitize_reference(reference).unwrap();
        let (bwt, _, _) = build_bwt_and_sa_samples(&clean, 1).unwrap();
        let bounded = position.min(bwt.len());
        bwt[..bounded].iter().filter(|&&ch| ch == base).count() as u32
    }

    #[test]
    fn sanitize_maps_iupac_ambiguity_codes_to_n() {
        // Real references carry IUPAC degeneracy codes; the index must ingest
        // them (mapped to N) rather than aborting. Genuinely invalid bytes still
        // error loudly.
        let clean = sanitize_reference(b"ACGTRYSWKMryswkmBDHV").unwrap();
        assert_eq!(&clean, b"ACGTNNNNNNNNNNNNNNNN");
        // A stock-reference-shaped sequence with ambiguity codes builds an index.
        let index = BlockedFMIndex::build(b"ACGTRYSWKMACGTACGTAC", 4)
            .expect("index build must succeed over IUPAC codes");
        let _ = index;
        // A non-sequence byte is still rejected.
        assert!(sanitize_reference(b"ACGT@CGT").is_err());
    }

    #[test]
    fn sa_at_recovers_reference_position() {
        let reference = b"ACGTACGT";
        let mut aligner = BWTAligner::new(reference).expect("aligner should initialize");
        let result = aligner
            .align_read(b"ACGT")
            .expect("alignment should succeed");
        assert!(result.has_candidates());

        let index = BlockedFMIndex::build(reference, 4).expect("index build should succeed");
        let position = index.sa_at(result.interval.lower as usize);
        assert!(position + 4 <= reference.len());
        assert_eq!(&reference[position..position + 4], b"ACGT");
    }

    #[test]
    fn owned_backing_matches_generic_ops() {
        use crate::genomics::fm_backing::{self, BwtBacking};

        let reference = b"ACGTNACGTACGTACGTNNACG";
        let index = BlockedFMIndex::build(reference, 5).expect("build");

        // The generic ops, run over the owned backing through the BwtBacking
        // trait, match a naive count over the raw BWT — i.e. the trait surface
        // feeds the algorithm correct data (not merely "equals the public method",
        // which now routes through the same generic op and so would be tautological).
        let clean = sanitize_reference(reference).unwrap();
        let (bwt, _, _) = build_bwt_and_sa_samples(&clean, 1).unwrap();
        for symbol in [
            FmSymbol::Sentinel,
            FmSymbol::Base(BaseCode::A),
            FmSymbol::Base(BaseCode::C),
            FmSymbol::Base(BaseCode::G),
            FmSymbol::Base(BaseCode::T),
            FmSymbol::Base(BaseCode::N),
        ] {
            // The raw BWT stores the sentinel as `$` and real `N`s as `N`, which
            // matches the FM-index's symbol semantics (sentinel counted separately
            // from `N`), so a byte count is the ground truth.
            let byte = match symbol {
                FmSymbol::Sentinel => b'$',
                FmSymbol::Base(BaseCode::A) => b'A',
                FmSymbol::Base(BaseCode::C) => b'C',
                FmSymbol::Base(BaseCode::G) => b'G',
                FmSymbol::Base(BaseCode::T) => b'T',
                FmSymbol::Base(BaseCode::N) => b'N',
            };
            for pos in 0..=index.len() {
                let naive = bwt[..pos.min(bwt.len())]
                    .iter()
                    .filter(|&&c| c == byte)
                    .count() as u32;
                assert_eq!(
                    fm_backing::rank(&index, symbol, pos),
                    naive,
                    "generic rank vs naive @ {pos} for {symbol:?}"
                );
                assert_eq!(
                    index.rank(symbol, pos),
                    naive,
                    "public rank vs naive @ {pos} for {symbol:?}"
                );
            }
        }

        // c_table by value (trait) equals the inherent c_table by reference.
        assert_eq!(BwtBacking::c_table(&index), *index.c_table());

        // sa_at via the generic op matches the public method over the interval.
        let interval = index.backward_search(b"ACGT");
        for bwt_idx in (interval.lower as usize)..(interval.upper as usize) {
            assert_eq!(fm_backing::sa_at(&index, bwt_idx), index.sa_at(bwt_idx));
        }
    }

    #[test]
    fn sampled_sa_is_sparse_not_dense() {
        // A reference long enough that a dense (one-slot-per-position) sample
        // array would be ~bwt_len; the compact structure stores ~bwt_len/rate.
        let reference = vec![b'A'; 4096];
        let index = BlockedFMIndex::build(&reference, 64).expect("build");
        let sampled = index.sampled();
        assert_eq!(sampled.len(), reference.len() + 1); // covers every BWT position
                                                        // Far fewer stored values than positions (sampled at `rate`).
        assert!(
            sampled.num_samples() <= sampled.len() / sampled.rate() + 1,
            "expected ~len/rate samples, got {} for len {} rate {}",
            sampled.num_samples(),
            sampled.len(),
            sampled.rate(),
        );
        assert!(sampled.num_samples() * 4 < sampled.len());
    }

    #[test]
    fn serialization_accessors_expose_backing() {
        let reference = b"ACGTNACGTACGTACGT";
        let index = BlockedFMIndex::build(reference, 6).expect("build");

        // Boundaries length is num_blocks + 1 (one per block + terminal).
        assert_eq!(index.boundaries().len(), index.num_blocks() + 1);

        // Every block exposes its 2-bit BWT, occ structure, and span.
        for block in index.blocks() {
            assert_eq!(block.end() - block.start(), block.bwt().len());
            // Each of the 5 occ rank bitvectors has ceil(n/64) words.
            let n = block.bwt().len();
            let expected_words = n.div_ceil(64);
            for bv in block.occ().bitvectors() {
                assert_eq!(bv.len(), expected_words);
            }
            // Each of the 5 occ superblock arrays has ceil(n/stride) + 1 entries.
            let expected_sb = n.div_ceil(block.occ().stride()) + 1;
            for sb in block.occ().superblocks() {
                assert_eq!(sb.len(), expected_sb);
            }
            let _ = block.sentinel_offset();
        }

        // The sampled SA exposes marks/superblocks/values.
        let s = index.sampled();
        assert_eq!(s.marks().len(), s.len().div_ceil(64));
        assert_eq!(s.values().len(), s.num_samples());
        assert!(!s.superblocks().is_empty());
    }
}
