//! Borrowed, zero-copy query view over a memory-mapped index.
//!
//! [`FmIndexView`] reads the FM-index sections directly from the mmap as
//! `&[u64]`/`&[u32]` (via a checked `slice::align_to`) and implements
//! [`crate::genomics::fm_backing::BwtBacking`]. Because the FM query algorithm is
//! generic over that trait, the view runs the *same* code as the in-RAM
//! [`crate::genomics::BlockedFMIndex`] — it can differ only in how it fetches
//! bytes, never in logic. The zero-copy reinterpretation is value-correct only on
//! little-endian hosts (the file is little-endian; [`FmIndexView::new`] rejects
//! big-endian hosts).

use crate::core::{ContigSet, Locus};
use crate::genomics::fm_backing::{self, BwtBacking};
use crate::genomics::index::format::{SectionEntry, SectionKind};
use crate::genomics::index::io::IndexIoError;
use crate::genomics::{popcount_range, BaseCode, FMInterval, FmSymbol, RANK_STRIDE};

/// A borrowed FM-index over memory-mapped bytes.
#[derive(Debug)]
pub struct FmIndexView<'a> {
    bwt_len: usize,
    block_size: usize,
    num_blocks: usize,
    sentinel_pos: usize,
    sample_rate: usize,
    c_table: [u32; 6],
    /// `6 * (num_blocks + 1)` u32: `[c0..c4, sentinel]` per boundary entry.
    boundaries: &'a [u32],
    /// `num_blocks` u64 record offsets, relative to the start of `blocks`.
    block_dir: &'a [u64],
    /// The whole `Blocks` section payload (records are sliced on demand).
    blocks: &'a [u8],
    sampled: SampledView<'a>,
}

#[derive(Debug)]
struct SampledView<'a> {
    bwt_len: usize,
    marks: &'a [u64],
    superblocks: &'a [u32],
    values: &'a [u32],
}

/// A parsed block record (computed on demand from the directory; holds only
/// borrowed slices into the mmap).
struct BlockView<'a> {
    start: usize,
    end: usize,
    sentinel_offset: Option<usize>,
    stride: usize,
    bwt_data: &'a [u64],
    bwt_amb: &'a [u64],
    occ_bitvecs: [&'a [u64]; 5],
    occ_superblocks: [&'a [u32]; 5],
}

impl<'a> FmIndexView<'a> {
    /// Build a view from the mmap `bytes` and the parsed `sections`. Validates
    /// section presence, lengths, and the little-endian host requirement.
    pub(crate) fn new(bytes: &'a [u8], sections: &[SectionEntry]) -> Result<Self, IndexIoError> {
        if cfg!(target_endian = "big") {
            return Err(IndexIoError::Invalid(
                "zero-copy index view requires a little-endian host".to_string(),
            ));
        }

        let fm_meta = section_bytes(bytes, sections, SectionKind::FmMeta)?;
        let mut o = 0usize;
        let block_size = read_u64(fm_meta, &mut o)? as usize;
        let bwt_len = read_u64(fm_meta, &mut o)? as usize;
        let sentinel_pos = read_u64(fm_meta, &mut o)? as usize;
        let sample_rate = read_u64(fm_meta, &mut o)? as usize;
        let num_blocks = read_u64(fm_meta, &mut o)? as usize;
        let mut c_table = [0u32; 6];
        for slot in &mut c_table {
            *slot = read_u32(fm_meta, &mut o)?;
        }

        let boundaries_bytes = section_bytes(bytes, sections, SectionKind::Boundaries)?;
        let boundaries = as_u32_slice(boundaries_bytes)?;
        if boundaries.len() != 6 * (num_blocks + 1) {
            return Err(IndexIoError::Invalid(
                "boundaries length mismatch".to_string(),
            ));
        }

        let blocks = section_bytes(bytes, sections, SectionKind::Blocks)?;
        let dir_bytes = num_blocks
            .checked_mul(8)
            .ok_or_else(|| IndexIoError::Invalid("block directory overflow".to_string()))?;
        if dir_bytes > blocks.len() {
            return Err(IndexIoError::Invalid(
                "block directory out of bounds".to_string(),
            ));
        }
        let block_dir = as_u64_slice(&blocks[..dir_bytes])?;
        validate_block_records(blocks, block_dir)?;

        // SaSamples.
        let sa = section_bytes(bytes, sections, SectionKind::SaSamples)?;
        let mut so = 0usize;
        let rate = read_u64(sa, &mut so)? as usize;
        if rate != sample_rate {
            return Err(IndexIoError::Invalid(format!(
                "SaSamples rate {rate} != FmMeta sa_sample_rate {sample_rate}"
            )));
        }
        let s_bwt_len = read_u64(sa, &mut so)? as usize;
        let marks_words = read_u64(sa, &mut so)? as usize;
        let sb_len = read_u64(sa, &mut so)? as usize;
        let values_len = read_u64(sa, &mut so)? as usize;
        let marks = as_u64_slice(slice_exact(sa, &mut so, marks_words * 8)?)?;
        let superblocks = as_u32_slice(slice_exact(sa, &mut so, sb_len * 4)?)?;
        let values = as_u32_slice(slice_exact(sa, &mut so, values_len * 4)?)?;

        Ok(Self {
            bwt_len,
            block_size,
            num_blocks,
            sentinel_pos,
            sample_rate,
            c_table,
            boundaries,
            block_dir,
            blocks,
            sampled: SampledView {
                bwt_len: s_bwt_len,
                marks,
                superblocks,
                values,
            },
        })
    }

    /// Exact FM-index backward search (delegates to the shared generic op).
    pub fn backward_search(&self, pattern: &[u8]) -> FMInterval {
        fm_backing::backward_search(self, pattern)
    }

    /// Suffix array value at BWT index `index` (delegates to the shared generic op).
    pub fn sa_at(&self, index: usize) -> usize {
        fm_backing::sa_at(self, index)
    }

    /// Locate up to `max_hits` 0-based reference positions for `interval`.
    pub fn locate_interval(&self, interval: FMInterval, max_hits: usize) -> Vec<u32> {
        fm_backing::locate_interval(self, interval, max_hits)
    }

    /// Parse block record `block_idx` from the directory (borrowed slices only).
    fn block(&self, block_idx: usize) -> BlockView<'a> {
        let rec = self.block_dir[block_idx] as usize;
        let s = &self.blocks[rec..];
        let mut o = 0usize;
        let start = read_u64_infallible(s, &mut o) as usize;
        let end = read_u64_infallible(s, &mut o) as usize;
        let sentinel_raw = read_i64_infallible(s, &mut o);
        let stride = read_u64_infallible(s, &mut o) as usize;
        let bwt_data_words = read_u64_infallible(s, &mut o) as usize;
        let bwt_amb_words = read_u64_infallible(s, &mut o) as usize;
        let occ_bitvec_words = read_u64_infallible(s, &mut o) as usize;
        let occ_superblock_len = read_u64_infallible(s, &mut o) as usize;

        let bwt_data = as_u64_slice(&s[o..o + bwt_data_words * 8]).expect("aligned");
        o += bwt_data_words * 8;
        let bwt_amb = as_u64_slice(&s[o..o + bwt_amb_words * 8]).expect("aligned");
        o += bwt_amb_words * 8;

        let mut occ_bitvecs: [&[u64]; 5] = [&[]; 5];
        for slot in &mut occ_bitvecs {
            *slot = as_u64_slice(&s[o..o + occ_bitvec_words * 8]).expect("aligned");
            o += occ_bitvec_words * 8;
        }
        let mut occ_superblocks: [&[u32]; 5] = [&[]; 5];
        for slot in &mut occ_superblocks {
            *slot = as_u32_slice(&s[o..o + occ_superblock_len * 4]).expect("aligned");
            o += occ_superblock_len * 4;
        }

        BlockView {
            start,
            end,
            sentinel_offset: if sentinel_raw < 0 {
                None
            } else {
                Some(sentinel_raw as usize)
            },
            stride,
            bwt_data,
            bwt_amb,
            occ_bitvecs,
            occ_superblocks,
        }
    }
}

impl BwtBacking for FmIndexView<'_> {
    fn bwt_len(&self) -> usize {
        self.bwt_len
    }
    fn block_size(&self) -> usize {
        self.block_size
    }
    fn num_blocks(&self) -> usize {
        self.num_blocks
    }
    fn sentinel_pos(&self) -> usize {
        self.sentinel_pos
    }
    fn sample_rate(&self) -> usize {
        self.sample_rate
    }
    fn c_table(&self) -> [u32; 6] {
        self.c_table
    }
    fn boundary_base(&self, block_idx: usize, base_index: usize) -> u32 {
        self.boundaries[block_idx * 6 + base_index]
    }
    fn boundary_sentinel(&self, block_idx: usize) -> u32 {
        self.boundaries[block_idx * 6 + 5]
    }
    fn block_rank(&self, block_idx: usize, symbol: FmSymbol, within: usize) -> u32 {
        let block = self.block(block_idx);
        let n = block.end - block.start;
        let bounded = within.min(n);
        match symbol {
            FmSymbol::Sentinel => match block.sentinel_offset {
                Some(off) if off < bounded => 1,
                _ => 0,
            },
            FmSymbol::Base(code) => {
                // Mirror RankSelectIndex::rank over the borrowed occ slices.
                let bi = code.index();
                let sb = bounded / block.stride;
                let within_start = sb * block.stride;
                let mut count = block.occ_superblocks[bi][sb]
                    + popcount_range(block.occ_bitvecs[bi], within_start, bounded);
                // Mirror BWTBlock::rank_symbol's N correction: the sentinel is
                // stored as `N` in the block BWT but is not a real `N`.
                if code == BaseCode::N {
                    if let Some(off) = block.sentinel_offset {
                        if off < bounded {
                            count = count.saturating_sub(1);
                        }
                    }
                }
                count
            }
        }
    }
    fn block_symbol(&self, block_idx: usize, within: usize) -> FmSymbol {
        let block = self.block(block_idx);
        // Mirror CompressedDNA::base_at: ambiguity (1 bit/base) marks `N`, else
        // a 2-bit code (32 bases per u64 word).
        let amb_word = block.bwt_amb[within / 64];
        if amb_word & (1u64 << (within % 64)) != 0 {
            return FmSymbol::Base(BaseCode::N);
        }
        let data_word = block.bwt_data[within / 32];
        let code = ((data_word >> ((within % 32) * 2)) & 0b11) as u8;
        let base = match code {
            0 => BaseCode::A,
            1 => BaseCode::C,
            2 => BaseCode::G,
            _ => BaseCode::T,
        };
        FmSymbol::Base(base)
    }
    fn sampled_at(&self, index: usize) -> Option<u32> {
        debug_assert!(index < self.sampled.bwt_len);
        let word = self.sampled.marks[index / 64];
        if word & (1u64 << (index % 64)) == 0 {
            return None;
        }
        let sb = index / RANK_STRIDE;
        let rank = self.sampled.superblocks[sb]
            + popcount_range(self.sampled.marks, sb * RANK_STRIDE, index);
        Some(self.sampled.values[rank as usize])
    }
}

/// An FM-index view paired with the contig set: exact-match queries to `Locus`es.
#[derive(Debug)]
pub struct GenomeIndexView<'a> {
    fm: FmIndexView<'a>,
    contigs: &'a ContigSet,
}

impl<'a> GenomeIndexView<'a> {
    pub(crate) fn new(fm: FmIndexView<'a>, contigs: &'a ContigSet) -> Self {
        Self { fm, contigs }
    }

    /// The borrowed FM-index.
    pub fn fm(&self) -> &FmIndexView<'a> {
        &self.fm
    }

    /// The contig set.
    pub fn contigs(&self) -> &ContigSet {
        self.contigs
    }

    /// Locate exact occurrences of `pattern`, returning up to `max_hits` `Locus`es
    /// with boundary-straddling hits removed, sorted by `(contig, pos)` — the same
    /// surface (and result) as `GenomeIndex::locate_exact`.
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

/// Validate that every block record (at the directory's offsets) lies fully
/// within the `Blocks` section, using checked arithmetic so a corrupt word-count
/// can never overflow into an in-bounds-but-wrong slice. Lets `block()` use
/// infallible slicing on a file that passed `open()`.
fn validate_block_records(blocks: &[u8], block_dir: &[u64]) -> Result<(), IndexIoError> {
    const HEADER: usize = 64; // 8 u64/i64 fields
    for &rec_off in block_dir {
        let rec = rec_off as usize;
        let header_end = rec
            .checked_add(HEADER)
            .filter(|&e| e <= blocks.len())
            .ok_or_else(|| {
                IndexIoError::Invalid("block record header out of bounds".to_string())
            })?;
        // Re-read the four word-count fields (header layout: start,end,sentinel,
        // stride at [rec..rec+32]; then the 4 counts at [rec+32..rec+64]).
        let mut o = rec + 32;
        let bwt_data_words = read_u64(blocks, &mut o)? as usize;
        let bwt_amb_words = read_u64(blocks, &mut o)? as usize;
        let occ_bitvec_words = read_u64(blocks, &mut o)? as usize;
        let occ_superblock_len = read_u64(blocks, &mut o)? as usize;
        debug_assert_eq!(o, header_end);
        // Bytes block() touches after the header: bwt_data + bwt_amb + 5×bitvec
        // (u64) + 5×superblock (u32). (totals/pad are not read by block().)
        let payload = bwt_data_words
            .checked_mul(8)
            .and_then(|x| bwt_amb_words.checked_mul(8).and_then(|y| x.checked_add(y)))
            .and_then(|x| {
                occ_bitvec_words
                    .checked_mul(8)
                    .and_then(|y| y.checked_mul(5))
                    .and_then(|y| x.checked_add(y))
            })
            .and_then(|x| {
                occ_superblock_len
                    .checked_mul(4)
                    .and_then(|y| y.checked_mul(5))
                    .and_then(|y| x.checked_add(y))
            })
            .ok_or_else(|| IndexIoError::Invalid("block record size overflow".to_string()))?;
        let rec_end = header_end
            .checked_add(payload)
            .filter(|&e| e <= blocks.len())
            .ok_or_else(|| {
                IndexIoError::Invalid("block record payload out of bounds".to_string())
            })?;
        let _ = rec_end;
    }
    Ok(())
}

/// The bytes of section `kind`.
fn section_bytes<'a>(
    bytes: &'a [u8],
    sections: &[SectionEntry],
    kind: SectionKind,
) -> Result<&'a [u8], IndexIoError> {
    let s = sections
        .iter()
        .find(|s| s.kind == kind)
        .ok_or_else(|| IndexIoError::Invalid(format!("missing section {kind:?}")))?;
    let start = s.offset as usize;
    let end = start + s.bytes as usize;
    bytes
        .get(start..end)
        .ok_or_else(|| IndexIoError::Invalid(format!("section {kind:?} out of bounds")))
}

/// Reinterpret a length-8k byte slice (8-aligned start) as `&[u64]`.
fn as_u64_slice(bytes: &[u8]) -> Result<&[u64], IndexIoError> {
    // SAFETY: any byte pattern is a valid `u64`; we require an empty prefix
    // (start is 8-aligned) and empty suffix (length is a multiple of 8). Values
    // are correct on little-endian hosts (checked in `FmIndexView::new`).
    let (prefix, mid, suffix) = unsafe { bytes.align_to::<u64>() };
    if !prefix.is_empty() || !suffix.is_empty() {
        return Err(IndexIoError::Invalid("misaligned u64 array".to_string()));
    }
    Ok(mid)
}

/// Reinterpret a length-4k byte slice (4-aligned start) as `&[u32]`.
fn as_u32_slice(bytes: &[u8]) -> Result<&[u32], IndexIoError> {
    // SAFETY: as `as_u64_slice`, for `u32` (4-aligned, length a multiple of 4).
    let (prefix, mid, suffix) = unsafe { bytes.align_to::<u32>() };
    if !prefix.is_empty() || !suffix.is_empty() {
        return Err(IndexIoError::Invalid("misaligned u32 array".to_string()));
    }
    Ok(mid)
}

/// Take the next `len` bytes from `bytes` starting at `*offset`, advancing it.
fn slice_exact<'a>(
    bytes: &'a [u8],
    offset: &mut usize,
    len: usize,
) -> Result<&'a [u8], IndexIoError> {
    let end = offset
        .checked_add(len)
        .ok_or_else(|| IndexIoError::Invalid("slice overflow".to_string()))?;
    let out = bytes
        .get(*offset..end)
        .ok_or_else(|| IndexIoError::Invalid("slice out of bounds".to_string()))?;
    *offset = end;
    Ok(out)
}

fn read_u32(bytes: &[u8], offset: &mut usize) -> Result<u32, IndexIoError> {
    let end = *offset + 4;
    let buf = bytes
        .get(*offset..end)
        .ok_or_else(|| IndexIoError::Invalid("unexpected EOF".to_string()))?;
    *offset = end;
    Ok(u32::from_le_bytes(buf.try_into().unwrap()))
}

fn read_u64(bytes: &[u8], offset: &mut usize) -> Result<u64, IndexIoError> {
    let end = *offset + 8;
    let buf = bytes
        .get(*offset..end)
        .ok_or_else(|| IndexIoError::Invalid("unexpected EOF".to_string()))?;
    *offset = end;
    Ok(u64::from_le_bytes(buf.try_into().unwrap()))
}

fn read_u64_infallible(bytes: &[u8], offset: &mut usize) -> u64 {
    let v = u64::from_le_bytes(bytes[*offset..*offset + 8].try_into().unwrap());
    *offset += 8;
    v
}

fn read_i64_infallible(bytes: &[u8], offset: &mut usize) -> i64 {
    let v = i64::from_le_bytes(bytes[*offset..*offset + 8].try_into().unwrap());
    *offset += 8;
    v
}
