//! Index reader/writer for Rosalind's on-disk reference index format.
//!
//! Single-file, little-endian, section-based (see `docs/index-format.md`). Every
//! multi-byte array is written as little-endian words at an 8-byte-aligned file
//! offset so the reader can reinterpret mmap bytes as `&[u64]`/`&[u32]` with a
//! checked `slice::align_to` (zero-copy, no rebuild). Serialization is
//! deterministic: fixed section order, fixed-width LE encoding, zero padding, no
//! timestamps ⇒ byte-identical across repeated builds of the same genome.

use std::fs::File;
use std::io::{self, BufWriter, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

use blake3::Hasher;
use thiserror::Error;

use crate::core::ContigSet;
use crate::genomics::index::format::{
    IndexHeader, IndexVersion, SectionEntry, SectionKind, ROSALIND_INDEX_MAGIC,
};
use crate::genomics::{CompressedDNA, GenomeIndex};
use crate::util::mmap::MmapReadOnly;

/// Errors for index IO.
#[derive(Debug, Error)]
pub enum IndexIoError {
    #[error("io error: {0}")]
    Io(#[from] io::Error),

    #[error("invalid index: {0}")]
    Invalid(String),
}

/// A loaded reference index: a memory-mapped byte buffer plus the parsed header,
/// section table, and contig set. The FM-index query surface is obtained as a
/// borrowed [`crate::genomics::FmIndexView`] via [`ReferenceIndex::view`].
#[derive(Debug)]
pub struct ReferenceIndex {
    /// Path to the index file.
    pub path: PathBuf,
    /// Parsed file header.
    pub header: IndexHeader,
    /// Parsed contig set (names, lengths, global offsets).
    pub(crate) contigs: ContigSet,
    /// Parsed section table.
    pub sections: Vec<SectionEntry>,
    pub(crate) mmap: MmapReadOnly,
}

impl ReferenceIndex {
    /// Return the raw bytes of the memory-mapped index file.
    pub fn bytes(&self) -> &[u8] {
        self.mmap.as_bytes()
    }

    /// The parsed contig set.
    pub fn contigs(&self) -> &ContigSet {
        &self.contigs
    }
}

/// Writes a new index file.
#[derive(Debug)]
pub struct IndexWriter {
    file: File,
}

impl IndexWriter {
    /// Create a new index file for writing (overwriting if it exists).
    pub fn create(path: impl AsRef<Path>) -> Result<Self, IndexIoError> {
        let file = File::create(path)?;
        Ok(Self { file })
    }

    /// Serialise a built [`GenomeIndex`] into the deterministic on-disk format.
    ///
    /// Sections are written in a fixed order, each padded to an 8-byte boundary:
    /// `Contigs`, `Reference2bit`, `FmMeta`, `Boundaries`, `Blocks`, `SaSamples`,
    /// then the section table; the header is rewritten last with the section
    /// table's location. `reference_blake3` is the BLAKE3 of the uppercased ASCII
    /// reference (the stable identity used by the B4 stale-index guard).
    pub fn write_genome_index(self, index: &GenomeIndex) -> Result<(), IndexIoError> {
        let fm = index.fm();
        let contigs = index.contigs();
        if contigs.is_empty() {
            return Err(IndexIoError::Invalid(
                "index must contain at least one contig".to_string(),
            ));
        }

        let reference_blake3 = {
            let mut hasher = Hasher::new();
            hasher.update(index.reference());
            *hasher.finalize().as_bytes()
        };
        let sa_sample_rate = u32::try_from(fm.sa_sample_rate())
            .map_err(|_| IndexIoError::Invalid("sa_sample_rate exceeds u32".to_string()))?;
        let contig_count = u32::try_from(contigs.len())
            .map_err(|_| IndexIoError::Invalid("too many contigs (exceeds u32)".to_string()))?;
        let mut header = IndexHeader::new_v1(contig_count, sa_sample_rate, reference_blake3);

        let mut w = BufWriter::new(self.file);

        // Reserve space for the header (rewritten at the end).
        w.write_all(&vec![0u8; IndexHeader::FIXED_SIZE])?;

        let mut sections: Vec<SectionEntry> = Vec::new();

        // --- Contigs (byte-copy read on load; alignment not required, but we pad
        // every section start to 8 uniformly). ---
        let off = pad_to_8(&mut w)?;
        for contig in contigs.iter() {
            let name = contig.name.as_bytes();
            let name_len = u32::try_from(name.len())
                .map_err(|_| IndexIoError::Invalid("contig name too long".to_string()))?;
            w.write_all(&name_len.to_le_bytes())?;
            w.write_all(name)?;
            w.write_all(&(contig.length as u64).to_le_bytes())?;
            w.write_all(&contig.global_offset.to_le_bytes())?;
        }
        push_section(&mut sections, SectionKind::Contigs, off, &mut w)?;

        // --- Reference2bit: len, data_words, amb_words, then data[u64], amb[u64]. ---
        let off = pad_to_8(&mut w)?;
        let reference_2bit = CompressedDNA::compress(index.reference())
            .map_err(|e| IndexIoError::Invalid(format!("reference compress: {e}")))?;
        w.write_all(&(reference_2bit.len() as u64).to_le_bytes())?;
        w.write_all(&(reference_2bit.words().len() as u64).to_le_bytes())?;
        w.write_all(&(reference_2bit.ambiguity().bits().len() as u64).to_le_bytes())?;
        write_u64_slice(&mut w, reference_2bit.words())?;
        write_u64_slice(&mut w, reference_2bit.ambiguity().bits())?;
        push_section(&mut sections, SectionKind::Reference2bit, off, &mut w)?;

        // --- FmMeta: 5 u64 then c_table[6 u32]. ---
        let off = pad_to_8(&mut w)?;
        w.write_all(&(fm.block_size() as u64).to_le_bytes())?;
        w.write_all(&(fm.len() as u64).to_le_bytes())?;
        w.write_all(&(fm.sentinel_position() as u64).to_le_bytes())?;
        w.write_all(&(fm.sa_sample_rate() as u64).to_le_bytes())?;
        w.write_all(&(fm.num_blocks() as u64).to_le_bytes())?;
        for v in *fm.c_table() {
            w.write_all(&v.to_le_bytes())?;
        }
        push_section(&mut sections, SectionKind::FmMeta, off, &mut w)?;

        // --- Boundaries: num_blocks+1 entries of [cumulative_counts; 5] + sentinel. ---
        let off = pad_to_8(&mut w)?;
        debug_assert_eq!(fm.boundaries().len(), fm.num_blocks() + 1);
        for boundary in fm.boundaries().iter() {
            for c in boundary.cumulative_counts {
                w.write_all(&c.to_le_bytes())?;
            }
            w.write_all(&boundary.sentinel_count.to_le_bytes())?;
        }
        push_section(&mut sections, SectionKind::Boundaries, off, &mut w)?;

        // --- Blocks: directory of num_blocks u64 record-offsets, then records. ---
        let off = pad_to_8(&mut w)?;
        write_blocks_section(&mut w, fm)?;
        push_section(&mut sections, SectionKind::Blocks, off, &mut w)?;

        // --- SaSamples: rate, bwt_len, marks_words, sb_len, values_len, then arrays. ---
        let off = pad_to_8(&mut w)?;
        let sampled = fm.sampled();
        w.write_all(&(sampled.rate() as u64).to_le_bytes())?;
        w.write_all(&(sampled.len() as u64).to_le_bytes())?;
        w.write_all(&(sampled.marks().len() as u64).to_le_bytes())?;
        w.write_all(&(sampled.superblocks().len() as u64).to_le_bytes())?;
        w.write_all(&(sampled.values().len() as u64).to_le_bytes())?;
        write_u64_slice(&mut w, sampled.marks())?;
        write_u32_slice(&mut w, sampled.superblocks())?;
        write_u32_slice(&mut w, sampled.values())?;
        push_section(&mut sections, SectionKind::SaSamples, off, &mut w)?;

        // --- Section table, then rewrite the header. ---
        let section_table_offset = pad_to_8(&mut w)?;
        write_section_table(&mut w, &sections)?;
        let section_table_end = w.stream_position()?;

        header.section_table_offset = section_table_offset;
        header.section_table_bytes = section_table_end - section_table_offset;

        w.seek(SeekFrom::Start(0))?;
        write_header(&mut w, &header)?;
        w.flush()?;
        Ok(())
    }
}

/// Pad `w` with zero bytes up to the next 8-byte boundary; return the (aligned)
/// current offset.
fn pad_to_8(w: &mut (impl Write + Seek)) -> Result<u64, IndexIoError> {
    let pos = w.stream_position()?;
    let rem = pos % 8;
    if rem == 0 {
        return Ok(pos);
    }
    let pad = (8 - rem) as usize;
    w.write_all(&[0u8; 8][..pad])?;
    Ok(pos + pad as u64)
}

/// Record a section spanning `[offset, current_position)`.
fn push_section(
    sections: &mut Vec<SectionEntry>,
    kind: SectionKind,
    offset: u64,
    w: &mut (impl Write + Seek),
) -> Result<(), IndexIoError> {
    let end = w.stream_position()?;
    sections.push(SectionEntry {
        kind,
        offset,
        bytes: end - offset,
    });
    Ok(())
}

fn write_u64_slice(w: &mut impl Write, words: &[u64]) -> Result<(), IndexIoError> {
    for &word in words {
        w.write_all(&word.to_le_bytes())?;
    }
    Ok(())
}

fn write_u32_slice(w: &mut impl Write, words: &[u32]) -> Result<(), IndexIoError> {
    for &word in words {
        w.write_all(&word.to_le_bytes())?;
    }
    Ok(())
}

/// Write the `Blocks` section: a `num_blocks × u64` directory of record offsets
/// (relative to the section start, each 8-aligned), then one self-describing
/// record per block. Two passes — sizes first (to lay out the directory), then
/// the records — so the writer streams without materialising the whole section.
fn write_blocks_section(
    w: &mut impl Write,
    fm: &crate::genomics::BlockedFMIndex,
) -> Result<(), IndexIoError> {
    let num_blocks = fm.num_blocks();
    let dir_bytes = (num_blocks as u64) * 8;

    // Pass 1: record offsets (relative to section start).
    let mut offsets = Vec::with_capacity(num_blocks);
    let mut cursor = dir_bytes;
    for block in fm.blocks() {
        offsets.push(cursor);
        cursor += block_record_len(block);
    }

    // Directory.
    for &o in &offsets {
        w.write_all(&o.to_le_bytes())?;
    }

    // Pass 2: records.
    for block in fm.blocks() {
        write_block_record(w, block)?;
    }
    Ok(())
}

/// Byte length of one block record (header + payload), padded to 8.
///
/// All five occ rank bitvectors share one length, and all five superblock arrays
/// share one length (`RankSelectIndex::build` allocates them uniformly). This is
/// load-bearing: the directory offsets computed here must match the bytes
/// `write_block_record` emits, which also use `[0].len()` for all five.
fn block_record_len(block: &crate::genomics::BWTBlock) -> u64 {
    let bwt_data_words = block.bwt().words().len() as u64;
    let bwt_amb_words = block.bwt().ambiguity().bits().len() as u64;
    let occ_bitvec_words = block.occ().bitvectors()[0].len() as u64;
    let occ_superblock_len = block.occ().superblocks()[0].len() as u64;
    debug_assert!(
        block
            .occ()
            .bitvectors()
            .iter()
            .all(|bv| bv.len() as u64 == occ_bitvec_words),
        "occ rank bitvectors must share one length (serialization invariant)"
    );
    debug_assert!(
        block
            .occ()
            .superblocks()
            .iter()
            .all(|sb| sb.len() as u64 == occ_superblock_len),
        "occ superblock arrays must share one length (serialization invariant)"
    );
    let body = 64 // 8 header u64/i64 fields
        + 8 * bwt_data_words
        + 8 * bwt_amb_words
        + 8 * 5 * occ_bitvec_words
        + 4 * 5 * occ_superblock_len
        + 4 * 5; // totals[5]
    (body + 7) / 8 * 8
}

/// Write one block record: 8 header fields (`start,end,sentinel_offset,stride,
/// bwt_data_words,bwt_amb_words,occ_bitvec_words,occ_superblock_len`), then
/// `bwt_data[u64]`, `bwt_amb[u64]`, 5×`occ_bitvec[u64]`, 5×`occ_superblock[u32]`,
/// `totals[u32;5]`, padded to 8.
fn write_block_record(
    w: &mut impl Write,
    block: &crate::genomics::BWTBlock,
) -> Result<(), IndexIoError> {
    let bwt = block.bwt();
    let occ = block.occ();
    let bwt_data_words = bwt.words().len() as u64;
    let bwt_amb_words = bwt.ambiguity().bits().len() as u64;
    let occ_bitvec_words = occ.bitvectors()[0].len() as u64;
    let occ_superblock_len = occ.superblocks()[0].len() as u64;

    let sentinel_offset: i64 = match block.sentinel_offset() {
        Some(o) => o as i64,
        None => -1,
    };

    w.write_all(&(block.start() as u64).to_le_bytes())?;
    w.write_all(&(block.end() as u64).to_le_bytes())?;
    w.write_all(&sentinel_offset.to_le_bytes())?;
    w.write_all(&(occ.stride() as u64).to_le_bytes())?;
    w.write_all(&bwt_data_words.to_le_bytes())?;
    w.write_all(&bwt_amb_words.to_le_bytes())?;
    w.write_all(&occ_bitvec_words.to_le_bytes())?;
    w.write_all(&occ_superblock_len.to_le_bytes())?;

    write_u64_slice(w, bwt.words())?;
    write_u64_slice(w, bwt.ambiguity().bits())?;
    for bv in occ.bitvectors() {
        write_u64_slice(w, bv)?;
    }
    for sb in occ.superblocks() {
        write_u32_slice(w, sb)?;
    }
    write_u32_slice(w, &occ.totals())?;

    // Pad the record to an 8-byte boundary (the u32 arrays may leave a 4-byte tail).
    let body = 64
        + 8 * bwt_data_words
        + 8 * bwt_amb_words
        + 8 * 5 * occ_bitvec_words
        + 4 * 5 * occ_superblock_len
        + 4 * 5;
    let pad = ((body + 7) / 8 * 8 - body) as usize;
    if pad != 0 {
        w.write_all(&[0u8; 8][..pad])?;
    }
    Ok(())
}

/// Reads an existing index file.
#[derive(Debug)]
pub struct IndexReader;

impl IndexReader {
    /// Open and memory-map an existing index file, parsing the header, section
    /// table, and contig set, and validating every section's extent + 8-alignment.
    /// The FM query surface is obtained via [`ReferenceIndex::view`].
    pub fn open(path: impl AsRef<Path>) -> Result<ReferenceIndex, IndexIoError> {
        let path = path.as_ref().to_path_buf();
        let file = File::open(&path)?;
        let mmap = MmapReadOnly::map(&file)?;
        let bytes = mmap.as_bytes();

        let header = read_header(bytes)?;
        validate_header(&header)?;

        let version = IndexVersion::from_u16(header.version).ok_or_else(|| {
            IndexIoError::Invalid(format!("unsupported index version {}", header.version))
        })?;
        if version != IndexVersion::V1 {
            return Err(IndexIoError::Invalid(format!(
                "unsupported index version {version:?}"
            )));
        }

        let sections = read_section_table(bytes, &header)?;
        validate_sections(bytes, &sections)?;
        let contigs = read_contigs(bytes, &header, &sections)?;

        Ok(ReferenceIndex {
            path,
            header,
            contigs,
            sections,
            mmap,
        })
    }
}

/// Validate that every section lies within the file and starts 8-aligned.
fn validate_sections(bytes: &[u8], sections: &[SectionEntry]) -> Result<(), IndexIoError> {
    for s in sections {
        if s.offset % 8 != 0 {
            return Err(IndexIoError::Invalid(format!(
                "section {:?} offset {} is not 8-aligned",
                s.kind, s.offset
            )));
        }
        let end = s
            .offset
            .checked_add(s.bytes)
            .ok_or_else(|| IndexIoError::Invalid("section extent overflow".to_string()))?;
        if end as usize > bytes.len() {
            return Err(IndexIoError::Invalid(format!(
                "section {:?} out of bounds",
                s.kind
            )));
        }
    }
    Ok(())
}

/// Parse the `Contigs` section into a `ContigSet`, recomputing global offsets via
/// `push` and asserting they match the stored offsets (an integrity check).
fn read_contigs(
    bytes: &[u8],
    header: &IndexHeader,
    sections: &[SectionEntry],
) -> Result<ContigSet, IndexIoError> {
    let section = sections
        .iter()
        .find(|s| s.kind == SectionKind::Contigs)
        .ok_or_else(|| IndexIoError::Invalid("missing contigs section".to_string()))?;
    let mut offset = section.offset as usize;
    // `validate_sections` already proved `offset + bytes <= file len`, so this
    // addition cannot overflow on the 64-bit targets Rosalind supports.
    let end = section.offset as usize + section.bytes as usize;

    let mut contigs = ContigSet::new();
    for _ in 0..header.contig_count {
        let name_len = read_u32(bytes, &mut offset)? as usize;
        if offset + name_len > end {
            return Err(IndexIoError::Invalid(
                "contig name out of bounds".to_string(),
            ));
        }
        let name = std::str::from_utf8(&bytes[offset..offset + name_len])
            .map_err(|_| IndexIoError::Invalid("contig name not valid utf-8".to_string()))?
            .to_string();
        offset += name_len;
        // The two u64 fields (length, global_offset) must lie within this section,
        // not spill into the next one (`read_u64` only bounds-checks against the
        // whole file). Reject a record that would cross the section boundary.
        if offset + 16 > end {
            return Err(IndexIoError::Invalid(
                "contig record crosses the section boundary".to_string(),
            ));
        }
        let length = read_u64(bytes, &mut offset)?;
        let stored_global = read_u64(bytes, &mut offset)?;

        let length_u32 = u32::try_from(length)
            .map_err(|_| IndexIoError::Invalid("contig length exceeds u32".to_string()))?;
        let id = contigs.push(name, length_u32);
        let recomputed = contigs.by_id(id).expect("just pushed").global_offset;
        if recomputed != stored_global {
            return Err(IndexIoError::Invalid(format!(
                "contig {id} global_offset {stored_global} != recomputed {recomputed}"
            )));
        }
    }
    Ok(contigs)
}

fn validate_header(header: &IndexHeader) -> Result<(), IndexIoError> {
    if header.magic != ROSALIND_INDEX_MAGIC {
        return Err(IndexIoError::Invalid("bad magic bytes".to_string()));
    }
    if header.endian != IndexHeader::ENDIAN_LITTLE {
        return Err(IndexIoError::Invalid(
            "index endian mismatch (expected little-endian)".to_string(),
        ));
    }
    if header.header_bytes as usize != IndexHeader::FIXED_SIZE {
        return Err(IndexIoError::Invalid(format!(
            "unexpected header size {}",
            header.header_bytes
        )));
    }
    if header.contig_count == 0 {
        return Err(IndexIoError::Invalid(
            "index must contain at least one contig".to_string(),
        ));
    }
    if header.section_table_offset == 0 {
        return Err(IndexIoError::Invalid(
            "missing section table offset".to_string(),
        ));
    }
    Ok(())
}

fn read_header(bytes: &[u8]) -> Result<IndexHeader, IndexIoError> {
    if bytes.len() < IndexHeader::FIXED_SIZE {
        return Err(IndexIoError::Invalid(
            "file too small for header".to_string(),
        ));
    }
    let mut offset = 0usize;
    let mut magic = [0u8; 8];
    magic.copy_from_slice(&bytes[offset..offset + 8]);
    offset += 8;

    let version = read_u16(bytes, &mut offset)?;
    let endian = bytes[offset];
    offset += 1;
    let reserved0 = bytes[offset];
    offset += 1;
    let flags = read_u32(bytes, &mut offset)?;
    let contig_count = read_u32(bytes, &mut offset)?;
    let sa_sample_rate = read_u32(bytes, &mut offset)?;
    let header_bytes = read_u64(bytes, &mut offset)?;
    let section_table_offset = read_u64(bytes, &mut offset)?;
    let section_table_bytes = read_u64(bytes, &mut offset)?;
    let mut reference_blake3 = [0u8; 32];
    reference_blake3.copy_from_slice(&bytes[offset..offset + 32]);

    Ok(IndexHeader {
        magic,
        version,
        endian,
        reserved0,
        flags,
        contig_count,
        sa_sample_rate,
        header_bytes,
        section_table_offset,
        section_table_bytes,
        reference_blake3,
    })
}

fn write_header(mut w: impl Write, header: &IndexHeader) -> Result<(), IndexIoError> {
    w.write_all(&header.magic)?;
    w.write_all(&header.version.to_le_bytes())?;
    w.write_all(&[header.endian])?;
    w.write_all(&[header.reserved0])?;
    w.write_all(&header.flags.to_le_bytes())?;
    w.write_all(&header.contig_count.to_le_bytes())?;
    w.write_all(&header.sa_sample_rate.to_le_bytes())?;
    w.write_all(&header.header_bytes.to_le_bytes())?;
    w.write_all(&header.section_table_offset.to_le_bytes())?;
    w.write_all(&header.section_table_bytes.to_le_bytes())?;
    w.write_all(&header.reference_blake3)?;
    Ok(())
}

fn write_section_table(mut w: impl Write, sections: &[SectionEntry]) -> Result<(), IndexIoError> {
    let count: u32 = sections
        .len()
        .try_into()
        .map_err(|_| IndexIoError::Invalid("too many sections".to_string()))?;
    w.write_all(&count.to_le_bytes())?;
    for section in sections {
        let kind_u32 = section.kind as u32;
        w.write_all(&kind_u32.to_le_bytes())?;
        w.write_all(&section.offset.to_le_bytes())?;
        w.write_all(&section.bytes.to_le_bytes())?;
    }
    Ok(())
}

fn read_section_table(
    bytes: &[u8],
    header: &IndexHeader,
) -> Result<Vec<SectionEntry>, IndexIoError> {
    let mut offset = header.section_table_offset as usize;
    let end = offset
        .checked_add(header.section_table_bytes as usize)
        .ok_or_else(|| IndexIoError::Invalid("section table overflow".to_string()))?;
    if end > bytes.len() {
        return Err(IndexIoError::Invalid(
            "section table out of bounds".to_string(),
        ));
    }
    let count = read_u32(bytes, &mut offset)? as usize;
    let mut sections = Vec::with_capacity(count);
    for _ in 0..count {
        let kind_raw = read_u32(bytes, &mut offset)?;
        let kind = SectionKind::from_u32(kind_raw)
            .ok_or_else(|| IndexIoError::Invalid(format!("unknown section kind {}", kind_raw)))?;
        let sec_offset = read_u64(bytes, &mut offset)?;
        let sec_bytes = read_u64(bytes, &mut offset)?;
        sections.push(SectionEntry {
            kind,
            offset: sec_offset,
            bytes: sec_bytes,
        });
    }
    Ok(sections)
}

fn read_u16(bytes: &[u8], offset: &mut usize) -> Result<u16, IndexIoError> {
    if *offset + 2 > bytes.len() {
        return Err(IndexIoError::Invalid("unexpected EOF".to_string()));
    }
    let mut buf = [0u8; 2];
    buf.copy_from_slice(&bytes[*offset..*offset + 2]);
    *offset += 2;
    Ok(u16::from_le_bytes(buf))
}

fn read_u32(bytes: &[u8], offset: &mut usize) -> Result<u32, IndexIoError> {
    if *offset + 4 > bytes.len() {
        return Err(IndexIoError::Invalid("unexpected EOF".to_string()));
    }
    let mut buf = [0u8; 4];
    buf.copy_from_slice(&bytes[*offset..*offset + 4]);
    *offset += 4;
    Ok(u32::from_le_bytes(buf))
}

fn read_u64(bytes: &[u8], offset: &mut usize) -> Result<u64, IndexIoError> {
    if *offset + 8 > bytes.len() {
        return Err(IndexIoError::Invalid("unexpected EOF".to_string()));
    }
    let mut buf = [0u8; 8];
    buf.copy_from_slice(&bytes[*offset..*offset + 8]);
    *offset += 8;
    Ok(u64::from_le_bytes(buf))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genomics::index::format::SectionKind;
    use crate::genomics::GenomeIndex;
    use std::env;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn temp_path(suffix: &str) -> PathBuf {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        env::temp_dir().join(format!("rosalind-b3b-{suffix}-{nanos}.idx"))
    }

    fn sample_index() -> GenomeIndex {
        GenomeIndex::from_named_sequences(&[
            ("chr1".to_string(), b"ACGTACGTNNACGTACG".to_vec()),
            ("chr2".to_string(), b"TTTTGGGGCCCCAAAANNNN".to_vec()),
        ])
        .expect("build")
    }

    #[test]
    fn serialized_index_is_byte_identical_across_builds() {
        let idx = sample_index();
        let p1 = temp_path("det1");
        let p2 = temp_path("det2");
        IndexWriter::create(&p1)
            .unwrap()
            .write_genome_index(&idx)
            .unwrap();
        IndexWriter::create(&p2)
            .unwrap()
            .write_genome_index(&idx)
            .unwrap();
        let a = std::fs::read(&p1).unwrap();
        let b = std::fs::read(&p2).unwrap();
        assert_eq!(
            a, b,
            "two serializations of the same genome must be byte-identical"
        );
        let _ = std::fs::remove_file(p1);
        let _ = std::fs::remove_file(p2);
    }

    #[test]
    fn open_parses_header_contigs_and_sections() {
        let idx = sample_index();
        let path = temp_path("open");
        IndexWriter::create(&path)
            .unwrap()
            .write_genome_index(&idx)
            .unwrap();

        let loaded = IndexReader::open(&path).expect("open");
        assert_eq!(loaded.header.magic, ROSALIND_INDEX_MAGIC);

        // Contigs round-trip with names, lengths, and global offsets.
        let contigs = loaded.contigs();
        assert_eq!(contigs.len(), 2);
        assert_eq!(contigs.by_name("chr1").unwrap().length, 17);
        assert_eq!(contigs.by_name("chr2").unwrap().global_offset, 17);

        // reference_blake3 matches the uppercased ASCII reference.
        let mut hasher = blake3::Hasher::new();
        hasher.update(idx.reference());
        assert_eq!(
            *hasher.finalize().as_bytes(),
            loaded.header.reference_blake3
        );

        // The FM sections are all present.
        for kind in [
            SectionKind::FmMeta,
            SectionKind::Boundaries,
            SectionKind::Blocks,
            SectionKind::SaSamples,
        ] {
            assert!(
                loaded.sections.iter().any(|s| s.kind == kind),
                "missing {kind:?}"
            );
        }
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn open_rejects_a_truncated_index() {
        let idx = sample_index();
        let path = temp_path("trunc");
        IndexWriter::create(&path)
            .unwrap()
            .write_genome_index(&idx)
            .unwrap();
        let mut bytes = std::fs::read(&path).unwrap();
        bytes.truncate(bytes.len() / 2);
        let truncated = temp_path("trunc-half");
        std::fs::write(&truncated, &bytes).unwrap();
        assert!(IndexReader::open(&truncated).is_err());
        let _ = std::fs::remove_file(path);
        let _ = std::fs::remove_file(truncated);
    }

    #[test]
    fn open_rejects_bad_magic() {
        let idx = sample_index();
        let path = temp_path("magic");
        IndexWriter::create(&path)
            .unwrap()
            .write_genome_index(&idx)
            .unwrap();
        let mut bytes = std::fs::read(&path).unwrap();
        bytes[0] = b'X';
        let corrupt = temp_path("magic-bad");
        std::fs::write(&corrupt, &bytes).unwrap();
        assert!(IndexReader::open(&corrupt).is_err());
        let _ = std::fs::remove_file(path);
        let _ = std::fs::remove_file(corrupt);
    }

    #[test]
    fn serialized_sections_are_present_aligned_and_in_bounds() {
        let idx = sample_index();
        let path = temp_path("sections");
        IndexWriter::create(&path)
            .unwrap()
            .write_genome_index(&idx)
            .unwrap();
        let bytes = std::fs::read(&path).unwrap();

        let header = read_header(&bytes).unwrap();
        validate_header(&header).unwrap();
        let sections = read_section_table(&bytes, &header).unwrap();

        for kind in [
            SectionKind::Contigs,
            SectionKind::Reference2bit,
            SectionKind::FmMeta,
            SectionKind::Boundaries,
            SectionKind::Blocks,
            SectionKind::SaSamples,
        ] {
            let s = sections
                .iter()
                .find(|s| s.kind == kind)
                .unwrap_or_else(|| panic!("missing section {kind:?}"));
            assert_eq!(s.offset % 8, 0, "{kind:?} not 8-aligned");
            let end = s.offset + s.bytes;
            assert!(end as usize <= bytes.len(), "{kind:?} out of bounds");
        }
        let _ = std::fs::remove_file(path);
    }
}
