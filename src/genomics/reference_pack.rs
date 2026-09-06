//! Deterministic, mmap-friendly analysis references (`.rref`).
//!
//! A reference pack deliberately contains no search index. It stores contig
//! metadata and 64-base blocks made from two 2-bit sequence words plus one
//! ambiguity word. FASTA construction is two-pass and line-buffer bounded, so
//! its memory use depends on contig metadata rather than reference length.

use std::collections::HashSet;
use std::fs::File;
use std::io::{self, BufRead, BufWriter, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use blake3::Hasher;
use thiserror::Error;

use crate::core::ContigSet;
use crate::genomics::{IndexReader, ReferenceIndex, ReferenceView};
use crate::io::decompress::open_input;
use crate::util::atomic::AtomicFile;
use crate::util::mmap::MmapReadOnly;

const MAGIC: [u8; 8] = *b"ROSRREF\0";
const VERSION: u16 = 1;
const ENDIAN_LITTLE: u8 = 1;
const HEADER_BYTES: u32 = 128;
const BLOCK_BYTES: u64 = 24;
const BASES_PER_BLOCK: u64 = 64;
const MAX_FASTA_LINE_BYTES: usize = 1 << 20;

/// Errors returned while building or opening an analysis reference pack.
#[derive(Debug, Error)]
pub enum ReferencePackError {
    /// Filesystem or stream failure.
    #[error("io error: {0}")]
    Io(#[from] io::Error),
    /// Malformed or unsupported `.rref` bytes.
    #[error("invalid reference pack: {0}")]
    Invalid(String),
    /// Malformed or unsupported FASTA input.
    #[error("invalid FASTA: {0}")]
    Fasta(String),
    /// Malformed or unsupported legacy index.
    #[error("legacy index error: {0}")]
    Index(#[from] crate::genomics::index::IndexIoError),
}

/// Bounded base access shared by borrowed legacy views and reference packs.
pub trait ReferenceSequence {
    /// Return the normalized base at a concatenated global coordinate.
    fn base_at(&self, global: usize) -> u8;

    /// Total reference length.
    fn len(&self) -> usize;

    /// Whether this reference has no bases.
    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Decode a bounded reference window into a caller-provided buffer.
    fn decode_window(&self, start: usize, end: usize, out: &mut Vec<u8>) {
        out.clear();
        let end = end.min(self.len());
        out.extend((start..end).map(|position| self.base_at(position)));
    }

    /// Decode a bounded reference window directly into shared storage.
    fn decode_window_arc(&self, start: usize, end: usize) -> Arc<[u8]> {
        let end = end.min(self.len());
        (start..end)
            .map(|position| self.base_at(position))
            .collect()
    }
}

impl ReferenceSequence for ReferenceView<'_> {
    fn base_at(&self, global: usize) -> u8 {
        ReferenceView::base_at(self, global)
    }

    fn len(&self) -> usize {
        ReferenceView::len(self)
    }

    fn decode_window(&self, start: usize, end: usize, out: &mut Vec<u8>) {
        ReferenceView::decode_window(self, start, end, out);
    }

    fn decode_window_arc(&self, start: usize, end: usize) -> Arc<[u8]> {
        ReferenceView::decode_window_arc(self, start, end)
    }
}

/// Common identity and coordinate access used by per-locus analyses.
///
/// Search-specific operations intentionally do not appear here. Both the new
/// `.rref` reader and the legacy `.idx` reader implement this interface.
pub trait ReferenceProvider: ReferenceSequence {
    /// Ordered contig metadata defining canonical coordinate order.
    fn contigs(&self) -> &ContigSet;
    /// BLAKE3 of the normalized (uppercase A/C/G/T/N) concatenated reference.
    fn source_reference_blake3(&self) -> [u8; 32];
}

impl ReferenceSequence for ReferenceIndex {
    fn len(&self) -> usize {
        self.contigs().total_length() as usize
    }

    fn base_at(&self, global: usize) -> u8 {
        self.reference_view()
            .expect("an opened index has a validated reference section")
            .base_at(global)
    }

    fn decode_window(&self, start: usize, end: usize, out: &mut Vec<u8>) {
        self.reference_view()
            .expect("an opened index has a validated reference section")
            .decode_window(start, end, out);
    }

    fn decode_window_arc(&self, start: usize, end: usize) -> Arc<[u8]> {
        self.reference_view()
            .expect("an opened index has a validated reference section")
            .decode_window_arc(start, end)
    }
}

impl ReferenceProvider for ReferenceIndex {
    fn contigs(&self) -> &ContigSet {
        self.contigs()
    }

    fn source_reference_blake3(&self) -> [u8; 32] {
        self.header.reference_blake3
    }
}

/// Header and integrity metadata exposed by `reference inspect`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReferencePackMetadata {
    /// On-disk `.rref` format version.
    pub format_version: u16,
    /// Concatenated reference length.
    pub total_bases: u64,
    /// Number of ordered contigs.
    pub contig_count: u32,
    /// BLAKE3 of normalized concatenated reference bases.
    pub source_reference_blake3: [u8; 32],
    /// BLAKE3 of all bytes following the fixed header.
    pub content_blake3: [u8; 32],
}

/// A validated, memory-mapped `.rref` analysis reference.
#[derive(Debug)]
pub struct ReferencePackReader {
    path: PathBuf,
    metadata: ReferencePackMetadata,
    contigs: ContigSet,
    blocks_offset: usize,
    mmap: MmapReadOnly,
}

/// An analysis reference opened from either the preferred `.rref` format or a
/// compatible legacy `.idx` artifact.
#[derive(Debug)]
pub enum AnalysisReference {
    /// Lightweight analysis reference.
    Pack(ReferencePackReader),
    /// Legacy search index used through its embedded reference section.
    LegacyIndex(ReferenceIndex),
}

impl AnalysisReference {
    /// Open by file magic, retaining legacy replay without conflating formats.
    pub fn open(path: impl AsRef<Path>) -> Result<Self, ReferencePackError> {
        let path = path.as_ref();
        let mut magic = [0u8; 8];
        let mut file = File::open(path)?;
        use std::io::Read;
        file.read_exact(&mut magic)?;
        if magic == MAGIC {
            Ok(Self::Pack(ReferencePackReader::open(path)?))
        } else {
            Ok(Self::LegacyIndex(IndexReader::open(path)?))
        }
    }

    /// Validated layout for the evidence pread provider. A None second offset
    /// identifies interleaved 24-byte/64-base rref blocks; Some identifies the
    /// separate 2-bit data and ambiguity arrays in a legacy index.
    pub(crate) fn file_window_layout(&self) -> Result<(u64, Option<u64>), ReferencePackError> {
        match self {
            Self::Pack(reference) => Ok((reference.blocks_offset as u64, None)),
            Self::LegacyIndex(reference) => {
                let (data, ambiguity) = reference.reference_file_offsets()?;
                Ok((data, Some(ambiguity)))
            }
        }
    }

    /// Whether compatibility mode is serving a legacy search index.
    pub fn is_legacy_index(&self) -> bool {
        matches!(self, Self::LegacyIndex(_))
    }

    /// Inspect only the magic prefix to choose the replay flag or artifact role.
    pub fn path_is_pack(path: impl AsRef<Path>) -> Result<bool, ReferencePackError> {
        let mut magic = [0u8; 8];
        let mut file = File::open(path)?;
        use std::io::Read;
        file.read_exact(&mut magic)?;
        Ok(magic == MAGIC)
    }
}

impl ReferenceSequence for AnalysisReference {
    fn base_at(&self, global: usize) -> u8 {
        match self {
            Self::Pack(reference) => reference.base_at(global),
            Self::LegacyIndex(reference) => reference.base_at(global),
        }
    }

    fn len(&self) -> usize {
        match self {
            Self::Pack(reference) => reference.len(),
            Self::LegacyIndex(reference) => reference.len(),
        }
    }

    fn decode_window(&self, start: usize, end: usize, out: &mut Vec<u8>) {
        match self {
            Self::Pack(reference) => reference.decode_window(start, end, out),
            Self::LegacyIndex(reference) => reference.decode_window(start, end, out),
        }
    }

    fn decode_window_arc(&self, start: usize, end: usize) -> Arc<[u8]> {
        match self {
            Self::Pack(reference) => reference.decode_window_arc(start, end),
            Self::LegacyIndex(reference) => reference.decode_window_arc(start, end),
        }
    }
}

impl ReferenceProvider for AnalysisReference {
    fn contigs(&self) -> &ContigSet {
        match self {
            Self::Pack(reference) => reference.contigs(),
            Self::LegacyIndex(reference) => reference.contigs(),
        }
    }

    fn source_reference_blake3(&self) -> [u8; 32] {
        match self {
            Self::Pack(reference) => reference.source_reference_blake3(),
            Self::LegacyIndex(reference) => reference.source_reference_blake3(),
        }
    }
}

impl ReferencePackReader {
    /// Open a reference pack and eagerly validate its layout and content hash.
    pub fn open(path: impl AsRef<Path>) -> Result<Self, ReferencePackError> {
        let path = path.as_ref().to_path_buf();
        let file = File::open(&path)?;
        let mmap = MmapReadOnly::map(&file)?;
        let bytes = mmap.as_bytes();
        if bytes.len() < HEADER_BYTES as usize {
            return Err(ReferencePackError::Invalid("truncated header".into()));
        }
        let mut offset = 0usize;
        if take(bytes, &mut offset, 8)? != MAGIC {
            return Err(ReferencePackError::Invalid("bad magic bytes".into()));
        }
        let format_version = read_u16(bytes, &mut offset)?;
        if format_version != VERSION {
            return Err(ReferencePackError::Invalid(format!(
                "unsupported format version {format_version}"
            )));
        }
        let endian = read_u8(bytes, &mut offset)?;
        let _reserved = read_u8(bytes, &mut offset)?;
        if endian != ENDIAN_LITTLE || cfg!(target_endian = "big") {
            return Err(ReferencePackError::Invalid(
                "reference pack requires a little-endian host".into(),
            ));
        }
        let header_bytes = read_u32(bytes, &mut offset)?;
        if header_bytes != HEADER_BYTES {
            return Err(ReferencePackError::Invalid(format!(
                "unexpected header size {header_bytes}"
            )));
        }
        let contig_count = read_u32(bytes, &mut offset)?;
        let _reserved = read_u32(bytes, &mut offset)?;
        let total_bases = read_u64(bytes, &mut offset)?;
        let contigs_offset = read_u64(bytes, &mut offset)? as usize;
        let contigs_bytes = read_u64(bytes, &mut offset)? as usize;
        let blocks_offset = read_u64(bytes, &mut offset)? as usize;
        let blocks_bytes = read_u64(bytes, &mut offset)? as usize;
        let source_reference_blake3 = read_hash(bytes, &mut offset)?;
        let content_blake3 = read_hash(bytes, &mut offset)?;

        if contig_count == 0 || total_bases == 0 {
            return Err(ReferencePackError::Invalid(
                "reference must contain at least one non-empty contig".into(),
            ));
        }
        if contigs_offset != HEADER_BYTES as usize
            || blocks_offset % 8 != 0
            || contigs_offset.checked_add(contigs_bytes) > Some(blocks_offset)
        {
            return Err(ReferencePackError::Invalid(
                "invalid or unaligned section offsets".into(),
            ));
        }
        let expected_blocks = total_bases.div_ceil(BASES_PER_BLOCK) * BLOCK_BYTES;
        if blocks_bytes as u64 != expected_blocks
            || blocks_offset.checked_add(blocks_bytes) != Some(bytes.len())
        {
            return Err(ReferencePackError::Invalid(
                "block section size does not match reference length".into(),
            ));
        }
        let content_start = contigs_offset;
        // Hash through bounded file I/O rather than touching every mmap page.
        // Sparse evidence queries should not make the complete reference resident
        // merely to validate it before their first planned reference window.
        let mut content_file = file.try_clone()?;
        content_file.seek(SeekFrom::Start(content_start as u64))?;
        let mut content_hasher = Hasher::new();
        let mut content_buffer = [0u8; 64 * 1024];
        loop {
            let count = std::io::Read::read(&mut content_file, &mut content_buffer)?;
            if count == 0 {
                break;
            }
            content_hasher.update(&content_buffer[..count]);
        }
        let actual_content = content_hasher.finalize();
        if actual_content.as_bytes() != &content_blake3 {
            return Err(ReferencePackError::Invalid(
                "content checksum mismatch".into(),
            ));
        }

        let mut contigs = ContigSet::new();
        let mut cursor = contigs_offset;
        let contigs_end = contigs_offset + contigs_bytes;
        for _ in 0..contig_count {
            let name_len = read_u32_bounded(bytes, &mut cursor, contigs_end)? as usize;
            let length = read_u32_bounded(bytes, &mut cursor, contigs_end)?;
            let global_offset = read_u64_bounded(bytes, &mut cursor, contigs_end)?;
            let name_bytes = take_bounded(bytes, &mut cursor, name_len, contigs_end)?;
            let name = std::str::from_utf8(name_bytes)
                .map_err(|_| ReferencePackError::Invalid("contig name is not UTF-8".into()))?;
            let id = contigs.push(name, length);
            if contigs.by_id(id).map(|contig| contig.global_offset) != Some(global_offset) {
                return Err(ReferencePackError::Invalid(
                    "contig global offset mismatch".into(),
                ));
            }
        }
        if cursor != contigs_end || contigs.total_length() != total_bases {
            return Err(ReferencePackError::Invalid(
                "contig table does not match declared size".into(),
            ));
        }

        Ok(Self {
            path,
            metadata: ReferencePackMetadata {
                format_version,
                total_bases,
                contig_count,
                source_reference_blake3,
                content_blake3,
            },
            contigs,
            blocks_offset,
            mmap,
        })
    }

    /// Path backing this memory map.
    pub fn path(&self) -> &Path {
        &self.path
    }

    /// Validated pack metadata and integrity hashes.
    pub fn metadata(&self) -> &ReferencePackMetadata {
        &self.metadata
    }

    fn word(&self, block: usize, within: usize) -> u64 {
        let offset = self.blocks_offset + block * BLOCK_BYTES as usize + within * 8;
        u64::from_le_bytes(
            self.mmap.as_bytes()[offset..offset + 8]
                .try_into()
                .expect("validated block extent"),
        )
    }
}

impl ReferenceSequence for ReferencePackReader {
    fn len(&self) -> usize {
        self.metadata.total_bases as usize
    }

    fn base_at(&self, global: usize) -> u8 {
        debug_assert!(global < self.metadata.total_bases as usize);
        let block = global / BASES_PER_BLOCK as usize;
        let position = global % BASES_PER_BLOCK as usize;
        let ambiguity = self.word(block, 2);
        if ambiguity & (1u64 << position) != 0 {
            return b'N';
        }
        let packed = self.word(block, usize::from(position >= 32));
        let shift = (position % 32) * 2;
        match (packed >> shift) & 0b11 {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        }
    }
}

impl ReferenceProvider for ReferencePackReader {
    fn contigs(&self) -> &ContigSet {
        &self.contigs
    }

    fn source_reference_blake3(&self) -> [u8; 32] {
        self.metadata.source_reference_blake3
    }
}

/// Streaming builder and legacy-index converter for `.rref` files.
#[derive(Debug, Default)]
pub struct ReferencePackBuilder;

impl ReferencePackBuilder {
    /// Build from FASTA using fixed-size input buffering and atomic publication.
    pub fn build(
        fasta: impl AsRef<Path>,
        output: impl AsRef<Path>,
        replace: bool,
    ) -> Result<ReferencePackMetadata, ReferencePackError> {
        let fasta = fasta.as_ref();
        if fasta.as_os_str() == "-" {
            return Err(ReferencePackError::Fasta(
                "reference build requires a seekable path because validation is two-pass".into(),
            ));
        }
        let (contigs, source_hash) = scan_fasta(fasta, |_| Ok(()))?;
        let output = output.as_ref();
        let mut atomic = AtomicFile::create(output)?;
        let metadata = write_pack(atomic.file_mut(), &contigs, source_hash, |emit| {
            let (_, second_hash) = scan_fasta(fasta, emit)?;
            if second_hash != source_hash {
                return Err(ReferencePackError::Fasta(
                    "FASTA changed between the metadata and sequence passes".into(),
                ));
            }
            Ok(())
        })?;
        atomic.commit(replace)?;
        Ok(metadata)
    }

    /// Convert the reference section of a legacy `.idx` without rebuilding its
    /// search structures.
    pub fn convert(
        index: impl AsRef<Path>,
        output: impl AsRef<Path>,
        replace: bool,
    ) -> Result<ReferencePackMetadata, ReferencePackError> {
        let index = IndexReader::open(index)?;
        let source_hash = index.header.reference_blake3;
        let output = output.as_ref();
        let mut atomic = AtomicFile::create(output)?;
        let metadata = write_pack(atomic.file_mut(), index.contigs(), source_hash, |emit| {
            let view = index.reference_view()?;
            for position in 0..view.len() {
                emit(view.base_at(position))?;
            }
            Ok(())
        })?;
        atomic.commit(replace)?;
        Ok(metadata)
    }
}

fn write_pack(
    file: &mut File,
    contigs: &ContigSet,
    source_hash: [u8; 32],
    produce: impl FnOnce(
        &mut dyn FnMut(u8) -> Result<(), ReferencePackError>,
    ) -> Result<(), ReferencePackError>,
) -> Result<ReferencePackMetadata, ReferencePackError> {
    let contig_table = encode_contigs(contigs)?;
    let contigs_offset = HEADER_BYTES as u64;
    let contigs_bytes = contig_table.len() as u64;
    let blocks_offset = align8(contigs_offset + contigs_bytes);
    let total_bases = contigs.total_length();
    let blocks_bytes = total_bases.div_ceil(BASES_PER_BLOCK) * BLOCK_BYTES;

    let mut writer = BufWriter::new(file);
    writer.write_all(&vec![0u8; HEADER_BYTES as usize])?;
    writer.write_all(&contig_table)?;
    writer.write_all(&vec![
        0u8;
        (blocks_offset - contigs_offset - contigs_bytes)
            as usize
    ])?;
    let mut hasher = Hasher::new();
    hasher.update(&contig_table);
    hasher.update(&vec![
        0u8;
        (blocks_offset - contigs_offset - contigs_bytes)
            as usize
    ]);

    let mut block = [0u64; 3];
    let mut position = 0u64;
    let mut emit = |base: u8| -> Result<(), ReferencePackError> {
        if position >= total_bases {
            return Err(ReferencePackError::Invalid(
                "sequence pass produced more bases than metadata pass".into(),
            ));
        }
        let within = (position % BASES_PER_BLOCK) as usize;
        let (code, ambiguous) = encode_base(base)?;
        if ambiguous {
            block[2] |= 1u64 << within;
        }
        block[usize::from(within >= 32)] |= (code as u64) << ((within % 32) * 2);
        position += 1;
        if within == 63 {
            write_block(&mut writer, &mut hasher, block)?;
            block = [0; 3];
        }
        Ok(())
    };
    produce(&mut emit)?;
    if position != total_bases {
        return Err(ReferencePackError::Invalid(
            "sequence pass produced fewer bases than metadata pass".into(),
        ));
    }
    if position % BASES_PER_BLOCK != 0 {
        write_block(&mut writer, &mut hasher, block)?;
    }
    writer.flush()?;
    let content_blake3 = *hasher.finalize().as_bytes();
    let metadata = ReferencePackMetadata {
        format_version: VERSION,
        total_bases,
        contig_count: contigs.len() as u32,
        source_reference_blake3: source_hash,
        content_blake3,
    };

    let file = writer.into_inner().map_err(|error| error.into_error())?;
    file.seek(SeekFrom::Start(0))?;
    write_header(
        file,
        &metadata,
        contigs_offset,
        contigs_bytes,
        blocks_offset,
        blocks_bytes,
    )?;
    Ok(metadata)
}

fn scan_fasta(
    path: &Path,
    mut emit: impl FnMut(u8) -> Result<(), ReferencePackError>,
) -> Result<(ContigSet, [u8; 32]), ReferencePackError> {
    let mut reader = open_input(path)?;
    let mut line = Vec::new();
    let mut contigs = ContigSet::new();
    let mut names = HashSet::new();
    let mut current_name: Option<String> = None;
    let mut current_len = 0u64;
    let mut hash = Hasher::new();

    while read_bounded_line(&mut *reader, &mut line)? {
        let trimmed = trim_ascii(&line);
        if trimmed.is_empty() {
            continue;
        }
        if let Some(header) = trimmed.strip_prefix(b">") {
            if let Some(name) = current_name.take() {
                push_contig(&mut contigs, name, current_len)?;
            }
            let name = header
                .split(|byte| byte.is_ascii_whitespace())
                .next()
                .filter(|name| !name.is_empty())
                .ok_or_else(|| ReferencePackError::Fasta("header has no contig name".into()))?;
            let name = std::str::from_utf8(name)
                .map_err(|_| ReferencePackError::Fasta("contig name is not UTF-8".into()))?
                .to_string();
            if !names.insert(name.clone()) {
                return Err(ReferencePackError::Fasta(format!(
                    "duplicate contig name '{name}'"
                )));
            }
            current_name = Some(name);
            current_len = 0;
            continue;
        }
        if current_name.is_none() {
            return Err(ReferencePackError::Fasta(
                "expected a header starting with '>' before sequence data".into(),
            ));
        }
        for &base in trimmed {
            let normalized = normalize_base(base)?;
            emit(normalized)?;
            hash.update(&[normalized]);
            current_len = current_len
                .checked_add(1)
                .ok_or_else(|| ReferencePackError::Fasta("reference length overflow".into()))?;
        }
    }
    if let Some(name) = current_name {
        push_contig(&mut contigs, name, current_len)?;
    }
    if contigs.is_empty() {
        return Err(ReferencePackError::Fasta(
            "reference contains no FASTA records".into(),
        ));
    }
    Ok((contigs, *hash.finalize().as_bytes()))
}

fn push_contig(
    contigs: &mut ContigSet,
    name: String,
    length: u64,
) -> Result<(), ReferencePackError> {
    if length == 0 {
        return Err(ReferencePackError::Fasta(format!(
            "FASTA record '{name}' has no sequence data"
        )));
    }
    let length = u32::try_from(length)
        .map_err(|_| ReferencePackError::Fasta(format!("contig '{name}' exceeds u32 length")))?;
    contigs.push(name, length);
    Ok(())
}

fn read_bounded_line(reader: &mut dyn BufRead, line: &mut Vec<u8>) -> io::Result<bool> {
    line.clear();
    loop {
        let available = reader.fill_buf()?;
        if available.is_empty() {
            return Ok(!line.is_empty());
        }
        let consumed = available
            .iter()
            .position(|byte| *byte == b'\n')
            .map_or(available.len(), |position| position + 1);
        if line.len().saturating_add(consumed) > MAX_FASTA_LINE_BYTES {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("FASTA line exceeds the {MAX_FASTA_LINE_BYTES}-byte build buffer"),
            ));
        }
        line.extend_from_slice(&available[..consumed]);
        let ended = available.get(consumed.wrapping_sub(1)) == Some(&b'\n');
        reader.consume(consumed);
        if ended {
            return Ok(true);
        }
    }
}

fn trim_ascii(mut bytes: &[u8]) -> &[u8] {
    while bytes.first().is_some_and(u8::is_ascii_whitespace) {
        bytes = &bytes[1..];
    }
    while bytes.last().is_some_and(u8::is_ascii_whitespace) {
        bytes = &bytes[..bytes.len() - 1];
    }
    bytes
}

fn normalize_base(base: u8) -> Result<u8, ReferencePackError> {
    match base.to_ascii_uppercase() {
        b'A' | b'C' | b'G' | b'T' | b'N' => Ok(base.to_ascii_uppercase()),
        b'U' => Ok(b'T'),
        other => Err(ReferencePackError::Fasta(format!(
            "unsupported nucleotide '{}'",
            char::from(other)
        ))),
    }
}

fn encode_base(base: u8) -> Result<(u8, bool), ReferencePackError> {
    match normalize_base(base)? {
        b'A' => Ok((0, false)),
        b'C' => Ok((1, false)),
        b'G' => Ok((2, false)),
        b'T' => Ok((3, false)),
        b'N' => Ok((0, true)),
        _ => unreachable!(),
    }
}

fn encode_contigs(contigs: &ContigSet) -> Result<Vec<u8>, ReferencePackError> {
    let mut bytes = Vec::new();
    for contig in contigs.iter() {
        let name = contig.name.as_bytes();
        let name_len = u32::try_from(name.len())
            .map_err(|_| ReferencePackError::Invalid("contig name exceeds u32".into()))?;
        bytes.extend_from_slice(&name_len.to_le_bytes());
        bytes.extend_from_slice(&contig.length.to_le_bytes());
        bytes.extend_from_slice(&contig.global_offset.to_le_bytes());
        bytes.extend_from_slice(name);
    }
    Ok(bytes)
}

fn write_block(
    writer: &mut impl Write,
    hasher: &mut Hasher,
    block: [u64; 3],
) -> Result<(), ReferencePackError> {
    for word in block {
        let bytes = word.to_le_bytes();
        writer.write_all(&bytes)?;
        hasher.update(&bytes);
    }
    Ok(())
}

fn write_header(
    writer: &mut impl Write,
    metadata: &ReferencePackMetadata,
    contigs_offset: u64,
    contigs_bytes: u64,
    blocks_offset: u64,
    blocks_bytes: u64,
) -> io::Result<()> {
    writer.write_all(&MAGIC)?;
    writer.write_all(&VERSION.to_le_bytes())?;
    writer.write_all(&[ENDIAN_LITTLE, 0])?;
    writer.write_all(&HEADER_BYTES.to_le_bytes())?;
    writer.write_all(&metadata.contig_count.to_le_bytes())?;
    writer.write_all(&0u32.to_le_bytes())?;
    writer.write_all(&metadata.total_bases.to_le_bytes())?;
    writer.write_all(&contigs_offset.to_le_bytes())?;
    writer.write_all(&contigs_bytes.to_le_bytes())?;
    writer.write_all(&blocks_offset.to_le_bytes())?;
    writer.write_all(&blocks_bytes.to_le_bytes())?;
    writer.write_all(&metadata.source_reference_blake3)?;
    writer.write_all(&metadata.content_blake3)?;
    Ok(())
}

fn align8(value: u64) -> u64 {
    value.div_ceil(8) * 8
}

fn take<'a>(
    bytes: &'a [u8],
    offset: &mut usize,
    len: usize,
) -> Result<&'a [u8], ReferencePackError> {
    take_bounded(bytes, offset, len, bytes.len())
}

fn take_bounded<'a>(
    bytes: &'a [u8],
    offset: &mut usize,
    len: usize,
    end: usize,
) -> Result<&'a [u8], ReferencePackError> {
    let next = offset
        .checked_add(len)
        .filter(|next| *next <= end && *next <= bytes.len())
        .ok_or_else(|| ReferencePackError::Invalid("unexpected EOF".into()))?;
    let result = &bytes[*offset..next];
    *offset = next;
    Ok(result)
}

fn read_u8(bytes: &[u8], offset: &mut usize) -> Result<u8, ReferencePackError> {
    Ok(take(bytes, offset, 1)?[0])
}

fn read_u16(bytes: &[u8], offset: &mut usize) -> Result<u16, ReferencePackError> {
    Ok(u16::from_le_bytes(
        take(bytes, offset, 2)?.try_into().unwrap(),
    ))
}

fn read_u32(bytes: &[u8], offset: &mut usize) -> Result<u32, ReferencePackError> {
    Ok(u32::from_le_bytes(
        take(bytes, offset, 4)?.try_into().unwrap(),
    ))
}

fn read_u64(bytes: &[u8], offset: &mut usize) -> Result<u64, ReferencePackError> {
    Ok(u64::from_le_bytes(
        take(bytes, offset, 8)?.try_into().unwrap(),
    ))
}

fn read_u32_bounded(
    bytes: &[u8],
    offset: &mut usize,
    end: usize,
) -> Result<u32, ReferencePackError> {
    Ok(u32::from_le_bytes(
        take_bounded(bytes, offset, 4, end)?.try_into().unwrap(),
    ))
}

fn read_u64_bounded(
    bytes: &[u8],
    offset: &mut usize,
    end: usize,
) -> Result<u64, ReferencePackError> {
    Ok(u64::from_le_bytes(
        take_bounded(bytes, offset, 8, end)?.try_into().unwrap(),
    ))
}

fn read_hash(bytes: &[u8], offset: &mut usize) -> Result<[u8; 32], ReferencePackError> {
    Ok(take(bytes, offset, 32)?.try_into().unwrap())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn paths(name: &str) -> (PathBuf, PathBuf) {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let root = std::env::temp_dir().join(format!("rosalind-rref-{name}-{nonce}"));
        std::fs::create_dir_all(&root).unwrap();
        (root.join("reference.fa"), root.join("reference.rref"))
    }

    #[test]
    fn build_is_deterministic_and_round_trips_contigs_and_ambiguity() {
        let (fasta, pack) = paths("roundtrip");
        std::fs::write(&fasta, b">chr1 description\nacgtn\n>weird:name-2\nTTuN\n").unwrap();
        let first = ReferencePackBuilder::build(&fasta, &pack, false).unwrap();
        let reader = ReferencePackReader::open(&pack).unwrap();
        assert_eq!(reader.contigs().len(), 2);
        assert_eq!(
            reader
                .contigs()
                .by_name("weird:name-2")
                .unwrap()
                .global_offset,
            5
        );
        let mut decoded = Vec::new();
        reader.decode_window(0, reader.len(), &mut decoded);
        assert_eq!(decoded, b"ACGTNTTTN");
        assert_eq!(reader.metadata(), &first);
        let bytes = std::fs::read(&pack).unwrap();
        ReferencePackBuilder::build(&fasta, &pack, true).unwrap();
        assert_eq!(std::fs::read(&pack).unwrap(), bytes);
        std::fs::remove_dir_all(fasta.parent().unwrap()).unwrap();
    }

    #[test]
    fn corruption_and_truncation_are_rejected_at_open() {
        let (fasta, pack) = paths("corrupt");
        std::fs::write(&fasta, b">chr1\nACGTN\n").unwrap();
        ReferencePackBuilder::build(&fasta, &pack, false).unwrap();
        let mut bytes = std::fs::read(&pack).unwrap();
        bytes[HEADER_BYTES as usize] ^= 1;
        let corrupt = pack.with_extension("corrupt.rref");
        std::fs::write(&corrupt, bytes).unwrap();
        assert!(ReferencePackReader::open(corrupt)
            .unwrap_err()
            .to_string()
            .contains("checksum"));
        let truncated = pack.with_extension("truncated.rref");
        let mut bytes = std::fs::read(&pack).unwrap();
        bytes.pop();
        std::fs::write(&truncated, bytes).unwrap();
        assert!(ReferencePackReader::open(truncated).is_err());
        std::fs::remove_dir_all(fasta.parent().unwrap()).unwrap();
    }

    #[test]
    fn invalid_fasta_is_rejected_without_publishing_destination() {
        let (fasta, pack) = paths("invalid");
        std::fs::write(&fasta, b">chr1\nACGX\n").unwrap();
        assert!(ReferencePackBuilder::build(&fasta, &pack, false).is_err());
        assert!(!pack.exists());
        std::fs::remove_dir_all(fasta.parent().unwrap()).unwrap();
    }
}
