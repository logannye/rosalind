//! Index reader/writer for Rosalind’s on-disk reference index format.

use std::fs::File;
use std::io::{self, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

use blake3::Hasher;
use thiserror::Error;

use crate::genomics::index::format::{
    IndexHeader, IndexVersion, SectionEntry, SectionKind, ROSALIND_INDEX_MAGIC,
};
use crate::util::mmap::MmapReadOnly;

/// Errors for index IO.
#[derive(Debug, Error)]
pub enum IndexIoError {
    #[error("io error: {0}")]
    Io(#[from] io::Error),

    #[error("invalid index: {0}")]
    Invalid(String),
}

/// Minimal contig metadata for the index.
#[derive(Debug, Clone)]
pub struct ContigInfo {
    /// Contig name (e.g., `chr1`).
    pub name: String,
    /// Contig length in bases.
    pub length: u64,
}

/// A loaded reference index.
///
/// For now this is a “container” around a memory-mapped byte buffer plus parsed
/// metadata. Subsequent todos (`sa-construction`, `fm-rank-structure`) will add
/// the actual FM-index sections and structured accessors.
#[derive(Debug)]
pub struct ReferenceIndex {
    /// Path to the index file.
    pub path: PathBuf,
    /// Parsed file header.
    pub header: IndexHeader,
    /// Parsed contig table.
    pub contigs: Vec<ContigInfo>,
    /// Parsed section table.
    pub sections: Vec<SectionEntry>,
    mmap: MmapReadOnly,
}

impl ReferenceIndex {
    /// Return the raw bytes of the memory-mapped index file.
    pub fn bytes(&self) -> &[u8] {
        self.mmap.as_bytes()
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

    /// Write a minimal v1 index.
    ///
    /// Current contents:
    /// - header
    /// - sections: contigs + reference payload + (empty) SA samples
    /// - section table at end
    ///
    /// Note: This is a scaffold. Later milestones will replace `reference_payload`
    /// with compressed reference storage and add FM-index sections.
    pub fn write_v1(
        mut self,
        contigs: &[ContigInfo],
        reference_payload: &[u8],
        sa_sample_rate: u32,
    ) -> Result<(), IndexIoError> {
        if contigs.is_empty() {
            return Err(IndexIoError::Invalid(
                "index must contain at least one contig".to_string(),
            ));
        }

        let mut hasher = Hasher::new();
        hasher.update(reference_payload);
        let reference_blake3 = *hasher.finalize().as_bytes();

        let mut header =
            IndexHeader::new_v1(contigs.len() as u32, sa_sample_rate, reference_blake3);

        // Reserve space for header.
        self.file.seek(SeekFrom::Start(0))?;
        self.file.write_all(&vec![0u8; IndexHeader::FIXED_SIZE])?;

        let mut sections: Vec<SectionEntry> = Vec::new();

        // Section: contigs
        let contigs_offset = self.file.stream_position()?;
        for contig in contigs {
            let name_bytes = contig.name.as_bytes();
            let name_len: u32 = name_bytes
                .len()
                .try_into()
                .map_err(|_| IndexIoError::Invalid("contig name too long".to_string()))?;
            self.file.write_all(&name_len.to_le_bytes())?;
            self.file.write_all(name_bytes)?;
            self.file.write_all(&contig.length.to_le_bytes())?;
        }
        let contigs_end = self.file.stream_position()?;
        sections.push(SectionEntry {
            kind: SectionKind::Contigs,
            offset: contigs_offset,
            bytes: contigs_end - contigs_offset,
        });

        // Section: reference payload (for now).
        let reference_offset = self.file.stream_position()?;
        self.file.write_all(reference_payload)?;
        let reference_end = self.file.stream_position()?;
        sections.push(SectionEntry {
            kind: SectionKind::Reference,
            offset: reference_offset,
            bytes: reference_end - reference_offset,
        });

        // Section: SA samples.
        //
        // Payload format (v1):
        // - u64 sample_count
        // - sample_count × u64 (suffix array positions)
        let sa_offset = self.file.stream_position()?;
        self.file.write_all(&0u64.to_le_bytes())?;
        let sa_end = self.file.stream_position()?;
        sections.push(SectionEntry {
            kind: SectionKind::SaSamples,
            offset: sa_offset,
            bytes: sa_end - sa_offset,
        });

        // Section table.
        let section_table_offset = self.file.stream_position()?;
        write_section_table(&mut self.file, &sections)?;
        let section_table_end = self.file.stream_position()?;

        self.file.flush()?;

        header.section_table_offset = section_table_offset;
        header.section_table_bytes = section_table_end - section_table_offset;

        // Write header.
        self.file.seek(SeekFrom::Start(0))?;
        write_header(&mut self.file, &header)?;

        Ok(())
    }
}

/// Reads an existing index file.
#[derive(Debug)]
pub struct IndexReader;

impl IndexReader {
    /// Open and memory-map an existing index file.
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
                "unsupported index version {:?}",
                version
            )));
        }

        let sections = read_section_table(bytes, &header)?;
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

fn read_contigs(
    bytes: &[u8],
    header: &IndexHeader,
    sections: &[SectionEntry],
) -> Result<Vec<ContigInfo>, IndexIoError> {
    let contig_section = sections
        .iter()
        .find(|s| s.kind == SectionKind::Contigs)
        .ok_or_else(|| IndexIoError::Invalid("missing contigs section".to_string()))?;
    let mut offset = contig_section.offset as usize;
    let end = offset
        .checked_add(contig_section.bytes as usize)
        .ok_or_else(|| IndexIoError::Invalid("contigs section overflow".to_string()))?;
    if end > bytes.len() {
        return Err(IndexIoError::Invalid(
            "contigs section out of bounds".to_string(),
        ));
    }

    let mut contigs = Vec::with_capacity(header.contig_count as usize);
    for _ in 0..header.contig_count {
        let name_len = read_u32(bytes, &mut offset)?;
        let name_len_usize: usize = name_len as usize;
        if offset + name_len_usize > end {
            return Err(IndexIoError::Invalid(
                "contig name out of bounds".to_string(),
            ));
        }
        let name = std::str::from_utf8(&bytes[offset..offset + name_len_usize])
            .map_err(|_| IndexIoError::Invalid("contig name not valid utf-8".to_string()))?
            .to_string();
        offset += name_len_usize;
        let length = read_u64(bytes, &mut offset)?;
        contigs.push(ContigInfo { name, length });
    }

    Ok(contigs)
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
    use std::env;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn temp_path(suffix: &str) -> PathBuf {
        let timestamp = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("time went backwards")
            .as_nanos();
        env::temp_dir().join(format!("rosalind-index-{suffix}-{timestamp}.idx"))
    }

    #[test]
    fn roundtrip_minimal_v1_index() {
        let path = temp_path("roundtrip");
        let contigs = vec![ContigInfo {
            name: "chr1".to_string(),
            length: 8,
        }];
        let reference_payload = b"ACGTACGT";
        IndexWriter::create(&path)
            .expect("create")
            .write_v1(&contigs, reference_payload, 32)
            .expect("write");

        let loaded = IndexReader::open(&path).expect("open");
        assert_eq!(loaded.header.magic, ROSALIND_INDEX_MAGIC);
        assert_eq!(loaded.contigs.len(), 1);
        assert_eq!(loaded.contigs[0].name, "chr1");
        assert_eq!(loaded.contigs[0].length, 8);
        assert!(loaded
            .sections
            .iter()
            .any(|s| s.kind == SectionKind::SaSamples));

        // Verify checksum matches the reference payload.
        let mut hasher = Hasher::new();
        hasher.update(reference_payload);
        assert_eq!(
            *hasher.finalize().as_bytes(),
            loaded.header.reference_blake3
        );

        let _ = std::fs::remove_file(path);
    }
}
