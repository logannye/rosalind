//! On-disk index format definitions.
//!
//! The format is a single-file, little-endian binary with a fixed header and a
//! sequence of sections. The header is designed to be read without parsing the
//! rest of the file and includes a BLAKE3 checksum of the reference sequence(s)
//! used to build the index.
//!
//! This is intentionally conservative:
//! - little-endian only (explicitly checked)
//! - versioned, with room for forward-compatible flags
//! - offsets/lengths are u64 to support large references

use std::fmt;

/// Magic bytes identifying a Rosalind index file.
pub const ROSALIND_INDEX_MAGIC: [u8; 8] = *b"ROSALIND";

/// Current on-disk index version.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IndexVersion {
    V1 = 1,
}

impl IndexVersion {
    pub fn from_u16(v: u16) -> Option<Self> {
        match v {
            1 => Some(IndexVersion::V1),
            _ => None,
        }
    }
}

/// File header (fixed-size).
///
/// All integers are **little-endian**.
#[derive(Clone, Copy)]
pub struct IndexHeader {
    /// File magic.
    pub magic: [u8; 8],
    /// Format version (see [`IndexVersion`]).
    pub version: u16,
    /// Endian marker (1 = little-endian).
    pub endian: u8,
    /// Reserved for future use.
    pub reserved0: u8,
    /// Format flags (currently 0).
    pub flags: u32,
    /// Number of contigs described by the index.
    pub contig_count: u32,
    /// Suffix-array sample rate for FM-index lookup (e.g., 32 or 64).
    pub sa_sample_rate: u32,
    /// Header size in bytes (fixed for v1).
    pub header_bytes: u64,
    /// Offset to the section table.
    pub section_table_offset: u64,
    /// Size of the section table in bytes.
    pub section_table_bytes: u64,
    /// BLAKE3 hash of the reference payload used to build the index.
    pub reference_blake3: [u8; 32],
}

impl fmt::Debug for IndexHeader {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("IndexHeader")
            .field("magic", &self.magic)
            .field("version", &self.version)
            .field("endian", &self.endian)
            .field("flags", &self.flags)
            .field("contig_count", &self.contig_count)
            .field("sa_sample_rate", &self.sa_sample_rate)
            .field("header_bytes", &self.header_bytes)
            .field("section_table_offset", &self.section_table_offset)
            .field("section_table_bytes", &self.section_table_bytes)
            .field("reference_blake3", &hex32(&self.reference_blake3))
            .finish()
    }
}

impl IndexHeader {
    /// Little-endian marker value.
    pub const ENDIAN_LITTLE: u8 = 1;
    /// Fixed header size in bytes for v1.
    pub const FIXED_SIZE: usize = 8  // magic
        + 2 // version
        + 1 // endian
        + 1 // reserved0
        + 4 // flags
        + 4 // contig_count
        + 4 // sa_sample_rate
        + 8 // header_bytes
        + 8 // section_table_offset
        + 8 // section_table_bytes
        + 32; // reference_blake3

    /// Construct a v1 header with placeholder section table fields (filled in by the writer).
    pub fn new_v1(contig_count: u32, sa_sample_rate: u32, reference_blake3: [u8; 32]) -> Self {
        Self {
            magic: ROSALIND_INDEX_MAGIC,
            version: IndexVersion::V1 as u16,
            endian: Self::ENDIAN_LITTLE,
            reserved0: 0,
            flags: 0,
            contig_count,
            sa_sample_rate: sa_sample_rate.max(1),
            header_bytes: Self::FIXED_SIZE as u64,
            section_table_offset: 0,
            section_table_bytes: 0,
            reference_blake3,
        }
    }
}

fn hex32(bytes: &[u8; 32]) -> String {
    const HEX: &[u8; 16] = b"0123456789abcdef";
    let mut out = vec![0u8; 64];
    for (i, b) in bytes.iter().enumerate() {
        out[2 * i] = HEX[(b >> 4) as usize];
        out[2 * i + 1] = HEX[(b & 0x0f) as usize];
    }
    String::from_utf8_lossy(&out).into_owned()
}

/// A section kind identifier used in the section table.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u32)]
pub enum SectionKind {
    /// Contig metadata table.
    Contigs = 1,
    /// Reference payload (v1 stores raw uppercase ASCII).
    Reference = 2,
    /// Sampled suffix array payload (may be empty in early versions).
    SaSamples = 3,
}

impl SectionKind {
    pub fn from_u32(v: u32) -> Option<Self> {
        match v {
            1 => Some(SectionKind::Contigs),
            2 => Some(SectionKind::Reference),
            3 => Some(SectionKind::SaSamples),
            _ => None,
        }
    }
}

/// One entry in the section table.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SectionEntry {
    pub kind: SectionKind,
    pub offset: u64,
    pub bytes: u64,
}

impl SectionEntry {
    pub const FIXED_SIZE: usize = 4 + 8 + 8;
}
