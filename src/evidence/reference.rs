use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::os::unix::fs::FileExt;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use super::EvidenceError;
use crate::core::ContigSet;
use crate::genomics::{AnalysisReference, ReferenceProvider};

#[derive(Debug)]
struct FaiEntry {
    offset: u64,
    line_bases: u64,
    line_bytes: u64,
}

/// Reference windows from existing .rref/.idx or an uncompressed local FASTA
/// with an existing .fai. FASTA reads do not materialize complete contigs.
#[derive(Debug, Clone)]
pub struct EvidenceReference {
    inner: Arc<ReferenceInner>,
}
#[derive(Debug)]
enum ReferenceInner {
    Packed {
        file: File,
        contigs: ContigSet,
        data_offset: u64,
        ambiguity_offset: Option<u64>,
    },
    Fasta {
        path: PathBuf,
        file: File,
        contigs: ContigSet,
        entries: Vec<FaiEntry>,
    },
    Unavailable(ContigSet),
}
impl EvidenceReference {
    /// Open and validate local inputs, retaining only indexed reference access and metadata.
    pub fn open(path: impl AsRef<Path>) -> Result<Self, EvidenceError> {
        Self::open_with_fai(path, None)
    }
    /// Open a FASTA whose content-addressed FAI was relocated independently.
    pub fn open_with_fai(
        path: impl AsRef<Path>,
        explicit_fai: Option<&Path>,
    ) -> Result<Self, EvidenceError> {
        let path = path.as_ref();
        let mut file = File::open(path)?;
        let mut prefix = [0u8; 8];
        let count = file.read(&mut prefix)?;
        if count > 0 && prefix[0] == b'>' {
            file.seek(SeekFrom::Start(0))?;
            let fai_path = explicit_fai
                .map(PathBuf::from)
                .unwrap_or_else(|| PathBuf::from(format!("{}.fai", path.display())));
            let index = File::open(&fai_path).map_err(|error| {
                EvidenceError::InvalidInput(format!(
                    "local FASTA requires existing {}: {error}",
                    fai_path.display()
                ))
            })?;
            let mut contigs = ContigSet::new();
            let mut entries = Vec::new();
            for (line_number, line) in BufReader::new(index).lines().enumerate() {
                let line = line?;
                let parts: Vec<_> = line.split('\t').collect();
                if parts.len() < 5 {
                    return Err(EvidenceError::InvalidInput(format!(
                        "invalid FAI line {}",
                        line_number + 1
                    )));
                }
                let number = |field: usize| {
                    parts[field].parse::<u64>().map_err(|_| {
                        EvidenceError::InvalidInput(format!(
                            "invalid FAI number on line {}",
                            line_number + 1
                        ))
                    })
                };
                let length = u32::try_from(number(1)?).map_err(|_| {
                    EvidenceError::InvalidInput("FASTA contig exceeds u32 coordinates".into())
                })?;
                let entry = FaiEntry {
                    offset: number(2)?,
                    line_bases: number(3)?,
                    line_bytes: number(4)?,
                };
                if length == 0
                    || entry.line_bases == 0
                    || entry.line_bytes < entry.line_bases
                    || contigs.by_name(parts[0]).is_some()
                {
                    return Err(EvidenceError::InvalidInput(
                        "invalid or duplicate FAI contig".into(),
                    ));
                }
                contigs.push(parts[0], length);
                entries.push(entry);
            }
            if contigs.is_empty() {
                return Err(EvidenceError::InvalidInput("empty FAI".into()));
            }
            Ok(Self {
                inner: Arc::new(ReferenceInner::Fasta {
                    path: path.to_path_buf(),
                    file,
                    contigs,
                    entries,
                }),
            })
        } else if prefix.starts_with(&[0x1f, 0x8b]) {
            Err(EvidenceError::InvalidInput(
                "evidence FASTA currently requires an uncompressed local FASTA with .fai".into(),
            ))
        } else {
            let reference = AnalysisReference::open(path)
                .map_err(|error| EvidenceError::InvalidInput(error.to_string()))?;
            let (data_offset, ambiguity_offset) = reference
                .file_window_layout()
                .map_err(|error| EvidenceError::InvalidInput(error.to_string()))?;
            let contigs = reference.contigs().clone();
            // Release the validated mmap before execution. Keeping it and
            // decoding successive windows through mmap would accumulate RSS
            // with the visited genome span on machines without memory pressure.
            drop(reference);
            Ok(Self {
                inner: Arc::new(ReferenceInner::Packed {
                    file,
                    contigs,
                    data_offset,
                    ambiguity_offset,
                }),
            })
        }
    }
    /// Construct reference metadata with N bases for coverage-only extraction.
    pub fn unavailable(contigs: ContigSet) -> Self {
        Self {
            inner: Arc::new(ReferenceInner::Unavailable(contigs)),
        }
    }
    /// Read the canonical reference contig dictionary.
    pub fn contigs(&self) -> &ContigSet {
        match self.inner.as_ref() {
            ReferenceInner::Packed { contigs, .. } => contigs,
            ReferenceInner::Fasta { contigs, .. } | ReferenceInner::Unavailable(contigs) => contigs,
        }
    }
    /// Return the local FASTA path when this provider uses indexed FASTA.
    pub fn fasta_path(&self) -> Option<&Path> {
        match self.inner.as_ref() {
            ReferenceInner::Fasta { path, .. } => Some(path),
            _ => None,
        }
    }
    /// Whether real reference sequence is available rather than N placeholders.
    pub fn has_sequence(&self) -> bool {
        !matches!(self.inner.as_ref(), ReferenceInner::Unavailable(_))
    }
    /// Fetch an exact half-open reference window using bounded storage.
    pub fn read_window(&self, contig: u32, start: u32, end: u32) -> Result<Vec<u8>, EvidenceError> {
        let metadata = self
            .contigs()
            .by_id(contig)
            .ok_or_else(|| EvidenceError::InvalidInput("unknown reference contig".into()))?;
        if start > end || end > metadata.length {
            return Err(EvidenceError::InvalidInput(
                "reference window out of bounds".into(),
            ));
        }
        let global_offset = metadata.global_offset;
        match self.inner.as_ref() {
            ReferenceInner::Packed {
                file,
                data_offset,
                ambiguity_offset,
                ..
            } => {
                let start = global_offset + u64::from(start);
                let end = global_offset + u64::from(end);
                read_packed_window(file, *data_offset, *ambiguity_offset, start, end)
            }
            ReferenceInner::Unavailable(_) => Ok(vec![b'N'; (end - start) as usize]),
            ReferenceInner::Fasta { file, entries, .. } => {
                let entry = &entries[contig as usize];
                let mut output = vec![0u8; (end - start) as usize];
                let mut position = start as u64;
                let mut written = 0usize;
                while position < end as u64 {
                    let line_position = position % entry.line_bases;
                    let n = (entry.line_bases - line_position).min(end as u64 - position) as usize;
                    let offset = entry
                        .offset
                        .checked_add(
                            (position / entry.line_bases)
                                .checked_mul(entry.line_bytes)
                                .ok_or(EvidenceError::CounterOverflow)?,
                        )
                        .and_then(|value| value.checked_add(line_position))
                        .ok_or(EvidenceError::CounterOverflow)?;
                    file.read_exact_at(&mut output[written..written + n], offset)?;
                    for base in &mut output[written..written + n] {
                        *base = match base.to_ascii_uppercase() {
                            b'A' => b'A',
                            b'C' => b'C',
                            b'G' => b'G',
                            b'T' => b'T',
                            b'N' | b'R' | b'Y' | b'S' | b'W' | b'K' | b'M' | b'B' | b'D' | b'H'
                            | b'V' => b'N',
                            _ => {
                                return Err(EvidenceError::InvalidInput(
                                    "FAI window does not match FASTA sequence layout".into(),
                                ))
                            }
                        };
                    }
                    written += n;
                    position += n as u64;
                }
                Ok(output)
            }
        }
    }
}

fn read_packed_window(
    file: &File,
    data_offset: u64,
    ambiguity_offset: Option<u64>,
    start: u64,
    end: u64,
) -> Result<Vec<u8>, EvidenceError> {
    if start == end {
        return Ok(Vec::new());
    }
    let word = |bytes: &[u8], offset: usize| {
        u64::from_le_bytes(bytes[offset..offset + 8].try_into().unwrap())
    };
    let mut output = Vec::with_capacity((end - start) as usize);
    if let Some(ambiguity_offset) = ambiguity_offset {
        let first_data = start / 32;
        let first_ambiguity = start / 64;
        let mut data = vec![0u8; ((end.div_ceil(32) - first_data) * 8) as usize];
        let mut ambiguity = vec![0u8; ((end.div_ceil(64) - first_ambiguity) * 8) as usize];
        file.read_exact_at(&mut data, data_offset + first_data * 8)?;
        file.read_exact_at(&mut ambiguity, ambiguity_offset + first_ambiguity * 8)?;
        for global in start..end {
            let ambiguous = word(&ambiguity, ((global / 64 - first_ambiguity) * 8) as usize)
                & (1u64 << (global % 64))
                != 0;
            let code =
                (word(&data, ((global / 32 - first_data) * 8) as usize) >> ((global % 32) * 2)) & 3;
            output.push(if ambiguous {
                b'N'
            } else {
                b"ACGT"[code as usize]
            });
        }
    } else {
        let first_block = start / 64;
        let mut packed = vec![0u8; ((end.div_ceil(64) - first_block) * 24) as usize];
        file.read_exact_at(&mut packed, data_offset + first_block * 24)?;
        for global in start..end {
            let block = ((global / 64 - first_block) * 24) as usize;
            let within = global % 64;
            let ambiguous = word(&packed, block + 16) & (1u64 << within) != 0;
            let code = (word(&packed, block + if within >= 32 { 8 } else { 0 })
                >> ((within % 32) * 2))
                & 3;
            output.push(if ambiguous {
                b'N'
            } else {
                b"ACGT"[code as usize]
            });
        }
    }
    Ok(output)
}
