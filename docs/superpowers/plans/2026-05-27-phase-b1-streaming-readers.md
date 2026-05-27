# Phase B1 — Streaming FASTA/FASTQ readers + transparent decompression Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Move FASTA/FASTQ parsing out of `main.rs` into a library-first `io/` layer of streaming,
multi-record readers that transparently read plain **or gzip/bgzf** input and accept `-` (stdin),
returning `core`-friendly owned records and typed `CoreError`s.

**Architecture:** Three new modules under `src/io/`: `decompress` (a magic-sniffing `BufRead` wrapper +
`open_input`/`open_output` that handle `-`), `fasta` (a multi-record streaming `FastaReader` iterator),
and `fastq` (a streaming `FastqReader` iterator). `main.rs` keeps thin adapter functions (`read_fasta`
single-record CLI policy, `read_fastq` collect, plus the existing pairing logic) that call the new
library readers — so the binary gains gzip + stdin support and the multi-record capability lives in the
library, while CLI behavior stays single-contig (multi-contig consumption is Phase B4).

**Tech Stack:** Rust 2021; `flate2` (new dep, default `miniz_oxide` pure-Rust backend) for gzip/bgzf
sequential decompression; `std::io::BufRead`; `thiserror` via the existing `core::CoreError`.

This is stage **B1** of the Phase B design spec
(`docs/superpowers/specs/2026-05-27-phase-b-genome-scale-design.md`, §9). It lands green and
independently; it does **not** make alignment/calling multi-contig (that is B2/B4).

---

## File structure

- `Cargo.toml` — add `flate2 = "1.0"`.
- `src/io/decompress.rs` — **Create.** `maybe_decompress<R: BufRead + 'static>(R) -> io::Result<Box<dyn BufRead>>`, `open_input<P: AsRef<Path>>(P)`, `open_output<P: AsRef<Path>>(P)`. No genomics types; pure plumbing.
- `src/io/fasta.rs` — **Create.** `FastaRecord { name, sequence }`, `FastaReader<R: BufRead>: Iterator<Item = Result<FastaRecord, CoreError>>`.
- `src/io/fastq.rs` — **Create.** `FastqRecord { name, sequence, qualities }`, `FastqReader<R: BufRead>: Iterator<Item = Result<FastqRecord, CoreError>>`.
- `src/io/mod.rs` — **Modify.** Add `pub mod decompress; pub mod fasta; pub mod fastq;`.
- `src/main.rs` — **Modify.** Delete the private `struct FastaRecord` / `struct FastqRecord`; import the `io::` ones; replace the bodies of `read_fasta` / `read_fastq` with adapters over the new readers; keep `read_fastq_pairs` / `normalize_read_name` / `FastqPair` / `ResolvedReads`.
- `tests/io_readers.rs` — **Create.** Public-API integration test: a gzipped multi-record FASTA read through `open_input` + `FastaReader` builds a 3-contig `ContigSet`.

---

## Task 1: `flate2` dependency + `io/decompress.rs`

**Files:**
- Modify: `Cargo.toml`
- Create: `src/io/decompress.rs`
- Modify: `src/io/mod.rs`

- [ ] **Step 1: Add the dependency.** In `Cargo.toml`, under `[dependencies]`, add this line immediately after the `blake3` line (keep the existing comment style):

```toml
# Transparent gzip/bgzf-sequential decompression for FASTA/FASTQ ingestion
flate2 = "1.0"
```

- [ ] **Step 2: Create the module with signatures + tests (bodies unimplemented).** Create `src/io/decompress.rs`:

```rust
//! Transparent decompression + stream plumbing for pipe-native IO.
//!
//! Input may be plain text or gzip-compressed; bgzf (the BAM/“block gzip”
//! format) is a series of concatenated gzip members, so a multi-member gzip
//! decoder reads it transparently for sequential access. Random-access bgzf
//! (virtual offsets, `.gzi`/`.csi`) is Phase D and is not handled here.

use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::Path;

use flate2::bufread::MultiGzDecoder;

/// The two-byte gzip magic prefix (`1f 8b`), shared by plain gzip and bgzf.
const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];

/// Wrap a buffered reader, transparently decompressing if its first two bytes
/// are the gzip magic. Plain input passes through unchanged.
pub fn maybe_decompress<R: BufRead + 'static>(reader: R) -> io::Result<Box<dyn BufRead>> {
    unimplemented!()
}

/// Open `path` (or `-` for stdin) as a transparently-decompressed buffered reader.
pub fn open_input<P: AsRef<Path>>(path: P) -> io::Result<Box<dyn BufRead>> {
    unimplemented!()
}

/// Open `path` (or `-` for stdout) as a writer.
pub fn open_output<P: AsRef<Path>>(path: P) -> io::Result<Box<dyn Write>> {
    unimplemented!()
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::io::{Cursor, Read, Write as _};

    fn gzip(bytes: &[u8]) -> Vec<u8> {
        let mut enc = GzEncoder::new(Vec::new(), Compression::default());
        enc.write_all(bytes).unwrap();
        enc.finish().unwrap()
    }

    #[test]
    fn plain_input_passes_through() {
        let mut out = String::new();
        maybe_decompress(Cursor::new(b"plain text".to_vec()))
            .unwrap()
            .read_to_string(&mut out)
            .unwrap();
        assert_eq!(out, "plain text");
    }

    #[test]
    fn gzip_input_is_transparently_decompressed() {
        let gz = gzip(b"hello world");
        let mut out = String::new();
        maybe_decompress(Cursor::new(gz))
            .unwrap()
            .read_to_string(&mut out)
            .unwrap();
        assert_eq!(out, "hello world");
    }

    #[test]
    fn empty_input_is_treated_as_plain_and_yields_nothing() {
        let mut out = String::new();
        maybe_decompress(Cursor::new(Vec::<u8>::new()))
            .unwrap()
            .read_to_string(&mut out)
            .unwrap();
        assert_eq!(out, "");
    }
}
```

- [ ] **Step 3: Register the module.** In `src/io/mod.rs`, add above `pub mod bam;`:

```rust
pub mod decompress;
```

- [ ] **Step 4: Run the tests to verify they fail.**

Run: `cargo test --lib io::decompress 2>&1 | tail -20`
Expected: tests run and FAIL — panics with `not implemented` from `unimplemented!()`.

- [ ] **Step 5: Implement the three functions.** Replace the three `unimplemented!()` bodies in `src/io/decompress.rs`:

```rust
pub fn maybe_decompress<R: BufRead + 'static>(reader: R) -> io::Result<Box<dyn BufRead>> {
    let mut reader = reader;
    let is_gzip = {
        let head = reader.fill_buf()?;
        head.len() >= 2 && head[0] == GZIP_MAGIC[0] && head[1] == GZIP_MAGIC[1]
    };
    if is_gzip {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(reader))))
    } else {
        Ok(Box::new(reader))
    }
}

pub fn open_input<P: AsRef<Path>>(path: P) -> io::Result<Box<dyn BufRead>> {
    let path = path.as_ref();
    if path.as_os_str() == "-" {
        maybe_decompress(BufReader::new(io::stdin()))
    } else {
        maybe_decompress(BufReader::new(File::open(path)?))
    }
}

pub fn open_output<P: AsRef<Path>>(path: P) -> io::Result<Box<dyn Write>> {
    let path = path.as_ref();
    if path.as_os_str() == "-" {
        Ok(Box::new(io::stdout()))
    } else {
        Ok(Box::new(File::create(path)?))
    }
}
```

- [ ] **Step 6: Run the tests to verify they pass.**

Run: `cargo test --lib io::decompress 2>&1 | tail -20`
Expected: `test result: ok. 3 passed`.

- [ ] **Step 7: Commit.**

```bash
git add Cargo.toml Cargo.lock src/io/decompress.rs src/io/mod.rs
git commit -m "feat(io/decompress): transparent gzip/bgzf + stdin/stdout plumbing (open_input/open_output)" \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: `io/fasta.rs` — streaming, multi-record FASTA reader

**Files:**
- Create: `src/io/fasta.rs`
- Modify: `src/io/mod.rs`

- [ ] **Step 1: Create the module with the type, signature, and tests (body unimplemented).** Create `src/io/fasta.rs`:

```rust
//! Streaming, multi-record FASTA reader. Library-first: owned records, typed
//! `CoreError`, no CLI (`anyhow`) or htslib types. One record per `>` header;
//! sequence bytes are uppercased and newline-stripped. Read-length/contig-count
//! agnostic — the multi-contig *consumer* is Phase B4; this reader already
//! yields every record.

use std::io::BufRead;

use crate::core::CoreError;

/// One FASTA record: a contig name and its uppercased sequence bytes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastaRecord {
    /// First whitespace-delimited token of the `>` header line.
    pub name: String,
    /// Sequence bytes, ASCII-uppercased, with line breaks removed.
    pub sequence: Vec<u8>,
}

/// A streaming reader over a FASTA source, yielding one [`FastaRecord`] per header.
#[derive(Debug)]
pub struct FastaReader<R: BufRead> {
    reader: R,
    /// The header (without `>`) of the record that the previous `next()` peeked.
    pending_header: Option<String>,
    line: String,
}

impl<R: BufRead> FastaReader<R> {
    /// Construct a reader over any buffered source.
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            pending_header: None,
            line: String::new(),
        }
    }
}

impl<R: BufRead> Iterator for FastaReader<R> {
    type Item = Result<FastaRecord, CoreError>;

    fn next(&mut self) -> Option<Self::Item> {
        unimplemented!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::ContigSet;
    use std::io::Cursor;

    fn read_all(input: &str) -> Vec<FastaRecord> {
        FastaReader::new(Cursor::new(input.as_bytes().to_vec()))
            .collect::<Result<Vec<_>, _>>()
            .expect("FASTA should parse")
    }

    #[test]
    fn reads_single_record_uppercased() {
        let recs = read_all(">chr1 some description\nacgt\nACGT\n");
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].name, "chr1");
        assert_eq!(recs[0].sequence, b"ACGTACGT");
    }

    #[test]
    fn reads_three_records_and_builds_a_three_contig_set() {
        let recs = read_all(">chr1\nACGT\n>chr2\nAACCGGTT\n>chr3\nGG\n");
        assert_eq!(recs.len(), 3);

        let mut contigs = ContigSet::new();
        for r in &recs {
            contigs.push(r.name.clone(), r.sequence.len() as u32);
        }
        assert_eq!(contigs.len(), 3);
        assert_eq!(contigs.by_name("chr2").unwrap().id, 1);
        // chr1(4) + chr2(8) → chr3 starts at global offset 12.
        assert_eq!(contigs.by_name("chr3").unwrap().global_offset, 12);
    }

    #[test]
    fn skips_blank_lines_between_records() {
        let recs = read_all("\n>chr1\nAC\n\nGT\n\n>chr2\nTT\n");
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].sequence, b"ACGT");
        assert_eq!(recs[1].sequence, b"TT");
    }

    #[test]
    fn missing_header_is_a_malformed_record_error() {
        let err = FastaReader::new(Cursor::new(b"ACGT\n".to_vec()))
            .next()
            .unwrap()
            .unwrap_err();
        assert!(matches!(err, CoreError::MalformedRecord(_)));
    }

    #[test]
    fn header_without_sequence_is_a_malformed_record_error() {
        let err = FastaReader::new(Cursor::new(b">chr1\n".to_vec()))
            .next()
            .unwrap()
            .unwrap_err();
        assert!(matches!(err, CoreError::MalformedRecord(_)));
    }

    #[test]
    fn empty_input_yields_no_records() {
        assert!(read_all("").is_empty());
    }
}
```

- [ ] **Step 2: Register the module.** In `src/io/mod.rs`, add (after `pub mod decompress;`):

```rust
pub mod fasta;
```

- [ ] **Step 3: Run the tests to verify they fail.**

Run: `cargo test --lib io::fasta 2>&1 | tail -20`
Expected: tests run and FAIL — `not implemented` from the `unimplemented!()` `next()`.

- [ ] **Step 4: Implement `next()`.** Replace the `unimplemented!()` in `impl Iterator for FastaReader`:

```rust
    fn next(&mut self) -> Option<Self::Item> {
        // Establish this record's header (from the peeked one, or by scanning).
        let header = match self.pending_header.take() {
            Some(h) => h,
            None => loop {
                self.line.clear();
                match self.reader.read_line(&mut self.line) {
                    Ok(0) => return None, // clean EOF: no more records
                    Ok(_) => {}
                    Err(e) => return Some(Err(CoreError::Io(e))),
                }
                let trimmed = self.line.trim();
                if trimmed.is_empty() {
                    continue;
                }
                match trimmed.strip_prefix('>') {
                    Some(rest) => break rest.to_string(),
                    None => {
                        return Some(Err(CoreError::MalformedRecord(format!(
                            "expected FASTA header starting with '>', found '{trimmed}'"
                        ))))
                    }
                }
            },
        };

        let name = match header.split_whitespace().next() {
            Some(n) => n.to_string(),
            None => {
                return Some(Err(CoreError::MalformedRecord(
                    "FASTA header has no name".to_string(),
                )))
            }
        };

        // Accumulate sequence lines until the next header or EOF.
        let mut sequence = Vec::new();
        loop {
            self.line.clear();
            match self.reader.read_line(&mut self.line) {
                Ok(0) => break, // EOF terminates the final record
                Ok(_) => {}
                Err(e) => return Some(Err(CoreError::Io(e))),
            }
            let trimmed = self.line.trim();
            if trimmed.is_empty() {
                continue;
            }
            if let Some(rest) = trimmed.strip_prefix('>') {
                self.pending_header = Some(rest.to_string());
                break;
            }
            sequence.extend(trimmed.bytes().map(|b| b.to_ascii_uppercase()));
        }

        if sequence.is_empty() {
            return Some(Err(CoreError::MalformedRecord(format!(
                "FASTA record '{name}' has no sequence data"
            ))));
        }

        Some(Ok(FastaRecord { name, sequence }))
    }
```

- [ ] **Step 5: Run the tests to verify they pass.**

Run: `cargo test --lib io::fasta 2>&1 | tail -20`
Expected: `test result: ok. 6 passed`.

- [ ] **Step 6: Commit.**

```bash
git add src/io/fasta.rs src/io/mod.rs
git commit -m "feat(io/fasta): streaming multi-record FASTA reader (library-first, typed errors)" \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 3: `io/fastq.rs` — streaming FASTQ reader

**Files:**
- Create: `src/io/fastq.rs`
- Modify: `src/io/mod.rs`

- [ ] **Step 1: Create the module with the type, signatures, and tests (bodies unimplemented).** Create `src/io/fastq.rs`:

```rust
//! Streaming FASTQ reader. Library-first: owned records, typed `CoreError`, no
//! CLI/htslib types. Four lines per record (`@name`, sequence, `+`, qualities);
//! sequence is uppercased; qualities are kept as raw ASCII (Phred+33) bytes.

use std::io::BufRead;

use crate::core::CoreError;

/// One FASTQ record.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastqRecord {
    /// First whitespace-delimited token of the `@` header line.
    pub name: String,
    /// Sequence bytes, ASCII-uppercased.
    pub sequence: Vec<u8>,
    /// Quality bytes as raw ASCII (Phred+33), same length as `sequence`.
    pub qualities: Vec<u8>,
}

/// A streaming reader over a FASTQ source, yielding one [`FastqRecord`] per 4 lines.
#[derive(Debug)]
pub struct FastqReader<R: BufRead> {
    reader: R,
    line: String,
}

impl<R: BufRead> FastqReader<R> {
    /// Construct a reader over any buffered source.
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            line: String::new(),
        }
    }

    /// Read the next non-blank line, or `None` at clean EOF.
    fn next_nonblank(&mut self) -> Result<Option<String>, CoreError> {
        unimplemented!()
    }

    /// Read the next line; a clean EOF here is a malformed (truncated) record.
    fn next_required(&mut self, ctx: &str) -> Result<String, CoreError> {
        unimplemented!()
    }
}

impl<R: BufRead> Iterator for FastqReader<R> {
    type Item = Result<FastqRecord, CoreError>;

    fn next(&mut self) -> Option<Self::Item> {
        unimplemented!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn read_all(input: &str) -> Result<Vec<FastqRecord>, CoreError> {
        FastqReader::new(Cursor::new(input.as_bytes().to_vec())).collect()
    }

    #[test]
    fn reads_two_records_and_trims_name_after_space() {
        let recs = read_all("@read1 1:N:0:CG\nacgt\n+\nIIII\n@read2\nTT\n+\n##\n").unwrap();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].name, "read1");
        assert_eq!(recs[0].sequence, b"ACGT");
        assert_eq!(recs[0].qualities, b"IIII");
        assert_eq!(recs[1].name, "read2");
    }

    #[test]
    fn missing_at_prefix_is_malformed() {
        let err = read_all("read1\nACGT\n+\nIIII\n").unwrap_err();
        assert!(matches!(err, CoreError::MalformedRecord(_)));
    }

    #[test]
    fn missing_plus_separator_is_malformed() {
        let err = read_all("@read1\nACGT\n-\nIIII\n").unwrap_err();
        assert!(matches!(err, CoreError::MalformedRecord(_)));
    }

    #[test]
    fn seq_qual_length_mismatch_is_malformed() {
        let err = read_all("@read1\nACGT\n+\nII\n").unwrap_err();
        assert!(matches!(err, CoreError::MalformedRecord(_)));
    }

    #[test]
    fn truncated_record_at_eof_is_malformed() {
        let err = read_all("@read1\nACGT\n").unwrap_err();
        assert!(matches!(err, CoreError::MalformedRecord(_)));
    }

    #[test]
    fn empty_input_yields_no_records() {
        assert!(read_all("").unwrap().is_empty());
    }
}
```

- [ ] **Step 2: Register the module.** In `src/io/mod.rs`, add (after `pub mod fasta;`):

```rust
pub mod fastq;
```

- [ ] **Step 3: Run the tests to verify they fail.**

Run: `cargo test --lib io::fastq 2>&1 | tail -20`
Expected: tests run and FAIL — `not implemented`.

- [ ] **Step 4: Implement the helpers and `next()`.** Replace the three `unimplemented!()` bodies:

```rust
    fn next_nonblank(&mut self) -> Result<Option<String>, CoreError> {
        loop {
            self.line.clear();
            let n = self.reader.read_line(&mut self.line)?;
            if n == 0 {
                return Ok(None);
            }
            let trimmed = self.line.trim();
            if !trimmed.is_empty() {
                return Ok(Some(trimmed.to_string()));
            }
        }
    }

    fn next_required(&mut self, ctx: &str) -> Result<String, CoreError> {
        self.line.clear();
        let n = self.reader.read_line(&mut self.line)?;
        if n == 0 {
            return Err(CoreError::MalformedRecord(format!(
                "unexpected end of FASTQ while reading {ctx}"
            )));
        }
        Ok(self.line.trim().to_string())
    }
```

```rust
    fn next(&mut self) -> Option<Self::Item> {
        let header = match self.next_nonblank() {
            Ok(Some(h)) => h,
            Ok(None) => return None, // clean EOF between records
            Err(e) => return Some(Err(e)),
        };
        let name = match header.strip_prefix('@') {
            Some(rest) => match rest.split_whitespace().next() {
                Some(n) => n.to_string(),
                None => {
                    return Some(Err(CoreError::MalformedRecord(
                        "FASTQ header has no read name".to_string(),
                    )))
                }
            },
            None => {
                return Some(Err(CoreError::MalformedRecord(format!(
                    "expected FASTQ header starting with '@', found '{header}'"
                ))))
            }
        };

        let sequence = match self.next_required("sequence") {
            Ok(s) => s.to_ascii_uppercase().into_bytes(),
            Err(e) => return Some(Err(e)),
        };
        let plus = match self.next_required("'+' separator") {
            Ok(s) => s,
            Err(e) => return Some(Err(e)),
        };
        if !plus.starts_with('+') {
            return Some(Err(CoreError::MalformedRecord(format!(
                "expected '+' separator for read '{name}', found '{plus}'"
            ))));
        }
        let qualities = match self.next_required("qualities") {
            Ok(s) => s.into_bytes(),
            Err(e) => return Some(Err(e)),
        };
        if sequence.len() != qualities.len() {
            return Some(Err(CoreError::MalformedRecord(format!(
                "sequence/quality length mismatch for read '{name}' ({} vs {})",
                sequence.len(),
                qualities.len()
            ))));
        }

        Some(Ok(FastqRecord {
            name,
            sequence,
            qualities,
        }))
    }
```

- [ ] **Step 5: Run the tests to verify they pass.**

Run: `cargo test --lib io::fastq 2>&1 | tail -20`
Expected: `test result: ok. 6 passed`.

- [ ] **Step 6: Commit.**

```bash
git add src/io/fastq.rs src/io/mod.rs
git commit -m "feat(io/fastq): streaming FASTQ reader (library-first, typed errors)" \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 4: Rewire `main.rs` onto the library readers (parsing leaves `main.rs`)

**Files:**
- Modify: `src/main.rs`

The goal: delete the private record structs and the hand-rolled parsing; route `read_fasta`/`read_fastq`
through the new library readers (gaining gzip + `-` support); keep CLI behavior single-contig and keep
the pairing logic. The existing `main.rs` tests must still pass unchanged.

- [ ] **Step 1: Import the library record types; delete the private structs.** At the top of
`src/main.rs`, with the other `use rosalind::...` imports, add:

```rust
use rosalind::io::decompress::open_input;
use rosalind::io::fasta::{FastaReader, FastaRecord};
use rosalind::io::fastq::{FastqReader, FastqRecord};
```

Then delete the two private struct definitions (currently around lines 152–162):

```rust
struct FastaRecord {
    name: String,
    sequence: Vec<u8>,
}

#[derive(Debug, Clone)]
struct FastqRecord {
    name: String,
    sequence: Vec<u8>,
    qualities: Vec<u8>,
}
```

(Leave `struct FastqPair`, `enum ResolvedReads`, and everything else in place — they now refer to the
imported `FastqRecord`, which has identical fields.)

- [ ] **Step 2: Replace the body of `read_fasta`.** Replace the entire `fn read_fasta(...) { ... }`
(currently lines ~850–892) with this adapter:

```rust
/// Read a reference FASTA (plain or gzip; `-` = stdin). Phase B1 keeps the
/// single-contig CLI policy: only the first record is used; additional records
/// are warned about (multi-contig consumption is Phase B4). The streaming
/// parser itself lives in `io::fasta`.
fn read_fasta(path: &PathBuf) -> Result<FastaRecord> {
    let reader =
        open_input(path).with_context(|| format!("failed to open {}", path.display()))?;
    let mut records = FastaReader::new(reader);
    let first = records
        .next()
        .ok_or_else(|| anyhow!("FASTA file {} is missing a record", path.display()))?
        .with_context(|| format!("failed to parse FASTA {}", path.display()))?;
    if records.next().is_some() {
        eprintln!(
            "warning: only the first FASTA record is currently used; ignoring the rest \
             (multi-contig lands in Phase B4)"
        );
    }
    Ok(first)
}
```

- [ ] **Step 3: Replace the body of `read_fastq`.** Replace the entire `fn read_fastq(...) { ... }`
(currently lines ~894–956) with this adapter:

```rust
/// Read a FASTQ file (plain or gzip; `-` = stdin) into a vector of records.
/// The streaming parser lives in `io::fastq`.
fn read_fastq(path: &PathBuf) -> Result<Vec<FastqRecord>> {
    let reader = open_input(path)
        .with_context(|| format!("failed to open FASTQ file {}", path.display()))?;
    FastqReader::new(reader)
        .collect::<std::result::Result<Vec<_>, _>>()
        .with_context(|| format!("failed to parse FASTQ {}", path.display()))
}
```

- [ ] **Step 4: Build and clear any newly-unused imports.**

Run: `cargo build 2>&1 | grep -iE 'error|warning'`
Expected: no errors. If the compiler warns that `std::io::BufReader` (or `File`, `bail`) is now unused
*because it was only used by the deleted parsing*, remove just those names from their `use` lines until
the build is warning-clean. Do **not** remove names still used elsewhere in `main.rs`.

- [ ] **Step 5: Run the existing reader tests + full suite to confirm no behavior change.**

Run: `cargo test 2>&1 | tail -30`
Expected: full suite green, including the preserved `fasta_parser_extracts_primary_name`,
`fastq_parser_trims_after_space`, and `align_reads_with_fm_index` (which now construct the imported
`FastqRecord` — identical fields). Report the totals.

- [ ] **Step 6: Smoke-test gzip transparency from the CLI.**

```bash
cargo build --release
gzip -c examples/data/reads.fastq > /tmp/reads.fastq.gz 2>/dev/null || \
  (python3 scripts/generate_toy_data.py /tmp/illumina_toy && gzip -c /tmp/illumina_toy/reads_R1.fastq > /tmp/reads.fastq.gz)
# Align reading the GZIPPED reads; should behave identically to the uncompressed run.
./target/release/rosalind align \
  --reference examples/data/ref.fa \
  --reads /tmp/reads.fastq.gz \
  --format sam | head -5
```
Expected: SAM header + alignment lines (no decompression error). Report the output.

- [ ] **Step 7: Commit.**

```bash
git add src/main.rs
git commit -m "refactor(main): route FASTA/FASTQ reading through io:: readers (adds gzip + stdin)" \
  -m "Deletes the private FastaRecord/FastqRecord structs and the hand-rolled parsers; read_fasta/read_fastq are now thin CLI adapters over io::fasta/io::fastq via io::decompress::open_input. CLI stays single-contig (multi-contig is Phase B4); existing reader tests are preserved." \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 5: Public-API integration test (gz multi-record → `ContigSet`)

**Files:**
- Create: `tests/io_readers.rs`

This is the B1 definition-of-done check, run by `cargo test` in CI: a gzipped, multi-record FASTA, read
through the public `open_input` + `FastaReader` API, builds the expected multi-contig `ContigSet`.

- [ ] **Step 1: Write the integration test.** Create `tests/io_readers.rs`:

```rust
//! Phase B1 DoD: the public reader API reads gzipped, multi-record FASTA and
//! produces a multi-contig ContigSet.

use std::io::{Cursor, Write};

use flate2::write::GzEncoder;
use flate2::Compression;
use rosalind::core::ContigSet;
use rosalind::io::decompress::maybe_decompress;
use rosalind::io::fasta::FastaReader;

fn gzip(bytes: &[u8]) -> Vec<u8> {
    let mut enc = GzEncoder::new(Vec::new(), Compression::default());
    enc.write_all(bytes).unwrap();
    enc.finish().unwrap()
}

#[test]
fn gzipped_multi_record_fasta_builds_a_contig_set() {
    let fasta = b">chr1 first\nACGTACGT\n>chr2\nAACC\n>chr3\nGGGGTTTT\n";
    let gz = gzip(fasta);

    let reader = maybe_decompress(Cursor::new(gz)).unwrap();
    let records: Vec<_> = FastaReader::new(reader)
        .collect::<Result<Vec<_>, _>>()
        .expect("gz FASTA should parse");

    assert_eq!(records.len(), 3);

    let mut contigs = ContigSet::new();
    for r in &records {
        contigs.push(r.name.clone(), r.sequence.len() as u32);
    }
    assert_eq!(contigs.len(), 3);
    assert_eq!(contigs.by_name("chr1").unwrap().global_offset, 0);
    assert_eq!(contigs.by_name("chr2").unwrap().global_offset, 8);
    assert_eq!(contigs.by_name("chr3").unwrap().global_offset, 12);
    assert_eq!(contigs.total_length(), 20);
}
```

- [ ] **Step 2: Run the integration test.**

Run: `cargo test --test io_readers 2>&1 | tail -20`
Expected: `test result: ok. 1 passed`.

- [ ] **Step 3: Confirm the whole suite is still green.**

Run: `cargo test 2>&1 | tail -15` then `cargo fmt --all -- --check`
Expected: full suite green; formatting clean.

- [ ] **Step 4: Commit.**

```bash
git add tests/io_readers.rs
git commit -m "test(io): gzipped multi-record FASTA → multi-contig ContigSet (B1 DoD)" \
  -m "Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Final verification (before opening the B1 PR)

- `cargo test` — full suite green (report totals).
- `cargo build 2>&1 | grep -iE 'error|warning'` — clean.
- `cargo fmt --all -- --check` — clean.
- `grep -n "fn read_to_string\|\.lines()\|starts_with('>')\|starts_with('@')" src/main.rs` — should be
  empty (no FASTA/FASTQ *parsing* remains in `main.rs`; only the thin adapters that call `io::`).
- CLI smoke: `align` reads a `.gz` reads file identically to the uncompressed file (Task 4 Step 6).

## Self-Review

- **Spec coverage (B1 row of §9):** streaming multi-record readers in `io/` ✔ (Tasks 2, 3); gz/bgzf
  transparent ✔ (Task 1, `MultiGzDecoder` handles concatenated members); `-`/stdin ✔ (`open_input`);
  library-first / no htslib/anyhow in `io/` signatures ✔ (readers return `CoreError`); parsing leaves
  `main.rs` ✔ (Task 4 + the grep gate); 3-record FASTA → 3-contig `ContigSet` ✔ (Task 2 + Task 5);
  malformed-input typed errors ✔ (Tasks 2, 3).
- **Out of scope confirmed:** alignment/calling stay single-contig (Task 4 keeps the first-record
  policy + warning); no FM-index, persistence, or `@SQ`/`##contig` changes here.
- **Type consistency:** `FastaRecord { name, sequence }` and `FastqRecord { name, sequence, qualities }`
  match the field names `main.rs` already uses (`.name`, `.sequence`, `.qualities`), so `FastqPair`,
  `ResolvedReads`, `read_fastq_pairs`, `normalize_read_name`, and `align_reads_with_fm_index` compile
  unchanged. `open_input`/`maybe_decompress`/`open_output` names are used consistently across tasks.
- **No placeholders:** every step ships complete code or an exact command + expected output.
