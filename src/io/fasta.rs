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
