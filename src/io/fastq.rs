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

    /// Read the next line; a clean EOF here is a malformed (truncated) record.
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
}

impl<R: BufRead> Iterator for FastqReader<R> {
    type Item = Result<FastqRecord, CoreError>;

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

    #[test]
    fn quality_line_starting_with_at_is_not_a_header() {
        // Classic FASTQ trap: a quality line may legitimately begin with '@'.
        // The parser reads 4 lines per record positionally, so this must NOT be
        // mistaken for the next record's header.
        let recs = read_all("@r1\nACGT\n+\n@III\n").unwrap();
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].qualities, b"@III");
    }
}
