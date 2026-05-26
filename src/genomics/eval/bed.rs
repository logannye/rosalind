use std::collections::BTreeMap;
use std::ops::Range;

use thiserror::Error;

/// Errors while parsing a BED file.
#[derive(Debug, Error)]
pub enum BedParseError {
    /// A line could not be parsed.
    #[error("invalid BED line {line}: {msg}")]
    InvalidLine {
        /// 1-based line number in the input.
        line: usize,
        /// Human-readable parse error.
        msg: String,
    },
}

/// A simple per-contig interval index for BED ranges.
///
/// Assumes 0-based half-open intervals `[start, end)`.
#[derive(Debug, Clone)]
pub struct BedIndex {
    intervals: BTreeMap<String, Vec<Range<u32>>>,
}

impl BedIndex {
    /// Create an empty BED index (matches nothing).
    pub fn empty() -> Self {
        Self {
            intervals: BTreeMap::new(),
        }
    }

    /// Parse and index a BED file given as a string.
    pub fn from_str(contents: &str) -> Result<Self, BedParseError> {
        let mut per_chrom: BTreeMap<String, Vec<Range<u32>>> = BTreeMap::new();

        for (i, line) in contents.lines().enumerate() {
            let line_no = i + 1;
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            let parts: Vec<&str> = trimmed.split_whitespace().collect();
            if parts.len() < 3 {
                return Err(BedParseError::InvalidLine {
                    line: line_no,
                    msg: "expected at least 3 columns".to_string(),
                });
            }
            let chrom = parts[0].to_string();
            let start: u32 = parts[1].parse().map_err(|_| BedParseError::InvalidLine {
                line: line_no,
                msg: "start is not a u32".to_string(),
            })?;
            let end: u32 = parts[2].parse().map_err(|_| BedParseError::InvalidLine {
                line: line_no,
                msg: "end is not a u32".to_string(),
            })?;
            if end < start {
                return Err(BedParseError::InvalidLine {
                    line: line_no,
                    msg: "end < start".to_string(),
                });
            }
            per_chrom.entry(chrom).or_default().push(start..end);
        }

        // Sort and merge overlapping/adjacent intervals deterministically.
        for intervals in per_chrom.values_mut() {
            intervals.sort_by(|a, b| a.start.cmp(&b.start).then_with(|| a.end.cmp(&b.end)));
            let mut merged: Vec<Range<u32>> = Vec::new();
            for r in intervals.drain(..) {
                match merged.last_mut() {
                    None => merged.push(r),
                    Some(last) => {
                        if r.start <= last.end {
                            last.end = last.end.max(r.end);
                        } else {
                            merged.push(r);
                        }
                    }
                }
            }
            *intervals = merged;
        }

        Ok(Self { intervals: per_chrom })
    }

    /// Returns true if `(chrom, pos0)` is inside any interval.
    pub fn contains(&self, chrom: &str, pos0: u32) -> bool {
        let Some(v) = self.intervals.get(chrom) else {
            return false;
        };
        // Binary search by range start.
        let idx = match v.binary_search_by(|r| {
            if pos0 < r.start {
                std::cmp::Ordering::Greater
            } else if pos0 >= r.end {
                std::cmp::Ordering::Less
            } else {
                std::cmp::Ordering::Equal
            }
        }) {
            Ok(i) => i,
            Err(i) => i,
        };
        v.get(idx)
            .is_some_and(|r| pos0 >= r.start && pos0 < r.end)
    }
}


