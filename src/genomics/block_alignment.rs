use crate::genomics::{BaseCode, BlockedFMIndex, FmSymbol};
use thiserror::Error;

/// Summary of an FM-index interval `[lower, upper)` in 0-based coordinates.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FMInterval {
    /// Lower bound (inclusive) of the FM interval.
    pub lower: u32,
    /// Upper bound (exclusive) of the FM interval.
    pub upper: u32,
}

impl FMInterval {
    /// Create an interval covering the entire BWT.
    pub fn full(length: usize) -> Self {
        Self {
            lower: 0,
            upper: length as u32,
        }
    }

    /// Interval width (number of matching suffixes).
    pub fn width(&self) -> u32 {
        self.upper.saturating_sub(self.lower)
    }

    /// Returns true when no matches remain.
    pub fn is_empty(&self) -> bool {
        self.width() == 0
    }
}

/// Result emitted by block-level alignment.
#[derive(Debug, Clone)]
pub struct BlockAlignmentSummary {
    /// Identifier of the FM-index block that produced this summary.
    pub block_id: usize,
    /// Range of block identifiers covered by the accumulated summary.
    pub span: (usize, usize),
    /// Number of read bases processed while producing this summary.
    pub processed_bases: usize,
    /// Number of mismatches or early termination events encountered.
    pub mismatches: usize,
    /// Current FM interval after the processed bases.
    pub interval: FMInterval,
    /// Whether the search exhausted all candidates.
    pub exhausted: bool,
    /// Heuristic alignment score for prioritisation.
    pub score: f32,
}

impl BlockAlignmentSummary {
    /// Create an empty summary for a block (used when the segment is empty).
    pub fn empty(block_id: usize, interval: FMInterval) -> Self {
        Self {
            block_id,
            span: (block_id, block_id),
            processed_bases: 0,
            mismatches: 0,
            interval,
            exhausted: false,
            score: 0.0,
        }
    }

    /// Whether the backward search still has candidates remaining.
    pub fn has_candidates(&self) -> bool {
        !self.interval.is_empty()
    }

    /// Merge this summary with another from an adjacent block.
    ///
    /// The right-hand summary is assumed to represent a continuation of the
    /// backward search. The merged summary preserves the newest FM interval
    /// while aggregating metrics.
    pub fn merge_with(&self, rhs: &Self) -> Self {
        Self {
            block_id: rhs.block_id,
            span: (self.span.0.min(rhs.span.0), self.span.1.max(rhs.span.1)),
            processed_bases: self.processed_bases + rhs.processed_bases,
            mismatches: self.mismatches + rhs.mismatches,
            interval: rhs.interval,
            exhausted: self.exhausted || rhs.exhausted,
            score: self.score + rhs.score,
        }
    }
}

/// Errors surfaced during backward search.
#[derive(Debug, Error)]
pub enum AlignmentError {
    /// Input read contained an unsupported base.
    #[error("unsupported base '{base}' at offset {offset}")]
    UnsupportedBase {
        /// Offending base.
        base: char,
        /// Offset within the read.
        offset: usize,
    },

    /// Workspace buffer was too small for the provided read segment.
    #[error("workspace too small: required {required}, provided {provided}")]
    WorkspaceTooSmall {
        /// Bytes required to process the segment.
        required: usize,
        /// Bytes that were actually supplied.
        provided: usize,
    },
}

/// Perform backward search for a read segment using the FM-index, returning a
/// compact block summary.
///
/// The workspace buffer must have capacity ≥ `read_segment.len()`; it is reused
/// to hold an uppercase copy of the segment and avoid allocations.
pub fn align_within_block(
    read_segment: &[u8],
    index: &BlockedFMIndex,
    block_id: usize,
    workspace: &mut [u8],
) -> Result<BlockAlignmentSummary, AlignmentError> {
    if workspace.len() < read_segment.len() {
        return Err(AlignmentError::WorkspaceTooSmall {
            required: read_segment.len(),
            provided: workspace.len(),
        });
    }

    // Copy into workspace in reverse order (tail-first) so we can iterate once.
    for (dst, &src) in workspace.iter_mut().zip(read_segment.iter().rev()) {
        *dst = src.to_ascii_uppercase();
    }

    let mut interval = FMInterval::full(index.len());
    let mut mismatches = 0usize;
    let mut processed = 0usize;
    let mut exhausted = false;

    for (_offset, &symbol) in workspace.iter().take(read_segment.len()).enumerate() {
        let base_code = match BaseCode::from_ascii(symbol) {
            Some(code) => code,
            None => {
                mismatches += 1;
                continue;
            }
        };
        let fm_symbol = FmSymbol::Base(base_code);
        let c_row = index.c_table()[fm_symbol.order()];
        let new_lower = c_row + index.rank(fm_symbol, interval.lower as usize);
        let new_upper = c_row + index.rank(fm_symbol, interval.upper as usize);

        processed += 1;
        interval = FMInterval {
            lower: new_lower,
            upper: new_upper,
        };

        if interval.is_empty() {
            mismatches += 1;
            exhausted = true;
            break;
        }
    }

    let score = (processed as f32) - (mismatches as f32 * 1.5);

    Ok(BlockAlignmentSummary {
        block_id,
        span: (block_id, block_id),
        processed_bases: processed,
        mismatches,
        interval,
        exhausted,
        score,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genomics::BlockedFMIndex;

    #[test]
    fn backward_search_produces_expected_interval() {
        let reference = b"ACGTCGTA";
        let index = BlockedFMIndex::build(reference, 4).unwrap();
        let mut workspace = vec![0u8; reference.len()];

        let summary =
            align_within_block(b"CGTA", &index, 0, &mut workspace).expect("alignment should work");
        assert!(summary.has_candidates());
        assert_eq!(summary.span, (0, 0));
        assert_eq!(summary.processed_bases, 4);
        assert_eq!(summary.mismatches, 0);
        assert!(summary.interval.width() >= 1);
    }

    #[test]
    fn unsupported_base_records_mismatch() {
        let reference = b"AAAA";
        let index = BlockedFMIndex::build(reference, 2).unwrap();
        let mut workspace = vec![0u8; 4];
        let summary =
            align_within_block(b"A?A", &index, 0, &mut workspace).expect("alignment should work");
        assert_eq!(summary.mismatches, 1);
    }

    #[test]
    fn merge_combines_metrics_and_intervals() {
        let reference = b"ACGTACGT";
        let index = BlockedFMIndex::build(reference, 4).unwrap();
        let mut workspace = vec![0u8; 4];

        let left = align_within_block(b"AC", &index, 0, &mut workspace).unwrap();
        let right = align_within_block(b"GT", &index, 1, &mut workspace).unwrap();

        let merged = left.merge_with(&right);
        assert_eq!(
            merged.processed_bases,
            left.processed_bases + right.processed_bases
        );
        assert_eq!(merged.mismatches, left.mismatches + right.mismatches);
        assert_eq!(merged.span.0, left.span.0);
        assert_eq!(merged.span.1, right.span.1);
        assert_eq!(merged.interval, right.interval);
    }
}
