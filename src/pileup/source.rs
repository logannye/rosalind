//! Read sources feed coordinate-sorted reads into the pileup engine. `SliceSource`
//! is an in-memory source (used by the Rust API and tests); a BAM-backed source
//! lands when the calling path is migrated.

use crate::core::{AlignedRead, CoreError};

/// A source of coordinate-sorted aligned reads for the pileup engine.
///
/// Implementations MUST yield reads sorted ascending by `(contig, pos)`.
pub trait ReadSource {
    /// Return the next read, or `None` when exhausted.
    fn next_read(&mut self) -> Result<Option<AlignedRead>, CoreError>;
}

/// An in-memory read source. Sorts the provided reads by `(contig, pos)` on
/// construction so callers need not pre-sort.
#[derive(Debug)]
pub struct SliceSource {
    reads: std::vec::IntoIter<AlignedRead>,
}

impl SliceSource {
    /// Build a source from reads (sorted by `(contig, pos)` here).
    pub fn new(mut reads: Vec<AlignedRead>) -> Self {
        reads.sort_by(|a, b| a.contig.cmp(&b.contig).then(a.pos.0.cmp(&b.pos.0)));
        Self {
            reads: reads.into_iter(),
        }
    }
}

impl ReadSource for SliceSource {
    fn next_read(&mut self) -> Result<Option<AlignedRead>, CoreError> {
        Ok(self.reads.next())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::{CigarOp, CigarOpKind, Position, SamFlags};
    use std::sync::Arc;

    fn read(contig: u32, pos: u32) -> AlignedRead {
        AlignedRead {
            contig,
            pos: Position(pos),
            mapq: 60,
            flags: SamFlags::default(),
            cigar: vec![CigarOp::new(CigarOpKind::Match, 4)],
            seq: Arc::from(b"ACGT".to_vec().into_boxed_slice()),
            qual: Arc::from(vec![30u8; 4].into_boxed_slice()),
        }
    }

    #[test]
    fn slice_source_yields_reads_sorted_by_contig_then_pos() {
        let mut src = SliceSource::new(vec![read(1, 50), read(0, 200), read(0, 100)]);
        let mut order = Vec::new();
        while let Some(r) = src.next_read().unwrap() {
            order.push((r.contig, r.pos.0));
        }
        assert_eq!(order, vec![(0, 100), (0, 200), (1, 50)]);
    }
}
