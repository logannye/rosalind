//! Contig-aware genomic coordinates. Per-contig positions are `u32`; the
//! concatenated-genome offset used by the multi-contig FM-index (Phase B) is
//! 64-bit-safe, so no single coordinate space is capped at 4.29 Gbp.

use std::sync::Arc;

/// A reference contig (chromosome / sequence) with a stable id.
#[derive(Debug, Clone)]
pub struct Contig {
    /// Dense 0-based id (index into the `ContigSet`).
    pub id: u32,
    /// Contig name as it appears in the reference / SAM `@SQ`.
    pub name: Arc<str>,
    /// Length in bases.
    pub length: u32,
    /// 0-based offset of this contig within the concatenated genome
    /// (64-bit-safe; used by the multi-contig index in Phase B).
    pub global_offset: u64,
}

/// A 0-based position within a single contig.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Position(
    /// 0-based offset within a contig.
    pub u32,
);

/// A canonical genomic coordinate: a contig id plus a position.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Locus {
    /// Contig id (index into the originating `ContigSet`).
    pub contig: u32,
    /// 0-based position within the contig.
    pub pos: Position,
}

/// The ordered set of reference contigs: name<->id lookup and global offsets.
#[derive(Debug, Clone, Default)]
pub struct ContigSet {
    contigs: Vec<Contig>,
}

impl ContigSet {
    /// Create an empty contig set.
    pub fn new() -> Self {
        Self {
            contigs: Vec::new(),
        }
    }

    /// Append a contig, assigning the next id and the running global offset.
    /// Returns the new contig's id.
    pub fn push(&mut self, name: impl Into<Arc<str>>, length: u32) -> u32 {
        let id = u32::try_from(self.contigs.len()).expect("contig count exceeds u32::MAX");
        let global_offset = self
            .contigs
            .last()
            .map_or(0, |c| c.global_offset + c.length as u64);
        self.contigs.push(Contig {
            id,
            name: name.into(),
            length,
            global_offset,
        });
        id
    }

    /// Look up a contig by id.
    pub fn by_id(&self, id: u32) -> Option<&Contig> {
        self.contigs.get(id as usize)
    }

    /// Look up a contig by name.
    ///
    /// O(n) in the number of contigs; for repeated lookups, build a name-indexed
    /// map from [`ContigSet::iter`].
    pub fn by_name(&self, name: &str) -> Option<&Contig> {
        self.contigs.iter().find(|c| c.name.as_ref() == name)
    }

    /// Number of contigs.
    pub fn len(&self) -> usize {
        self.contigs.len()
    }

    /// Whether the set is empty.
    pub fn is_empty(&self) -> bool {
        self.contigs.is_empty()
    }

    /// Iterate contigs in id order.
    pub fn iter(&self) -> impl Iterator<Item = &Contig> {
        self.contigs.iter()
    }

    /// Total length of the concatenated genome (64-bit-safe).
    pub fn total_length(&self) -> u64 {
        self.contigs
            .last()
            .map_or(0, |c| c.global_offset + c.length as u64)
    }

    /// Resolve a 0-based offset in the concatenated genome to a [`Locus`].
    /// Returns `None` if `global` lies at or beyond the end of the last contig.
    /// O(log n) in the number of contigs.
    pub fn resolve(&self, global: u64) -> Option<Locus> {
        // INVARIANT: `self.contigs` is sorted by ascending `global_offset`
        // (maintained by `push`); `partition_point` relies on this. A future
        // mutator must preserve it or this binary search breaks.
        // `partition_point` gives the count of contigs whose offset is <= global;
        // the owning contig is the one just before that boundary.
        let after = self.contigs.partition_point(|c| c.global_offset <= global);
        let contig = self.contigs.get(after.checked_sub(1)?)?;
        let pos = global - contig.global_offset;
        if pos < contig.length as u64 {
            Some(Locus {
                contig: contig.id,
                pos: Position(pos as u32),
            })
        } else {
            None
        }
    }

    /// Resolve a `len`-byte span starting at `global`, returning its start
    /// [`Locus`] only if the entire span `[global, global + len)` lies within a
    /// single contig (the bwa "bns" boundary rule). Returns `None` if the span
    /// crosses a contig boundary or runs past the end of the genome.
    pub fn resolve_span(&self, global: u64, len: u32) -> Option<Locus> {
        let start = self.resolve(global)?;
        let contig = self.by_id(start.contig)?;
        let end = global.checked_add(len as u64)?;
        if end <= contig.global_offset + contig.length as u64 {
            Some(start)
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn contig_set_assigns_ids_and_global_offsets() {
        let mut set = ContigSet::new();
        let chr1 = set.push("chr1", 1_000);
        let chr2 = set.push("chr2", 500);

        assert_eq!(chr1, 0);
        assert_eq!(chr2, 1);
        assert_eq!(set.by_id(chr1).unwrap().global_offset, 0);
        assert_eq!(set.by_id(chr2).unwrap().global_offset, 1_000);
        assert_eq!(set.by_name("chr2").unwrap().id, 1);
        assert_eq!(set.total_length(), 1_500);
    }

    #[test]
    fn loci_order_by_contig_then_position() {
        let a = Locus {
            contig: 0,
            pos: Position(100),
        };
        let b = Locus {
            contig: 0,
            pos: Position(200),
        };
        let c = Locus {
            contig: 1,
            pos: Position(0),
        };
        assert!(a < b && b < c);
    }

    #[test]
    fn by_id_out_of_range_is_none() {
        let mut set = ContigSet::new();
        set.push("chr1", 1_000);
        assert!(set.by_id(99).is_none());
    }

    #[test]
    fn empty_contig_set_is_empty_with_zero_total_length() {
        let set = ContigSet::new();
        assert!(set.is_empty());
        assert_eq!(set.len(), 0);
        assert_eq!(set.total_length(), 0);
    }

    #[test]
    fn resolve_maps_global_offset_to_locus() {
        let mut set = ContigSet::new();
        set.push("chr1", 10); // global 0..10
        set.push("chr2", 5); // global 10..15
        assert_eq!(
            set.resolve(0),
            Some(Locus {
                contig: 0,
                pos: Position(0)
            })
        );
        assert_eq!(
            set.resolve(9),
            Some(Locus {
                contig: 0,
                pos: Position(9)
            })
        );
        assert_eq!(
            set.resolve(10),
            Some(Locus {
                contig: 1,
                pos: Position(0)
            })
        );
        assert_eq!(
            set.resolve(14),
            Some(Locus {
                contig: 1,
                pos: Position(4)
            })
        );
    }

    #[test]
    fn resolve_out_of_range_is_none() {
        let mut set = ContigSet::new();
        set.push("chr1", 10);
        assert_eq!(set.resolve(10), None); // exactly past chr1's end (the only contig)
        assert_eq!(set.resolve(100), None);
        assert_eq!(ContigSet::new().resolve(0), None); // empty set
    }

    #[test]
    fn resolve_span_within_a_contig_returns_start() {
        let mut set = ContigSet::new();
        set.push("chr1", 10);
        set.push("chr2", 5);
        assert_eq!(
            set.resolve_span(8, 2),
            Some(Locus {
                contig: 0,
                pos: Position(8)
            })
        ); // [8,10) within chr1
        assert_eq!(
            set.resolve_span(10, 5),
            Some(Locus {
                contig: 1,
                pos: Position(0)
            })
        ); // [10,15) within chr2
    }

    #[test]
    fn resolve_span_crossing_a_boundary_is_rejected() {
        let mut set = ContigSet::new();
        set.push("chr1", 10);
        set.push("chr2", 5);
        assert_eq!(set.resolve_span(8, 4), None); // [8,12) crosses chr1->chr2
        assert_eq!(set.resolve_span(12, 4), None); // [12,16) runs past genome end
    }
}
