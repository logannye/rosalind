//! The per-position pileup column and its observations — the public substrate
//! type produced by `PileupEngine` and consumed by callers and plugins.

use crate::core::Locus;

/// One base observation at a pileup position. Only callable (A/C/G/T) bases are
/// recorded; `allele` is the 0..=3 index (A=0, C=1, G=2, T=3).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Obs {
    /// Allele index 0..=3 (A/C/G/T).
    pub allele: u8,
    /// Base quality (Phred) of the observed base.
    pub base_qual: u8,
    /// Mapping quality of the read this observation came from.
    pub mapq: u8,
    /// Whether the read maps to the reverse strand (for strand-bias use).
    pub reverse: bool,
}

/// All callable observations stacked at a single reference position.
#[derive(Debug, Clone, PartialEq)]
pub struct PileupColumn {
    /// The reference coordinate of this column.
    pub locus: Locus,
    /// The reference base at this locus (uppercase ASCII; `b'N'` if unknown).
    pub ref_base: u8,
    /// Callable observations, in deterministic (active-read insertion) order.
    pub obs: Vec<Obs>,
}

impl PileupColumn {
    /// Number of callable observations.
    pub fn depth(&self) -> u32 {
        self.obs.len() as u32
    }

    /// Per-allele observation counts, indexed `[A, C, G, T]`.
    pub fn allele_counts(&self) -> [u32; 4] {
        let mut counts = [0u32; 4];
        for o in &self.obs {
            counts[o.allele as usize] += 1;
        }
        counts
    }

    /// Per-allele, per-strand counts: `[allele][0 = forward, 1 = reverse]`.
    pub fn strand_counts(&self) -> [[u32; 2]; 4] {
        let mut counts = [[0u32; 2]; 4];
        for o in &self.obs {
            counts[o.allele as usize][o.reverse as usize] += 1;
        }
        counts
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::Position;

    fn obs(allele: u8, reverse: bool) -> Obs {
        Obs {
            allele,
            base_qual: 30,
            mapq: 60,
            reverse,
        }
    }

    #[test]
    fn depth_allele_and_strand_counts() {
        let col = PileupColumn {
            locus: Locus {
                contig: 0,
                pos: Position(100),
            },
            ref_base: b'A',
            obs: vec![obs(0, false), obs(0, true), obs(1, false)],
        };
        assert_eq!(col.depth(), 3);
        assert_eq!(col.allele_counts(), [2, 1, 0, 0]);
        // [allele][0=fwd,1=rev]: A has 1 fwd + 1 rev, C has 1 fwd.
        let sc = col.strand_counts();
        assert_eq!(sc[0], [1, 1]);
        assert_eq!(sc[1], [1, 0]);
    }
}
