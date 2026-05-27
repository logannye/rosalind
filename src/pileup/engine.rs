//! The streaming pileup engine: CIGAR-aware, read-filtered, strand-aware,
//! bounded-memory. Yields one `PileupColumn` per covered reference position.
//!
//! (The engine's imports are added in Task 4, alongside the engine itself.)

/// Read-level filters applied as reads enter the pileup.
#[derive(Debug, Clone)]
pub struct PileupParams {
    /// Minimum mapping quality; reads below this are skipped.
    pub min_mapq: u8,
    /// Minimum base quality; observations below this are dropped.
    pub min_base_qual: u8,
    /// Skip secondary alignments (SAM flag 0x100).
    pub skip_secondary: bool,
    /// Skip supplementary alignments (SAM flag 0x800).
    pub skip_supplementary: bool,
    /// Skip PCR/optical duplicates (SAM flag 0x400).
    pub skip_duplicate: bool,
}

impl Default for PileupParams {
    fn default() -> Self {
        Self {
            min_mapq: 0,
            min_base_qual: 0,
            skip_secondary: true,
            skip_supplementary: true,
            skip_duplicate: true,
        }
    }
}

/// Counts of reads skipped during pileup, by reason (surfaced to callers/CLI).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SkipCounts {
    /// Reads with the unmapped flag set.
    pub unmapped: u64,
    /// Reads on a contig other than the target.
    pub wrong_contig: u64,
    /// Secondary alignments skipped.
    pub secondary: u64,
    /// Supplementary alignments skipped.
    pub supplementary: u64,
    /// Duplicate reads skipped.
    pub duplicate: u64,
    /// Reads below the MAPQ threshold.
    pub low_mapq: u64,
}

impl SkipCounts {
    /// Total reads skipped across all reasons.
    pub fn total(&self) -> u64 {
        self.unmapped
            + self.wrong_contig
            + self.secondary
            + self.supplementary
            + self.duplicate
            + self.low_mapq
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_params_skip_noise_and_keep_quality_open() {
        let p = PileupParams::default();
        assert!(p.skip_secondary && p.skip_supplementary && p.skip_duplicate);
        assert_eq!(p.min_mapq, 0);
        assert_eq!(p.min_base_qual, 0);
    }

    #[test]
    fn skip_counts_total_sums_all_reasons() {
        let s = SkipCounts {
            unmapped: 1,
            wrong_contig: 2,
            secondary: 3,
            supplementary: 4,
            duplicate: 5,
            low_mapq: 6,
        };
        assert_eq!(s.total(), 21);
    }
}
