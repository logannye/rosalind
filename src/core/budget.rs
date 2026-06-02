//! Memory-as-a-contract primitives. Streaming stages report a `WorkingSet`
//! bound so a run can be checked against a `MemoryBudget` *before* it starts
//! (the foundation for `rosalind plan`).

/// Per-base cost of a pileup read's reference→read-offset projection map (one
/// `HashMap<u32, usize>` entry ≈ 16 bytes). Shared by the realized accountant
/// (`PileupEngine::current_working_set`) and the `rosalind plan` estimator so the
/// two cannot drift.
pub const PILEUP_MAP_BYTES_PER_BASE: u64 = 16;
/// Per-base cost of a read's `seq` + `qual` byte buffers (1 byte each).
pub const PILEUP_SEQQUAL_BYTES_PER_BASE: u64 = 2;
/// Fixed per-active-read overhead (handles + struct).
pub const PILEUP_PER_READ_OVERHEAD: u64 = 64;
/// Fixed per-engine overhead.
pub const PILEUP_ENGINE_OVERHEAD: u64 = 256;

/// A declared cap on a streaming stage's working set.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MemoryBudget {
    /// Maximum permitted working-set size, in bytes.
    pub bytes: u64,
}

impl MemoryBudget {
    /// A budget expressed in mebibytes.
    pub fn from_mb(mb: u64) -> Self {
        Self {
            bytes: mb.saturating_mul(1024 * 1024),
        }
    }

    /// An effectively unbounded budget.
    pub fn unlimited() -> Self {
        Self { bytes: u64::MAX }
    }

    /// Whether a working set of `working_set_bytes` is permitted.
    pub fn admits(self, working_set_bytes: u64) -> bool {
        working_set_bytes <= self.bytes
    }
}

/// A reported working-set bound for a streaming stage.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct WorkingSet {
    /// Estimated peak working-set size, in bytes.
    pub bytes: u64,
}

impl WorkingSet {
    /// Whether this working set fits within `budget`.
    pub fn fits(self, budget: MemoryBudget) -> bool {
        budget.admits(self.bytes)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn budget_admits_within_and_rejects_beyond() {
        let budget = MemoryBudget::from_mb(2);
        assert!(budget.admits(1_000_000));
        assert!(!budget.admits(3 * 1024 * 1024));
        assert!(budget.admits(2 * 1024 * 1024)); // boundary: working set == budget is admitted
    }

    #[test]
    fn unlimited_admits_everything_and_working_set_checks_fit() {
        assert!(MemoryBudget::unlimited().admits(u64::MAX));
        let ws = WorkingSet { bytes: 512 };
        assert!(ws.fits(MemoryBudget::from_mb(1)));
        assert!(!WorkingSet {
            bytes: 5 * 1024 * 1024
        }
        .fits(MemoryBudget::from_mb(1)));
    }
}
