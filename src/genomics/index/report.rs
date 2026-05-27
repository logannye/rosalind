//! Pure, testable presentation + planning helpers for `rosalind index`.
//!
//! Kept out of `main.rs` (a thin CLI handler) so the build working-set estimate,
//! the build receipt, and the budget plan line are unit-tested without spawning a
//! process. Nothing here enforces anything — the `MemoryBudget` plan line is
//! record-only (honor-or-refuse is Phase C).

use crate::core::{MemoryBudget, WorkingSet};

/// A coarse, **record-only** estimate of the peak working set of building a
/// `GenomeIndex` over a reference of `reference_len` bases.
///
/// The build is dominated by SA-IS over the `u32` text (text + suffix array +
/// workspace) plus the in-RAM index structures — roughly **12 bytes per base**.
/// This is an intentionally coarse upper-ish model for the budget *seam*; precise
/// accounting is Phase C and the build cost itself is what Phase D reduces. It is
/// not a guarantee.
pub fn estimate_build_working_set(reference_len: u64) -> WorkingSet {
    const BYTES_PER_BASE: u64 = 12;
    const BASE_OVERHEAD: u64 = 1 << 20; // 1 MiB floor for short references
    WorkingSet {
        bytes: reference_len
            .saturating_mul(BYTES_PER_BASE)
            .saturating_add(BASE_OVERHEAD),
    }
}

/// The deterministic build receipt for a persisted index. Per-run fields (e.g.
/// realized RSS) are intentionally excluded — the caller prints those separately.
#[derive(Debug, Clone)]
pub struct IndexBuildReport {
    /// Path the index was written to.
    pub index_path: String,
    /// Per-contig `(name, length)` in id order.
    pub contigs: Vec<(String, u32)>,
    /// Total reference length in bases.
    pub total_bp: u64,
    /// BLAKE3 of the (uppercased ASCII) reference.
    pub reference_blake3: [u8; 32],
    /// On-disk size of the index file, in bytes.
    pub index_bytes: u64,
}

impl IndexBuildReport {
    /// Render the receipt as a deterministic multi-line string (trailing newline).
    pub fn render(&self) -> String {
        let mut out = String::new();
        out.push_str(&format!("index: {}\n", self.index_path));
        out.push_str(&format!(
            "contigs: {} ({} bp total)\n",
            self.contigs.len(),
            self.total_bp
        ));
        for (name, length) in &self.contigs {
            out.push_str(&format!("  {name}\t{length}\n"));
        }
        out.push_str(&format!(
            "reference_blake3: {}\n",
            hex32(&self.reference_blake3)
        ));
        out.push_str(&format!("index_bytes: {}\n", self.index_bytes));
        out
    }
}

/// Render the record-only budget plan line: estimated build peak vs the declared
/// budget, tagged `[OK]` (estimate fits) or `[OVER]` (estimate exceeds — the build
/// proceeds anyway; enforcement is Phase C).
pub fn render_plan_line(estimate: WorkingSet, budget: MemoryBudget) -> String {
    let verdict = if estimate.fits(budget) { "OK" } else { "OVER" };
    format!(
        "plan: est. build peak ~{} MiB / budget {} MiB  [{}]",
        estimate.bytes / (1 << 20),
        budget.bytes / (1 << 20),
        verdict
    )
}

/// Lowercase hex of a 32-byte digest.
fn hex32(bytes: &[u8; 32]) -> String {
    const HEX: &[u8; 16] = b"0123456789abcdef";
    let mut out = vec![0u8; 64];
    for (i, b) in bytes.iter().enumerate() {
        out[2 * i] = HEX[(b >> 4) as usize];
        out[2 * i + 1] = HEX[(b & 0x0f) as usize];
    }
    String::from_utf8(out).expect("hex digits are valid ASCII")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn estimate_grows_with_length_and_does_not_overflow() {
        let small = estimate_build_working_set(1_000).bytes;
        let large = estimate_build_working_set(1_000_000).bytes;
        assert!(large > small, "estimate must grow with reference length");
        assert!(large >= 12_000_000, "≈12 bytes/base");
        let _ = estimate_build_working_set(u64::MAX); // must not panic/overflow
    }

    #[test]
    fn receipt_render_is_deterministic_and_contains_fields() {
        let report = IndexBuildReport {
            index_path: "ref.idx".to_string(),
            contigs: vec![("chr1".to_string(), 100), ("chr2".to_string(), 50)],
            total_bp: 150,
            reference_blake3: [0xab; 32],
            index_bytes: 4096,
        };
        let a = report.render();
        assert_eq!(a, report.render(), "render must be deterministic");
        assert!(a.contains("index: ref.idx"));
        assert!(a.contains("contigs: 2 (150 bp total)"));
        assert!(a.contains("  chr1\t100"));
        assert!(a.contains("  chr2\t50"));
        assert!(a.contains(&format!("reference_blake3: {}", "ab".repeat(32))));
        assert!(a.contains("index_bytes: 4096"));
    }

    #[test]
    fn plan_line_reports_ok_and_over() {
        let budget = MemoryBudget::from_mb(100);
        let under = WorkingSet {
            bytes: 50 * (1 << 20),
        };
        let over = WorkingSet {
            bytes: 200 * (1 << 20),
        };
        assert!(render_plan_line(under, budget).ends_with("[OK]"));
        assert!(render_plan_line(over, budget).ends_with("[OVER]"));
        assert!(render_plan_line(under, budget).contains("budget 100 MiB"));
    }
}
