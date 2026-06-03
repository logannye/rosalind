//! Pure, testable presentation + planning helpers for `rosalind index`.
//!
//! Kept out of `main.rs` (a thin CLI handler) so the build memory model, the build
//! receipt, and the budget plan line are unit-tested without spawning a process.
//! Nothing here enforces anything — the `MemoryBudget` plan line is record-only
//! (honor-or-refuse is Phase D).

use crate::core::{MemoryBudget, WorkingSet};

/// A code-grounded model of the **peak simultaneously-live n-scale memory** of the
/// SA-IS index build, broken down by the arrays `sais_impl`/`induce_sort`/the
/// FM-index actually allocate. Element sizes are the real Rust types: the SA-IS
/// input text and the suffix array are `u32`/`i32` (4 B); the LMS index vectors
/// are `usize` (8 B); `types` is one byte. `#LMS ≤ n/2` is modeled at n/2. The
/// recursion (`sais_impl(&reduced, …)`) keeps the parent's text/types/lms arrays
/// live while a ~n/2-size child runs, so its overhead is modeled as a geometric
/// tail (~1× the parent's live n-scale arrays). This is a peak-SET estimate, not a
/// sum of every allocation ever made — D0 MEASURES whether it captures the
/// realized peak (the gate); it is not tuned to the gate.
///
/// This is the model `plan --reference` and the `index` build receipt SHARE (anti-drift).
/// D0 MEASURED 41 B/base realized (constant across E. coli + yeast); this model is ~45 B/base,
/// so it ENVELOPES the realized peak (model/realized ≈ 1.08). It is a code-grounded model
/// validated by D0 — NOT a proven worst-case bound; the realized peak in the build receipt is
/// the backstop.
#[derive(Debug, Clone)]
pub struct BuildMemoryModel {
    /// Per-component `(name, bytes)` of the modeled peak set.
    pub components: Vec<(String, u64)>,
    /// Sum of the components.
    pub total_bytes: u64,
}

impl BuildMemoryModel {
    /// Model the peak n-scale build memory for a reference of `n` bases.
    pub fn from_reference_len(n: u64) -> Self {
        const SA_ELEM: u64 = 4; // u32 / i32 suffix-array + text element
        const USIZE: u64 = 8; // LMS index vectors are Vec<usize>
        let lms = n / 2; // #LMS ≤ n/2 (upper-ish)
        let comps: Vec<(&str, u64)> = vec![
            ("text(u32)", n.saturating_mul(SA_ELEM)),
            ("types(1B)", n),
            ("lms_positions(usize)", lms.saturating_mul(USIZE)),
            ("induce-sort suffix array(i32)", n.saturating_mul(SA_ELEM)),
            ("lms_in_sa_order(usize)", lms.saturating_mul(USIZE)),
            ("lms_name(u32)", n.saturating_mul(SA_ELEM)),
            ("reduced string(u32)", lms.saturating_mul(SA_ELEM)),
            ("returned suffix array(u32)", n.saturating_mul(SA_ELEM)),
            ("fm-index (bwt + rank + C-table)", n.saturating_mul(3)),
        ];
        let level0: u64 = comps
            .iter()
            .fold(0u64, |acc, (_, b)| acc.saturating_add(*b));
        // Recursion: parent text/types/lms_positions/lms_name/reduced stay live
        // (~text 4n + types 1n + lms_positions 4n + lms_name 4n + reduced 2n = 15n)
        // while a ~n/2 child runs; the geometric tail ≈ that parent-live amount.
        let recursion = n.saturating_mul(15);
        let mut components: Vec<(String, u64)> =
            comps.into_iter().map(|(s, b)| (s.to_string(), b)).collect();
        components.push(("recursion (geometric tail)".to_string(), recursion));
        let total_bytes = level0.saturating_add(recursion);
        Self {
            components,
            total_bytes,
        }
    }

    /// The modeled peak as a `WorkingSet` for the budget plan line. This is the model
    /// `plan --reference`, the `index` plan line, and the build receipt all share, so
    /// the up-front prediction and the realized accounting cannot drift.
    pub fn working_set(&self) -> WorkingSet {
        WorkingSet {
            bytes: self.total_bytes,
        }
    }

    /// Render the model as a deterministic multi-line breakdown (bytes + per-base).
    pub fn render(&self, total_bp: u64) -> String {
        let mut out = String::from("build memory model (n-scale peak set):\n");
        let denom = total_bp.max(1);
        for (name, bytes) in &self.components {
            out.push_str(&format!(
                "  {name:<34} {:>6} MiB  ({} B/base)\n",
                bytes / (1 << 20),
                bytes / denom
            ));
        }
        out.push_str(&format!(
            "  {:-<34} {:>6} MiB  ({} B/base)\n",
            "total ",
            self.total_bytes / (1 << 20),
            self.total_bytes / denom
        ));
        out
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
    fn model_envelopes_the_d0_measured_realized_peak() {
        // D0 measured 41 B/base realized; the model is ~45 B/base, so it envelopes the
        // realized peak (the soundness margin) without ballooning. Pinned against formula drift.
        let n = 100_000_000u64; // 100 Mbp
        let total = BuildMemoryModel::from_reference_len(n).total_bytes;
        assert!(
            total >= 41 * n,
            "model {} B/base must envelope the D0-measured 41 B/base realized peak",
            total / n
        );
        assert!(
            total <= 50 * n,
            "model {} B/base is implausibly large (formula drift?)",
            total / n
        );
        assert_eq!(
            BuildMemoryModel::from_reference_len(n).working_set().bytes,
            total
        );
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

    #[test]
    fn build_model_breaks_down_and_sums() {
        let m = BuildMemoryModel::from_reference_len(1_000_000);
        let sum: u64 = m.components.iter().map(|(_, b)| b).sum();
        assert_eq!(sum, m.total_bytes, "components must sum to total");
        let names: Vec<&str> = m.components.iter().map(|(n, _)| n.as_str()).collect();
        for needed in ["text(u32)", "suffix array", "lms_name"] {
            assert!(
                names.iter().any(|n| n.contains(needed)),
                "missing component {needed}: {names:?}"
            );
        }
        // Honest sanity: SA-IS over a u32 text is many bytes/base, not a handful.
        assert!(
            m.total_bytes >= 20_000_000,
            "model must be ≥20 B/base (~{} B/base)",
            m.total_bytes / 1_000_000
        );
    }

    #[test]
    fn build_model_is_monotonic_and_overflow_safe() {
        assert!(
            BuildMemoryModel::from_reference_len(2_000_000).total_bytes
                > BuildMemoryModel::from_reference_len(1_000_000).total_bytes
        );
        let _ = BuildMemoryModel::from_reference_len(u64::MAX); // must not panic
    }
}
