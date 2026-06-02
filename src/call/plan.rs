//! Pure, testable planning helpers for `rosalind plan` and `variants --enforce`.
//!
//! The estimator predicts the **data-dependent working set** of a bounded
//! whole-genome germline call from the index header plus the declared depth cap
//! and an assumed max read length. It shares its cost constants with the realized
//! accountant (`PileupEngine::current_working_set`) so a passing plan and the
//! realized receipt cannot silently diverge. The *predicted peak RSS* adds a
//! process baseline the caller measures at runtime (`peak_rss_bytes()`), so the
//! prediction is comparable to the realized `peak_rss` the post-run check uses.

use crate::core::{
    MemoryBudget, WorkingSet, PILEUP_ENGINE_OVERHEAD, PILEUP_MAP_BYTES_PER_BASE,
    PILEUP_PER_READ_OVERHEAD, PILEUP_SEQQUAL_BYTES_PER_BASE,
};

/// Estimate the peak streaming working set of a whole-genome germline call: the
/// largest contig's reference (decoded, 1×) + the depth-capped active read set +
/// the fixed engine overhead. A true upper bound when actual reads do not exceed
/// `max_read_len` and depth is capped at `max_depth` (both enforced at runtime —
/// `max_read_len` is the one assumption, with the post-run check as backstop).
pub fn estimate_variants_working_set(
    largest_contig_len: u64,
    max_depth: u32,
    max_read_len: u32,
) -> WorkingSet {
    let per_read = (max_read_len as u64)
        .saturating_mul(PILEUP_MAP_BYTES_PER_BASE + PILEUP_SEQQUAL_BYTES_PER_BASE)
        .saturating_add(PILEUP_PER_READ_OVERHEAD);
    let active = (max_depth as u64).saturating_mul(per_read);
    WorkingSet {
        bytes: largest_contig_len
            .saturating_add(active)
            .saturating_add(PILEUP_ENGINE_OVERHEAD),
    }
}

/// Predicted peak process RSS = a measured process baseline + the estimated
/// working set. Comparable to the realized `peak_rss` the post-run check uses.
pub fn predicted_peak_rss_bytes(
    largest_contig_len: u64,
    max_depth: u32,
    max_read_len: u32,
    baseline_rss_bytes: u64,
) -> u64 {
    baseline_rss_bytes.saturating_add(
        estimate_variants_working_set(largest_contig_len, max_depth, max_read_len).bytes,
    )
}

/// Render the `rosalind plan --index` breakdown: a measured baseline + the
/// working-set components → predicted peak, tagged `[FITS]`/`[REFUSE]` against the
/// budget (or `[no budget]` when none is declared). Deterministic given inputs.
pub fn render_variants_plan(
    largest_contig_len: u64,
    max_depth: u32,
    max_read_len: u32,
    baseline_rss_bytes: u64,
    budget_mb: Option<u64>,
) -> String {
    const MIB: u64 = 1 << 20;
    let per_read = (max_read_len as u64)
        .saturating_mul(PILEUP_MAP_BYTES_PER_BASE + PILEUP_SEQQUAL_BYTES_PER_BASE)
        .saturating_add(PILEUP_PER_READ_OVERHEAD);
    let active = (max_depth as u64).saturating_mul(per_read);
    let predicted = predicted_peak_rss_bytes(
        largest_contig_len,
        max_depth,
        max_read_len,
        baseline_rss_bytes,
    );
    let verdict = match budget_mb {
        Some(mb) => {
            if MemoryBudget::from_mb(mb).admits(predicted) {
                format!("/ budget {mb} MiB  [FITS]")
            } else {
                format!("/ budget {mb} MiB  [REFUSE]")
            }
        }
        None => "[no budget]".to_string(),
    };
    format!(
        "plan: predicted peak RSS (upper bound)\n  \
         process baseline (measured):        {} MiB\n  \
         reference decode (largest contig):  {} MiB\n  \
         active set @ max-depth {}:           {} MiB\n  \
         engine overhead:                    {} MiB\n  \
         -------------------------------------------------\n  \
         predicted peak: ~{} MiB {}\n",
        baseline_rss_bytes / MIB,
        largest_contig_len / MIB,
        max_depth,
        active / MIB,
        PILEUP_ENGINE_OVERHEAD / MIB,
        predicted / MIB,
        verdict,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn estimate_grows_with_inputs_and_does_not_overflow() {
        let small = estimate_variants_working_set(1_000, 100, 150).bytes;
        let deeper = estimate_variants_working_set(1_000, 2_000, 150).bytes;
        let bigger_ref = estimate_variants_working_set(1_000_000, 100, 150).bytes;
        assert!(deeper > small, "more depth → larger working set");
        assert!(bigger_ref > small, "larger contig → larger working set");
        // Active term for D=100, L=150: 100 * (150*18 + 64) = 100 * 2764 = 276_400.
        assert_eq!(small, 1_000 + 276_400 + 256);
        let _ = estimate_variants_working_set(u64::MAX, u32::MAX, u32::MAX); // no panic
    }

    #[test]
    fn predicted_peak_is_baseline_plus_working_set() {
        let ws = estimate_variants_working_set(248_000_000, 1_000, 250).bytes;
        let predicted = predicted_peak_rss_bytes(248_000_000, 1_000, 250, 50_000_000);
        assert_eq!(predicted, 50_000_000 + ws);
    }

    #[test]
    fn render_reports_fits_and_refuse() {
        // Tiny working set; generous budget → FITS.
        let fits = render_variants_plan(1_000, 100, 150, 1_000_000, Some(4096));
        assert!(
            fits.contains("[FITS]"),
            "generous budget should fit: {fits}"
        );
        // 248 MiB contig + baseline 50 MiB ≫ 64 MiB budget → REFUSE.
        let refuse = render_variants_plan(248 * (1 << 20), 1000, 250, 50 * (1 << 20), Some(64));
        assert!(
            refuse.contains("[REFUSE]"),
            "tight budget should refuse: {refuse}"
        );
        // No budget → advisory.
        assert!(render_variants_plan(1_000, 100, 150, 0, None).contains("[no budget]"));
    }
}
