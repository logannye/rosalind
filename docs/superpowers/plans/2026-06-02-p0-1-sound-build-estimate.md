# P0.1 — Sound Build-Feasibility Prediction: Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans. Steps use checkbox (`- [ ]`) syntax. Run INLINE.

**Goal:** Make `plan --reference` and the `index --memory-budget-mb` plan line predict the build peak from the code-grounded `BuildMemoryModel` (~45 B/base, envelopes the D0-measured 41 B/base) instead of the stale 12 B/base estimate — so the prediction is conservative (predicted ≥ realized), not 3.4× optimistic.

**Architecture:** Add a `BuildMemoryModel::working_set()` helper; point both plan sites at it; render the breakdown in `plan --reference`; delete the stale `estimate_build_working_set`; pin the soundness margin with a test.

**Tech Stack:** Rust; `src/genomics/index/report.rs` + `src/main.rs`.

**Spec:** `docs/superpowers/specs/2026-06-02-p0-1-sound-build-estimate-design.md`

---

### Task 1: Point the build plan at `BuildMemoryModel` + delete the stale estimate

**Files:** `src/genomics/index/report.rs`, `src/genomics/index/mod.rs`, `src/genomics/mod.rs`, `src/main.rs`, `tests/plan_enforce.rs`

- [ ] **Step 1: Strengthen the integration test (failing)** — update `plan_reference_reports_build_estimate` in `tests/plan_enforce.rs` to assert the new breakdown + the ~45 B/base total (the model's `render()` prints B/base per line, so a toy FASTA suffices):

```rust
#[test]
fn plan_reference_reports_build_estimate() {
    let (dir, _idx) = build_index();
    let fa = dir.join("ref.fa");
    let out = Command::new(bin())
        .args(["plan", "--reference"])
        .arg(&fa)
        .args(["--budget-mb", "4096"])
        .output()
        .unwrap();
    assert!(out.status.success());
    let stdout = String::from_utf8_lossy(&out.stdout);
    // The model breakdown is rendered, and the build peak is now the code-grounded
    // ~45 B/base model (envelopes the D0-measured 41 B/base), not the stale 12 B/base.
    assert!(
        stdout.contains("build memory model"),
        "missing model breakdown: {stdout}"
    );
    assert!(
        stdout.contains("45 B/base"),
        "expected the ~45 B/base model total, not the old 12 B/base: {stdout}"
    );
    assert!(stdout.contains("plan:"), "missing build plan verdict line: {stdout}");
    std::fs::remove_dir_all(&dir).ok();
}
```

- [ ] **Step 2: Run to verify it fails**

Run: `cargo test --test plan_enforce plan_reference_reports_build_estimate -- --nocapture 2>&1 | grep -E "test result|FAILED|missing"`
Expected: FAIL — today `plan --reference` prints the 12 B/base one-liner with no breakdown.

- [ ] **Step 3: Add `working_set()` + the soundness test; delete the stale estimator** (`src/genomics/index/report.rs`)

Add the helper inside `impl BuildMemoryModel` (after `from_reference_len`):

```rust
    /// The modeled peak as a `WorkingSet` for the budget plan line. This is the model
    /// `plan --reference`, the `index` plan line, and the build receipt all share, so
    /// the up-front prediction and the realized accounting cannot drift.
    pub fn working_set(&self) -> WorkingSet {
        WorkingSet {
            bytes: self.total_bytes,
        }
    }
```

Extend the doc-comment on `BuildMemoryModel` (the `struct` doc) with a final sentence:

```rust
/// … (existing text) …
///
/// This is the model `plan --reference` and the `index` build receipt SHARE (anti-drift).
/// D0 MEASURED 41 B/base realized (constant across E. coli + yeast); this model is ~45 B/base,
/// so it ENVELOPES the realized peak (model/realized ≈ 1.08). It is a code-grounded model
/// validated by D0 — NOT a proven worst-case bound; the realized peak in the build receipt is
/// the backstop.
```

Delete the stale `estimate_build_working_set` fn entirely (the `pub fn estimate_build_working_set(reference_len: u64) -> WorkingSet { … }` block, ~lines 10–26 including its doc-comment).

Delete its unit test `estimate_grows_with_length_and_does_not_overflow` (the whole `#[test] fn estimate_grows_with_length_and_does_not_overflow() { … }`).

Add a dedicated soundness-envelope test (leave the existing `build_model_breaks_down_and_sums` as-is — its `≥ 20 B/base` sanity check still passes at 45 B/base):

```rust
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
    }
```

- [ ] **Step 4: Drop the stale estimator from the re-exports** — remove `estimate_build_working_set` from:
  - `src/genomics/index/mod.rs` (the `pub use report::{… estimate_build_working_set …}` line),
  - `src/genomics/mod.rs` (the `pub use index::{… estimate_build_working_set …}` line),
  - `src/main.rs` (the `use rosalind::genomics::{… estimate_build_working_set …}` import).

- [ ] **Step 5: Point `run_plan --reference` at the model + render the breakdown** (`src/main.rs`, the `--reference` else-branch)

Replace:

```rust
        let estimate = estimate_build_working_set(total_bp);
        match budget_mb {
            Some(mb) => println!("{}", render_plan_line(estimate, MemoryBudget::from_mb(mb))),
            None => println!(
                "plan: est. build peak ~{} MiB (advisory; build is O(reference)) [no budget]",
                estimate.bytes / (1 << 20)
            ),
        }
```

with:

```rust
        // Predict the build peak from the code-grounded BuildMemoryModel (~45 B/base,
        // envelopes the D0-measured 41 B/base realized) — NOT the old 12 B/base estimate,
        // which under-predicted ~3.4x. Advisory: the build is still O(reference) until D1a.
        let model = BuildMemoryModel::from_reference_len(total_bp);
        print!("{}", model.render(total_bp));
        match budget_mb {
            Some(mb) => println!(
                "{}",
                render_plan_line(model.working_set(), MemoryBudget::from_mb(mb))
            ),
            None => println!(
                "plan: est. build peak ~{} MiB (advisory; build is O(reference)) [no budget]",
                model.total_bytes / (1 << 20)
            ),
        }
```

- [ ] **Step 6: Point the `run_index` plan line at the same model** (`src/main.rs`, ~lines 595–634)

Hoist the model so the plan line and the receipt share one instance. After `let total_bp: u64 = records.iter().map(...).sum();` add:

```rust
    let model = BuildMemoryModel::from_reference_len(total_bp);
```

Change the plan line (in `if let Some(mb) = memory_budget_mb { … }`) from:

```rust
        let estimate = estimate_build_working_set(total_bp);
        eprintln!("{}", render_plan_line(estimate, MemoryBudget::from_mb(mb)));
```

to:

```rust
        eprintln!("{}", render_plan_line(model.working_set(), MemoryBudget::from_mb(mb)));
```

Remove the now-duplicate `let model = BuildMemoryModel::from_reference_len(total_bp);` later in the receipt block (it now reuses the hoisted `model`).

- [ ] **Step 7: Run the test + full suite + gates**

Run: `cargo test --test plan_enforce plan_reference_reports_build_estimate -- --nocapture 2>&1 | grep -E "test result|FAILED"` → PASS.
Run: `cargo test 2>&1 | grep -iE "FAILED" || echo "no failures"` → `no failures`.
Run: `cargo clippy --all-targets -- -D warnings 2>&1 | tail -2` → clean (the gate stays green).
Run: `cargo build --release 2>&1 | grep -i warning || echo "no warnings"; cargo fmt && cargo fmt --check && echo "fmt clean"`.

- [ ] **Step 8: Commit**

```bash
git add src/genomics/index/report.rs src/genomics/index/mod.rs src/genomics/mod.rs src/main.rs tests/plan_enforce.rs
git commit -m "fix(plan): sound build-feasibility prediction — BuildMemoryModel (~45 B/base) not the stale 12 B/base"
```

---

## Final verification

- [ ] `cargo test` green; `cargo clippy --all-targets -- -D warnings` clean; `cargo build --release` 0 warnings; `cargo fmt --check` clean.
- [ ] Smoke: `cargo run --release -- plan --reference examples/data/illumina_toy/reference.fa --budget-mb 1` prints the model breakdown + a `[OK]`/`[OVER]` verdict; the total line shows ~45 B/base. `index --memory-budget-mb` and the build receipt now agree (no within-command 12-vs-45 contradiction).
- [ ] CI green on GitHub.
