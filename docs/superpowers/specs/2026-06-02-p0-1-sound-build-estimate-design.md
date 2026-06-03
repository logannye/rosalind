# P0.1 — Sound Build-Feasibility Prediction (design)

**Status:** Approved design — 2026-06-02. Phase 0 of the implementation roadmap
([`docs/ROADMAP.md`](../../ROADMAP.md) §6, P0.1). The do-first soundness foundation: a `plan` that
says FITS then OOMs would destroy the one thing the contract sells.

---

## 1. Problem

`rosalind plan --reference` (and the `rosalind index --memory-budget-mb` plan line) predict the index
build's peak memory at a **stale, dangerously optimistic 12 B/base** (`estimate_build_working_set`,
`genomics/index/report.rs:19`), while:

- the D0 build-memory probe **measured 41 B/base realized** (constant across E. coli + yeast =
  n-scale), and
- the code-grounded `BuildMemoryModel` (`report.rs:46`) sums the real SA-IS array sizes to **~45 B/base**
  and is already used by `rosalind index`'s build receipt (`main.rs:634`).

So today the build-feasibility prediction **under-estimates the realized peak by ~3.4×**, and a single
command (`index --memory-budget-mb`) prints *two different numbers*: a 12 B/base plan line and a 45
B/base receipt model. A `plan` that says FITS at 12 B/base then OOMs at 41 B/base is the exact
"predict-before-you-commit" failure that destroys the contract's value — and it is the gating
prerequisite for any build-enforcement (`index --enforce`, Phase D1a).

## 2. Goals / non-goals

**Goals**
- Make `plan --reference` and the `index --memory-budget-mb` plan line predict the build peak from the
  **`BuildMemoryModel` (~45 B/base)**, so the prediction is a **conservative upper-ish bound** (envelopes
  the D0-measured 41 B/base by ~10%) instead of under-predicting 3.4×.
- One model: `plan`, the `index` plan line, and the `index` build receipt all use `BuildMemoryModel`
  (anti-drift, exactly as the *calling* path shares its cost constants between predictor and accountant).
- `plan --reference` renders the model's **per-component breakdown** (transparency: *why* the peak is
  what it is).
- Delete the now-unused stale estimator.

**Non-goals**
- Making the build itself bounded — that is Phase D1a (the external-memory build). This increment makes
  the *prediction* sound; `plan --reference` stays honestly labeled **advisory** (build is still
  O(reference)).
- Claiming a *proven worst-case* bound. `BuildMemoryModel` is a code-grounded model validated by D0
  (attribution ~1.08 = model/realized), not a proof. The realized post-build peak (the receipt) remains
  the backstop; the doc-comment + output stay honest about this.
- Changing the `BuildMemoryModel` formula, the build itself, or `index --enforce` (not yet wired).

## 3. Design

### 3.1 Swap the estimator at the two plan sites

Both sites currently call `estimate_build_working_set(total_bp)` (12 B/base). Replace with the model:

- **`run_plan` `--reference` path** (`main.rs:~732`): build `let model = BuildMemoryModel::from_reference_len(total_bp);` and predict from `model.total_bytes`. Print the **breakdown** (`model.render(total_bp)`) followed by the verdict line (`render_plan_line(WorkingSet { bytes: model.total_bytes }, budget)`), or the `[no budget]` advisory line. Keep the **"advisory; build is O(reference)"** honesty.
- **`run_index --memory-budget-mb` plan line** (`main.rs:~599`): use the same `BuildMemoryModel` total for `render_plan_line` (one line; the full breakdown already prints in the build receipt below it, which `run_index` computes at `:634`).

`render_plan_line` is unchanged — it takes a `WorkingSet`; we hand it `WorkingSet { bytes: model.total_bytes }`.

### 3.2 Delete the stale estimator

`estimate_build_working_set` becomes unused after §3.1. Remove:
- the fn + its doc-comment (`report.rs:10–26`),
- its unit test `estimate_grows_with_length_and_does_not_overflow` (`report.rs:167–174`),
- its re-exports: `genomics/index/mod.rs:16`, `genomics/mod.rs:36`, and the `main.rs:11` import.

(The 12 B/base number was explicitly "not a guarantee"; nothing else depends on it.)

### 3.3 Soundness test + honest doc

- **Test (in `report.rs`):** assert the model **envelopes the D0-measured realized peak** — for a
  representative `n` (e.g. 100 Mbp), `BuildMemoryModel::from_reference_len(n).total_bytes ≥ 41 * n`
  (the measured B/base) **and** `≤ 50 * n` (a sanity ceiling so the model can't silently balloon). This
  documents the soundness margin (45 ≥ 41) and pins it against regressions in the formula.
- **Doc-comment (on `BuildMemoryModel` or `from_reference_len`):** note that this is the model `plan`
  and the build receipt share, that D0 measured 41 B/base realized (model/realized ≈ 1.08), and that it
  is a code-grounded *model* validated by D0, **not** a proven worst-case bound — the realized peak in
  the build receipt is the backstop.

## 4. File-by-file change list

| File | Change |
|---|---|
| `src/main.rs` | `run_plan --reference`: predict from `BuildMemoryModel` + render the breakdown; `run_index`: plan line uses the model total; drop the `estimate_build_working_set` import. |
| `src/genomics/index/report.rs` | Delete `estimate_build_working_set` + its test; extend the `BuildMemoryModel` doc-comment (the D0 41 B/base / model-not-proof note); add the envelope soundness test. |
| `src/genomics/index/mod.rs`, `src/genomics/mod.rs` | Drop `estimate_build_working_set` from the re-exports. |

## 5. Testing

1. **Unit (`report.rs`):** the model envelopes the D0-measured peak (`total_bytes ≥ 41·n`, `≤ 50·n`);
   the existing model breakdown/monotonic/overflow tests stay green.
2. **Integration (`tests/index_cli.rs` or `plan_enforce.rs`):** `plan --reference <fa> --budget-mb N`
   now prints the model **breakdown** (`render()` shows per-line **B/base**, which is `total_bytes/total_bp
   = 45` at *any* reference size, so a small toy FASTA suffices). Assert the output contains the
   breakdown header (`build memory model`) and a total line showing **~45 B/base** (≫ the old 12 B/base
   one-liner with no breakdown), plus an `[OK]` verdict under a generous budget and `[OVER]` under a tiny
   one. Use a toy FASTA via the existing index fixtures.
3. **Regression:** full suite green; the `index --memory-budget-mb` plan line and the build receipt now
   agree (no within-command 12-vs-45 contradiction); clippy `-D warnings` + MSRV jobs stay green.

## 6. References

- `src/genomics/index/report.rs` — `estimate_build_working_set` (12 B/base, to delete) + `BuildMemoryModel` (~45 B/base, to use).
- `src/main.rs:599` (index plan line), `:634` (index receipt model), `:732` (`plan --reference`).
- `docs/findings/2026-06-02-d0-build-memory-probe.md` — the 41 B/base measurement + the ~1.08 attribution.
- `docs/ROADMAP.md` §6 Phase 0 — this increment.
