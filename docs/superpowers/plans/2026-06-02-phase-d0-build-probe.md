# Phase D0 — Measure-First Build Probe Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans (inline, chosen for this work). Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Instrument the index build to emit a build receipt (realized peak RSS + a code-grounded n-scale SA-IS memory accountant), run it on E. coli + yeast, and record a pre-registered verdict that gates the rest of Phase D — without changing the construction algorithm.

**Architecture:** Add a structured `BuildMemoryModel` to `genomics/index/report.rs` that enumerates the real n-scale arrays SA-IS + the FM-index allocate (using true element sizes + a geometric recursion term); have `rosalind index` emit a build receipt (realized peak via `util/rss.rs` + the model breakdown + attribution ratio + bytes/base); a probe script builds E. coli + yeast and captures both; a findings doc records the CONFIRM/NULL verdict. Plus a focused `OPEN_PROBLEMS.md` reframing (√t = framing, native-EM = mechanism).

**Tech Stack:** Rust 1.72, `cargo test`/`fmt`/`build`, bash + curl + python3. No new deps. Branch `rosalind/phase-d0-build-probe` off merged `main`.

**Spec:** [`docs/superpowers/specs/2026-06-02-phase-d0-build-probe-design.md`](../specs/2026-06-02-phase-d0-build-probe-design.md).

**Honesty guardrail (do not tune to the gate):** the accountant must be a *code-grounded enumeration* of the arrays `sais_impl`/`induce_sort`/the FM-index actually allocate, with their true byte sizes — NOT a coefficient reverse-engineered to clear ≥70%. If the honest enumeration lands well below the realized peak, that is a real signal (investigate / NULL), not a number to fudge. The most robust CONFIRM evidence is **bytes/base ≈ constant across E. coli and yeast** (the peak scales with n, not a fixed overhead).

---

## File Structure

- **Modify** `src/genomics/index/report.rs` — add `BuildMemoryModel` (the n-scale accountant) + render; keep `estimate_build_working_set`/`render_plan_line` (still used by `run_plan` + the `index --memory-budget-mb` line). Tests in-file.
- **Modify** `src/main.rs` — `run_index` emits the structured build receipt (realized peak + model + attribution + bytes/base).
- **Create** `scripts/build_memory_probe.sh` — build E. coli + yeast, capture both receipts + the verdict.
- **Modify** `.gitignore` — ignore `results/` (already ignored by the Move-#5 `/results/` rule; verify).
- **Create** `docs/findings/2026-06-02-d0-build-memory-probe.md` — the recorded experiment + verdict.
- **Modify** `docs/OPEN_PROBLEMS.md` — the focused √t-as-framing reframing.

---

## Task 1: `BuildMemoryModel` — the code-grounded n-scale accountant

**Files:**
- Modify: `src/genomics/index/report.rs`

- [ ] **Step 1: Write the failing tests.** Append to the `tests` module in `src/genomics/index/report.rs`:

```rust
    #[test]
    fn build_model_breaks_down_and_sums() {
        let m = BuildMemoryModel::from_reference_len(1_000_000);
        // Components sum to the total.
        let sum: u64 = m.components.iter().map(|(_, b)| b).sum();
        assert_eq!(sum, m.total_bytes, "components must sum to total");
        // The major SA-IS arrays are present and n-scale.
        let names: Vec<&str> = m.components.iter().map(|(n, _)| n.as_str()).collect();
        for needed in ["text(u32)", "suffix array", "lms_name"] {
            assert!(names.iter().any(|n| n.contains(needed)), "missing component {needed}: {names:?}");
        }
        // Honest sanity: SA-IS over a u32 text is many bytes/base, not a handful.
        assert!(m.total_bytes >= 20_000_000, "model must be ≥20 B/base (~{} B/base)", m.total_bytes / 1_000_000);
    }

    #[test]
    fn build_model_is_monotonic_and_overflow_safe() {
        assert!(BuildMemoryModel::from_reference_len(2_000_000).total_bytes
            > BuildMemoryModel::from_reference_len(1_000_000).total_bytes);
        let _ = BuildMemoryModel::from_reference_len(u64::MAX); // must not panic
    }
```

- [ ] **Step 2: Run to verify failure**

Run: `cd ~/rosalind && cargo test -p rosalind --lib genomics::index::report::tests::build_model 2>&1 | tail -10`
Expected: FAIL — `BuildMemoryModel` not found.

- [ ] **Step 3: Implement `BuildMemoryModel`.** In `src/genomics/index/report.rs`, add (after `estimate_build_working_set`):

```rust
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
        // Real element sizes.
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
        // Single-level n-scale footprint.
        let level0: u64 = comps.iter().map(|(_, b)| *b).sum();
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
```

- [ ] **Step 4: Run to verify pass**

Run: `cd ~/rosalind && cargo test -p rosalind --lib genomics::index::report 2>&1 | tail -10`
Expected: PASS (the two new tests + the existing report tests). The model totals ~33 B/base (level0 ~18n + recursion 15n) — a code-grounded estimate, not tuned.

- [ ] **Step 5: Commit**

```bash
cd ~/rosalind && git add src/genomics/index/report.rs && git commit -m "feat(index): BuildMemoryModel — code-grounded n-scale SA-IS build accountant (D0)"
```

---

## Task 2: `rosalind index` emits the build receipt

**Files:**
- Modify: `src/main.rs` (`run_index`)

- [ ] **Step 1: Replace the bare peak-RSS line with a structured receipt.** In `run_index` (`src/main.rs`, ~line 439), the build report import + the final lines are:

```rust
    print!("{}", report.render());

    // Realized peak RSS (per-run, informational) → stderr.
    eprintln!("build peak RSS: {} MiB", peak_rss_bytes() / (1 << 20));
    Ok(())
}
```

Replace those final lines (from `// Realized peak RSS` through `Ok(())`) with:

```rust
    // Build receipt: realized peak RSS vs the modeled n-scale SA-IS build memory —
    // the D0 measure-first probe. The realized peak is machine-dependent; the
    // breakdown + attribution ratio are the analysis payload.
    let peak = peak_rss_bytes();
    let model = BuildMemoryModel::from_reference_len(total_bp);
    let denom = total_bp.max(1);
    eprintln!(
        "build: realized peak RSS {} MiB ({} B/base) over {} bp",
        peak / (1 << 20),
        peak / denom,
        total_bp
    );
    eprint!("{}", model.render(total_bp));
    let ratio = if peak > 0 {
        model.total_bytes as f64 / peak as f64
    } else {
        0.0
    };
    eprintln!(
        "build: model/realized attribution = {:.2} [{}]",
        ratio,
        if ratio >= 0.70 { "CONFIRM ≥0.70" } else { "below 0.70" }
    );
    Ok(())
}
```

And add `BuildMemoryModel` to the `report` import. Find the existing `use` of the report items in `run_index` (the `estimate_build_working_set, render_plan_line, IndexBuildReport` import — at the top of `main.rs` or local to the file) and add `BuildMemoryModel`. (Check: `grep -n "estimate_build_working_set\|IndexBuildReport\|use rosalind::genomics::index" src/main.rs` and extend the matching `use` to include `BuildMemoryModel`.)

- [ ] **Step 2: Build**

Run: `cd ~/rosalind && cargo build 2>&1 | tail -5`
Expected: success, 0 warnings.

- [ ] **Step 3: Smoke-test on the bundled toy** (the receipt renders + the build still works)

Run: `cd ~/rosalind && cargo run -q -- index --reference examples/data/illumina_toy/reference.fa --output /tmp/d0toy.idx 2>&1 | tail -20`
Expected: the build report (contigs/total_bp/blake3/index_bytes) on stdout, then the build-memory model breakdown + `build: realized peak RSS … ` + the attribution line on stderr. (On a tiny toy the realized peak is dominated by the binary baseline, so the ratio will be tiny/`below 0.70` — that's expected for a toy; the gate is assessed on real genomes in Task 4.)

- [ ] **Step 4: Run the existing index CLI gates**

Run: `cd ~/rosalind && cargo test --test index_cli 2>&1 | tail -10`
Expected: PASS (the build receipt is additive; the deterministic stdout build report is unchanged).

- [ ] **Step 5: Commit**

```bash
cd ~/rosalind && git add src/main.rs && git commit -m "feat(cli): index emits a build receipt (realized peak vs modeled n-scale memory) (D0)"
```

---

## Task 3: The probe script

**Files:**
- Create: `scripts/build_memory_probe.sh`
- Modify: `.gitignore` (verify `results/` already ignored)

- [ ] **Step 1: Write the script.** Create `scripts/build_memory_probe.sh`:

```bash
#!/usr/bin/env bash
# D0 measure-first probe: build the index for a small and a ~3x-larger REAL genome,
# capture each build receipt (realized peak RSS + modeled n-scale memory + bytes/base),
# and surface the pre-registered gate verdict. Build-only — no reads/alignment.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RESULTS="$ROOT/results/build-probe"
mkdir -p "$RESULTS"
BIN="$ROOT/target/release/rosalind"
[ -x "$BIN" ] || (cd "$ROOT" && cargo build --release)

# E. coli (gate genome) — reuse the Move-#5 cached reference if present.
ECOLI="$ROOT/results/flagship-ecoli/ecoli.fa"
if [ ! -s "$ECOLI" ]; then
  ECOLI="$RESULTS/ecoli.fa"
  if [ ! -s "$ECOLI" ]; then
    curl -fsSL "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz" -o "$RESULTS/ecoli.fa.gz"
    gunzip -f "$RESULTS/ecoli.fa.gz"
  fi
fi

# S. cerevisiae R64 (~12.1 Mbp, ~3x; multi-contig — build handles it).
YEAST="$RESULTS/yeast.fa"
if [ ! -s "$YEAST" ]; then
  curl -fsSL "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz" -o "$RESULTS/yeast.fa.gz"
  gunzip -f "$RESULTS/yeast.fa.gz"
fi

SUMMARY="$RESULTS/SUMMARY.txt"
: > "$SUMMARY"
probe() {
  local name="$1" ref="$2"
  echo "================ $name ================" | tee -a "$SUMMARY"
  "$BIN" index --reference "$ref" --output "$RESULTS/$name.idx" 2>&1 | tee -a "$SUMMARY"
  echo | tee -a "$SUMMARY"
}
probe "ecoli" "$ECOLI"
probe "yeast" "$YEAST"

echo ">> D0 gate (assess in the findings doc): realized peak within ±25% of ~185 MiB on E. coli" | tee -a "$SUMMARY"
echo ">> AND ≥70% model/realized attribution; bytes/base ~constant across E. coli & yeast = n-scale." | tee -a "$SUMMARY"
echo ">> Summary: $SUMMARY"
```

- [ ] **Step 2: Make executable + confirm gitignore.**

```bash
cd ~/rosalind && chmod +x scripts/build_memory_probe.sh && git check-ignore results/build-probe/ecoli.fa && echo "(results ignored ✓)"
```
(If `results/` is not ignored, append `/results/` to `.gitignore` — but the Move-#5 rule should already cover it.)

- [ ] **Step 3: Commit**

```bash
cd ~/rosalind && git add scripts/build_memory_probe.sh && git commit -m "feat(scripts): build_memory_probe.sh — D0 build-memory probe on real genomes (D0)"
```

---

## Task 4: Run the probe (produce the verdict numbers)

**Files:** none (produces `results/build-probe/`, gitignored)

- [ ] **Step 1: Run it.** (E. coli build is ~seconds; yeast ~12 Mbp build is the longer one — generous timeout. Downloads are cached after first run.)

Run: `cd ~/rosalind && bash scripts/build_memory_probe.sh 2>&1 | tail -40`
Expected: a build report + build-memory model + realized peak + attribution for both E. coli and yeast.

- [ ] **Step 2: Read the captured numbers + compute bytes/base scaling.**

Run: `cd ~/rosalind && grep -E "realized peak RSS|attribution|^=====" results/build-probe/SUMMARY.txt`
Expected: the realized peak RSS + B/base + attribution ratio for each genome. Note for the findings doc: (a) E. coli realized peak vs ~185 MiB (±25%?), (b) the attribution ratio (≥0.70?), (c) E. coli B/base vs yeast B/base (≈ constant → n-scale → CONFIRM the build is SA-workspace-bound).

(No commit — `results/` is gitignored.)

---

## Task 5: Findings doc + `OPEN_PROBLEMS.md` reframing

**Files:**
- Create: `docs/findings/2026-06-02-d0-build-memory-probe.md`
- Modify: `docs/OPEN_PROBLEMS.md`

- [ ] **Step 1: Write the findings doc** with the **actual captured numbers** from Task 4. Sections:
  1. **Result + verdict (one line):** "E. coli build realized peak `<X>` MiB (`<Y>` B/base); yeast `<X2>` MiB (`<Y2>` B/base); model/realized attribution `<r>`. Verdict: **CONFIRM/NULL**."
  2. **The pre-registered gate** (verbatim from the spec) and whether each criterion is met (±25% of 185 MiB; ≥70% attribution; bytes/base constant across the two genomes).
  3. **Table:** genome, bp, realized peak RSS, B/base, model total, attribution ratio.
  4. **Interpretation:** if CONFIRM — the build is intermediate-state-bound (the SA/text/workspace arrays), the baseline curve point (b=n, full-RAM) is set, **D1 (blocked external-memory SA construction) is greenlit**. If NULL — what dominates instead, and the re-scope.
  5. **Reproduce:** `bash scripts/build_memory_probe.sh`. Honest note: realized peak is machine-dependent; the model is a code-grounded peak-set estimate.

- [ ] **Step 2: Reframe `OPEN_PROBLEMS.md`.** Add a short, dated note at the top of the §3.2 "beachhead" / §5-D section (do not rewrite the whole doc) stating the accepted pivot: the √t/Cook–Mertz machinery is **theoretical framing** for "recomputation along a √-shaped space/time curve," **not** the literal construction kernel (SA-IS produces a permutation with no low-degree-extension structure; Cook–Mertz is super-polynomial with no systems realization); Phase D's mechanism is a **native budget-tunable external-memory SA/BWT constructor** honoring a declared `MemoryBudget` along a measured curve, wrapped in the `plan`/`--enforce`/`verify` contract; the differentiation is the **contract + curve + verifiable receipt**, not a new complexity bound. Link the D0 findings doc.

- [ ] **Step 3: Commit**

```bash
cd ~/rosalind && git add docs/findings/2026-06-02-d0-build-memory-probe.md docs/OPEN_PROBLEMS.md && git commit -m "docs(findings): D0 build-memory probe verdict + OPEN_PROBLEMS √t-as-framing reframing (D0)"
```

---

## Task 6: Final verification

**Files:** none

- [ ] **Step 1: Format + zero-warning builds.**

Run: `cd ~/rosalind && cargo fmt --all && cargo fmt --all -- --check && echo FMT_CLEAN; cargo build 2>&1 | grep -c "warning:"; cargo build --release 2>&1 | grep -c "warning:"`
Expected: `FMT_CLEAN`; `0`; `0`.

- [ ] **Step 2: Full suite + tree clean of regenerated data.**

Run: `cd ~/rosalind && cargo test 2>&1 | grep -E "FAILED|panicked|[1-9][0-9]* failed" || echo "no failures"; git status --short`
Expected: no failures; no `results/` tracked.

- [ ] **Step 3: Commit any fmt fixups**

```bash
cd ~/rosalind && git add -A && git commit -m "style: rustfmt fixups (D0)" || true
```

---

## Self-Review notes

- **Spec coverage:** §2 accountant → Task 1; §3 build receipt → Task 2; §4 probe script → Task 3; §4 run → Task 4; §5 findings + §6 OPEN_PROBLEMS reframing → Task 5; §7 testing → Tasks 1/2/6.
- **Placeholder scan:** `<X>/<Y>/<r>` in Task 5 are real numbers captured at run time (Task 4) for the findings doc — not plan placeholders. All code steps are complete.
- **Type consistency:** `BuildMemoryModel::from_reference_len(n) -> Self` with `.components: Vec<(String,u64)>`, `.total_bytes: u64`, `.render(total_bp)` — defined in Task 1, consumed in Task 2. The accountant constants (4 B SA elem, 8 B usize, 1 B types, n/2 LMS) are stated once and reused.
- **Honesty:** the model is a code-grounded enumeration (Task 1 comment + the guardrail at top), not tuned to the gate; the robust CONFIRM signal is bytes/base-constant-across-genomes (Task 4/5). A NULL is an acceptable, recorded outcome.
