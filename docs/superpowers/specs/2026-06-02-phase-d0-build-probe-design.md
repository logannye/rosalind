# Phase D0 — measure-first build-memory probe + gate (design)

**Status:** Spec for review — 2026-06-02. The first increment of **Phase D** (sublinear-space index
construction), after a scoping workflow (`wlb1l1924`) **pivoted the kernel**: the literal √t/Cook–Mertz
combiner is retired as the construction mechanism (it is stubbed *and* discarded on the correctness path;
SA-IS produces a permutation with no low-degree-extension structure to compress; Cook–Mertz is
super-polynomial with no systems realization). Phase D's real target is a **native, polynomial-time,
budget-tunable external-memory SA/BWT constructor** along a *measured* space/time curve, wrapped in the
shipped `plan`/`--enforce`/`verify` contract; **√t/Cook–Mertz is kept as honest theoretical framing only.**
D0 is the cheap, decisive **measure-first gate** that confirms the build is SA-workspace-bound before any
external-memory work begins. Built on a fresh branch off the merged `main`.

## 1. Goal + the precise claim

Instrument the existing index build to emit a **build receipt** (realized peak RSS + a logical accountant
that attributes the peak to the simultaneously-live n-scale SA-IS arrays), run it on a small and a
~3× larger real genome, and record a **pre-registered verdict** that gates the rest of Phase D:

> **CONFIRM** — on E. coli (4.64 Mbp): realized peak RSS within **±25%** of the audited ~185 MiB **AND**
> ≥**70%** of the realized peak attributable to the modeled SA/text/workspace arrays → the build is
> intermediate-state-bound (OPEN_PROBLEMS §3.1), the external-memory path (D1) is **greenlit**.
>
> **NULL** — peak is dominated by something else (allocator slack, the FM-index rank bitvectors, the
> persistence buffers), or attribution < 70% → **stop and re-scope**: a block-tunable SA build would not
> reduce the realized peak; this is a publishable honest-NULL that redirects Phase D and saves months.

D0 changes **no construction algorithm** — it is instrumentation + measurement + the recorded verdict
(plus a small docs reframing, §6). It also ships a real *build* receipt the contract was missing on the
index path.

## 2. The logical build-memory accountant (`genomics/index/report.rs`)

The current `estimate_build_working_set` is a coarse "~12 bytes/base" scalar. Replace/augment it with a
structured, pure, unit-tested model derived from the **actual `sais_impl` allocations**
(`src/genomics/suffix_array.rs`), so the model and the kernel cannot silently drift (the `plan.rs`
discipline). The dominant simultaneously-live n-scale terms (n = reference length in bases):

| Component | Modeled bytes | Source |
|---|---|---|
| `u32` text (SA-IS input alphabet) | 4n | `sais_u32(text: &[u32], …)` |
| `types` (`Vec<SuffixType>`, 1 byte/elem) | 1n | `classify_types` |
| induce-sort SA (`Vec<i32>`) | 4n | `induce_sort` → `vec![-1i32; n]` |
| `lms_name` (`Vec<u32>`) | 4n | `sais_impl` line ~223 |
| reduced string + recursion (geometric, ≤ n/2 first level) | ~2n | `reduced` + `sais_impl(&reduced, …)` |
| returned SA (`Vec<u32>`) | 4n | `sais_impl` return |
| FM-index structures (BWT + rank/C-table) | ~2n | `BlockedFMIndex` build |

`BuildMemoryModel { total_bytes, components: Vec<(name, bytes)> }` with a `from_reference_len(n)`
constructor and a `render()` for the receipt. (The model is an estimate of the *peak set*, not a sum of
all allocations ever made; D0 measures whether it captures the realized peak — that is the experiment.)

## 3. The build receipt (`rosalind index`)

`run_index` (`src/main.rs`) currently prints `build peak RSS: {} MiB`. Extend it to emit a structured
build receipt — to stderr and into the existing index build report (`IndexBuildReport::render` /
`report.rs`):

- realized **peak RSS** (`util::rss::peak_rss_bytes()`, captured after the build) and realized **bytes/base**
  (peak ÷ total_bp);
- the **`BuildMemoryModel` breakdown** (per-component bytes + total) and the modeled bytes/base;
- the **attribution ratio** = modeled total ÷ realized peak, with a `[CONFIRM ≥0.70]` / `[below 0.70]` tag.

The realized peak is machine-dependent; the breakdown + attribution ratio are the analysis payload.

## 4. The probe script (`scripts/build_memory_probe.sh`)

A reproducible bash script (curl + the release binary):

1. **E. coli** (gate genome) — reuse the cached `results/flagship-ecoli/ecoli.fa` if present (from Move #5),
   else fetch GCF_000005845.2 (4.64 Mbp). `rosalind index` → capture the build receipt.
2. **A ~3× larger reference** — *S. cerevisiae* R64 (NCBI RefSeq GCF_000146045.2, ~12.1 Mbp, multi-contig;
   the index handles multi-contig — build only, **no reads/align**). `rosalind index` → capture the receipt.
3. Print both receipts + the **bytes/base for each** (is it ~constant across the two → the build is n-scale)
   and the **gate verdict** for E. coli. Results → `results/build-probe/` (gitignored).

## 5. The findings doc (`docs/findings/2026-06-02-d0-build-memory-probe.md`)

The recorded experiment: realized peak RSS + bytes/base for E. coli and yeast, the component attribution,
and the **pre-registered CONFIRM/NULL verdict** from §1 — which is the D1 go/no-go. Honest notes: realized
peak is machine-dependent; the model is a peak-set estimate; the verdict is what matters. If CONFIRM,
state the greenlight for D1 (blocked external-memory SA construction) and the baseline curve point
(b = n, full-RAM) every later increment must beat.

## 6. Docs reframing (the pivot, made honest)

A focused edit to `docs/OPEN_PROBLEMS.md` (§3.2 / §5-D) reflecting the accepted pivot: the √t/Cook–Mertz
machinery is **theoretical framing** for "recomputation along a √-shaped curve," **not** the literal
construction kernel; Phase D's mechanism is a **native budget-tunable external-memory SA/BWT constructor**;
the differentiation is the **contract + measured curve + verifiable receipt**, *not* a new complexity bound.
(One honest paragraph; no rewrite of the whole doc.)

## 7. Testing

- Unit tests for `BuildMemoryModel`: monotonic in n, the per-component breakdown sums to the total, no
  overflow at u32::MAX-scale n, and the constants match the `sais_impl` allocation sizes.
- `cargo test --test index_cli` (the existing build CLI gates) stays green; the build receipt is additive.
- The probe script + findings doc **are the experiment** (not a permanent unit test — realized RSS is
  machine-dependent). The script asserts the build succeeds; the verdict is transcribed into the doc.
- Gates unchanged: `cargo fmt --check`, `cargo build` 0 warnings (debug + release), full suite.

## 8. Non-goals

No construction-algorithm change; no blocked/external-memory build (that is D1, gated on D0's CONFIRM); no
u64 widening (D4); no `sa_sample_rate` knob (it shrinks *query*, not *build* peak — explicitly out, per the
guardrail); no new complexity-bound claim; no Cook–Mertz combiner work (retired). MSRV 1.72; no new deps.

## 9. Branch

`rosalind/phase-d0-build-probe`, off the merged `main`. Its own PR when done (the user decides merge).
