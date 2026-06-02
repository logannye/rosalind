# Rosalind — Open Problems & Research Thesis

**Status:** Living research-direction document — 2026-05-27. Audience: builders and researchers
who fork, extend, or build on Rosalind. Companion to the engineering roadmap in
[`docs/superpowers/specs/2026-05-26-rosalind-target-architecture.md`](superpowers/specs/2026-05-26-rosalind-target-architecture.md).

---

## 1. The thesis

Rosalind does not try to out-align bwa-mem2/minimap2 or out-call DeepVariant — that race is a
non-goal. Its differentiated, defensible contribution is a property almost no production genomics
tool offers:

> **Memory is a declared, predictable, never-refusing, verifiable contract** — across the entire
> genomics lifecycle (index construction, alignment, pileup, calling, query), with a **continuous
> space/time tradeoff** the builder selects by stating a RAM budget.

You tell Rosalind the RAM you have; it tells you — *before you commit* — whether the job fits and
how long it will take; it then honors that ceiling, sliding along a space/time curve to do so; and
when the budget is below the comfortable in-RAM regime it **degrades gracefully** (more
recomputation, then spill to disk) rather than refusing. Every run emits a receipt recording the
*realized* peak working set, and CI enforces the bound. This is the capability that matters to the
people who actually need Rosalind: **edge, field, clinical, and large-/meta-/pan-genome** settings
where RAM is fixed, swap is death, and an unpredictable OOM at hour 20 of a build is unacceptable.

## 2. The mechanism: a space/time *curve*, not an operating point

The enabling machinery is the square-root-space evaluation framework Rosalind is built over
(Williams 2025, *Simulating Time with Square-Root Space*; Cook–Mertz 2024 tree evaluation). The
space bound is **~√t** — the square root of the running time, up to lower-order (sub-polynomial)
factors. (The "O(√t)" in the source headers denotes the same bound.)

The key reframing — and what makes this a *systems* contribution rather than a theory citation —
is that this is **not a single low-space operating point.** Block-respecting simulation with a
tunable block size `b` (and checkpoint density) yields a **continuous curve**: from `b = t`
(standard, fast, O(t) space) through `b = √t` (minimal space ~√t, more
recomputation) and onward toward external-memory spill. **A declared RAM budget simply selects the
point on that curve.** The theory layer's job is to *be the knob.*

## 3. The open problem

The Williams/Cook–Mertz result is a **space** upper bound that, in the general simulation, **buys
space by spending time** (recomputation). So the real research question is not "can we hit
~√t space" — it is:

> **For which genomics computations can we realize the full space/time curve with *tolerable time
> overhead*?** I.e., where is the trade a small polynomial / modest constant factor rather than a
> blow-up?

### 3.1 The decisive lens: intermediate-state-bound vs. input-bound

The simulation shrinks the **working state** of a time-*t* computation. It does **nothing** to
shrink the memory needed to *hold the inputs or outputs.* This single fact tells you where the
framework can and cannot help:

| Computation | Dominant memory | Framework applies? |
|---|---|---|
| **FM-index / suffix-array construction** | the suffix array + recursion **workspace** (≈4–8·n), not the reference (≈n/4 packed) — **intermediate state** | **Yes — sweet spot** |
| Long-sequence/graph **DP** (alignment) | the DP matrix — **intermediate state** (though Hirschberg/banding already help the common cases) | Partially |
| **Joint / cohort calling** | holding *N* samples' data — **input** (O(N)) | No — beat by streaming samples, not √t |
| **Pangenome-graph** analysis | the graph/index itself — **input** | No — the graph dominates; √t doesn't shrink it |

This is *why* construction is the beachhead, and the honest test any future target must pass:
**is the memory bottleneck the computation's working state, or its inputs?** Input-bound problems
(cohort calling, pangenomes) are real and important, but they need *orthogonal* techniques
(streaming, succinct/compressed inputs) — they are framed here as separate open problems, not as
things this framework will magically solve.

### 3.2 The beachhead: sublinear-space full-text index construction

> **Pivot — Phase D mechanism (2026-06-02; see [`findings/2026-06-02-d0-build-memory-probe.md`](findings/2026-06-02-d0-build-memory-probe.md)).**
> A scoping study established that the literal √t / Cook–Mertz machinery is **theoretical framing**, not the
> construction *kernel*: SA-IS produces a permutation (random-access induced sorting) whose memory is the
> suffix array itself — there is no low-degree-extension structure to compress, the in-repo combiner is a
> stub whose output is discarded on the correctness path, and Cook–Mertz is super-polynomial with no systems
> realization. So Phase D's **mechanism** is a **native, polynomial-time, budget-tunable external-memory
> SA/BWT constructor** (pSAscan/eSAIS/Big-BWT-style: SA-IS each budget-sized block in RAM → exact disk-backed
> merge) that honors a declared `MemoryBudget` along a **measured** space/time curve, **degrading to disk
> rather than refusing**, wrapped in the shipped `plan`/`--enforce`/`verify` contract. The √t bound (Williams
> 2025; Cook–Mertz 2024) is kept as the asymptotic *justification* for "recomputation along a √-shaped
> curve" — **not** as the runtime kernel, and Rosalind makes **no new space-complexity claim** (prior
> succinct-construction work already hits low-space fixed points). The defensible contribution is the
> **declared-budget contract + measured curve + verifiable receipt over a real memory-bound build.** The D0
> probe (above link) confirmed the build is intermediate-state-bound by the SA workspace (41 B/base, constant
> across E. coli and yeast), greenlighting this path.

Building the FM-index/BWT is Rosalind's own headline memory limitation ("the FM-index is built in
memory at the start of each run; memory proportional to the reference"). It is the ideal first
target because it is:

1. **Intermediate-state-bound** — the SA-IS suffix array + recursion workspace dominate; the input
   (reference) is small. Exactly the framework's sweet spot.
2. **One-time and offline** — so the time-for-space trade is *most tolerable*: build slowly, even on
   the constrained device itself; distribute the artifact; mmap it forever.
3. **On the critical path** — it unblocks indexing references that do not fit in RAM today (T2T,
   wheat ~16 Gbp, conifers, metagenomic catalogs, pangenome linearizations) on commodity hardware,
   and completes Rosalind's bounded-memory story *end to end* (build *and* use within a declared
   budget).
4. **A clean, composable primitive** — useful to a builder immediately, not a giant end-to-end
   subsystem they must adopt wholesale.

**Open question (D):** characterize and minimize the time overhead of budget-tunable, sublinear-
space full-text index construction; identify the regimes where it is practical.

### 3.3 The broader landscape (framed, not yet committed)

- **DP-heavy steps under budget** — long-read/graph alignment DP as a tunable space/time computation.
- **Input-bound frontiers** — space-bounded joint/cohort calling and pangenome-graph analysis;
  these require streaming + succinct-structure techniques, *complementary* to (not solved by) the
  √t machinery. Honest open problems.

## 4. The contract design (what "declared memory" means concretely)

A single declared budget governs the whole system, not just one stage:

- **Declare → predict.** `rosalind plan --budget 4G …` reports *feasibility*, the *predicted peak
  working set*, the *estimated time*, and the *curve* ("at 8 GB → ~3× faster") — so you commit with
  your eyes open. No surprise OOM, no surprise multi-day run.
- **Honor or refuse — then degrade rather than refuse.** Every streaming stage exposes a provable
  working-set bound and honors the ceiling; below the in-RAM threshold it slides further down the
  curve (more recomputation → external-memory spill) and **still completes**, slower. A 2 GB field
  device builds a human index; it does not get "won't fit."
- **Lifecycle-wide.** The same budget bounds **construction, mmap-query (tunable SA-sampling),
  alignment, pileup, calling, and sort.** State your RAM once.
- **Verifiable.** The run receipt records the *realized* peak working set (content-addressed,
  alongside tool version + input/param/output hashes); CI enforces the envelope on growing inputs;
  `rosalind verify` re-checks a receipt without re-running.
- **Deterministic regardless of the point chosen.** The same inputs + budget produce byte-identical
  artifacts; parallelism is thread-count-invariant. The space/time point never changes the *output*.

## 5. Roadmap (B → E), re-derived around the thesis

Prioritized for the builders who need the unique capability first; performance/breadth follow. The
practical pillars (B, C) stand alone — a great engine even if the research bet (D) runs long — and
D lands the space-complexity headline on top.

- **B — Genome-scale kernel (finish).** Zero-copy persisted lean multi-contig index (build-once →
  mmap), wired consumers, whole-genome pileup, pipe-native. Lay the `MemoryBudget` / space-accounting
  hooks. The theory layer (`algebra`/`ledger`/`machine`/`tree`/`space`) is **kept and repurposed**,
  not feature-gated away.
- **C — Memory as a verifiable contract.** Declared `MemoryBudget` threaded through every streaming
  stage + mmap-query; `rosalind plan` (envelope + curve); honor-or-refuse + graceful degradation;
  real RSS gates in CI; deterministic (thread-count-invariant) execution + receipts + `rosalind
  verify`. → "declare your RAM → predictable, honored, verifiable, end-to-end" for
  align/call/query/sort.
- **D — Sublinear-space index construction (beachhead research contribution).** The √t-family knob:
  build the index under the declared budget across the full curve (in-RAM-fast → checkpointed
  ~√n → external-memory spill), with the time overhead characterized; makes the theory
  layer load-bearing; erases the O(reference) build-RAM caveat; completes the contract for the build
  step.
- **E — Reach + substrate.** Fast deterministic-parallel aligner + germline indels + richer read QC
  (broad usability); the Python/tensor **pileup-feature substrate** + pluggable models (the
  ML-researcher/builder surface); rigorous eval + calibration validation + the open-problems
  benchmark suite.

(Deterministic-parallel *execution* is part of the contract in **C**; aligner-*algorithm* quality and
micro-optimization are **E**. The existing exact-match aligner is sufficient to exercise C and D on
real data — the point there is memory, not speed.)

## 6. What this means for you

- **Edge / field / clinical builders** — declare your device's RAM; get a checkable guarantee it
  fits, a time estimate, byte-identical + auditable outputs, and a build that degrades rather than
  dies.
- **Large / meta / pangenome researchers** — index references that do not fit in RAM today, on
  commodity hardware.
- **ML-genomics builders** — (Phase E) a fast, deterministic, memory-bounded, Python/tensor-
  accessible per-locus feature stream over a whole genome to compute on.
- **Algorithms / systems researchers** — a working testbed for space-bounded computation on a real,
  memory-bound domain, and a concrete open problem (§3) with a clean first target.

## 7. References

- R. R. Williams. *Simulating Time with Square-Root Space.* 2025.
- J. Cook, I. Mertz. *Tree Evaluation in Sublinear Space* (catalytic / tree-evaluation lineage). 2024.
- The "O(√t)" label in the source headers denotes the same ~√t bound (square-root space, up to lower-order factors).
