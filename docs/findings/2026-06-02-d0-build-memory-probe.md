# D0 — measure-first build-memory probe (the Phase-D gate)

**2026-06-02.** Before committing the Phase-D external-memory construction arc, a cheap, decisive
measure-first probe: is the index *build* peak actually dominated by the SA-IS workspace (so a
budget-tunable / blocked build can reduce it), or by something a block-tunable build wouldn't touch?

**Verdict: CONFIRM.** On the real *E. coli* K-12 MG1655 chromosome (4,641,652 bp), the realized build peak
is **182 MiB** — within **1.6%** of the audited ~185 MiB (gate: ±25%) — at **41 bytes/base**. On
*S. cerevisiae* R64 (12,157,105 bp, 17 contigs) it is **482 MiB**, also **41 bytes/base**. The bytes/base is
**identical across a 2.6× size jump**, and a code-grounded model of the n-scale SA-IS arrays envelopes the
realized peak (model/realized = 1.09 / 1.08). The build is **intermediate-state-bound** by the
suffix-array workspace, exactly as `docs/OPEN_PROBLEMS.md` §3.1 claimed → **D1 (blocked external-memory SA
construction) is greenlit**, with the baseline curve point set at **b = n (full-RAM) = 41 B/base**.

## The pre-registered gate (and the result)

| Criterion | Threshold | Result | |
|---|---|---|---|
| E. coli realized peak vs audited ~185 MiB | within ±25% | **182 MiB (−1.6%)** | ✅ |
| model/realized attribution (peak ⊆ modeled SA arrays) | ≥ 0.70 | **1.09** (model ≥ realized) | ✅ |
| bytes/base constant across genomes (peak is n-scale, not fixed overhead) | ≈ constant | **41 → 41** (E. coli → yeast) | ✅ |

All three met. The peak is the SA construction workspace; it scales with reference length, not with a
fixed baseline — so reducing the *simultaneously-live* SA workspace (the D1+ external-memory/blocked plan)
will reduce the realized build peak.

## The numbers

| Genome | bp | realized peak RSS | B/base | model total | model B/base | attribution |
|---|---|---|---|---|---|---|
| *E. coli* K-12 MG1655 (NC_000913.3) | 4,641,652 | 182 MiB | 41 | 199 MiB | 45 | 1.09 |
| *S. cerevisiae* R64 (17 contigs) | 12,157,105 | 482 MiB | 41 | 521 MiB | 44 | 1.08 |

The code-grounded `BuildMemoryModel` (from the real `sais_impl` allocations) breaks the ~45 B/base down as:
`u32` text 4 + `types` 1 + `lms_positions` (usize) 4 + induce-sort SA (i32) 4 + `lms_in_sa_order` (usize) 4 +
`lms_name` 4 + reduced string 2 + returned SA 4 + FM-index 3 + recursion (geometric tail) 15. The realized
41 B/base sitting just under the 45 B/base model confirms these arrays are the peak set (not, say, allocator
slack, the FM-index rank bitvectors, or the persistence buffers — any of which would have produced a
realized peak well above the n-scale model, i.e. a NULL).

## What this greenlights (D1+)

The realized peak is the **simultaneously-live** SA-IS workspace at full RAM (b = n). The D1 plan — SA-IS
each budget-sized **block** of the text in RAM, then an exact disk-backed merge into the full SA/BWT — keeps
only one block's workspace live at a time, so a declared `MemoryBudget` selects the block size and the
realized peak drops along a measured space/time curve (~pSAscan's O(n²/M) shape), degrading to disk rather
than refusing, wrapped in the shipped `plan`/`--enforce`/`verify` contract. D0 is the baseline every later
increment must beat.

## Honest notes

- The realized `peak_rss` is **machine-dependent**; the model is a code-grounded *peak-set* estimate, not a
  byte-exact predictor. What is robust is the **bytes/base constancy** (n-scale) and that the realized peak
  sits **within** the modeled SA-workspace envelope on two genomes.
- This is **not** a new complexity-bound claim and **not** the √t simulation made real — see the Phase-D
  reframing in [`../OPEN_PROBLEMS.md`](../OPEN_PROBLEMS.md). The differentiation is the contract + measured
  curve + verifiable receipt over a real memory-bound build.

## Reproduce

```bash
bash scripts/build_memory_probe.sh   # builds E. coli + yeast, prints both build receipts + the gate
```

Build-only (no reads/alignment). Downloads cached after the first run; regenerated data lands in
`results/build-probe/` (gitignored). The build receipt is emitted by every `rosalind index` run.
