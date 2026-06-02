# Phase D1 — blocked external-memory SA/BWT construction (design notes, pre-brainstorm)

**Status:** DESIGN NOTES for the morning brainstorm — 2026-06-02 (overnight). **Not a finished spec.** D0
CONFIRMED (build is SA-workspace-bound, 41 B/base constant across E. coli & yeast — see
`docs/findings/2026-06-02-d0-build-memory-probe.md`, PR #23), greenlighting D1. D1 is the genuinely hard,
correctness-critical, research-grade increment of Phase D; it carries real **design forks** that want
Daisy/Logan's call before implementation (the way the √t pivot did). These notes lay out the decomposition,
the algorithm options, and the forks so the morning brainstorm → spec → plan → inline execution is fast.
**These were prepared overnight while the implementation was deliberately deferred to supervised execution.**

## 1. Goal (recap)

A **native, polynomial-time, budget-tunable external-memory SA/BWT constructor**: index a reference *larger
than RAM* under a declared `MemoryBudget`, along a **measured** space/time curve (more recomputation/spill
→ less RAM, more time), **degrading to disk rather than refusing**, producing a **byte-identical** index to
the in-RAM `sais_u32` build, wrapped in the shipped `plan`/`--enforce`/`verify` contract. D0 set the baseline
curve point: **b = n (full RAM) = 41 B/base**. The win is reducing the *simultaneously-live* SA workspace by
processing the text in budget-sized blocks and spilling.

## 2. Algorithm options (FORK — needs the morning call)

| Option | Knob | Pros | Cons |
|---|---|---|---|
| **pSAscan-style** (block SA-IS + gap-array merge; Kärkkäinen-Kempa-Puglisi) — *recommended* | block size = f(budget) | reuses the tested `sais_u32` per block; the canonical external-memory SA; clean budget→peak | the **gap-array merge** (backward-search "matching" of block suffixes vs the processed tail) is the hard, bug-prone crux |
| **DC3 / skew** (difference-cover) | sample param v | cleanest built-in √-knob (O(vn) time / O(n/√v) space, *linear work*) | discards the tested SA-IS kernel; recursive-skew EM spill is awkward; revalidation risk on the most load-bearing primitive |
| **Big-BWT / prefix-free parsing** | window/modulus | excellent for *repetitive* inputs (pangenomes) | not a general budget knob for non-repetitive refs; a different artifact — better as the D5 second mode |

**Recommendation:** pSAscan-style (it reuses our verified SA-IS, is the canonical EM-SA, and maps cleanly
onto the contract). DC3 only if we decide the gap-merge risk is worse than re-validating a skew kernel.

## 3. Proposed decomposition (each sub-step independently verifiable + gated)

- **D1a — byte-identical test harness + blocked-construction API.** A `genomics::ext_sa` module with
  `build_sa_blocked(text, max_symbol, block_size) -> Vec<u32>` and a **property test** asserting it equals
  `sais_u32` for all block sizes on the existing SA-IS test inputs + random small inputs. The *first*
  implementation may be a correctness-reference (even a simple block-sort + k-way-merge) purely to lock the
  API + the byte-identical gate; it need not yet reduce peak or be fast. **Gate:** byte-identical to
  `sais_u32`, all block sizes. *(This is the one piece that is safe to build unsupervised — pure correctness
  with a total-order gate — but it is mostly scaffolding; the real value is D1b/c.)*
- **D1b — disk-spill → real peak reduction.** Spill per-block partial results + the output SA to a `mmap`'d
  temp file; stream the merge; keep only one block's workspace + merge front in RAM. **Gate:** byte-identical
  index AND `rosalind index` build receipt shows realized peak ≈ f(block_size) **< the D0 baseline** (41
  B/base) on E. coli. This is where the curve becomes real.
- **D1c — efficient per-block construction (the gap-array merge).** Replace the reference merge with SA-IS
  per block + the pSAscan gap-array merge so the *time* is tolerable at genome scale (the D0 baseline runs
  in ~seconds; the blocked build must stay within a small factor). **Gate:** byte-identical AND wall-time
  within a characterized factor of the b=n baseline on E. coli + yeast.
- **D2 — contract wiring.** block size = f(`MemoryBudget`); extend `report.rs`/`plan` so
  `rosalind plan --reference --budget` predicts the build peak+time curve (sharing constants with the
  realized accountant, the `plan.rs` discipline); `index --enforce`.
- **D3 — characterize the measured curve** (peak vs wall-time at b = n, n/2, √n) on E. coli + a larger
  genome; the honest answer to OPEN_PROBLEMS §3 ("tolerable time overhead?").
- **D4 — u64 widening** (break `MAX_GENOME_LEN` 4.29 Gbp for wheat/pangenome): SA element type throughout
  `suffix_array.rs` (induce_sort `i32` + `-1` sentinel → `i64`), `SampledSuffixArray.values`, the on-disk
  format value arrays + a **format version bump**, the C-table, every cumulative `u32`, + every persistence
  test. Mechanical but wide; easy to under-scope.
- **D5 — Big-BWT/PFP second mode** for repetitive pangenome inputs.

## 4. Cross-cutting gates + risks

- **Byte-identical is the non-negotiable safety net** at every step (vs `sais_u32` on test inputs; vs the
  monolithic index on E. coli). It is what makes a wrong merge *fail loud* instead of shipping silently.
- **Determinism + thread-count-invariance** must survive spilling (a stated contract guarantee) — temp-file
  layout + merge order must be deterministic.
- **The hard part is the gap-array merge** (D1c). pSAscan/eSAIS are real systems, not refactors — budget
  L–XL, and this is exactly where supervised execution + checkpoints matter.
- **Time may be intolerable at the low-RAM end for huge genomes** (pSAscan reports ~170 h for 200 GiB at
  3.5 GiB RAM); "tolerable" is the open question D3 answers — the honest outcome may be "practical down to a
  budget floor."
- **Subagent shared-tree hazard** (if ever parallelized): never let subagents run `cargo fix`/mutating git;
  verify clean tree + commit-stat at each task boundary.

## 5. Forks for the morning brainstorm

1. **Algorithm:** pSAscan gap-array (recommended) vs DC3 vs straight-to-Big-BWT.
2. **D1a framing:** build the throwaway correctness-reference first (locks the API + gate, safe) — or go
   straight for the efficient pSAscan path and accept slower convergence?
3. **Spill format:** raw `mmap`'d temp files vs a structured on-disk run format; determinism guarantees.
4. **Block-size → budget policy:** how `--memory-budget-mb` maps to block size + the predicted curve.
5. **Scope of the FIRST shippable D1 PR:** D1a alone (correctness skeleton + harness), or D1a+D1b (the first
   real measured peak reduction)?

These notes are committed for review only; nothing here is executed. D1 implementation should begin from a
branch off the merged D0 (it builds on D0's `BuildMemoryModel` + build receipt).
