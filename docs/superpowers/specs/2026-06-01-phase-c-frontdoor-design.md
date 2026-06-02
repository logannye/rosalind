# Move #4 — The Front Door (positioning + discoverability) (design)

**Status:** Spec for review — 2026-06-01. Strategy Move #4 (the 2026-06-01 strategy synthesis; companion
to `docs/OPEN_PROBLEMS.md`): convert the *watching* fork wave into *building* forkers by making the shipped
breakthrough **discoverable** and the repo's front door **point at it**. Built on a branch stacked on
`rosalind/phase-c-contract` (the docs reference the contract verbs `plan`/`--enforce`/`verify`, which land
with Phase C). Pure positioning / docs / re-exports / examples — **no code-behavior changes.**

## 1. Why (the bottleneck)

The analysis found the community is in *watching* mode — every GitHub fork is byte-identical to or behind
`main` (zero commits on top). The repo's own front door actively misdirects them: `src/lib.rs:1` opens on
*"O(√t) Space Simulation via Height Compression"* with a `TuringMachine` usage example (so `cargo doc`
lands on the wrong product); the crate-root re-exports are **only** the √t theory types; **all** four
`examples/*.rs` are √t simulation demos; the toy generator emits a single `>chrToy` (no runnable
multi-contig `variants --index` demo); and the README "Extend" section routes builders to the
**non-bounded** `GenomicPlugin` path. This move fixes the front door so a drive-by forker's first five
minutes land on the bounded-memory contract — the genuinely unique, now-shipped capability.

## 2. Resolved positioning decisions (2026-06-01 brainstorm)

- **√t framing: contract-first, √t as an honest research footer.** README + `lib.rs` lead with the
  bounded-memory contract / genomics product (the shipped, tested capability). √t appears as a clearly
  labeled "Research direction (Phase D)" section — framed as future / not-yet-load-bearing, **no
  overclaiming** a layer that is currently a stub. (Anchored on Williams' peer-reviewed O(√(t log t)); the
  withdrawn arXiv 2508.14831 is never cited.)
- **Plugin path: demote + label, do NOT remove.** Route the front door to the bounded substrate (the
  `PileupColumn` iterator + `ReadSource`/`call_germline_whole_genome`) as THE way to build bounded
  analytics, with a runnable cookbook example. Clearly label the `GenomicPlugin` trait + `framework/` + the
  Python RNA-seq demo as **legacy / non-bounded** (still works; does NOT inherit the memory contract). No
  `#[deprecated]`, no code removal.

## 3. The seven deliverables

### 3.1 Crate-root re-exports (`src/lib.rs`)
Add a curated **genomics product surface** so `use rosalind::{…}` lands on the engine, not the theory
layer. The existing √t re-exports are KEPT (regrouped under a `// Research layer (√t)` comment). New
re-exports (all already reachable via module paths today — this is the convenience + the signal of "this
is the product"):
- substrate: `PileupEngine, PileupColumn, Obs, ReadSource, SliceSource, StreamingBamSource, PileupParams`
- calling: `call_germline_whole_genome, call_germline_region_streaming, GermlineCall, GermlineParams`
- contract: `MemoryBudget, WorkingSet, estimate_variants_working_set, predicted_peak_rss_bytes`
- index + receipt: `GenomeIndex, IndexReader, ReferenceView, RunManifest`

(Curation rule: re-export what a builder *composes on* — the substrate, the bounded drive, the contract
types, the index reader, the receipt. Do NOT re-export internal/legacy types.)

### 3.2 `src/lib.rs` top-level rustdoc rewrite
Replace the `//! # O(√t) Space Simulation via Height Compression` opener (lines 1–26) with a
genomics-contract lead: what Rosalind is (deterministic, low-memory genomics engine; declare your RAM →
`plan` → `--enforce` → `verify`; bounded whole-genome calling) and a **runnable doctest** that builds a
`PileupEngine` over a tiny in-memory `SliceSource` and pulls a `PileupColumn` (covered by
`cargo test --doc`). Add a clearly-labeled `## Research direction (Phase D)` section that honestly
describes the √t space-bounded-construction ambition as future work (not yet load-bearing). Also fix the
false docline at `src/pileup/mod.rs` that claims plugins build on the streaming engine — state plainly
that the bounded contract applies to the germline/pileup path and the plugin/framework lineage is
non-bounded.

### 3.3 `CONTRACT.md` (new, repo root)
The authoritative contract document, linked from the README:
- the verbs: **declare** a budget → **`rosalind plan`** (predict before committing) → **`--enforce`**
  (honor: refuse exit 3 / fail exit 4, never a silent OOM-kill) → **`rosalind verify`** (re-check the
  receipt without re-running);
- the honest brand line: *"never silently OOM-kills you — it fits, or it tells you up front, and proves
  the realized peak with a receipt"* (explicitly NOT "never refuses" — graceful degrade/spill is Phase D);
- an **Extend** section routing builders to the `PileupColumn` iterator substrate (pointing at
  `examples/custom_pileup_analytics.rs`) and labeling the `GenomicPlugin`/`framework/`/Python-RNA-seq
  lineage as legacy / non-bounded;
- honest scope: the contract is the **germline `variants --index`** path; somatic is region-bounded;
  index build is O(reference) (Phase D).

### 3.4 README rewrite (`README.md`)
- Lead with the one-command contract story:
  `rosalind plan` → `rosalind variants --index --enforce` → `rosalind verify`.
- Rewrite the **Extend** section per §2 (substrate primary + cookbook; plugin lineage labeled legacy).
- Replace any "never refuses" / over-claim language with the honest brand from §3.3.
- Add the multi-contig demo (§3.5) to the runnable examples.
- Keep √t as a short "Research direction (Phase D)" footer, linking `docs/OPEN_PROBLEMS.md`.
- Link `CONTRACT.md`.

### 3.5 Multi-contig demo fixture (`scripts/generate_toy_data.py`)
Extend the generator to emit a small **2–3 contig** reference (it writes a single `>chrToy` today) so the
flagship path runs out of the box: `index → sort → plan → variants --index --enforce → verify`. Keep it
deterministic (seeded) and tiny. Update the `SHA256SUMS`/manifest emission accordingly. (A regenerated
`examples/data/` multi-contig fixture may be committed so the README commands run without invoking Python.)

### 3.6 Substrate cookbook example (`examples/custom_pileup_analytics.rs`, new)
A **non-caller** consumer that computes a per-locus metric (e.g. coverage + a simple QC count) directly
over the `PileupColumn` iterator from a `SliceSource` — demonstrating the substrate as a platform pattern
that inherits bounded memory + determinism for free, with no variant calling involved. Must compile and
run via `cargo run --example custom_pileup_analytics`.

### 3.7 issue #3 rewrite (GitHub — confirm-first)
Rewrite the public roadmap issue (currently says "separate out the theory layer as a demo," contradicting
the load-bearing-√t thesis) to the contract framing: "memory as a declared, predicted, honored, verifiable
contract — beachheaded on the bounded whole-genome caller — with √t as the FUTURE space/time knob (Phase
D)." **Outward-facing: draft the new body, show it for approval, and only post on explicit go-ahead.**

## 4. Testing

- `cargo test --doc` — the new `lib.rs` doctest builds a `PileupEngine` and pulls a column (proves the
  documented product surface actually compiles + runs).
- `cargo run --example custom_pileup_analytics` — the cookbook compiles and runs.
- A smoke test (subprocess, `tests/`) drives the **multi-contig** fixture through
  `index → sort → variants --index → (plan / verify)` and asserts success + a multi-contig VCF — so the
  README's headline commands are proven, not just claimed.
- README / `CONTRACT.md` are prose, but every command shown must be copy-pasteable-correct (verified by
  the smoke test + manual run).
- Gates unchanged: `cargo fmt --all -- --check`, `cargo build` 0 warnings (debug + release), full suite.

## 5. Non-goals

No code-behavior changes (this is positioning/docs/re-exports/examples only); no `#[deprecated]` markers or
removal of the plugin/framework lineage; the PyO3 zero-copy NumPy pileup-channel stream (the larger Phase-E
substrate wedge) is OUT; the flagship real-genome artifact is **Move #5** (separate). MSRV 1.72 preserved;
no new dependencies.

## 6. Branch / sequencing

Built on `rosalind/phase-c-frontdoor`, stacked on `rosalind/phase-c-contract` (PR #21). The README/CONTRACT
reference the contract verbs, so this merges *after* (or together with) Phase C. When Phase C lands on
`main`, this rebases trivially.
