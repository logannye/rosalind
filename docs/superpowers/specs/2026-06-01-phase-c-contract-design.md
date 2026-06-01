# Phase C — memory as a verifiable contract (variants path) (design)

**Status:** Spec for review — 2026-06-01. The Phase-C sub-stage that turns Rosalind's
already-real *memory receipt* into a *contract* on the bounded whole-genome variants path:
**declare → predict → honor-or-refuse → verify.** Builds directly on the Phase-B4 caller
(`rosalind variants --index`, merged — PR #19). Under the contract thesis in
[`docs/OPEN_PROBLEMS.md`](../../OPEN_PROBLEMS.md) and the strategy reframe
(the contract is the shippable breakthrough; √t is the future Phase-D knob).

## 1. The capability we are shooting for

> **You declare a RAM budget; `rosalind plan` tells you *before you commit a byte* whether your
> whole-genome germline call will fit; the run *honors* that budget — fitting cleanly or refusing
> cleanly, never silently OOM-killed mid-job; and `rosalind verify` re-checks a deterministic,
> BLAKE3-stamped receipt proving the realized peak landed inside your budget and the VCF came from
> exactly these inputs.**

Today (post-B4) the engine already *measures* and *records* peak memory — `peak_rss_bytes()` is a real
`getrusage` signal (`util/rss.rs`), and the manifest carries `peak_rss_bytes` + `max_working_set_bytes`
as BLAKE3-stamped canonical JSON (`provenance/mod.rs`, written at `main.rs:1083-1094`). But the budget
check is **record-only** (`main.rs:1109-1120` computes `budget.admits(peak_rss)` then prints
*"record-only, run completed"*), there is **no pre-run prediction** for the streaming path, the
**stdout path writes no manifest** (`main.rs:1097-1101`), and the reported `max_working_set_bytes` is a
**real under-count** of true peak. Phase C closes exactly those four gaps. No incumbent ships this
four-property contract; it is unique *without* √t (which extends the same contract down to index
construction in Phase D).

## 2. The contract math — the four working-set terms

For the bounded whole-genome variants drive (`call/whole_genome.rs`), true peak working set is the sum
of four terms. Phase C makes each one either *exactly known* or *hard-bounded*, so a prediction can be a
true upper bound and the realized high-water can be measured exactly.

| Term | Bytes | Knowable a priori? | How Phase C bounds it |
|---|---|---|---|
| **Reference** (current contig) | `largest_contig_len` | ✅ from index header | decode-then-**move** into the `Arc` — eliminates the `buf`+`Arc` double-copy at `whole_genome.rs:63-64` (**2× → 1×**) |
| **Active read set** | `D × max_read_len × per_base_cost` | ⚠️ depth **capped**, read-len **assumed** | **hard cap at `--max-depth D`** (deterministic); read length is the one residual assumption |
| **Resident VCF rows** | O(1) (writer buffer) | ✅ constant | **stream each call to the writer as produced** — no genome-wide row Vec; the only resident "rows" memory is the `BufWriter`'s fixed buffer |
| **Fixed overhead** | const | ✅ | constant |

Two numbers fall out:

- **Predicted** (`rosalind plan`, pre-run, from index header + declared `D` + assumed max-read-len `L`):
  a true upper bound *modulo* the read-length assumption.
- **Realized** (the receipt's high-water): **exact** — sampled from a *corrected* accountant that
  includes all four terms.

**The honest guarantee.** Depth is hard-capped (guaranteed), rows are streamed (guaranteed small),
reference is exactly the largest contig (1× after the move fix). The only residual assumption is max
read length — and `--enforce` **fails loud post-run if realized > budget** (a clean non-zero exit, never
a silent overrun). So the brand is **"fits-or-tells-you up front, and proves the realized peak with a
receipt"** — explicitly **not** "never refuses" (graceful degrade/spill needs the Phase-D ladder and is
out of scope here; without it `--enforce` can only refuse cleanly).

## 3. Scope

**In (the variants `--index` path only):**
- **C1 — soundness fix.** Make `max_working_set_bytes` a true conservative upper bound; deterministic
  `--max-depth` cap; per-call streaming VCF flush; decode-then-move reference (2×→1×). *Output-preserving
  by default* (cap defaults off in C1).
- **C2 — `rosalind plan` + `--enforce`.** A pre-run streaming estimator; the `plan` subcommand;
  honor-or-refuse enforcement (refuse-up-front-on-predicted / fail-loud-post-run-on-realized). The
  `--max-depth` default (1000) lands here.
- **C3 — `rosalind verify` + receipt-on-stdout + CI contract gate.** A parser for the canonical manifest;
  the `verify` subcommand; always-persist-a-receipt; the CI contract suite.

**Out (deferred, by design):**
- **Index-build enforcement.** The build is O(reference); `rosalind index --memory-budget-mb` /
  `plan --reference` stay **advisory/record-only** (enforcement waits for Phase-D sublinear construction).
- **Somatic.** `call/pipeline.rs` collects both column streams into `Vec`s — region-bounded, not
  whole-genome-bounded. A region-bound refactor is out of scope; the bounded contract is scoped to the
  germline `--index` path and the docs say so.
- **Graceful degrade / external-memory spill.** The Phase-D √t/spill ladder. `--enforce` here refuses
  cleanly; it does not degrade.
- **On-demand `ref_base`** (shrinking the reference term to O(pileup window)) → Phase D/E.
- **Deterministic multithreading / thread-invariance.** The engine is single-threaded; we do not claim
  thread-invariance (it would be vacuously true) until a parallel path + gate exist (Phase E).

## 4. Decomposition (one spec → three green sub-stages, one PR each)

Mirrors the B3/B4 pattern. Each sub-stage lands green (`cargo test`, `cargo fmt --check`, 0 warnings)
with its own plan + PR.

```
Phase C — memory as a verifiable contract (variants path)
├─ C1  soundness fix              [the hard correctness core]
├─ C2  rosalind plan + --enforce
└─ C3  rosalind verify + stdout receipt + CI contract gate
```

## 5. C1 — the soundness fix (the hard correctness core)

The prerequisite for everything: until the modeled number is a true upper bound, any `plan` lies and the
brand inverts on first contact. C1 is **output-preserving by default** (the cap defaults to `None`); it
changes only the *accounting*, the *reference copy count*, and the *row-flush plumbing*.

### 5.1 Corrected working-set accountant (`pileup/engine.rs:135`)
`PileupEngine::current_working_set()` today sums only `ref_to_read.len()*16 + 64` per active read + 256.
It omits (a) the engine's own `self.reference` bytes and (b) each active read's `seq`/`qual` bytes. Fix it
to count what is actually resident:
- `self.reference.len()` (the decoded contig the engine holds),
- per active read: `ref_to_read.len()*16` (map) **+ `seq.len()` + `qual.len()`** (the byte buffers) + a
  small per-read constant,
- fixed engine overhead.

This is the single load-bearing correctness change. (It only ever *grows* the reported number — it
cannot newly pass a budget it failed before.)

### 5.2 Deterministic depth cap (`pileup/engine.rs`)
Add `PileupParams.max_depth: Option<u32>` (**default `None` in C1** — no behavior change; C2 sets the CLI
default). When `Some(D)` and the active set already holds `D` reads, **drop the incoming read** and count
it under a new `SkipCounts.over_max_depth`. Because the source is `(contig,pos)`-sorted (and `SliceSource`
sorts on construction), "the first `D` reads to overlap a position" is deterministic ⇒ capped output is
deterministic and order-independent. Output changes *only* at sites deeper than `D` (the artifact
pileups) — standard caller behavior. Applied in `advance_to`/`ingest`.

### 5.3 Incremental VCF writer (`io/vcf.rs`)
Split the monolithic `write_germline_vcf(writer, contigs, sample, &rows)` into:
- `write_germline_header(writer, contigs, sample)` — header + `##contig` lines (once),
- `write_germline_row(writer, &row)` — one record.

`write_germline_vcf` is retained as a thin wrapper (header + loop) for back-compat / the single-contig
path. Rows are produced in `(contig,pos)` order, so streaming them is **byte-identical** to the batch
write — the `golden_vcf` + `determinism` gates stay green.

### 5.4 Streaming sink + true high-water (`call/whole_genome.rs:47`)
Change the drive from returning `(Vec<rows>, WorkingSet)` to driving a **row sink** and returning the
*true* high-water working set:
```rust
pub fn call_germline_whole_genome<S: ReadSource>(
    source: S,
    ref_view: &ReferenceView,
    contigs: &ContigSet,
    pileup_params: PileupParams,
    germline_params: &GermlineParams,
    on_row: &mut dyn FnMut((Locus, u8, GermlineCall)) -> Result<(), CoreError>,
) -> Result<WorkingSet, CoreError>
```
- **decode-then-move:** `decode_window(.., &mut buf)` then `Arc::<[u8]>::from(std::mem::take(&mut buf))`
  (moves the allocation; no second copy). Reference term 2× → 1×.
- each emitted `(Locus, ref_base, call)` is passed straight to `on_row` (no genome-wide `rows` Vec).
- the returned `WorkingSet` is `max over contigs of (engine peak working set + resident-row bytes for the
  in-flight write batch)` — a sound upper bound on realized peak.

**API-shape decision:** a **sink callback**, not an `Iterator`. It is bounded, trivial to retrofit, and a
clean substrate primitive in its own right (a consumer passes *any* sink: stream-to-VCF, collect,
featurize-into-tensors). An `Iterator` is more idiomatic but fights the borrow checker across the
per-contig reference-decode / contig-switch state machine; a `.calls()` iterator adapter over this sink is
an easy future addition for the front-door cookbook.

`main.rs::run_variants_index` (`main.rs:1044-1101`) calls `write_germline_header` once, then passes a sink
that `write_germline_row`s into the (file or stdout) writer; the returned `WorkingSet` feeds the receipt.

### 5.5 C1 tests (the proof)
- **Real-RSS soundness gate** (`tests/`): on a generated growing input with an explicit `--max-depth D`,
  assert realized `peak_rss` is flat as the BAM grows **and** the realized working-set ≤ the modeled
  bound.
- **Depth-cap determinism**: capped output byte-identical across shuffled input order; sites > `D` are the
  only ones changed; `over_max_depth` counted.
- **Accountant completeness**: `current_working_set` includes reference + seq/qual (unit test on a known
  active set).
- **VCF byte-identity**: header-split + streamed rows == the old batch write (extend `golden_vcf`).

## 6. C2 — `rosalind plan` + `--enforce`

### 6.1 Streaming estimator (`src/call/plan.rs`, new — pure, unit-tested)
A sibling to `genomics/index/report.rs::estimate_build_working_set`, for the call path:
```rust
pub fn estimate_variants_working_set(largest_contig_len: u64, max_depth: u32, max_read_len: u32) -> WorkingSet
```
`predicted = largest_contig_len (ref 1×) + D·L·per_base_cost (active) + writer_buf_const + fixed`, using the
*same* per-base/per-read constants as the corrected C1 accountant (single source of truth — a shared
`const`/helper so the estimate and the realized accountant cannot drift). Needs **only** the index header
(`contigs.iter().map(|c| c.length).max()`) plus declared `D`/`L` — it never opens the BAM. Plus a pure
renderer `render_variants_plan(estimate_breakdown, budget) -> String` producing the multi-line
`[FITS]`/`[REFUSE]` breakdown.

### 6.2 `rosalind plan` subcommand (new `Commands::Plan`)
- `plan --index <idx> [--max-depth D] [--max-read-len L] [--budget-mb B]` → predicts the **variants** peak
  (the flagship): reference / active / rows / fixed → predicted peak / budget → `[FITS]`/`[REFUSE]`.
- `plan --reference <fa> [--budget-mb B]` → reuses `estimate_build_working_set` for the **build**
  (advisory; build is O(reference) and Phase-D enforces). `--index` XOR `--reference`.
- Output is pure-rendered + deterministic (testable, like `render_plan_line`); exit 0 always (planning is
  advisory; refusal is `variants --enforce`).

### 6.3 `--enforce` on `variants --index` (new flag on `Commands::Variants`)
Replaces the record-only verdict (`main.rs:1109-1120`) when set; both checks require `--memory-budget-mb`:
- **Pre-run:** `estimate_variants_working_set(largest_contig, D, L)` > budget → **refuse**, exit **3**,
  actionable stderr ("declared B MiB; predicted peak P MiB @ max-depth D / max-read-len L; raise
  --memory-budget-mb, lower --max-depth, or drop --enforce"). No work performed.
- **Post-run:** realized `peak_rss` > budget → write the VCF **and** the manifest (recording the
  violation), then **fail loud**, exit **4**. (The artifact + the proof-of-overrun are preserved; the
  exit code signals the violated contract.)
- **Without `--enforce`:** today's record-only behavior, exit 0 (back-compat).

A small `ContractStatus` / exit-code helper centralizes the codes (`0` ok, `3` predicted-over-refused,
`4` realized-over-completed) and the honest messages.

### 6.4 Depth-cap default lands here
C2 sets the CLI `--max-depth` **default to 1000** (always applied ⇒ the engine is bounded by default);
`--max-depth 0` = uncapped opt-out (and then `plan` cannot promise a hard bound — say so). `--max-read-len`
defaults to a short-read value (e.g. 250) with a flag to raise it for long reads. This is the one
deliberate default-output change in Phase C; it is documented (changes calls only at >1000× artifact
sites) and motivated by the contract surface.

### 6.5 C2 tests
`estimate_variants_working_set` monotonic + no-overflow; the breakdown renderer `[FITS]`/`[REFUSE]`;
`plan --index` / `plan --reference` subprocess tests; `--enforce` refuses up front (exit 3) on a tiny
budget with `--max-depth`; `--enforce` fails post-run (exit 4) on a budget below realized peak;
record-only path unchanged without `--enforce`.

## 7. C3 — `rosalind verify` + receipt-on-stdout + CI contract gate

### 7.1 Persist a receipt on stdout runs (`main.rs:1097-1101`)
Today the stdout branch writes no manifest, so "every run emits a receipt" is false on the default path.
Fix: stdout output → write a `rosalind.variants.manifest.json` sidecar in the cwd + announce it on stderr;
`--manifest <path>` redirects (works for both stdout and file output). Document the cwd-sidecar behavior.

### 7.2 Self-describing manifest params
Extend the params written at `main.rs:1083-1094` with `memory_budget_mb` (if declared), `contract_verdict`
(`within`/`over`/`unset`), `enforced` (bool), `max_depth`, `max_read_len`. A receipt then carries
everything `verify` needs without the user re-supplying flags.

### 7.3 Canonical-manifest parser (`provenance/mod.rs`)
Add `RunManifest::from_canonical_json(&str) -> Result<RunManifest, ManifestError>` — a small parser for the
**fixed** canonical shape we emit (not a general JSON parser), keeping deps lean (no `serde_json` in a
public surface) and matching the hand-rolled writer. Guarded by a `serialize → parse → serialize == identity`
property test over varied manifests.

### 7.4 `rosalind verify` subcommand (new `Commands::Verify`)
`verify --manifest <path> [--budget-mb B]` — read-and-assert, **no re-run**:
- parse the manifest; re-hash each listed input/output file (BLAKE3, streamed via `blake3_file`) and
  compare to the recorded digest → detects drift / proves "this VCF came from exactly these inputs";
- re-check recorded `peak_rss_bytes` ≤ recorded `memory_budget_mb` (or the supplied `--budget-mb`);
- exit 0 if all hold; non-zero with a per-check report otherwise. Missing files / hash mismatch / over
  budget each produce a clear, distinct message.

### 7.5 CI contract suite
Extend `tests/rss_budget.rs` into a subprocess contract suite (style of `tests/index_cli.rs`):
1. predicted envelope ≥ realized peak on a generated growing input;
2. realized working-set flat as the BAM grows (bounded-by-coverage);
3. `--enforce` refuses up front (exit 3);
4. `--enforce` fails post-run on a tiny budget (exit 4);
5. `verify` round-trips a real run's manifest and catches a tampered output hash.

Scope = add these to the existing suite (which CI already runs). A full CI-matrix overhaul (macOS runners,
clippy gate) stays Phase-E; clippy remains a Phase-E call under MSRV 1.72.

## 8. Full CLI surface after Phase C

```
rosalind plan (--index <idx> | --reference <fa>) [--max-depth D] [--max-read-len L] [--budget-mb B]
rosalind variants --index <idx> --alignments <sorted.bam>
    [--mapq-threshold N] [--quality-threshold Q]
    [--max-depth D] [--max-read-len L]
    [--memory-budget-mb M] [--enforce] [--manifest <path>] [-o out.vcf]
rosalind verify --manifest <path> [--budget-mb B]
```
`variants --reference` (single-contig) is unchanged; `--enforce`/`--max-depth`/`--max-read-len`/`--manifest`
are additive and back-compatible (absent ⇒ today's behavior, except the `--max-depth` default of 1000 from
C2).

## 9. Determinism & honesty constraints (cross-cutting)

- **Byte-identical output preserved.** The writer split + per-call streaming must not change VCF bytes
  (gated by `golden_vcf` + `determinism`). The depth cap is deterministic; capped output is reproducible
  and order-independent.
- **The estimate and the realized accountant share one set of constants** (§6.1) so a passing `plan` and
  the realized receipt cannot silently diverge.
- **Branding is honest in every message and doc:** "fits-or-tells-you / never silently OOM-kills," never
  "never refuses." The README/issue-3/`lib.rs` positioning fixes are tracked separately (front-door move),
  not in this spec.

## 10. Per-sub-stage gates

Each of C1/C2/C3: `cargo test` green (full suite), `cargo fmt --all -- --check`, `cargo build` 0 warnings
(debug **and** release), no htslib types in any new public signature, and the sub-stage's own gates
(C1 §5.5, C2 §6.5, C3 §7.5). Final whole-branch review per the established workflow.

## 11. Non-goals (so the promises stay airtight)

Index-build enforcement; somatic whole-genome bounding; graceful degrade / spill; on-demand `ref_base`;
deterministic multithreading / thread-invariance claims; out-accuracy-ing GATK/DeepVariant (the demo is a
memory/reproducibility contract win on SNVs, never an accuracy head-to-head).

## 12. Resolved decisions (2026-06-01 brainstorm)

- **Working-set model:** declared depth ceiling + deterministic cap (realized ≤ predicted by construction
  for the coverage term; read length is the residual assumption with the realized post-run check as the
  backstop).
- **Scope:** full contract on the variants path, decomposed C1/C2/C3; index build stays record-only;
  somatic out of scope.
- **Row flush:** stream each call to the writer (resident rows = O(batch)), via a **sink callback** on
  `call_germline_whole_genome` (not a per-contig batch, not an iterator).
- **Exit codes:** `3` = predicted over budget (refused, no work); `4` = realized over budget (completed,
  contract violated).
- **`--max-depth` default = 1000** (lands in C2; `0` = uncapped); `--max-read-len` default short-read (250).
- **Stdout runs** write a cwd manifest sidecar by default; `--manifest <path>` redirects.
- **Manifest parsing** for `verify` = a small hand-parser for the fixed canonical shape (no `serde_json`),
  round-trip property-tested.
