# Design: `emit_all_positions` — bounded, reference-complete per-locus output

**Status:** Approved design — 2026-06-09. Audience: the implementer. Companion to
[`CONTRACT.md`](../../../CONTRACT.md) (the bounded `PileupColumn` kernel) and the ColumnKit SDK.

---

## 1. The goal, in one sentence

Let a per-locus analyzer report **zero-coverage loci** (depth 0) instead of silently dropping
them — without changing the memory bound — retiring a standing honesty liability: the README's
flagship ColumnKit coverage example omits every uncovered base on real sparse data.

## 2. Why this, why now

`PileupEngine::next` (`src/pileup/engine.rs:397`) already walks the region **position by position**
(`while self.pos < self.region.end { … self.pos += 1; }`) and calls `build_column()` at **every**
position. At a coverage gap the column is built with empty observations, and line 407 —
`if !column.obs.is_empty()` — **suppresses** it (the `while` loop continues). So every uncovered
locus is visited and constructed, then dropped.

The shipped coverage example (`examples/columnkit_coverage.rs`) only *looks* correct because its toy
reads happen to tile most of the reference: reads at positions 0 and 4 (length 8) cover 0–11 of a
16 bp reference, so **positions 12–15 are silently omitted**. A real coverage/QC track must report
those as depth 0. For an honesty-branded repo, a flagship example that is silently wrong on sparse
input is a credibility liability.

## 3. The change

### 3.1 `PileupParams::emit_all_positions: bool` (default `false`)

Add the field to `PileupParams` (`engine.rs:45`) and `emit_all_positions: false` to its `Default`
(`engine.rs:67`). Behavior is unchanged for every existing caller (`PileupParams::default()` and any
struct literal that uses `..Default::default()`); literals that name every field are not used in the
codebase, so this is non-breaking.

### 3.2 `engine.rs:407` honors the flag

```rust
        if self.params.emit_all_positions || !column.obs.is_empty() {
            return Some(Ok(column));
        }
```

`self.params` is the engine's stored `PileupParams` (`engine.rs:153`). A gap position then emits a
valid `PileupColumn` with the correct `ref_base` (decoded from the reference at that position) and
`depth() == 0`.

### 3.3 Fix the flagship example

In `examples/columnkit_coverage.rs`, pass `emit_all_positions: true` so the coverage track reports
positions 12–15 as depth 0:

```rust
        PileupParams {
            emit_all_positions: true,
            ..Default::default()
        },
```

(plus a one-line comment that this makes the track reference-complete, and an `eprintln!`/comment
noting the uncovered tail now reports depth 0).

## 4. The memory bound is provably preserved

`emit_all_positions` never touches `advance_to` or `self.active` (the held read set); the column at
each position is built **before** the line-407 check regardless of the flag. The flag only decides
whether the already-built column is **returned** or **skipped**. Therefore `current_working_set()`
is byte-identical with the flag on or off, at every position.

**Honest scope note:** emit-all is a *volume* tradeoff — a whole-genome emit-all stream is large
even though peak memory stays bounded by coverage — so it is a **region/panel-oriented opt-in**, not
a whole-genome default. (Hence default `false`.)

## 5. Testing

In `src/pileup/engine.rs`'s `#[cfg(test)] mod tests` (mirroring the existing position-list tests):

- **Emits a column at every position when set (NEW test):** over a region with a gap (reads leaving
  an uncovered position), `emit_all_positions: true` yields a column at **every** position; the gap
  column has `depth() == 0` and the correct `ref_base`. Where the existing default test asserts a
  gapped position list (e.g. `vec![0, 1, 3, 4]`), the new emit-all variant asserts the complete list
  (`vec![0, 1, 2, 3, 4]`).
- **Default behavior unchanged:** the existing suppression tests are **left as-is** and stay green by
  construction (`emit_all_positions: false` still suppresses the gap position).
- **Bound preserved:** drive the same sparse fixture through two engines (flag on vs off) and assert
  the **peak `current_working_set().bytes`** is identical — the provable claim, test-pinned.

**Gates:** `cargo test`, `cargo clippy --all-targets -- -D warnings`, `cargo fmt --check`, and
`cargo run --example columnkit_coverage` (the example compiles and now prints depth-0 rows for the
uncovered tail).

## 6. Scope / non-goals

- **In scope:** the `PileupParams` flag, the one-line `engine.rs` change, the example fix, tests.
- **No CLI flag** (`features --all-positions` / a `variants` flag) in this PR — a separate decision;
  exposing genome-wide emit-all is a volume footgun to gate carefully.
- Not wired into `variants`/`gvcf` (empty rows in a VCF are wrong); the flag is for per-locus
  analytics (ColumnKit / `features`).
- No change to `run_bounded_whole_genome` — it already threads `params` through to the per-contig
  engine, so setting the flag in the params it receives is sufficient.
