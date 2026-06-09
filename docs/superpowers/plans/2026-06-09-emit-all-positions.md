# emit_all_positions Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a `PileupParams::emit_all_positions` flag (default `false`) so a per-locus analyzer reports zero-coverage loci, and fix the flagship ColumnKit coverage example to use it — without changing the memory bound.

**Architecture:** The engine already walks the region position-by-position and builds a column at every position; line `engine.rs:407` suppresses the empty ones. The flag makes that suppression conditional (one line). It never touches `advance_to`/`self.active`, so `current_working_set()` is byte-identical with the flag on or off. The whole-genome driver already threads `params` through, so no driver change is needed.

**Tech Stack:** Rust 1.83. The bounded `PileupColumn` kernel (`src/pileup/engine.rs`), the ColumnKit example (`examples/columnkit_coverage.rs`). Tests: `#[cfg(test)] mod tests` in `engine.rs`.

**Spec:** `docs/superpowers/specs/2026-06-09-emit-all-positions-design.md`.

**Verified:** all 11 `PileupParams { … }` literals in the tree use `..PileupParams::default()` (or `..Default::default()`) spread, so adding a defaulted field is non-breaking.

---

## File Structure

| File | Change | Responsibility |
|---|---|---|
| `src/pileup/engine.rs` | Modify `PileupParams` (line 45) + its `Default` (line 67) + the suppression at line 407; add 2 unit tests + an `engine_emit_all` helper in `#[cfg(test)] mod tests` | The flag + the conditional emit |
| `examples/columnkit_coverage.rs` | Modify the `run_bounded_whole_genome` call (≈ line 88) to pass `emit_all_positions: true` | Make the flagship coverage track reference-complete |

---

## Task 1: the `emit_all_positions` flag

**Files:**
- Modify: `src/pileup/engine.rs` (`PileupParams`, its `Default`, line 407; tests in `mod tests`)

- [ ] **Step 1: Write the failing tests**

In `src/pileup/engine.rs`, inside `#[cfg(test)] mod tests` (after the `engine` helper at line 446), add an `engine_emit_all` helper and two tests:

```rust
    fn engine_emit_all(reads: Vec<AlignedRead>, reference: &[u8]) -> PileupEngine<SliceSource> {
        PileupEngine::new(
            SliceSource::new(reads),
            Arc::from(reference.to_vec().into_boxed_slice()),
            0,
            0..reference.len() as u32,
            PileupParams {
                emit_all_positions: true,
                ..PileupParams::default()
            },
        )
    }

    #[test]
    fn emit_all_positions_reports_gaps_at_depth_zero() {
        // Same fixture as `deletion_leaves_a_reference_gap_with_no_observation`:
        // 2M1D2M at ref 0 over an 8 bp reference. Observed: 0,1,3,4. Gaps: 2 (deleted),
        // 5,6,7 (uncovered). With emit_all_positions EVERY position emits — gaps depth 0.
        let reference = b"AAAAAAAA";
        let read = AlignedRead {
            contig: 0,
            pos: Position(0),
            mapq: 60,
            flags: SamFlags::default(),
            cigar: vec![
                CigarOp::new(CigarOpKind::Match, 2),
                CigarOp::new(CigarOpKind::Deletion, 1),
                CigarOp::new(CigarOpKind::Match, 2),
            ],
            seq: Arc::from(b"GGGG".to_vec().into_boxed_slice()),
            qual: Arc::from(vec![30u8; 4].into_boxed_slice()),
        };
        let cols = columns(engine_emit_all(vec![read], reference));
        let positions: Vec<u32> = cols.iter().map(|c| c.locus.pos.0).collect();
        assert_eq!(positions, vec![0, 1, 2, 3, 4, 5, 6, 7], "every position emits");
        for p in [2u32, 5, 6, 7] {
            let col = cols.iter().find(|c| c.locus.pos.0 == p).unwrap();
            assert_eq!(col.depth(), 0, "gap position {p} must be depth 0");
        }
    }

    #[test]
    fn emit_all_positions_preserves_the_working_set_bound() {
        // The flag only decides whether an already-built column is returned; it never
        // touches the active set, so the peak working set is identical on/off.
        let reference = b"AAAAAAAA";
        let reads = || vec![mread(0, b"GG", false), mread(3, b"GG", false)];
        let peak = |mut e: PileupEngine<SliceSource>| -> u64 {
            let mut p = e.current_working_set().bytes;
            while let Some(c) = e.next() {
                c.expect("column");
                p = p.max(e.current_working_set().bytes);
            }
            p
        };
        let off = peak(engine(reads(), reference));
        let on = peak(engine_emit_all(reads(), reference));
        assert_eq!(off, on, "emit_all_positions must not change the working-set bound");
    }
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `cargo test -p rosalind-bio --lib pileup::engine::tests::emit_all 2>&1 | tail -20`
Expected: FAIL — a compile error, `struct PileupParams has no field named emit_all_positions` (the field does not exist yet).

- [ ] **Step 3: Add the field, its default, and honor it at line 407**

In `src/pileup/engine.rs`, add the field to `PileupParams` (after `max_read_len`, line 64):

```rust
    /// When `true`, emit a column at EVERY reference position in the region, including
    /// zero-coverage loci (depth 0) — for reference-complete per-locus analytics
    /// (coverage/QC). Default `false` (covered loci only). The working-set bound is
    /// unchanged: the empty column is already built at each position; this only decides
    /// whether it is returned. A volume tradeoff — a region/panel opt-in, not a
    /// whole-genome default.
    pub emit_all_positions: bool,
```

Add the default to the `Default` impl (after `max_read_len: None,`, line 76):

```rust
            emit_all_positions: false,
```

Make the suppression at line 407 conditional:

```rust
            if self.params.emit_all_positions || !column.obs.is_empty() {
                return Some(Ok(column));
            }
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `cargo test -p rosalind-bio --lib pileup::engine::tests::emit_all 2>&1 | tail -12`
Expected: PASS — `emit_all_positions_reports_gaps_at_depth_zero` and `emit_all_positions_preserves_the_working_set_bound` both pass.

- [ ] **Step 5: Confirm default behavior is unchanged**

Run: `cargo test -p rosalind-bio --lib pileup::engine 2>&1 | tail -6`
Expected: PASS — all existing engine tests (including `deletion_leaves_a_reference_gap_with_no_observation`, which asserts `vec![0, 1, 3, 4]` under default params) stay green.

- [ ] **Step 6: Commit**

```bash
git add src/pileup/engine.rs
git commit -m "feat(pileup): emit_all_positions — bounded reference-complete output

PileupParams.emit_all_positions (default false) makes engine.next() emit a
column at every reference position, including zero-coverage loci (depth 0), for
reference-complete per-locus analytics. The working-set bound is unchanged: the
empty column is already built at each position; the flag only decides whether it
is returned (test-pinned). Non-breaking (all literals use ..Default).

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: fix the flagship coverage example

**Files:**
- Modify: `examples/columnkit_coverage.rs` (the `run_bounded_whole_genome` call, ≈ line 88)

- [ ] **Step 1: Set `emit_all_positions: true` in the example**

In `examples/columnkit_coverage.rs`, replace the `run_bounded_whole_genome` call's `PileupParams::default()` argument:

```rust
    // emit_all_positions makes the coverage track REFERENCE-COMPLETE: every locus
    // emits, so uncovered positions (here ref 12–15, which no read reaches) report
    // depth 0 instead of being silently dropped. The memory bound is unchanged.
    let (ws, _skips) = run_bounded_whole_genome(
        &mut analyzer,
        SliceSource::new(reads),
        &ref_view,
        contigs,
        PileupParams {
            emit_all_positions: true,
            ..PileupParams::default()
        },
        &mut out,
    )?;
```

- [ ] **Step 2: Run the example and confirm the zero-depth tail is reported**

Run: `cargo run --quiet --example columnkit_coverage 2>/dev/null | grep -cE '\t(13|14|15|16)\t0$'`
Expected: `4` — the 16 bp reference's uncovered tail (1-based positions 13–16, which no read reaches) now reports depth 0 (before this change, the example printed nothing for them).

- [ ] **Step 3: Commit**

```bash
git add examples/columnkit_coverage.rs
git commit -m "fix(example): coverage track reports zero-coverage loci (emit_all_positions)

The flagship ColumnKit coverage example silently dropped every uncovered base on
sparse data (it only worked because the toy reads tiled most of the reference).
emit_all_positions makes it reference-complete — uncovered loci report depth 0.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 3: full gate

**Files:** none (verification only)

- [ ] **Step 1: Full suite + lint + format**

Run:
```bash
cargo test
cargo clippy --all-targets -- -D warnings
cargo fmt --check
```
Expected: PASS on all three. (If `cargo fmt --check` reports diffs, run `cargo fmt` and amend the relevant commit.)

- [ ] **Step 2: No commit** (verification only).

---

## Self-Review

**1. Spec coverage:**
- §3.1 the field + `Default` → **Task 1 Step 3**.
- §3.2 honor at line 407 → **Task 1 Step 3**.
- §3.3 fix the example → **Task 2**.
- §4 bound preserved → **Task 1 Step 1/4** (`emit_all_positions_preserves_the_working_set_bound`).
- §5 tests: new emit-all gap test (Task 1), default-unchanged (Task 1 Step 5), bound-preserved (Task 1), example runs + prints depth-0 tail (Task 2 Step 2), gates (Task 3).
- §6 non-goals: no CLI flag, not wired into variants/gvcf, no driver change — respected (only `PileupParams` + `engine.rs:407` + the example are touched).

**2. Placeholder scan:** No TBD/TODO. Every step has literal code or an exact command + expected output.

**3. Type consistency:** `emit_all_positions` (the field) is defined in Task 1 Step 3 and used identically in the `engine_emit_all` helper, both tests, and the example (Task 2). The test helpers (`engine`, `columns`, `mread`, `PileupEngine::new`, `current_working_set`, `col.depth()`, `c.locus.pos.0`) all match the existing `mod tests` API. The example's `run_bounded_whole_genome` signature is unchanged (still takes `params: PileupParams` by value).
