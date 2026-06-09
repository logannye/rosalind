# predicted_working_set Claim + verify Soundness Re-check Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the memory contract's core inequality — the prediction upper-bounds reality — a content-addressed, offline-auditable receipt property: record the deterministic `predicted_working_set_bytes` in the claim, and have `verify` re-check `predicted ≥ realized` offline (hard under `--enforce`, advisory otherwise).

**Architecture:** Part 1 adds one claim param in `run_variants_index`'s receipt block, computed by the existing pure `estimate_variants_working_set` (no baseline → cross-machine stable; not in `MEASUREMENT_KEYS` → stays in the claim). Part 2 adds a four-arm match to `verify_receipt` that re-derives the inequality from the recorded claim + measurement, failing only an enforced violation. Pure-`std`; no schema bump.

**Tech Stack:** Rust 1.83, the `rosalind-receipt` leaf crate (hand-rolled canonical JSON, `blake3`). Tests: `#[cfg(test)]` unit tests in the receipt crate + a CLI integration test in `tests/plan_enforce.rs` (`env!("CARGO_BIN_EXE_rosalind")`).

**Spec:** `docs/superpowers/specs/2026-06-09-predicted-working-set-claim-design.md`.

---

## File Structure

| File | Change | Responsibility |
|---|---|---|
| `src/main.rs` | Modify `run_variants_index` receipt block (≈ line 2260, after the `max_working_set_bytes` insert) | Record `predicted_working_set_bytes` as a claim param |
| `crates/receipt/src/lib.rs` | Modify `verify_receipt` (≈ line 755, after the `max_working_set ≤ peak_rss` check) + add unit tests in `#[cfg(test)] mod tests` | The offline soundness re-check |
| `tests/plan_enforce.rs` | Append one integration test | Prove the field lands in the claim, is deterministic, and `verify` passes |

**Key facts (verified in the code):**
- `estimate_variants_working_set(largest_contig_len, max_depth, max_read_len) -> WorkingSet` is a pure fn at `src/call/plan.rs:23`; call it as `rosalind::call::plan::estimate_variants_working_set(...)` (same path style as the existing `predicted_peak_rss_bytes` use). `.bytes` is the `u64`.
- In `run_variants_index`, `largest`, `max_depth`, `max_read_len` are all in scope at the receipt block.
- `MEASUREMENT_KEYS` (`lib.rs:100`) does NOT include `predicted_working_set_bytes`, so it stays in the claim through `finalize()`.
- In `verify_receipt`: `recorded_ws` (= `max_working_set_bytes`, `Option<u64>`) is already parsed (≈ line 725); `parse_num(key, &mut problems)` is the strict in-scope closure; `notes`/`problems` are the accumulators; `manifest.params.get("enforce")` reads the enforce flag (a claim param, set only under `--enforce`).

---

## Task 1: record `predicted_working_set_bytes` in the variants claim

**Files:**
- Modify: `src/main.rs` (`run_variants_index` receipt block, after `max_working_set_bytes`)
- Test: `tests/plan_enforce.rs` (append)

- [ ] **Step 1: Write the failing integration test**

Append to `tests/plan_enforce.rs` (reuses the in-file `bin()`, `unique_dir`, `run`, and `build_big_contig_fixture` helpers):

```rust
#[test]
fn variants_enforce_records_predicted_working_set_in_the_claim() {
    use rosalind::provenance::RunManifest;

    let (dir, idx, bam) = build_big_contig_fixture(1024 * 1024);

    // Run the enforced call twice → two receipts. The deterministic working-set
    // prediction must be identical (that determinism is what makes it a claim).
    let run_once = |vcf: &std::path::Path| {
        let out = Command::new(bin())
            .args(["variants", "--index"])
            .arg(&idx)
            .arg("--alignments")
            .arg(&bam)
            .args(["--max-depth", "1000", "--max-read-len", "250"])
            .args(["--memory-budget-mb", "512", "--enforce", "-o"])
            .arg(vcf)
            .output()
            .unwrap();
        assert!(out.status.success(), "enforced run failed: {out:?}");
        let text = std::fs::read_to_string(format!("{}.manifest.json", vcf.display())).unwrap();
        RunManifest::from_canonical_json(&text).unwrap()
    };

    let a = run_once(&dir.join("a.vcf"));
    let b = run_once(&dir.join("b.vcf"));

    // It is a CLAIM field (in params), not a measurement.
    let pred_a = a
        .params
        .get("predicted_working_set_bytes")
        .expect("predicted_working_set_bytes must be a claim param");
    assert!(
        !a.measurements.contains_key("predicted_working_set_bytes"),
        "predicted_working_set_bytes must live in the claim, not the measurement block"
    );
    // Deterministic across runs (the cross-machine claim property).
    assert_eq!(
        pred_a,
        b.params.get("predicted_working_set_bytes").unwrap(),
        "the predicted working set must be identical across runs"
    );

    // The receipt verifies (an enforced run satisfies predicted >= realized).
    let v = Command::new(bin())
        .args(["verify", "--manifest"])
        .arg(dir.join("a.vcf.manifest.json"))
        .output()
        .unwrap();
    assert!(
        v.status.success(),
        "verify must pass on the enforced receipt: {}",
        String::from_utf8_lossy(&v.stderr)
    );

    std::fs::remove_dir_all(&dir).ok();
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cargo test --test plan_enforce variants_enforce_records_predicted_working_set_in_the_claim`
Expected: FAIL — `predicted_working_set_bytes must be a claim param` panics (the field is not recorded yet).

- [ ] **Step 3: Record the claim param**

In `src/main.rs`, in `run_variants_index`'s receipt block, find the `max_working_set_bytes` insert and add the new param immediately after it:

```rust
        manifest.params.insert(
            "max_working_set_bytes".to_string(),
            max_ws.bytes.to_string(),
        );
        // The deterministic working-set PREDICTION (index header + declared caps, no
        // baseline) — a CLAIM field (not in MEASUREMENT_KEYS), so it is cross-machine
        // stable and `verify` can re-check `predicted >= realized` offline.
        manifest.params.insert(
            "predicted_working_set_bytes".to_string(),
            rosalind::call::plan::estimate_variants_working_set(largest, max_depth, max_read_len)
                .bytes
                .to_string(),
        );
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `cargo test --test plan_enforce variants_enforce_records_predicted_working_set_in_the_claim`
Expected: PASS — the field is in `params` (not `measurements`), identical across runs, and `verify` exits 0.

- [ ] **Step 5: Commit**

```bash
git add src/main.rs tests/plan_enforce.rs
git commit -m "feat(contract): record predicted_working_set_bytes in the variants claim

The deterministic working-set prediction (estimate_variants_working_set: index
header + declared caps, no baseline) now lands in the content-addressed claim
(not MEASUREMENT_KEYS), making the contract's prediction cross-machine stable
and offline-auditable. Additive; schema unchanged.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: the offline soundness re-check in `verify_receipt`

**Files:**
- Modify: `crates/receipt/src/lib.rs` (`verify_receipt`, after the `max_working_set ≤ peak_rss` check; unit tests in `#[cfg(test)] mod tests`)

- [ ] **Step 1: Write the failing unit tests**

In `crates/receipt/src/lib.rs`, inside `#[cfg(test)] mod tests` (next to `verify_receipt_passes_a_clean_receipt`), add a helper and four tests:

```rust
    /// Build + finalize a receipt exercising the working-set soundness check.
    fn soundness_receipt(predicted: Option<&str>, realized: Option<&str>, enforced: bool) -> String {
        let mut m = RunManifest::new("variants");
        if let Some(p) = predicted {
            m.params
                .insert("predicted_working_set_bytes".to_string(), p.to_string());
        }
        if let Some(r) = realized {
            m.params
                .insert("max_working_set_bytes".to_string(), r.to_string());
        }
        if enforced {
            m.params.insert("enforce".to_string(), "true".to_string());
        }
        m.finalize();
        m.to_canonical_json()
    }

    #[test]
    fn soundness_ok_when_enforced_and_prediction_bounds_realized() {
        let report = verify_receipt(
            &soundness_receipt(Some("2000000"), Some("1000000"), true),
            &VerifyOpts::default(),
        );
        assert!(report.ok, "{:?}", report.problems);
        assert!(
            report.notes.iter().any(|n| n.contains(">= realized")),
            "{:?}",
            report.notes
        );
    }

    #[test]
    fn soundness_fails_when_enforced_and_prediction_underbounds_realized() {
        let report = verify_receipt(
            &soundness_receipt(Some("500000"), Some("1000000"), true),
            &VerifyOpts::default(),
        );
        assert!(!report.ok, "an enforced under-prediction must fail verify");
        assert!(
            report.problems.iter().any(|p| p.contains("unsound prediction")),
            "{:?}",
            report.problems
        );
    }

    #[test]
    fn soundness_is_advisory_when_not_enforced() {
        let report = verify_receipt(
            &soundness_receipt(Some("500000"), Some("1000000"), false),
            &VerifyOpts::default(),
        );
        assert!(
            report.ok,
            "a record-only under-prediction must NOT fail: {:?}",
            report.problems
        );
        assert!(
            report.notes.iter().any(|n| n.contains("advisory")),
            "{:?}",
            report.notes
        );
    }

    #[test]
    fn soundness_skipped_when_no_prediction_recorded() {
        let report = verify_receipt(
            &soundness_receipt(None, Some("1000000"), true),
            &VerifyOpts::default(),
        );
        assert!(report.ok, "a pre-change receipt must verify: {:?}", report.problems);
        assert!(
            report
                .notes
                .iter()
                .any(|n| n.contains("no predicted_working_set_bytes")),
            "{:?}",
            report.notes
        );
    }
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `cargo test -p rosalind-receipt soundness_`
Expected: FAIL — `soundness_fails_when_enforced_and_prediction_underbounds_realized` passes verify (no check yet) so its `assert!(!report.ok)` fails; the note-based assertions also fail (notes absent). The check is not implemented.

- [ ] **Step 3: Implement the check**

In `crates/receipt/src/lib.rs::verify_receipt`, immediately after the `max_working_set ≤ peak_rss` consistency block (the `if let (Some(ws), Some(peak)) = (recorded_ws, recorded_peak) { ... }`), insert:

```rust
    // Offline soundness re-check: the deterministic working-set PREDICTION must upper-bound the
    // realized working set. Both are machine-independent (shared cost constants), so it is a true
    // cross-machine claim — HARD under --enforce (the prediction is a guaranteed upper bound there:
    // max_read_len is enforced at ingest, depth is always capped), advisory otherwise (a longer read
    // on a record-only run can legitimately exceed the assumed max_read_len).
    let predicted_ws = parse_num("predicted_working_set_bytes", &mut problems);
    let enforced = manifest.params.get("enforce").map(String::as_str) == Some("true");
    match (predicted_ws, recorded_ws) {
        (Some(pred), Some(real)) if pred >= real => notes.push(format!(
            "predicted working set {} MiB >= realized {} MiB",
            pred / (1 << 20),
            real / (1 << 20)
        )),
        (Some(pred), Some(real)) if enforced => problems.push(format!(
            "unsound prediction: predicted working set {pred} bytes < realized {real} bytes (enforced run)"
        )),
        (Some(pred), Some(real)) => notes.push(format!(
            "predicted working set {} MiB < realized {} MiB (advisory; run was not enforced)",
            pred / (1 << 20),
            real / (1 << 20)
        )),
        _ => notes.push(
            "no predicted_working_set_bytes to re-check (a pre-soundness-check receipt)".to_string(),
        ),
    }
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `cargo test -p rosalind-receipt soundness_`
Expected: PASS — all four `soundness_*` tests pass (OK+note / hard-fail+problem / advisory+note / skip+note).

- [ ] **Step 5: Commit**

```bash
git add crates/receipt/src/lib.rs
git commit -m "feat(verify): offline re-check that predicted >= realized working set

verify_receipt now re-derives the contract's core inequality from the recorded
claim (predicted_working_set_bytes) + measurement (max_working_set_bytes): a
hard FAIL (exit 5) on an enforced under-prediction, an advisory note on a
record-only run, and a clean skip for pre-change receipts. Tamper-protected from
both sides (claim self-hash + measurement self-hash).

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 3: full gate + end-to-end smoke

**Files:** none (verification only)

- [ ] **Step 1: Full suite + lint + format**

Run:
```bash
cargo test
cargo clippy --all-targets -- -D warnings
cargo fmt --check
```
Expected: PASS on all three. (If `cargo fmt --check` reports diffs, run `cargo fmt` and amend the last commit.)

- [ ] **Step 2: End-to-end smoke — verify surfaces the soundness line**

Run (a real enforced call, then read the receipt + run verify):
```bash
BIN=target/debug/rosalind; D=$(mktemp -d)
printf '>chr1\n' > "$D/ref.fa"; head -c 1000000 /dev/zero | tr '\0' 'A' >> "$D/ref.fa"; echo >> "$D/ref.fa"
printf '@r0\nAAAAAAAAAAAAAAAAAAAA\n+\nIIIIIIIIIIIIIIIIIIII\n' > "$D/reads.fq"
$BIN index --reference "$D/ref.fa" --output "$D/ref.idx" >/dev/null 2>&1
$BIN align --reference "$D/ref.fa" --reads "$D/reads.fq" --format bam --output "$D/raw.bam" >/dev/null 2>&1
$BIN sort --input "$D/raw.bam" --output "$D/s.bam" >/dev/null 2>&1
$BIN variants --index "$D/ref.idx" --alignments "$D/s.bam" --memory-budget-mb 512 --enforce -o "$D/c.vcf" >/dev/null 2>&1
grep -o '"predicted_working_set_bytes":"[0-9]*"' "$D/c.vcf.manifest.json"
$BIN verify --manifest "$D/c.vcf.manifest.json"
rm -rf "$D"
```
Expected: the `grep` prints a `predicted_working_set_bytes` line (the field is in the claim), and `verify` prints a `predicted working set … >= realized …` line and exits 0.

- [ ] **Step 3: No commit** (verification only; nothing changed).

---

## Self-Review

**1. Spec coverage:**
- §4 record `predicted_working_set_bytes` in the claim (not in `MEASUREMENT_KEYS`, deterministic, no baseline) → **Task 1**.
- §5 the four-arm `verify_receipt` check (OK note / hard-if-enforced / advisory / skip), reusing `recorded_ws`/`parse_num`, `enforce` from claim params, exit 5 via `problems` → **Task 2**.
- §6 soundness argument → encoded as the enforced-vs-advisory split (Task 2 Step 3 + the two enforced/non-enforced unit tests).
- §7 tests: unit OK/hard/advisory/skip → Task 2 Step 1; integration claim-location + determinism + verify-passes → Task 1 Step 1; gates → Task 3.
- §8 non-goals (no peak-RSS hard check, no schema bump, don't move `max_working_set_bytes`, no `features`) → respected: only an additive claim param + a working-set-level check; nothing touches `predicted_peak_rss_bytes`, `MEASUREMENT_KEYS`, schema, or `features`.

**2. Placeholder scan:** No TBD/TODO/"handle edge cases". Every step has literal code or an exact command + expected result.

**3. Type consistency:** `predicted_working_set_bytes` (claim param, string) is written in Task 1 and read in Task 2's check + tests + Task 1's integration test — same key everywhere. `estimate_variants_working_set(largest, max_depth, max_read_len).bytes` matches the real signature (`plan.rs:23`). `parse_num`, `recorded_ws`, `notes`, `problems`, `manifest.params` match the in-scope names in `verify_receipt`. The unit-test helper `soundness_receipt` builds via the documented `RunManifest` API (`new`/`params.insert`/`finalize`/`to_canonical_json`) used by the adjacent tests. MiB rendering uses `/ (1 << 20)` to match the file's existing style.
