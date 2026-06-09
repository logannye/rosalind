# Design: `predicted_working_set` as a claim + an offline soundness re-check in `verify`

**Status:** Approved design — 2026-06-09. Companion to [`CONTRACT.md`](../../../CONTRACT.md) and
[`docs/ROADMAP.md`](../../ROADMAP.md) (Phase-0 soundness spirit). Audience: the implementer.

---

## 1. The goal, in one sentence

Make the memory contract's core inequality — **the prediction upper-bounds reality** — a
**content-addressed, offline-auditable** property of the receipt, instead of only a runtime
governor side-effect.

## 2. Why this, why now

The headline claim is *"predict your peak before you commit a byte."* Today `verify_receipt`
(`crates/receipt/src/lib.rs`) re-checks `max_working_set ≤ peak_rss` and verdict-vs-budget
consistency, but it **never re-checks the prediction against reality.** Soundness is asserted only
*live* by the runtime governor and *in-process* by `tests/plan_enforce.rs` (which already proves
`predicted ≥ realized`). The gap: a third party holding only the receipt **cannot re-prove** the
prediction was sound — the deterministic half of the headline number isn't even in the
content-addressed claim.

Two facts in the code make this a clean, S-effort change:

1. `estimate_variants_working_set(largest_contig_len, max_depth, max_read_len) -> WorkingSet`
   (`src/call/plan.rs:23`) is **already a pure, deterministic function** of the index header +
   declared caps, sharing cost constants with the realized accountant
   (`PileupEngine::current_working_set`). It has **no baseline term**, so it is machine-independent.
2. `MEASUREMENT_KEYS` (`crates/receipt/src/lib.rs:100`) =
   `{peak_rss_bytes, max_working_set_bytes, predicted_peak_rss_bytes, baseline_rss_bytes,
   rss_residual_bytes, governor, contract_verdict}`. A param **not** in that list stays in the
   **claim** through `finalize()` — hash-protected and cross-machine stable.

So the deterministic working-set prediction can become a claim field with no schema change, and
`verify` can re-derive the inequality offline.

## 3. Scope

**In scope (this PR):**
- Record `predicted_working_set_bytes` as a **claim** param in the `variants --index` receipt.
- Add an offline re-check to `verify_receipt`: the prediction must upper-bound the realized working
  set — **hard fail under `--enforce`, advisory note otherwise.**

**Out of scope (named):**
- `features --index` parity — it shares the estimator and is an identical ~1-line follow-on; kept
  out to keep this PR focused (immediate fast-follow).
- Promoting `max_working_set_bytes` from a measurement into the claim — it is data-deterministic,
  but reclassifying it changes `content_hash()` and the reproduce/chain comparison surface; a
  separate, riskier change.
- Any hard check at the **peak-RSS** level — `predicted_peak_rss_bytes` folds in the measured
  baseline + a fixed allocator-slack margin, so a hard `predicted_peak ≥ realized_peak` check would
  false-fail; it stays a recorded advisory measurement.
- No schema-version bump (the new claim field is additive).

## 4. Part 1 — record `predicted_working_set_bytes` in the claim

In `run_variants_index`'s receipt block (`src/main.rs`, alongside the existing
`predicted_peak_rss_bytes` / `max_working_set_bytes` inserts), add one param. `largest`,
`max_depth`, and `max_read_len` are already in scope (they feed `predicted_peak`):

```rust
manifest.params.insert(
    "predicted_working_set_bytes".to_string(),
    rosalind::call::plan::estimate_variants_working_set(largest, max_depth, max_read_len)
        .bytes
        .to_string(),
);
```

- Deterministic (index header + declared caps), **no baseline** → cross-machine stable.
- **Not** in `MEASUREMENT_KEYS`, so `finalize()` keeps it in the claim → covered by `manifest_blake3`,
  part of `content_hash()`. The baseline-folded `predicted_peak_rss_bytes` stays a measurement.

## 5. Part 2 — the offline check in `verify_receipt`

In `crates/receipt/src/lib.rs::verify_receipt`, after the existing `max_working_set ≤ peak_rss`
consistency check, add (reusing the already-parsed `recorded_ws` = `max_working_set_bytes` and the
strict `parse_num` closure that flags a present-but-unparseable field):

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
        pred >> 20,
        real >> 20
    )),
    (Some(pred), Some(real)) if enforced => problems.push(format!(
        "unsound prediction: predicted working set {pred} bytes < realized {real} bytes (enforced run)"
    )),
    (Some(pred), Some(real)) => notes.push(format!(
        "predicted working set {} MiB < realized {} MiB (advisory; run was not enforced)",
        pred >> 20,
        real >> 20
    )),
    _ => notes.push(
        "no predicted_working_set_bytes to re-check (a pre-soundness-check receipt)".to_string(),
    ),
}
```

- **Tamper-protected both ways:** the prediction lives in the claim (`manifest_blake3` catches
  edits); the realized `max_working_set_bytes` lives in the measurement block (`measurement_blake3`
  catches edits). A forger cannot lower either to fake a pass.
- **`enforce`** is a claim param (recorded only when `--enforce`, via
  `cmd.flag_if(enforce, "--enforce")`), so the hard/advisory split is itself hash-protected.
- **Exit code:** a violation appends to `problems`, so `verify` exits **5** through its existing
  path. No new exit code.
- **Back-compat:** a pre-change receipt has no `predicted_working_set_bytes` → `parse_num` returns
  `None` → the `_` arm records a skip note, never a false failure (the pre-1.2 self-hash-skip
  pattern).

## 6. Soundness argument (why hard-under-enforce is correct)

Under `--enforce`: `max_read_len` is checked at ingest (a longer read aborts the run, exit 4), depth
is always capped at `max_depth` by the pileup, and the contig term is bounded by the largest contig.
The estimator and the realized accountant use the **same** per-base/per-read constants
(`PILEUP_MAP_BYTES_PER_BASE`, `PILEUP_SEQQUAL_BYTES_PER_BASE`, `PILEUP_PER_READ_OVERHEAD`,
`PILEUP_ENGINE_OVERHEAD`). Therefore `estimate_variants_working_set(...) ≥ realized working set` is a
**true invariant** for an enforced run, and a violation is a genuine estimator/accountant bug worth
failing on. Without `--enforce`, `max_read_len` is only an *assumption* (the run was record-only), so
the same comparison is reported as advisory — failing it would punish a correct record-only run.

## 7. Testing

**Unit (`crates/receipt/src/lib.rs`, `verify_receipt`)** — build the manifest, then `finalize()` it
so the claim self-hash is valid, then verify. *(The violation cases must be tested here, not via an
integration receipt: editing a finalized receipt to violate the inequality would trip the self-hash
first.)*
- enforced + `predicted ≥ realized` → `ok`, the `>= realized` note present.
- enforced + `predicted < realized` → **`ok == false`**, an `unsound prediction` problem present.
- non-enforced + `predicted < realized` → `ok == true`, an `advisory` note (no problem).
- no `predicted_working_set_bytes` → `ok == true`, the skip note (no false problem).

**Integration (`tests/plan_enforce.rs`)** — through the real CLI:
- `variants --index … --enforce -o calls.vcf` writes a receipt whose **claim** contains
  `predicted_working_set_bytes`, and `verify` exits 0.
- The recorded `predicted_working_set_bytes` is **byte-identical across two independent runs** (its
  determinism is exactly what makes it a cross-machine claim).

**Gates:** `cargo test`, `cargo clippy --all-targets -- -D warnings`, `cargo fmt --check` all clean.

## 8. Non-goals / guardrails (do not drift)

- No hard check at the peak-RSS level (baseline-folded → false-fail landmine).
- No schema bump; the new claim field is additive and pre-change receipts skip the check cleanly.
- Do not move `max_working_set_bytes` into the claim in this PR.
- `features` parity is a separate follow-on, not this PR.
