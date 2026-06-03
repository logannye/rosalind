# Sprint 1.2 — The Unforgeable Receipt (design)

**Status:** Approved design — 2026-06-02. Increment 1.2 of the engineering roadmap
([`docs/ROADMAP.md`](../../ROADMAP.md), Sprint 1). Builds on Sprint 1.1 (the runtime governor).

---

## 1. Problem

The run receipt is canonical JSON with BLAKE3 hashes of the inputs/outputs and the recorded
memory numbers, but **the manifest itself is unauthenticated**. `run_verify` already cross-checks
internal consistency (its own comment, `main.rs:926`, says *"The manifest has no self-hash yet"*),
but a hand-edited receipt whose fields are mutually consistent passes. Two gaps:

1. **No self-hash** — corruption or a casual edit of any field (peak, budget, verdict, a content
   hash) is not detected as long as the edited fields stay internally consistent.
2. **No schema version** — there is no version token on the receipt/feature schema, so a future
   columnar (Parquet) egress or any consumer cannot tell which schema produced a given artifact.
   Stamping one is cheap now and expensive to retrofit after a schema freezes.

This increment closes both, entirely on the existing canonical-JSON machinery.

## 2. Goals / non-goals

**Goals**
- A **self-hash** (`manifest_blake3`) over the canonical JSON, re-derived and checked by `verify` —
  upgrading the receipt from *corruption-proof* to **tamper-evident**.
- A **`schema_version`** token stamped into every receipt.
- Backward compatibility: a pre-1.2 receipt (no `manifest_blake3`) still verifies, with a note.

**Non-goals (explicitly deferred)**
- Cryptographic signatures / tamper-*proof* receipts (the "verify-attest" follow-up; a self-hash is
  tamper-*evident*, not tamper-proof — a determined forger can recompute it).
- Removing the `verify` budget-from-manifest fallback (decided: **keep it** — the self-hash makes the
  recorded budget tamper-evident; an auditor wanting tamper-proof supplies `--budget-mb` externally).
- Changing the `FEATURE_HEADER` TSV bytes (the schema version lives in the manifest, not the table —
  see §3.2).
- Arrow/Parquet egress (Sprint 3) — `schema_version` is laid down now so it is *ready* for it.

## 3. Design

### 3.1 Self-hash (`provenance::RunManifest`)

The self-hash is a `params` entry `manifest_blake3` — no change to the canonical-JSON top-level shape
or the hand-parser. It is computed over the canonical JSON **with that field excluded**, so `verify`
can re-derive it deterministically.

```rust
/// Current receipt/feature schema version. Bump on any breaking schema change.
pub const MANIFEST_SCHEMA_VERSION: u32 = 1;

impl RunManifest {
    /// BLAKE3 hex of the canonical JSON with the self-hash field excluded — the
    /// content this manifest commits to. Deterministic; `verify` re-derives it.
    pub fn content_hash(&self) -> String {
        let mut m = self.clone();
        m.params.remove("manifest_blake3");
        blake3_hex(m.to_canonical_json().as_bytes())
    }

    /// Stamp the schema version + the self-hash. Call LAST, immediately before
    /// serialization, so the hash covers every other field (including the version).
    pub fn finalize(&mut self) {
        self.params
            .insert("schema_version".to_string(), MANIFEST_SCHEMA_VERSION.to_string());
        let h = self.content_hash();
        self.params.insert("manifest_blake3".to_string(), h);
    }

    /// `Some(true)`/`Some(false)` if a self-hash is recorded and matches / mismatches;
    /// `None` if none is recorded (a pre-1.2 receipt).
    pub fn self_hash_ok(&self) -> Option<bool> {
        self.params
            .get("manifest_blake3")
            .map(|recorded| *recorded == self.content_hash())
    }
}
```

`RunManifest` already derives `Clone`, so `content_hash` cloning is free of new bounds. The order in
`finalize` matters: `schema_version` is inserted first so the self-hash commits to it; `content_hash`
removes `manifest_blake3` defensively (so `finalize` is idempotent and `self_hash_ok` is symmetric).

### 3.2 `schema_version` placement — manifest only

`schema_version` is recorded **only in the manifest `params`**, never in the `FEATURE_HEADER` TSV
bytes. Rationale: the feature table's byte-for-byte reproducibility is a shipped guarantee (the golden
feature test pins the exact header); injecting a version line would break it and is itself a breaking
feature-schema change. The provenance record is the correct, non-breaking home for a schema version —
it is exactly where a consumer (or the future Parquet egress) reads "which schema is this".

### 3.3 Wiring at the four receipt-write sites

Call `manifest.finalize()` once, immediately before serialization, at each site:

| Site | `main.rs` (approx) | Serialization |
|---|---|---|
| `run_somatic` | `:1148`/`:1168` | `write_manifest(&out, &manifest)` |
| `run_variants` (`--reference`) | `:1515`/`:1538` | `write_manifest(&path, &manifest)` |
| `run_features` | `:1767`/`:1838` | `std::fs::write(dest, manifest.to_canonical_json())` |
| `run_variants_index` | `:2146`/`:2218` | `std::fs::write(dest, manifest.to_canonical_json())` |

`finalize` takes `&mut self`, so each `let mut manifest = …` binding (already `mut`) calls
`manifest.finalize();` just before the write. `write_manifest` and `to_canonical_json` are unchanged.

### 3.4 `verify` (`run_verify`)

After parsing, add a self-hash check to the `problems` accumulator, alongside the existing checks:

```rust
match manifest.self_hash_ok() {
    Some(true) => {}
    Some(false) => problems.push(
        "manifest_blake3 mismatch: the receipt was modified after it was written".into(),
    ),
    None => println!("verify: note — no manifest_blake3 (a pre-1.2 receipt); skipping self-hash"),
}
```

The budget-from-manifest fallback (`main.rs:895`) stays. The `verify: OK` / exit-5 logic is unchanged
otherwise — a self-hash mismatch simply joins the existing `problems`, so it fails the same way.

### 3.5 Determinism & compatibility

- **Determinism preserved.** `schema_version` is constant; `manifest_blake3` is a deterministic hash
  of deterministic content. On one machine, identical inputs → byte-identical manifest *including* the
  self-hash. (Across machines the manifest already differs by `peak_rss_bytes`; the self-hash differs
  too — it authenticates *this* receipt, not cross-machine equality, exactly as intended.)
- **Backward compatible.** A v0.1.0 receipt has no `manifest_blake3` → `self_hash_ok()` is `None` →
  `verify` notes it and proceeds. No existing receipt breaks.
- **Existing tests.** `to_canonical_json` is unchanged, so the provenance unit tests that build a
  manifest by hand (no `finalize`) are unaffected. Integration tests assert via `.contains(…)`, so the
  two new params don't break them.

## 4. File-by-file change list

| File | Change |
|---|---|
| `src/provenance/mod.rs` | `MANIFEST_SCHEMA_VERSION` const; `content_hash`, `finalize`, `self_hash_ok` on `RunManifest`; unit tests. |
| `src/main.rs` | `manifest.finalize();` before each of the 4 receipt writes; the self-hash check in `run_verify`. |
| `tests/plan_enforce.rs` | Integration: a receipt carries `manifest_blake3` + `schema_version`; `verify` passes untampered; `verify` FAILS on a tampered field whose file hashes still match. |
| `CONTRACT.md` | The "Verify" section: note the receipt is now self-hashing (tamper-evident). |

## 5. Testing

1. **Unit (`provenance`)** — `content_hash` is stable across re-serialization; `finalize` makes
   `self_hash_ok() == Some(true)`; flipping any param after `finalize` makes it `Some(false)`; a
   manifest without the field returns `None`; `finalize` is idempotent (second call re-stamps to the
   same value).
2. **Integration (`tests/plan_enforce.rs`)** — a real `variants --index` receipt contains
   `manifest_blake3` + `schema_version=1`; `verify` passes; then hand-edit a *consistent* field (e.g.
   bump `memory_budget_mb` and `peak_rss` together so the consistency checks still pass) and assert
   `verify` now **exits 5** on the self-hash mismatch — the gap the consistency checks alone cannot
   close.
3. **Regression** — full suite green; the existing `verify_passes_on_an_untampered_run…` and
   `verify_rejects_an_internally_inconsistent_manifest` tests still pass (the latter now ALSO trips the
   self-hash, but still exits 5 — assert it stays 5).

## 6. References

- `src/provenance/mod.rs` — `RunManifest`, `to_canonical_json`, `from_canonical_json`, `blake3_hex`.
- `src/main.rs:870` (`run_verify`), `:1148`/`:1515`/`:1767`/`:2146` (the four receipt-write sites).
- `docs/ROADMAP.md` Sprint 1 (QW-3) — this increment.
