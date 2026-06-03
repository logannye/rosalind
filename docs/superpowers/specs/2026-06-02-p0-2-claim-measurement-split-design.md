# P0.2 — Receipt claim/measurement split (cross-machine-stable receipts)

**Status:** design
**Date:** 2026-06-02
**Roadmap:** `docs/ROADMAP.md` §6 Phase 0, item P0.2
**Predecessor:** Sprint 1.2 (tamper-evident self-hashing receipt), P0.1 (sound build estimate)

## The defect

`RunManifest::content_hash()` (`src/provenance/mod.rs`) clones the manifest, removes only
`manifest_blake3` from `params`, and hashes the **entire** canonical JSON — which still
includes the machine-dependent measured fields the run records in `params`:

- `peak_rss_bytes` (realized resident high-water on *this* machine)
- `max_working_set_bytes` (realized pileup working set)
- `predicted_peak_rss_bytes` (a prediction *anchored to this machine's measured baseline RSS*)
- `baseline_rss_bytes` (measured process baseline at start)
- `rss_residual_bytes` (realized peak − predicted)
- `governor` (runtime governor status)
- `contract_verdict` (derived from the measured peak vs the declared budget)

Consequence: **the same correct run on two machines produces two different receipt hashes.**
The receipt was supposed to be the reproducibility anchor — "running this produces an output
with this exact hash" — but the anchor moves whenever the cost of the run moves. That makes
`manifest_blake3` unusable as a cross-machine identity, which in turn blocks every downstream
proof that needs a stable receipt identity: chaining (P3), `reproduce` (P3), and cohort
receipts (P3). The hash must commit to *what was computed*, not *what it cost here*.

## The lossless-split principle

Partition the receipt into two parts with different reproducibility semantics:

- **Claim** — deterministic, portable, identical for the same logical run on any machine:
  `inputs`, `outputs` (content hashes), `params` (declared knobs: budget, thresholds,
  deterministic output counts), `subcommand`, `tool_version`, `schema_version`.
- **Measurement** — machine-/run-dependent observed cost: the seven fields above.

`content_hash()` (the claim hash, recorded as `manifest_blake3`) is computed over the **claim
only**, excluding the measurement block. Two machines that run the same logical computation
now produce the **same** `manifest_blake3`.

### Why the split must be lossless (the subtle part)

Before this change, the single self-hash covered `peak_rss_bytes` (it lived in the claimed
`params`), so the receipt was tamper-evident *against an edit to the measured peak* — e.g.
someone lowering `peak_rss_bytes` to make a run look like it fit a budget. If we simply drop
the measurement from the hash, **we lose that tamper-evidence** — a strict regression of the
Sprint 1.2 guarantee.

So the split adds a **second hash**: `measurement_blake3`, a BLAKE3 over the canonical
measurement block (excluding `measurement_blake3` itself). It is **not** cross-machine stable
(it hashes machine-dependent numbers) and therefore lives **inside** the measurement block
(putting it in the claim would re-break claim stability). Result — two-layer integrity:

| Layer | Hash | Covers | Reproducible across machines? | Catches |
|---|---|---|---|---|
| Claim | `manifest_blake3` (in `params`) | inputs, outputs, params, subcommand, versions | **Yes** | edits to declared inputs/outputs/params |
| Measurement | `measurement_blake3` (in `measurements`) | the 7 measured fields | No (local attestation) | edits to the measured cost on this machine |

The measurement is an **on-machine attestation**: tamper-evident locally, cross-checked for
internal consistency by `verify`, but inherently non-portable. Cryptographic *signing* of the
measurement by the running machine is out of scope here (a later attestation phase). This
spec's job is the split + the portable claim hash, **without regressing** the measured-field
tamper-evidence Sprint 1.2 shipped.

## Design

### Data model (`src/provenance/mod.rs`)

Add a field to `RunManifest`:

```rust
pub measurements: BTreeMap<String, String>,
```

Initialized empty by `RunManifest::new`. The measurement keys are a single audited policy
list — the source of truth for the partition:

```rust
pub const MEASUREMENT_KEYS: &[&str] = &[
    "peak_rss_bytes",
    "max_working_set_bytes",
    "predicted_peak_rss_bytes",
    "baseline_rss_bytes",
    "rss_residual_bytes",
    "governor",
    "contract_verdict",
];
```

Bump `MANIFEST_SCHEMA_VERSION` `1 → 2`.

### Centralized partition (Option B — relocation in `finalize`)

Callers keep inserting fields wherever is natural today; `finalize()` enforces the partition
centrally so no call site can leave a measured field in the claim:

```rust
pub fn finalize(&mut self) {
    // 1. Partition: relocate machine-dependent fields OUT of the claim.
    for key in MEASUREMENT_KEYS {
        if let Some(v) = self.params.remove(*key) {
            self.measurements.insert((*key).to_string(), v);
        }
    }
    // 2. Local measurement attestation (only when there is a measurement).
    if !self.measurements.is_empty() {
        let mh = self.measurement_hash();
        self.measurements.insert("measurement_blake3".to_string(), mh);
    }
    // 3. Stamp the schema version into the claim, then the claim self-hash LAST.
    self.params.insert("schema_version".to_string(), MANIFEST_SCHEMA_VERSION.to_string());
    let h = self.content_hash();
    self.params.insert("manifest_blake3".to_string(), h);
}
```

This keeps the existing manifest-building sites in `main.rs` unchanged (lower diff risk); the
policy lives in one list + one method. A `record_measurement(key, val)` helper is also exposed
for call sites/tests that prefer to be explicit.

### Serialization

`measurements` is a `string → string` object, serialized identically to `params`. Canonical
key order stays alphabetical, so `measurements` sits between `inputs` and `outputs`:

```
{"inputs":[…],"measurements":{…},"outputs":[…],"params":{…},"subcommand":"…","tool_version":"…"}
```

The `measurements` key is **emitted only when non-empty**. This keeps:
- pre-v2 receipts (no measurements) and v2 receipts with no measurements byte-identical, and
- the existing exact-JSON unit test (`canonical_json_has_sorted_keys_and_is_exact`, empty
  measurements) passing unchanged.

Refactor serialization into one helper parameterized by whether to include measurements:

- `to_canonical_json()` — full receipt (claim + measurements). Used for `write_manifest`.
- `to_canonical_claim_json()` — claim only (never emits measurements). Used by `content_hash`.

`content_hash()` removes `manifest_blake3` from a clone's `params`, then hashes
`to_canonical_claim_json()`. `measurement_hash()` hashes the canonical measurement object with
`measurement_blake3` removed.

### Hashing methods

```rust
pub fn content_hash(&self) -> String;        // claim hash (manifest_blake3) — cross-machine stable
pub fn measurement_hash(&self) -> String;    // measurement block hash (measurement_blake3) — local
pub fn self_hash_ok(&self) -> Option<bool>;          // claim: Some(match) / None if absent
pub fn measurement_hash_ok(&self) -> Option<bool>;   // measurement: Some(match) / None if absent
```

### Parser (back-compat)

`parse_manifest` accepts an **optional** `measurements` key after `inputs`: parse the next key
string and branch — if `"measurements"`, consume the object then expect `outputs`; if
`"outputs"`, treat measurements as empty. This parses both v1 (no key) and v2 receipts.

### Backward compatibility (graceful pre-v2 degradation)

A pre-v2 receipt put `peak_rss_bytes` in `params` and had no `measurements` key. On parse,
`measurements` is empty. Because `to_canonical_claim_json()` of a parsed v1 receipt (empty
measurements, peak still in params) is **byte-identical** to the old full canonical form,
`content_hash()` reproduces the old v1 hash exactly — so a v1 receipt's `self_hash_ok()`
still returns `Some(true)`. `measurement_hash_ok()` returns `None` (no measurement block);
`verify` notes and skips it. No version branching is needed in `content_hash` — the partition
is correct by construction for both shapes.

### `verify` (`run_verify` in `src/main.rs`)

- Read measured/declared fields via a new `manifest.get_recorded(key)` that checks
  `measurements` then `params` — so `verify` reads both v2 receipts (measured fields in
  `measurements`) and pre-v2 receipts (everything in `params`) with one call.
- Add a `measurement_hash_ok()` check beside `self_hash_ok()`: `Some(false)` → a
  "measurement_blake3 mismatch" problem; `None` → noted and skipped (pre-v2).
- Keep the existing internal-consistency cross-checks (ws ≤ peak; verdict vs peak vs budget) —
  they remain the semantic guard for the measurement and now read via `get_recorded`.
- Update the stale "the manifest has no self-hash yet" comment to the two-layer model.

## What does NOT move

Deterministic counts stay in the claim (same input → same value on any machine):
`feature_rows`, `over_max_depth`, `reads_skipped_total`, `io_rss_overhead_assumed_bytes` (a
compile-time constant), and all declared knobs (`min_qual`, `max_depth`, `max_read_len`,
`mapq_threshold`, `enforced`, `memory_budget_mb`, `region_start`, `somatic_snv_only`).

## Testing

New `provenance` unit tests:
1. **Cross-machine claim stability** — two manifests, identical claim, different measured
   values → identical `content_hash()`; both `self_hash_ok() == Some(true)`.
2. **Claim excludes / full includes** — `to_canonical_claim_json()` omits `peak_rss_bytes`;
   `to_canonical_json()` contains it.
3. **Measurement tamper-evidence** — finalize, edit a value in `measurements` →
   `measurement_hash_ok() == Some(false)` while `self_hash_ok()` stays `Some(true)`.
4. **v1 back-compat** — a v1-shaped receipt (peak in params, no measurements, schema 1) still
   `self_hash_ok() == Some(true)`; round-trips with no `measurements` key; parses to empty
   measurements; `measurement_hash_ok() == None`.
5. **v2 round-trip** — a finalized v2 receipt round-trips through `from_canonical_json` (the
   optional-measurements parser), preserving both maps and both hashes.

Integration-test updates (`tests/`): point measured-field reads at `.measurements`
(`plan_enforce.rs`: `max_working_set_bytes`, `predicted_peak_rss_bytes`, `peak_rss_bytes`,
`rss_residual_bytes`, `governor`; `gvcf.rs`: `peak_rss_bytes`); bump the `schema_version`
assertion `"1"→"2"` and assert `measurement_blake3` is present; **add** an end-to-end test that
tampers a measured field while keeping the receipt internally consistent and asserts
`verify` fails with `measurement_blake3 mismatch` (the new lossless property, mirroring the
existing claim-field tamper test). Substring assertions over raw JSON are unaffected (the
substrings still appear, now inside `measurements`).

## Out of scope

Cryptographic signing/attestation of the measurement by the running machine; build-identity
fields (`code_git_sha`, `rustc_version`, lockfile hash) — that is P0.3.

## Risks

- **Missing a measured key** → it stays in the claim and re-breaks cross-machine stability.
  Mitigated by `MEASUREMENT_KEYS` being the single audited list + a unit test asserting no
  measured key affects the claim hash.
- **Parser regression on the optional key** → covered by the v1 and v2 round-trip tests.
