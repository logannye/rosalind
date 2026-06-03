# P0.2 — Receipt claim/measurement split — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the reproducibility receipt's self-hash cross-machine stable by hashing only the
deterministic *claim*, while preserving tamper-evidence for the machine-dependent *measurement*
via a second, local hash.

**Architecture:** Add a `measurements: BTreeMap` block to `RunManifest`. `finalize()`
relocates the seven machine-dependent keys (`MEASUREMENT_KEYS`) out of `params` into
`measurements`, stamps a measurement hash, then stamps the claim self-hash over the
claim-only canonical form. The parser optionally consumes the new key; `verify` reads both
maps for back-compat. Bump schema version 1 → 2.

**Tech Stack:** Rust, BLAKE3 (`blake3` crate), hand-written canonical-JSON serializer/parser.

**Reference:** spec `docs/superpowers/specs/2026-06-02-p0-2-claim-measurement-split-design.md`.

---

### Task 1: Provenance unit tests (RED) — the split's core properties

**Files:**
- Modify: `src/provenance/mod.rs` (append to `mod tests`)

- [ ] **Step 1: Write the failing tests** — append to `mod tests`:

```rust
    #[test]
    fn claim_hash_is_stable_across_machine_dependent_measurements() {
        // Same logical run on two machines: identical claim, different measured cost.
        // The claim self-hash must match; the measurement is excluded from it.
        let mk = |peak: &str, ws: &str| {
            let mut m = RunManifest::new("variants");
            m.tool_version = "0.1.0".to_string();
            m.inputs.push(FileHash {
                path: "ref.fa".to_string(),
                blake3: "aa".to_string(),
            });
            m.outputs.push(FileHash {
                path: "out.vcf".to_string(),
                blake3: "bb".to_string(),
            });
            m.params.insert("min_qual".to_string(), "30".to_string());
            m.record_measurement("peak_rss_bytes", peak);
            m.record_measurement("max_working_set_bytes", ws);
            m.finalize();
            m
        };
        let a = mk("1000000", "4096");
        let b = mk("9999999", "8192");
        assert_eq!(
            a.content_hash(),
            b.content_hash(),
            "measured cost must not change the claim hash"
        );
        assert_eq!(a.self_hash_ok(), Some(true));
        assert_eq!(b.self_hash_ok(), Some(true));
        // Differing measurements DO change the measurement hash.
        assert_ne!(
            a.measurements.get("measurement_blake3"),
            b.measurements.get("measurement_blake3")
        );
    }

    #[test]
    fn claim_excludes_but_full_form_includes_measurements() {
        let mut m = RunManifest::new("variants");
        m.tool_version = "0.1.0".to_string();
        m.record_measurement("peak_rss_bytes", "123");
        m.finalize();
        assert!(
            !m.to_canonical_claim_json().contains("peak_rss_bytes"),
            "claim form must not carry the measurement"
        );
        assert!(
            m.to_canonical_json().contains("peak_rss_bytes"),
            "full form must record the measurement"
        );
    }

    #[test]
    fn editing_a_measurement_breaks_only_the_measurement_hash() {
        let mut m = RunManifest::new("variants");
        m.tool_version = "0.1.0".to_string();
        m.record_measurement("peak_rss_bytes", "123");
        m.finalize();
        assert_eq!(m.self_hash_ok(), Some(true));
        assert_eq!(m.measurement_hash_ok(), Some(true));
        // Lower the recorded peak WITHOUT re-finalizing (a tampered receipt).
        m.measurements
            .insert("peak_rss_bytes".to_string(), "1".to_string());
        assert_eq!(
            m.self_hash_ok(),
            Some(true),
            "claim hash is unaffected by the measurement edit"
        );
        assert_eq!(
            m.measurement_hash_ok(),
            Some(false),
            "measurement hash must catch the edit"
        );
    }

    #[test]
    fn finalize_relocates_measured_keys_out_of_the_claim() {
        let mut m = RunManifest::new("variants");
        m.tool_version = "0.1.0".to_string();
        // Insert measured fields the legacy way (into params); finalize must relocate.
        m.params
            .insert("peak_rss_bytes".to_string(), "555".to_string());
        m.params
            .insert("contract_verdict".to_string(), "within".to_string());
        m.params.insert("min_qual".to_string(), "30".to_string());
        m.finalize();
        for k in ["peak_rss_bytes", "contract_verdict"] {
            assert!(!m.params.contains_key(k), "{k} must leave the claim");
            assert!(m.measurements.contains_key(k), "{k} must enter measurements");
        }
        assert!(m.params.contains_key("min_qual"), "claim params stay put");
    }

    #[test]
    fn pre_v2_receipt_with_measurements_in_params_still_verifies() {
        // A pre-v2 receipt: peak in params, no measurements key, schema 1, self-hash
        // over the params-inclusive claim. It must still self-verify (graceful degrade).
        let mut m = RunManifest::new("variants");
        m.tool_version = "0.1.0".to_string();
        m.params
            .insert("peak_rss_bytes".to_string(), "123".to_string());
        m.params
            .insert("schema_version".to_string(), "1".to_string());
        let h = m.content_hash();
        m.params.insert("manifest_blake3".to_string(), h);

        assert_eq!(m.self_hash_ok(), Some(true));
        assert_eq!(m.measurement_hash_ok(), None, "no measurement block in v1");
        let json = m.to_canonical_json();
        assert!(!json.contains("\"measurements\""), "v1 emits no measurements key");
        let parsed = RunManifest::from_canonical_json(&json).expect("parse v1");
        assert!(parsed.measurements.is_empty());
        assert_eq!(parsed.self_hash_ok(), Some(true));
    }

    #[test]
    fn v2_receipt_round_trips_through_the_parser() {
        let mut m = RunManifest::new("features");
        m.tool_version = "9.9.9".to_string();
        m.inputs.push(FileHash {
            path: "a.idx".to_string(),
            blake3: "aa".to_string(),
        });
        m.outputs.push(FileHash {
            path: "out.tsv".to_string(),
            blake3: "bb".to_string(),
        });
        m.params
            .insert("feature_rows".to_string(), "42".to_string());
        m.record_measurement("peak_rss_bytes", "1000");
        m.record_measurement("governor", "enforced");
        m.finalize();

        let json = m.to_canonical_json();
        let parsed = RunManifest::from_canonical_json(&json).expect("parse v2");
        assert_eq!(parsed.to_canonical_json(), json, "round-trip is the identity");
        assert_eq!(
            parsed.measurements.get("peak_rss_bytes").map(String::as_str),
            Some("1000")
        );
        assert_eq!(parsed.params.get("feature_rows").map(String::as_str), Some("42"));
        assert_eq!(parsed.self_hash_ok(), Some(true));
        assert_eq!(parsed.measurement_hash_ok(), Some(true));
    }
```

- [ ] **Step 2: Run tests to verify they fail to compile**

Run: `cargo test -p rosalind --lib provenance 2>&1 | head -40`
Expected: compile errors — `record_measurement`, `measurements`, `to_canonical_claim_json`,
`measurement_hash_ok` not found.

---

### Task 2: Provenance implementation (GREEN)

**Files:**
- Modify: `src/provenance/mod.rs`

- [ ] **Step 1: Bump the schema version + add the policy list**

Replace the version const:

```rust
/// Current receipt/feature schema version. Bump on any breaking schema change.
/// v2: the receipt is split into a deterministic *claim* and a machine-dependent
/// *measurement* block; the self-hash (`manifest_blake3`) covers the claim only.
pub const MANIFEST_SCHEMA_VERSION: u32 = 2;

/// Keys whose values are machine-/run-dependent measurements, not part of the
/// deterministic claim. `finalize` relocates these out of `params` into the
/// `measurements` block so the claim hash is identical for the same logical run on
/// any machine. The single audited source of truth for the claim/measurement split.
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

- [ ] **Step 2: Add the `measurements` field + init in `new`**

In `struct RunManifest`, after `outputs`:

```rust
    /// Machine-/run-dependent measured cost (peak RSS, working set, verdict, …),
    /// excluded from the claim hash. Carries its own `measurement_blake3`.
    pub measurements: BTreeMap<String, String>,
```

In `RunManifest::new`, add `measurements: BTreeMap::new(),` to the struct literal.

- [ ] **Step 3: Add `record_measurement` + replace serialization with a shared helper**

Add to `impl RunManifest` (near `new`):

```rust
    /// Record a machine-/run-dependent measurement (excluded from the claim hash).
    pub fn record_measurement(&mut self, key: impl Into<String>, value: impl Into<String>) {
        self.measurements.insert(key.into(), value.into());
    }
```

Replace `to_canonical_json` with the claim/full pair backed by one helper:

```rust
    /// Serialize to canonical JSON: keys sorted, arrays sorted by path, no
    /// timestamps. Includes the measurement block (when non-empty). Used on disk.
    pub fn to_canonical_json(&self) -> String {
        self.push_canonical(true)
    }

    /// The claim-only canonical JSON: never emits the measurement block. This is the
    /// portion the self-hash commits to, so the hash is identical for the same
    /// logical run on any machine.
    pub fn to_canonical_claim_json(&self) -> String {
        self.push_canonical(false)
    }

    fn push_canonical(&self, include_measurements: bool) -> String {
        let mut out = String::new();
        out.push('{');
        out.push_str("\"inputs\":");
        push_file_hashes(&mut out, &self.inputs);
        if include_measurements && !self.measurements.is_empty() {
            out.push_str(",\"measurements\":");
            push_string_map(&mut out, &self.measurements);
        }
        out.push_str(",\"outputs\":");
        push_file_hashes(&mut out, &self.outputs);
        out.push_str(",\"params\":");
        push_string_map(&mut out, &self.params);
        out.push_str(",\"subcommand\":\"");
        out.push_str(&json_escape(&self.subcommand));
        out.push_str("\",\"tool_version\":\"");
        out.push_str(&json_escape(&self.tool_version));
        out.push_str("\"}");
        out
    }
```

Add the shared string-map serializer near `push_file_hashes`:

```rust
/// Render a `{"k":"v",…}` object, entries in the map's (sorted) key order.
fn push_string_map(out: &mut String, map: &BTreeMap<String, String>) {
    out.push('{');
    for (i, (k, v)) in map.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        out.push('"');
        out.push_str(&json_escape(k));
        out.push_str("\":\"");
        out.push_str(&json_escape(v));
        out.push('"');
    }
    out.push('}');
}
```

- [ ] **Step 4: Claim + measurement hashing**

Replace `content_hash` and add the measurement methods:

```rust
    /// BLAKE3 hex of the **claim** canonical JSON with the self-hash excluded — the
    /// content this manifest commits to, identical on any machine. `verify` re-derives it.
    pub fn content_hash(&self) -> String {
        let mut m = self.clone();
        m.params.remove("manifest_blake3");
        blake3_hex(m.to_canonical_claim_json().as_bytes())
    }

    /// BLAKE3 hex of the measurement block with `measurement_blake3` excluded — a
    /// LOCAL attestation of the measured cost. Not cross-machine stable by design.
    pub fn measurement_hash(&self) -> String {
        let mut m = self.measurements.clone();
        m.remove("measurement_blake3");
        let mut s = String::new();
        push_string_map(&mut s, &m);
        blake3_hex(s.as_bytes())
    }

    /// `Some(match)` if a measurement self-hash is recorded; `None` if absent
    /// (no measurement block, or a pre-v2 receipt).
    pub fn measurement_hash_ok(&self) -> Option<bool> {
        self.measurements
            .get("measurement_blake3")
            .map(|recorded| *recorded == self.measurement_hash())
    }

    /// Look up a recorded value by key, checking `measurements` then `params`. Lets
    /// `verify` read v2 receipts (measured fields in `measurements`) and pre-v2
    /// receipts (everything in `params`) uniformly.
    pub fn get_recorded(&self, key: &str) -> Option<&String> {
        self.measurements.get(key).or_else(|| self.params.get(key))
    }
```

- [ ] **Step 5: `finalize` — relocate, then stamp both hashes**

Replace `finalize`:

```rust
    /// Seal the receipt: partition measured fields out of the claim, stamp the
    /// measurement hash, the schema version, then the claim self-hash LAST (so it
    /// covers every other claim field, including the version). Idempotent.
    pub fn finalize(&mut self) {
        // 1. Partition: relocate machine-dependent fields OUT of the claim.
        for key in MEASUREMENT_KEYS {
            if let Some(v) = self.params.remove(*key) {
                self.measurements.insert((*key).to_string(), v);
            }
        }
        // 2. Local measurement attestation (only when a measurement exists).
        if !self.measurements.is_empty() {
            let mh = self.measurement_hash();
            self.measurements
                .insert("measurement_blake3".to_string(), mh);
        }
        // 3. Stamp the version into the claim, then the claim self-hash last.
        self.params.insert(
            "schema_version".to_string(),
            MANIFEST_SCHEMA_VERSION.to_string(),
        );
        let h = self.content_hash();
        self.params.insert("manifest_blake3".to_string(), h);
    }
```

(`self_hash_ok` is unchanged — it already compares `manifest_blake3` to `content_hash()`.)

- [ ] **Step 6: Parser — optionally consume `measurements`**

In `parse_manifest`, replace the `inputs → outputs` section. After parsing `inputs` and its
trailing comma, branch on the next key:

```rust
    fn parse_manifest(&mut self) -> Result<RunManifest, ManifestError> {
        self.expect(b'{')?;
        self.expect_key("inputs")?;
        let inputs = self.parse_file_array()?;
        self.expect(b',')?;
        // `measurements` is optional (absent in pre-v2 receipts and empty-measurement runs).
        let key = self.parse_string()?;
        self.expect(b':')?;
        let (measurements, outputs) = if key == "measurements" {
            let m = self.parse_params()?;
            self.expect(b',')?;
            self.expect_key("outputs")?;
            (m, self.parse_file_array()?)
        } else if key == "outputs" {
            (BTreeMap::new(), self.parse_file_array()?)
        } else {
            return Err(self.err(&format!(
                "expected \"measurements\" or \"outputs\", got \"{key}\""
            )));
        };
        self.expect(b',')?;
        self.expect_key("params")?;
        let params = self.parse_params()?;
        self.expect(b',')?;
        self.expect_key("subcommand")?;
        let subcommand = self.parse_string()?;
        self.expect(b',')?;
        self.expect_key("tool_version")?;
        let tool_version = self.parse_string()?;
        self.expect(b'}')?;
        Ok(RunManifest {
            tool_version,
            subcommand,
            inputs,
            params,
            outputs,
            measurements,
        })
    }
```

- [ ] **Step 7: Update the module doc + `content_hash`/`finalize` references**

Update the module header comment (top of file) to describe the claim/measurement split in one
sentence. Update the `finalize_is_idempotent` / `tampering_*` test comments only if they
reference "no self-hash" (they don't — leave them).

- [ ] **Step 8: Run the provenance tests**

Run: `cargo test -p rosalind --lib provenance 2>&1 | tail -25`
Expected: all green, including the new Task-1 tests and the pre-existing ones
(`canonical_json_has_sorted_keys_and_is_exact`, `finalize_*`, `tampering_*`, `parse_*`).

- [ ] **Step 9: `cargo fmt` + commit**

```bash
cargo fmt
git add src/provenance/mod.rs docs/superpowers/specs/2026-06-02-p0-2-claim-measurement-split-design.md docs/superpowers/plans/2026-06-02-p0-2-claim-measurement-split.md
git commit -m "feat(provenance): split receipt into deterministic claim + machine-dependent measurement"
```

---

### Task 3: `verify` reads both maps + checks the measurement hash

**Files:**
- Modify: `src/main.rs` (`run_verify`, ~lines 885-1006)

- [ ] **Step 1: Route lookups through `get_recorded`**

In `run_verify`: change the budget fallback (`manifest.params.get("memory_budget_mb")` →
`manifest.get_recorded("memory_budget_mb")`), the peak read
(`manifest.params.get("peak_rss_bytes")` → `manifest.get_recorded("peak_rss_bytes")`), the
`recorded_u64` closure (`manifest.params.get(k)` → `manifest.get_recorded(k)`), and the
verdict read (`manifest.params.get("contract_verdict")` → `manifest.get_recorded("contract_verdict")`).

- [ ] **Step 2: Add the measurement-hash check + reframe the stale comment**

Replace the "Internal-consistency cross-checks" comment block's stale "no self-hash yet"
wording with the two-layer model, and after the existing `self_hash_ok()` match add:

```rust
    // Measurement attestation: catches an edit to a measured field (e.g. lowering
    // peak_rss_bytes to fake a fit) that the claim hash cannot see by design.
    match manifest.measurement_hash_ok() {
        Some(true) => {}
        Some(false) => problems.push(
            "measurement_blake3 mismatch: a measured field was modified after the run".to_string(),
        ),
        None => {}
    }
```

- [ ] **Step 3: Build + run the verify-related integration tests**

Run: `cargo test -p rosalind --test plan_enforce verify 2>&1 | tail -25`
Expected: `verify_rejects_a_self_consistent_but_tampered_receipt` and
`verify_rejects_an_internally_inconsistent_manifest` pass; others compile.

- [ ] **Step 4: `cargo fmt` + commit**

```bash
cargo fmt
git add src/main.rs
git commit -m "feat(verify): read measurement block + check measurement_blake3 (claim/measurement split)"
```

---

### Task 4: Update integration tests to the new partition

**Files:**
- Modify: `tests/plan_enforce.rs`, `tests/gvcf.rs`

- [ ] **Step 1: Point measured-field reads at `.measurements`**

In `tests/plan_enforce.rs`:
- `estimator_upper_bounds_the_realized_working_set`: `m.params.get("max_working_set_bytes")`
  → `m.measurements.get("max_working_set_bytes")`.
- `predicted_peak_rss_upper_bounds_realized_peak`: `m.params.get("predicted_peak_rss_bytes")`
  and `m.params.get("peak_rss_bytes")` → `m.measurements.get(...)`.
- `receipt_records_residual_and_governor_fields_on_a_fitting_run`:
  `m.params.get("governor")` → `m.measurements.get("governor")`.
- the residual/assumed test (~711-721): `predicted_peak_rss_bytes`, `peak_rss_bytes`,
  `rss_residual_bytes` → `m.measurements.get(...)`; **keep** `io_rss_overhead_assumed_bytes`
  on `m.params` (a constant — it stays in the claim).

In `tests/gvcf.rs`: `rm.params.contains_key("peak_rss_bytes")` →
`rm.measurements.contains_key("peak_rss_bytes")`.

- [ ] **Step 2: Bump the schema-version assertion + assert the measurement hash**

In `tests/plan_enforce.rs::receipt_is_self_hashing_and_schema_versioned`, change
`text.contains("\"schema_version\":\"1\"")` → `"\"schema_version\":\"2\""`, and add:

```rust
    assert!(
        text.contains("\"measurement_blake3\":"),
        "no measurement self-hash: {text}"
    );
```

- [ ] **Step 3: Add the lossless-measurement-tamper end-to-end test**

Append to `tests/plan_enforce.rs`:

```rust
#[test]
fn verify_rejects_a_consistent_measurement_tamper_via_the_measurement_hash() {
    // Lower the recorded peak while keeping the receipt internally consistent (peak
    // still within budget, verdict still 'within'). The claim hash cannot see a
    // measurement edit — only the measurement hash catches it.
    let (dir, idx, bam) = build_sorted_bam_fixture();
    let vcf = dir.join("calls.vcf");
    let manifest = dir.join("calls.vcf.manifest.json");
    let out = Command::new(bin())
        .args(["variants", "--index"])
        .arg(&idx)
        .arg("--alignments")
        .arg(&bam)
        .args(["--memory-budget-mb", "4096", "--enforce", "-o"])
        .arg(&vcf)
        .output()
        .unwrap();
    assert!(out.status.success(), "run failed: {out:?}");

    let text = std::fs::read_to_string(&manifest).unwrap();
    let m = rosalind::provenance::RunManifest::from_canonical_json(&text).unwrap();
    let recorded_peak = m.measurements.get("peak_rss_bytes").unwrap().clone();
    // Replace the measured peak with a smaller, still-within-budget value.
    let tampered = text.replace(
        &format!("\"peak_rss_bytes\":\"{recorded_peak}\""),
        "\"peak_rss_bytes\":\"1\"",
    );
    assert_ne!(text, tampered, "the replace must have changed the peak");
    std::fs::write(&manifest, &tampered).unwrap();

    let v = Command::new(bin())
        .args(["verify", "--manifest"])
        .arg(&manifest)
        .output()
        .unwrap();
    assert_eq!(
        v.status.code(),
        Some(5),
        "a measurement tamper must fail verify: {v:?}"
    );
    let stderr = String::from_utf8_lossy(&v.stderr);
    assert!(
        stderr.contains("measurement_blake3 mismatch"),
        "expected a measurement-hash mismatch: {stderr}"
    );
    std::fs::remove_dir_all(&dir).ok();
}
```

- [ ] **Step 4: Run the full integration suites**

Run: `cargo test -p rosalind --test plan_enforce --test gvcf --test features 2>&1 | tail -30`
Expected: all green.

- [ ] **Step 5: `cargo fmt` + commit**

```bash
cargo fmt
git add tests/plan_enforce.rs tests/gvcf.rs
git commit -m "test: adapt receipt reads to the claim/measurement split + measurement-tamper e2e"
```

---

### Task 5: Full verification + clippy + MSRV

- [ ] **Step 1: Whole test suite**

Run: `cargo test 2>&1 | tail -30`
Expected: all green (lib + every integration test).

- [ ] **Step 2: Clippy (the CI gate)**

Run: `cargo clippy --all-targets -- -D warnings 2>&1 | tail -20`
Expected: no warnings.

- [ ] **Step 3: MSRV build (the CI gate)**

Run: `cargo +1.83 build 2>&1 | tail -10` (skip if 1.83 toolchain absent — CI covers it)
Expected: builds clean.

- [ ] **Step 4: Push + open PR + watch CI**

```bash
git push -u origin rosalind/p0-2-claim-measurement-split
gh pr create --title "feat(provenance): cross-machine-stable receipts — claim/measurement split (P0.2)" --body "<summary>"
gh pr checks --watch
```

Merge only after all 5 CI jobs are green.
