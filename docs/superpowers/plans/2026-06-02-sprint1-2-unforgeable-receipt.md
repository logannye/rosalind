# Sprint 1.2 — The Unforgeable Receipt: Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the run receipt tamper-evident (a `manifest_blake3` self-hash re-derived by `verify`) and schema-versioned (`schema_version`), on the existing canonical-JSON machinery.

**Architecture:** Add `content_hash`/`finalize`/`self_hash_ok` + a `MANIFEST_SCHEMA_VERSION` const to `provenance::RunManifest`. `finalize()` (called once before each of the 4 receipt writes) stamps the schema version then the self-hash over the canonical JSON with the hash field excluded. `verify` re-derives and checks it; a pre-1.2 receipt (no self-hash) verifies with a note.

**Tech Stack:** Rust, the existing `blake3_hex` + canonical-JSON serializer in `src/provenance/mod.rs`.

**Spec:** `docs/superpowers/specs/2026-06-02-sprint1-2-unforgeable-receipt-design.md`

---

### Task 1: Self-hash + schema-version primitives (`provenance::RunManifest`)

**Files:**
- Modify: `src/provenance/mod.rs`

- [ ] **Step 1: Write the failing unit tests** (append to the `#[cfg(test)] mod tests` in `src/provenance/mod.rs`)

```rust
    #[test]
    fn finalize_stamps_schema_version_and_a_matching_self_hash() {
        let mut m = RunManifest::new("variants");
        m.tool_version = "0.1.0".to_string();
        m.params.insert("peak_rss_bytes".to_string(), "123".to_string());
        m.finalize();
        assert_eq!(
            m.params.get("schema_version").map(String::as_str),
            Some(MANIFEST_SCHEMA_VERSION.to_string().as_str())
        );
        assert!(m.params.contains_key("manifest_blake3"));
        assert_eq!(m.self_hash_ok(), Some(true), "fresh finalize must verify");
    }

    #[test]
    fn tampering_any_field_breaks_the_self_hash() {
        let mut m = RunManifest::new("variants");
        m.tool_version = "0.1.0".to_string();
        m.params.insert("peak_rss_bytes".to_string(), "123".to_string());
        m.finalize();
        // Flip a field WITHOUT re-finalizing — the recorded hash no longer matches.
        m.params
            .insert("peak_rss_bytes".to_string(), "999".to_string());
        assert_eq!(m.self_hash_ok(), Some(false));
    }

    #[test]
    fn a_manifest_without_a_self_hash_returns_none() {
        let mut m = RunManifest::new("variants");
        m.tool_version = "0.1.0".to_string();
        assert_eq!(m.self_hash_ok(), None);
    }

    #[test]
    fn finalize_is_idempotent() {
        let mut m = RunManifest::new("variants");
        m.tool_version = "0.1.0".to_string();
        m.params.insert("peak_rss_bytes".to_string(), "123".to_string());
        m.finalize();
        let first = m.params.get("manifest_blake3").cloned();
        m.finalize();
        assert_eq!(m.params.get("manifest_blake3").cloned(), first);
        assert_eq!(m.self_hash_ok(), Some(true));
    }
```

- [ ] **Step 2: Run to verify it fails**

Run: `cargo test --lib provenance 2>&1 | grep -E "error|test result"`
Expected: FAIL — `MANIFEST_SCHEMA_VERSION`, `finalize`, `self_hash_ok` do not exist.

- [ ] **Step 3: Add the const + methods**

Add the const near the top of `src/provenance/mod.rs` (after the module doc-comment / imports, before `ManifestError`):

```rust
/// Current receipt/feature schema version. Bump on any breaking schema change.
pub const MANIFEST_SCHEMA_VERSION: u32 = 1;
```

Add the methods to the existing `impl RunManifest { … }` block (after `from_canonical_json`):

```rust
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
        self.params.insert(
            "schema_version".to_string(),
            MANIFEST_SCHEMA_VERSION.to_string(),
        );
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
```

- [ ] **Step 4: Run to verify it passes**

Run: `cargo test --lib provenance 2>&1 | grep -E "test result"`
Expected: PASS (the four new tests + the existing provenance tests).

- [ ] **Step 5: Verify no warnings, fmt, commit**

Run: `cargo build --release 2>&1 | grep -i warning || echo "no warnings"; cargo fmt && cargo fmt --check && echo "fmt clean"`

```bash
git add src/provenance/mod.rs
git commit -m "feat(provenance): self-hash + schema_version primitives (finalize/content_hash/self_hash_ok)"
```

---

### Task 2: Stamp at the 4 receipt sites + check in `verify`

**Files:**
- Modify: `src/main.rs` (4 receipt-write sites + `run_verify`)
- Modify: `CONTRACT.md`
- Test: `tests/plan_enforce.rs`

- [ ] **Step 1: Write the failing integration tests** (append to `tests/plan_enforce.rs`)

```rust
#[test]
fn receipt_is_self_hashing_and_schema_versioned() {
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
    let text = std::fs::read_to_string(&manifest).expect("manifest");
    assert!(text.contains("\"manifest_blake3\":"), "no self-hash: {text}");
    assert!(text.contains("\"schema_version\":\"1\""), "no schema_version: {text}");
    // An untampered receipt verifies.
    let v = Command::new(bin())
        .args(["verify", "--manifest"])
        .arg(&manifest)
        .output()
        .unwrap();
    assert!(v.status.success(), "untampered verify should pass: {v:?}");
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn verify_rejects_a_self_consistent_but_tampered_receipt() {
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
    // Tamper a field while keeping the receipt internally CONSISTENT: raise the
    // recorded budget (a 'within' verdict stays valid since peak << 8192 MiB), so the
    // verdict/consistency cross-checks still pass. Only the self-hash catches this.
    let text = std::fs::read_to_string(&manifest).unwrap();
    let tampered =
        text.replace("\"memory_budget_mb\":\"4096\"", "\"memory_budget_mb\":\"8192\"");
    assert_ne!(text, tampered, "the replace must have changed something");
    std::fs::write(&manifest, &tampered).unwrap();
    let v = Command::new(bin())
        .args(["verify", "--manifest"])
        .arg(&manifest)
        .output()
        .unwrap();
    assert_eq!(v.status.code(), Some(5), "tampered receipt must fail verify: {v:?}");
    let stderr = String::from_utf8_lossy(&v.stderr);
    assert!(
        stderr.contains("manifest_blake3 mismatch"),
        "expected a self-hash mismatch: {stderr}"
    );
    std::fs::remove_dir_all(&dir).ok();
}
```

- [ ] **Step 2: Run to verify they fail**

Run: `cargo test --test plan_enforce receipt_is_self_hashing verify_rejects_a_self_consistent -- --nocapture 2>&1 | grep -E "test result|FAILED"`
Expected: FAIL — no `manifest_blake3`/`schema_version` in the receipt yet, and `verify` does not check a self-hash.

- [ ] **Step 3: Call `finalize()` at the 4 receipt-write sites**

In `src/main.rs`, add `manifest.finalize();` immediately before the serialization at each site. The bindings are already `let mut manifest = RunManifest::new(…)`.

(a) `run_somatic` — before `let manifest_path = write_manifest(&output_vcf, &manifest)?;` (~`:1168`):

```rust
    manifest.finalize();
    let manifest_path = write_manifest(&output_vcf, &manifest)?;
```

(b) `run_variants` (`--reference`) — before `let manifest_path = write_manifest(&path, &manifest)?;` (~`:1538`):

```rust
            manifest.finalize();
            let manifest_path = write_manifest(&path, &manifest)?;
```

(c) `run_features` — before `std::fs::write(&dest, manifest.to_canonical_json())` (~`:1838`):

```rust
        manifest.finalize();
        std::fs::write(&dest, manifest.to_canonical_json())
```

(d) `run_variants_index` — before `std::fs::write(&dest, manifest.to_canonical_json())` (~`:2218`):

```rust
        manifest.finalize();
        std::fs::write(&dest, manifest.to_canonical_json())
```

> If two sites share identical surrounding text and an exact-match edit is ambiguous, include the
> distinguishing preceding line (the last `manifest.params.insert(…)` of that block) in the match.

- [ ] **Step 4: Add the self-hash check in `run_verify`**

In `src/main.rs` `run_verify`, after the internal-consistency cross-checks and before the
`if problems.is_empty()` block, add:

```rust
    // Self-hash: catches any post-write edit (even one that keeps the other fields
    // mutually consistent). A pre-1.2 receipt has no self-hash — note and skip.
    match manifest.self_hash_ok() {
        Some(true) => {}
        Some(false) => problems.push(
            "manifest_blake3 mismatch: the receipt was modified after it was written".to_string(),
        ),
        None => println!("verify: note — no manifest_blake3 (a pre-1.2 receipt); skipping self-hash"),
    }
```

- [ ] **Step 5: Run the new tests + the existing verify tests**

Run: `cargo test --test plan_enforce verify -- --nocapture 2>&1 | grep -E "test result|FAILED"`
Expected: PASS — the two new tests, plus the existing `verify_passes_on_an_untampered_run…` and
`verify_rejects_an_internally_inconsistent_manifest` (the latter now also trips the self-hash but
still exits 5).

- [ ] **Step 6: Document in CONTRACT.md**

In `CONTRACT.md`, in the `### 4. Verify` section, after the sentence describing what `verify` re-checks,
add:

```markdown
The receipt is **self-hashing**: a `manifest_blake3` over its own canonical JSON (and a `schema_version`)
is stamped at write time and re-derived by `verify`, so any post-write edit — even one that keeps the
other fields mutually consistent — is caught (exit 5). This is tamper-*evident*; a cryptographically
signed, tamper-*proof* receipt is a separate planned feature.
```

- [ ] **Step 7: Full suite + warnings + fmt, then commit**

Run: `cargo test 2>&1 | grep -cE "test result: ok"; cargo test 2>&1 | grep -iE "FAILED" || echo "no failures"; cargo build --release 2>&1 | grep -i warning || echo "no warnings"; cargo fmt && cargo fmt --check && echo "fmt clean"`
Expected: all sections ok, no failures, no warnings, fmt clean.

```bash
git add src/main.rs tests/plan_enforce.rs CONTRACT.md
git commit -m "feat(provenance): stamp the self-hash + schema_version on every receipt; verify checks it"
```

---

## Final verification

- [ ] **Full suite green:** `cargo test 2>&1 | tail -5`.
- [ ] **Zero warnings:** `cargo build --release 2>&1 | grep -i warning || echo "no warnings"` (lib + examples).
- [ ] **Smoke:** a real `variants --index … --enforce -o calls.vcf` receipt contains `manifest_blake3` + `schema_version`; `verify` passes; editing one byte of `calls.vcf.manifest.json` then `verify` exits 5 with "manifest_blake3 mismatch".
- [ ] CI green on GitHub before reporting done.
