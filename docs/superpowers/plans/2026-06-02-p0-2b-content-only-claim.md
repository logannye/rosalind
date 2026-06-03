# P0.2b — Content-only claim hash — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:executing-plans. Steps use `- [ ]`.

**Goal:** Make `manifest_blake3` a cross-machine content-address by hashing inputs/outputs as their
sorted `blake3` digests (dropping machine-specific paths) in the claim form, version-gated for
back-compat.

**Architecture:** A `FileRender` enum threads through the canonical serializer. The on-disk form
stays `{path, blake3}`; the claim form (for `schema_version >= 3`) becomes `["<blake3>",…]` sorted.
`finalize` bumps schema 2→3. `verify` unchanged (re-hashes from the on-disk paths).

**Reference:** spec `docs/superpowers/specs/2026-06-02-p0-2b-content-only-claim-design.md`.

---

### Task 1: Unit tests (RED)

**Files:** Modify `src/provenance/mod.rs` (append to `mod tests`).

- [ ] **Step 1: Write the failing tests**

```rust
    #[test]
    fn claim_hash_is_stable_across_machine_dependent_paths() {
        // Same content (blake3) at DIFFERENT paths on two machines → SAME claim hash.
        // The keystone of P0.2b: paths are dropped from the claim form.
        let mk = |idx_path: &str, out_path: &str| {
            let mut m = RunManifest::new("variants");
            m.tool_version = "0.1.0".to_string();
            m.inputs.push(FileHash {
                path: idx_path.to_string(),
                blake3: "aa".to_string(),
            });
            m.outputs.push(FileHash {
                path: out_path.to_string(),
                blake3: "bb".to_string(),
            });
            m.params.insert("min_qual".to_string(), "30".to_string());
            m.finalize();
            m
        };
        let a = mk("/home/alice/ref.idx", "/tmp/run-1/out.vcf");
        let b = mk("/data/ref.idx", "out.vcf");
        assert_eq!(
            a.content_hash(),
            b.content_hash(),
            "the claim hash must not depend on recorded paths"
        );
        assert_eq!(a.self_hash_ok(), Some(true));
        assert_eq!(b.self_hash_ok(), Some(true));
    }

    #[test]
    fn claim_hash_still_tracks_content() {
        // Different content (blake3) → different claim hash (the address is the content).
        let mk = |digest: &str| {
            let mut m = RunManifest::new("variants");
            m.tool_version = "0.1.0".to_string();
            m.inputs.push(FileHash {
                path: "ref.idx".to_string(),
                blake3: digest.to_string(),
            });
            m.finalize();
            m
        };
        assert_ne!(mk("aa").content_hash(), mk("bb").content_hash());
    }

    #[test]
    fn claim_drops_paths_but_the_on_disk_form_keeps_them() {
        let mut m = RunManifest::new("variants");
        m.tool_version = "0.1.0".to_string();
        m.inputs.push(FileHash {
            path: "/home/alice/secret/ref.idx".to_string(),
            blake3: "aa".to_string(),
        });
        m.finalize();
        assert!(
            !m.to_canonical_claim_json().contains("/home/alice"),
            "the claim form must not carry the recorded path"
        );
        assert!(
            m.to_canonical_json().contains("/home/alice"),
            "the on-disk form must keep the recorded path"
        );
        // The content digest is present in BOTH.
        assert!(m.to_canonical_claim_json().contains("aa"));
    }

    #[test]
    fn pre_p0_2b_receipt_with_paths_in_the_claim_still_verifies() {
        // A schema-2 receipt hashed paths INTO the claim. The version gate must
        // reproduce that path-inclusive form so it still self-verifies; a schema-3
        // receipt over the same files hashes differently (the form changed).
        let mk = |schema: &str| {
            let mut m = RunManifest::new("variants");
            m.tool_version = "0.1.0".to_string();
            m.inputs.push(FileHash {
                path: "ref.idx".to_string(),
                blake3: "aa".to_string(),
            });
            m.params
                .insert("schema_version".to_string(), schema.to_string());
            let h = m.content_hash();
            m.params.insert("manifest_blake3".to_string(), h);
            m
        };
        let v2 = mk("2");
        assert_eq!(v2.self_hash_ok(), Some(true), "schema-2 must self-verify");
        let v3 = mk("3");
        assert_eq!(v3.self_hash_ok(), Some(true), "schema-3 must self-verify");
        assert_ne!(
            v2.params.get("manifest_blake3"),
            v3.params.get("manifest_blake3"),
            "the claim form genuinely differs between schema 2 and 3"
        );
    }
```

- [ ] **Step 2: Confirm RED**

Run: `cargo test -p rosalind --lib provenance 2>&1 | tail -15`
Expected: `claim_hash_is_stable_across_machine_dependent_paths` FAILS (paths still in the claim →
`a` and `b` differ); `claim_drops_paths_...` FAILS; the back-compat test FAILS (no version gate yet,
so schema-2 and schema-3 hash identically).

---

### Task 2: Implementation (GREEN)

**Files:** Modify `src/provenance/mod.rs`.

- [ ] **Step 1: Bump the schema version**

```rust
/// Current receipt/feature schema version. Bump on any breaking schema change.
/// v2: claim/measurement split. v3: the claim hashes inputs/outputs by their sorted
/// `blake3` digests (paths dropped) so it is a cross-machine content-address.
pub const MANIFEST_SCHEMA_VERSION: u32 = 3;
```

- [ ] **Step 2: Add the `FileRender` enum + thread it through `push_canonical`**

Add above `RunManifest`:

```rust
/// How input/output file entries render in a canonical form. The on-disk receipt
/// keeps full `{path, blake3}`; the schema-3 claim drops the path and hashes only the
/// content digest, so the claim hash does not depend on where files live.
#[derive(Clone, Copy)]
enum FileRender {
    WithPath,
    ContentOnly,
}
```

Replace `push_canonical` to take a `FileRender` and dispatch the file arrays:

```rust
    fn push_canonical(&self, include_measurements: bool, files: FileRender) -> String {
        let mut out = String::new();
        out.push('{');
        out.push_str("\"inputs\":");
        push_files(&mut out, &self.inputs, files);
        if include_measurements && !self.measurements.is_empty() {
            out.push_str(",\"measurements\":");
            push_string_map(&mut out, &self.measurements);
        }
        out.push_str(",\"outputs\":");
        push_files(&mut out, &self.outputs, files);
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

- [ ] **Step 3: Update the two public canonical methods + add the version gate**

```rust
    pub fn to_canonical_json(&self) -> String {
        self.push_canonical(true, FileRender::WithPath)
    }

    pub fn to_canonical_claim_json(&self) -> String {
        self.push_canonical(false, self.claim_file_render())
    }

    /// Schema >= 3 → content-only claim (paths dropped, cross-machine). Older receipts
    /// hashed paths into the claim; reproduce their form so they still self-verify.
    fn claim_file_render(&self) -> FileRender {
        match self
            .params
            .get("schema_version")
            .and_then(|v| v.parse::<u32>().ok())
        {
            Some(v) if v >= 3 => FileRender::ContentOnly,
            _ => FileRender::WithPath,
        }
    }
```

- [ ] **Step 4: Add the `push_files` dispatcher + `push_blake3_list`**

Next to `push_file_hashes`:

```rust
/// Render a file array in the requested form: `[{"blake3","path"}]` sorted by path
/// (on-disk), or `["<blake3>",…]` sorted by blake3 (the content-only claim).
fn push_files(out: &mut String, files: &[FileHash], render: FileRender) {
    match render {
        FileRender::WithPath => push_file_hashes(out, files),
        FileRender::ContentOnly => push_blake3_list(out, files),
    }
}

/// Render `["<blake3>",…]`, sorted by digest (a content multiset — duplicates kept).
fn push_blake3_list(out: &mut String, files: &[FileHash]) {
    let mut digests: Vec<&str> = files.iter().map(|f| f.blake3.as_str()).collect();
    digests.sort_unstable();
    out.push('[');
    for (i, d) in digests.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        out.push('"');
        out.push_str(&json_escape(d));
        out.push('"');
    }
    out.push(']');
}
```

- [ ] **Step 5: Update doc comments for the dropped-paths reality**

Update the module "Scope note" and `content_hash`'s "(Recorded paths remain in the claim…)" note to
say paths are now dropped from the schema-3 claim. Update the P0.2-era line accordingly.

- [ ] **Step 6: Run provenance tests (GREEN)**

Run: `cargo test -p rosalind --lib provenance 2>&1 | tail -28`
Expected: all green (new P0.2b tests + every P0.2 test, including the v1/v2 back-compat ones).

- [ ] **Step 7: `cargo fmt` + commit**

```bash
cargo fmt
git add src/provenance/mod.rs docs/superpowers/specs/2026-06-02-p0-2b-content-only-claim-design.md docs/superpowers/plans/2026-06-02-p0-2b-content-only-claim.md
git commit -m "feat(provenance): content-only claim hash — cross-machine content-address (P0.2b)"
```

---

### Task 3: Integration assertion + full verification

**Files:** Modify `tests/plan_enforce.rs`.

- [ ] **Step 1: Bump the schema assertion**

In `receipt_is_self_hashing_and_schema_versioned`: `"\"schema_version\":\"2\""` → `"\"schema_version\":\"3\""`.

- [ ] **Step 2: Full suite + clippy + MSRV**

```bash
cargo test 2>&1 | grep -E "test result: FAILED|FAILED|^error"   # expect none
cargo clippy --all-targets -- -D warnings 2>&1 | tail -3        # clean
cargo +1.83 build 2>&1 | tail -3                                # MSRV (CI covers if absent)
```

- [ ] **Step 3: `cargo fmt` + commit**

```bash
cargo fmt
git add tests/plan_enforce.rs
git commit -m "test: receipt now schema 3 (content-only claim)"
```

- [ ] **Step 4: Push, PR, watch CI, merge on green**

```bash
git push -u origin rosalind/p0-2b-content-only-claim
gh pr create --title "feat(provenance): content-only claim hash — cross-machine content-address (P0.2b)" --body "<summary>"
gh pr checks --watch
```
