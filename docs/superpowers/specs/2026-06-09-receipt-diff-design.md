# Design: `rosalind diff <a> <b>` — a claim-level divergence localizer

**Status:** Approved design — 2026-06-09. Audience: the implementer. Companion to the receipt model
(`crates/receipt/src/lib.rs`), `reproduce` (`src/reproduce.rs`), and `chain verify`.

---

## 1. The goal, in one sentence

Turn divergence/tamper-evidence from *scary* ("a byte changed") into a one-line **localization**:
which hashed field differs between two receipts, bucketed by **causal role**.

## 2. Why this, why now

`verify` says a receipt was edited; `reproduce` says a re-run did/didn't match; `chain verify` walks
the provenance graph. But nothing answers *"how do these two receipts differ, and what's the likely
cause?"* — the natural question when a result diverges or two runs disagree. `reproduce::classify_outputs`
compares a receipt's recorded outputs against **freshly produced** ones (positional) — it never compares
**two receipts to each other**, so this is non-overlapping. Because Rosalind's claim is a typed,
content-addressed, code-identity-stamped object, the diff can do what a generic `jq`/textual diff can't:
**bucket the differences by causal salience** ("only `--alignments` changed; the output change is its
effect" vs "outputs differ with identical inputs/params/code → nondeterminism or corruption").

## 3. The core (`crates/receipt/src/diff.rs`)

Pure, `std`-only, wasm-clean (reuses `RunManifest`):

```rust
/// One differing scalar claim/measurement field. `None` = absent on that side.
pub struct FieldChange { pub key: String, pub a: Option<String>, pub b: Option<String> }

/// One differing input/output operand, labeled by its CLI flag (recovered from `command`).
pub struct OperandChange { pub flag: String, pub a: Option<String>, pub b: Option<String> }

pub struct ReceiptDiff {
    pub subcommand: Option<(String, String)>, // Some((a,b)) iff they differ
    pub inputs: Vec<OperandChange>,
    pub outputs: Vec<OperandChange>,
    pub code_identity: Vec<FieldChange>,
    pub science_params: Vec<FieldChange>,
    pub measurements: Vec<FieldChange>,
    /// `content_hash(a) == content_hash(b)` — the cross-machine claim addresses match.
    pub claims_identical: bool,
}

pub fn diff_receipts(a: &RunManifest, b: &RunManifest) -> ReceiptDiff;
```

**Bucketing rules:**
- `claims_identical = a.content_hash() == b.content_hash()`.
- **subcommand:** `Some((a,b))` iff `a.subcommand != b.subcommand`.
- **inputs / outputs:** a small `diff.rs`-local helper parses each receipt's `command` param for
  `<flag> @in:<h>` and `<flag> @out:<h>` operand pairs (the chain walker has a private `@in:`-only
  variant; `diff` needs both). Match by **flag** (robust to operand order);
  union the flags; emit an `OperandChange` for each flag whose hash differs or is absent on one side.
  Fallback when `command` is absent: positional, labeled `input[i]` / `output[i]`.
- **params** (the `params` BTreeMap): union the keys; for each key whose value differs (or is one-sided):
  - in **`BUILD_IDENTITY_KEYS`** → `code_identity`;
  - `manifest_blake3` (the self-hash, derived) or `command` (the recipe — its operands/opts are already
    surfaced) → **skip** (derived/redundant noise);
  - else → `science_params`.
- **measurements** (the `measurements` BTreeMap): union the keys; differing → `measurements`, skipping
  the derived `measurement_blake3`.

**New `pub const BUILD_IDENTITY_KEYS`** in `lib.rs` (next to `MEASUREMENT_KEYS`):
```rust
pub const BUILD_IDENTITY_KEYS: &[&str] =
    &["code_git_sha", "code_dirty", "rustc_version", "target_triple", "deps_lock_blake3"];
```
A unit test asserts `build_identity_pairs()`'s keys equal `BUILD_IDENTITY_KEYS` so they cannot drift.

## 4. The causal verdict + exit codes

A `ReceiptDiff::verdict() -> String` and `exit_code() -> i32`:
- **`claims_identical`** → `IDENTICAL claims`; if `measurements` is non-empty, append
  *"— only machine-dependent measurements differ (<keys>)"*. **Exit 0.**
- else, **a cause is present** (`subcommand` / `inputs` / `science_params` / `code_identity` non-empty)
  → *"claims DIFFER — cause: <bucket summary>"*; the `outputs` bucket is framed as the **effect**.
  **Exit 1.**
- else, **only `outputs` differ** (inputs+params+code identical) → *"claims DIFFER — outputs differ with
  identical inputs/params/code → nondeterminism or corruption"*. **Exit 1.**
- Read/parse failure (handled in the CLI) → **exit 2**.

(Because `content_hash()` excludes only `manifest_blake3`, `claims_identical == false` always implies
some non-derived field differs, so a cause/effect bucket always catches it.)

## 5. The CLI (`src/main.rs`)

```
rosalind diff <a.manifest.json> <b.manifest.json>          # human report
rosalind diff <a.manifest.json> <b.manifest.json> --json   # the structured ReceiptDiff
```
`run_diff` reads both files, parses via `RunManifest::from_canonical_json` (on error: stderr + exit 2),
calls `diff_receipts`, renders, and exits `0`/`1`. A new `Diff` clap variant + dispatch arm.
`--json` emits a compact hand-rolled object (like `ChainReport::to_json`) — no serde.

Human report shape:
```
diff: variants  vs  variants
CAUSE — input  --alignments  a1b2c3… → d4e5f6…
EFFECT — output -o            9f8e7d… → 1a2b3c…
(inputs[--index], all science params, code-identity: identical)
VERDICT: claims DIFFER — localized to 1 input (--alignments); the output change is its effect.
```

## 6. Testing

**Unit (`crates/receipt/src/diff.rs`):** build pairs of `RunManifest`s (via `new` + `params`/`inputs`/
`outputs` + `finalize`) and assert the bucketing + verdict + exit code:
- single input hash change → one `inputs` `OperandChange`, `exit 1`, verdict names the input.
- single science-param change (`max_depth`) → one `science_params` `FieldChange`.
- code-identity drift (`code_git_sha`) → one `code_identity` change, segregated from science params.
- outputs differ but inputs+params+code identical → the *nondeterminism/corruption* verdict.
- identical claims, only a measurement differs (`peak_rss_bytes`) → `claims_identical`, `exit 0`,
  the measurements note.
- `BUILD_IDENTITY_KEYS` matches `build_identity_pairs()` keys.

**Integration (`tests/diff.rs`):** produce two real receipts via the CLI (e.g. a `features` run, then the
same run with one input file's bytes flipped) and assert `rosalind diff` exits 1 and the report names the
changed input; and that diffing a receipt against itself exits 0 (`IDENTICAL`).

**Gates:** `cargo test`, `cargo clippy --all-targets -- -D warnings`, `cargo fmt --check`.

## 7. Non-goals / guardrails

- **Claim-level only.** It names the differing hashed field; it does **not** open a VCF/TSV or diff
  BAM bytes (a swamp).
- **No re-derivation.** It compares two existing receipts; re-running is `reproduce`'s job.
- **No overlap with `classify_outputs`** (recorded-vs-produced, positional) — `diff` is receipt-vs-receipt.
- Pure `std` in the receipt crate (wasm-portable for a future browser side-by-side); the CLI does the I/O.
