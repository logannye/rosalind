# Track D — `rosalind reproduce` + the reproduction certificate (design)

**Status:** Approved design — 2026-06-05. Author: Logan Nye. Companion to
[`docs/ROADMAP.md`](../../ROADMAP.md) (P3.2 + the reproducibility CI badge) and
[`CONTRACT.md`](../../../CONTRACT.md). Implements the standalone-roadmap critical-path track
(`project_rosalind_roadmap`): third-party byte re-derivation from a receipt, plus a
reproducibility/resource CI fence + self-hosted badge.

This is a **standalone** feature — no dependency on any other repository. It builds only on
Rosalind's already-shipped determinism + content-addressed receipt.

---

## 1. Goal

Make **"reproduce this exact result"** a command a stranger can run — and turn each successful
reproduction into a **portable, content-addressed certificate that chains to the original**, so
independent reproductions accumulate into a decentralized "reproducibility web" with no server.

The flagship demo this unlocks: a stranger takes a real `*.vcf` + its `*.manifest.json`, on a
*different machine*, runs one offline command, and gets a green **REPRODUCED** (byte-identical) —
with no GATK, no Docker, no Nextflow, no re-aligning — plus a `.repro.json` certificate that
countersigns the original claim. Flip one input byte → **DIVERGED**, naming the exact field. No
incumbent caller can do this (a non-deterministic caller reports DIVERGED on a *correct* run).

## 2. Non-goals (explicit)

- **Not** a correctness/biology claim. `reproduce` attests *re-derivation of output bytes from
  content-addressed inputs under a recorded code identity* — never that the calls are biologically
  right. The verdict text and docs say so.
- **Not** cryptographic authentication (yet). The certificate is tamper-**evident** (self-hashing
  BLAKE3), **signing-ready** for the later Ed25519 track (P3.3) — not tamper-**proof** today.
- **Not** BAM/bgzf reproduction in v1 (see §7 honesty scope).
- **Not** chain-of-artifacts traversal — that is Track C (`verify --chain`). `reproduce` operates
  on a single receipt; the two compose later but are independent here.
- **No cross-repo integration** of any kind.

## 3. Background — the current code (verified 2026-06-05)

- `src/provenance/mod.rs`: `RunManifest { tool_version, subcommand, inputs: Vec<FileHash>,
  params: BTreeMap<String,String>, outputs: Vec<FileHash>, measurements: BTreeMap<String,String> }`.
  `params` is the deterministic **claim**; `measurements` is machine-dependent and excluded from
  the claim hash (`MEASUREMENT_KEYS`, mod.rs:37). `content_hash()` (mod.rs:183) = BLAKE3 of the
  claim canonical JSON with `manifest_blake3` removed and (schema ≥ 3) paths dropped — a
  cross-machine content-address. `MANIFEST_SCHEMA_VERSION = 4` (mod.rs:31). `finalize()` (mod.rs:259)
  partitions measurements out, stamps `measurement_blake3`, `schema_version`, then `manifest_blake3`
  last. `from_canonical_json` (mod.rs:170) is a hand-parser round-tripping `to_canonical_json`.
  `write_manifest` (mod.rs:594) writes the sidecar.
- **The gap (confirmed):** the `subcommand` string is recorded, but the *args needed to replay it*
  are written ad hoc per call site and several are missing — `gvcf` (main.rs:1977, unrecorded),
  `chrom` (arg at main.rs:84, omitted in the receipt at ~main.rs:1599-1622), and the
  index-vs-reference **mode** (only inferable from input file type). A naive replay would run the
  **wrong** command and report a **false DIVERGED**.
- Receipt-writing sites are not centralized: `run_somatic` (~main.rs:1229-1250), `run_variants`
  (reference path, ~1599-1623), `run_features` (~1852-1923), `run_variants_index` (variants + gVCF,
  ~2233-2305). All call `RunManifest::new(...)` + ad-hoc `params` inserts + `finalize()` +
  `write_manifest`.
- `run_verify` (main.rs:890-1060) is the verb template: parse → re-derive → per-check report →
  exit 5 on drift. It is **not** factored into a reusable library function today.
- `action.yml` is the shipped composite Action (plan → variants --enforce → upload receipt).

## 4. Components & module layout (each one purpose, testable in isolation)

| New/changed | File | Responsibility |
|---|---|---|
| **new** | `src/provenance/command.rs` | `CommandCapture` — the single capture chokepoint: record a normalized, replayable invocation into the claim, and reconstruct an argv from it. |
| **new** | `src/provenance/repro.rs` | `ReproReceipt` — the reproduction certificate: type, self-hash, chaining, canonical JSON render/parse. |
| **new** | `src/provenance/badge.rs` | `badge` — emit a self-hosted shields-endpoint JSON + a static SVG. |
| **new** | `src/reproduce.rs` | the `reproduce` driver (fs + blake3 + process only; no htslib): locate inputs by content hash, re-exec, compare, mint the certificate. |
| **changed** | `src/provenance/mod.rs` | bump schema 4→5 (version-gated); add `verify_receipt` (extracted, see below); wire `CommandCapture`. |
| **changed** | `src/main.rs` | adopt `CommandCapture` at all four receipt sites; add `Reproduce` + `Badge` clap variants + thin dispatch; extract `run_verify`'s core into the library `verify_receipt`. |
| **changed** | `action.yml` (+ a new CI workflow) | the reproduce fence (cross-machine REPRODUCED + one-byte-flip negative) and badge publication. |

**Anti-drift extraction:** lift the receipt-checking core of `run_verify` (main.rs:890-1060) into a
library `verify_receipt(receipt_text, opts) -> VerifyReport`, consumed by both `verify` and
`reproduce` (and the future wasm verifier, Track B) so the three cannot diverge. `run_verify`
becomes a thin CLI shell over it.

## 5. Schema 5 — normalized command capture

A single builder, used at every receipt site, that is the same structure `reproduce` replays — so
recording and replay cannot drift.

```rust
// src/provenance/command.rs
let mut cmd = CommandCapture::new("variants");
cmd.input("--index", &index_path, &index_blake3);   // content-addressed input operand
cmd.input("--alignments", &bam_path, &bam_blake3);
cmd.opt("--mapq-threshold", mapq.to_string());      // defaults resolved explicitly
cmd.opt("--max-depth", max_depth.to_string());
cmd.opt("--max-read-len", max_read_len.to_string());
cmd.opt("--memory-budget-mb", budget.to_string());
cmd.flag_if(enforce, "--enforce");
cmd.flag_if(gvcf, "--gvcf");
cmd.output("-o", &out_path, &out_blake3);           // content-addressed output operand
cmd.record_into(&mut manifest);                      // the ONE place that writes the claim command
```

`record_into`:
1. **Populates `inputs[]` and `outputs[]`** from the `input()`/`output()` calls — so these stop
   being hand-maintained per site (drift eliminated; one source of truth).
2. **Writes `params["command"]`** — a deterministic, machine-independent normalized form where
   input/output operands are replaced by placeholders `@in:<blake3>` / `@out:<blake3>`. Example
   recorded value:
   `variants --index @in:<h1> --alignments @in:<h2> --mapq-threshold 20 --max-depth 1000 --max-read-len 250 --memory-budget-mb 256 --enforce -o @out:<h3>`
   Because operands are content hashes (not paths), `params["command"]` is hash-protected by the
   existing claim hash (editing the recorded command is caught) and stays a cross-machine
   content-address.
3. **Token ordering is normalized** (a fixed canonical order: positional/required flags, then
   options sorted by flag, then bare flags, then output) so the recorded command is stable
   run-to-run and the claim hash is meaningful.

**Relationship to discrete params (no drift):** `params["command"]` is the **authoritative replay
recipe**. To preserve human-readability and the existing back-compat `verify` behavior, the *same*
builder calls also write discrete claim params via a **mechanical** flag→key mapping (strip `--`,
dashes→underscores: `--max-depth` → `max_depth`, `--gvcf` → `gvcf`, etc.), so the discrete params
cannot drift from `command` — they are projections of the same calls, not a second hand-maintained
list. The index-vs-reference **mode** is captured implicitly by which input flag the recipe carries
(`--index` vs `--reference`) and surfaced as a discrete `mode` param for readability. Build-identity
(`code_git_sha`, …) and contract/measurement fields are untouched — still set by `finalize()` /
`record_measurement()`.

**Reconstruction** (`CommandCapture::to_argv(resolver)`): parse `params["command"]`, substitute
each `@in:<h>` with the resolver's content-located path and each `@out:<h>` with a caller-provided
temp path, yielding a `Vec<String>` argv for `current_exe()`. **No per-subcommand logic** lives in
`reproduce`; new subcommands/flags are captured automatically by using the builder.

**Schema bump 4 → 5**, version-gated exactly like the prior bumps (`claim_file_render` pattern,
mod.rs:131): schema-5 receipts carry `params["command"]`; schema ≤ 4 receipts have no `command`
and `reproduce` reports INCONCLUSIVE ("pre-schema-5 receipt: no recorded command to replay") rather
than guessing. `verify` continues to handle ≤ 4 unchanged.

## 6. `reproduce` — algorithm

`rosalind reproduce --manifest <m> --inputs <dir> [--no-attest] [-o <repro.json>]`

1. **Read + integrity-check** the receipt via `verify_receipt` (self-hash intact? schema known?).
   Tampered/malformed → exit 5 (reuse verify's code) with the per-check report.
2. **Schema/command gate:** schema ≥ 5 with a `command` param, and a supported output type
   (§7). Otherwise → INCONCLUSIVE (exit 7) with the reason.
3. **Content-locate inputs:** index every file under `--inputs` by BLAKE3 once
   (`hash → path`); bind each `@in:<h>` placeholder. Any unresolved input → INCONCLUSIVE
   (exit 7), naming the missing hash (it is "couldn't check," not "the result changed").
4. **Bind outputs** to fresh temp paths; **re-exec `current_exe()`** with the reconstructed argv
   (re-running enforces the recorded budget too, so the run's own contract still applies).
5. **Compare:** BLAKE3 each produced output vs the recorded `outputs[].blake3`.
6. **Resource-aware context:** read the re-run's realized peak from its temp receipt; report
   `peak X ≤/> declared Y (here: this machine)` — diagnostic, explicitly machine-local.
7. **Code context:** compare the current binary's build-identity (`code_git_sha`, `code_dirty`, …)
   to the receipt's; surface a mismatch prominently. A byte-match under a *different* code SHA is
   still REPRODUCED (a stronger result), with the difference noted.
8. **Verdict (on output bytes):** REPRODUCED (all supported outputs byte-identical) → exit 0;
   DIVERGED (any mismatch) → exit 6, with a **forensic diff** (which output; for text outputs, the
   first differing line / locus).
9. **Mint the certificate** (§8) unless `--no-attest`.

Verdicts/exit codes (Decision 2, **resolved**): `0` REPRODUCED · `6` DIVERGED · `7` INCONCLUSIVE
(couldn't run: pre-schema-5 / input not located / unsupported output) · `5` malformed/tampered
receipt (reused from `verify`). The 6-vs-7 split lets CI and humans distinguish "the result
*changed*" from "I *couldn't check*."

## 7. Honesty scope (the brand)

- **v1 covers the deterministic text outputs Rosalind emits: `variants → VCF` and `features → TSV`**
  — the flagship paths. Byte-equality there is airtight (canonical, deterministic rendering).
- **BAM/bgzf-producing verbs (`align`, `sort`) are out of scope for v1.** bgzf rests on a C zlib
  not captured by `deps_lock_blake3`, so `reproduce` reports them as INCONCLUSIVE — "not
  byte-comparable (bgzf/zlib) in v1" — rather than risking a false DIVERGED. (A future
  decompressed-content comparison is noted in §11.)
- The resource line is always labeled `(here: this machine)`; realized peak is machine-dependent.
- The certificate attests re-derivation, not correctness; its rendered verdict and the docs state
  this in one sentence.

## 8. The reproduction certificate (`.repro.json`)

A new lightweight, content-addressed type (`src/provenance/repro.rs`) — deliberately **not** a
`RunManifest`, to keep the boundary clean:

```
ReproReceipt {
  schema_version,
  parent_claim,                 // the original receipt's content_hash() — the cross-machine ID it chains to
  parent_subcommand,
  verdict,                      // REPRODUCED | DIVERGED
  outputs: [ { role, recorded_blake3, observed_blake3, match } ],
  reproducer_code: { git_sha, dirty, rustc, target_triple, deps_lock_blake3 },  // who/what reproduced
  resource_here: { peak_rss_bytes, declared_budget_mb, fit },   // machine-local, labeled
  chain_depth,                  // parent's depth + 1 (1 if the parent is an original run receipt)
  repro_blake3,                 // self-hash over the canonical form (tamper-evident; Ed25519-signable later)
}
```

- **Chaining / the web:** keyed on the parent's `content_hash` (the stable cross-machine ID). N
  independent reproductions of the same result mint N certificates pointing at one `parent_claim`
  → "+N independent confirmations," emergent with **no server**. A reproduce-of-a-`.repro.json`
  reads the parent's `chain_depth` and increments it.
- **Canonical render + self-hash** mirror `RunManifest`'s discipline (sorted keys, no timestamps,
  `repro_blake3` stamped last over the rest). A `from_canonical_json` parser round-trips it so a
  `.repro.json` is itself reproduce-able / chainable.
- **Write behavior (Decision 1, resolved):** by default write a `<manifest>.repro.json` sidecar
  next to the input `--manifest` (deterministic, non-cwd location — avoids the C3 "don't pollute a
  pipe user's cwd" footgun). `--no-attest` skips writing; `-o <path>` redirects.

## 9. CI fence + self-hosted badge

- **Fence (the substance):** a committed golden chain (a small input set + the receipt + the
  expected VCF) and a CI job that runs `reproduce` on the GitHub runner — a *different machine than
  the author's* — asserting **REPRODUCED**, plus a negative (a one-byte-flipped input → DIVERGED,
  exit 6). This is a regression fence on determinism itself.
- **Badge (cosmetic):** `rosalind badge --manifest <m> [--repro <r>] -o badge.{json,svg}` emits a
  shields.io-*endpoint-format* JSON **and** a static SVG — **self-hosted, no shields.io runtime
  dependency** — e.g. `● reproducible · fits 256 MiB`. The Action publishes it (artifact / committed
  to a badges branch).

## 10. Acceptance gates (spec §-gates; all must pass)

1. **Capture/replay roundtrip:** `CommandCapture` → `params["command"]` → `to_argv` reconstructs
   the original argv (operands substituted); property-tested over flags/opts/inputs/outputs.
2. **Single-source inputs/outputs:** `inputs[]`/`outputs[]` produced by `record_into` match the
   `input()`/`output()` calls; a consistency test asserts no drift between `command` operands and
   `inputs[]`/`outputs[]`.
3. **Real reproduce (flagship):** a `variants --index → VCF` run, then `reproduce` against
   content-located inputs → REPRODUCED, exit 0, and a `.repro.json` that self-hashes and whose
   `parent_claim == original.content_hash()`.
4. **Tamper → DIVERGED:** flip one output byte → DIVERGED, exit 6, forensic first-diff locus
   correct.
5. **Missing input → INCONCLUSIVE:** exit 7, names the absent hash; not DIVERGED.
6. **Back-compat:** a schema-4 receipt still `verify`s; `reproduce` on it → INCONCLUSIVE (no
   `command`), never a panic.
7. **BAM honesty:** a BAM-output verb → INCONCLUSIVE "not byte-comparable (bgzf/zlib)", not a false
   DIVERGED.
8. **Chaining:** `reproduce` of a `.repro.json` increments `chain_depth`; two certificates over the
   same parent share `parent_claim`.
9. **Cross-machine CI fence:** golden chain REPRODUCED on the runner; the one-byte-flip negative
   DIVERGED (exit 6). Determinism unaffected by the schema-5 change (existing `determinism` /
   `golden_vcf` suites stay green; goldens refreshed only if the *claim shape* legitimately changes,
   documented).
10. **Build hygiene:** `cargo fmt --check`, `cargo build` 0 warnings (debug + release), full suite
    green; `verify`/`--expect-code` unaffected by the schema bump.

## 11. Decisions resolved & out-of-scope/future

**Resolved:** Decision 1 → certificate default-writes a `<manifest>.repro.json` sidecar
(`--no-attest`/`-o`). Decision 2 → exit codes `0/6/7/5`.

**Future (not this track):** Ed25519-signing the `.repro.json` (Track A / P3.3) → counter-*signed*
reproduction; decompressed-content comparison for BAM; chain-aware `reproduce --chain` composing
with Track C; a public reproduction index/web (emergent from many certificates — no server built
here).

## 12. Risks

- **Schema-evolution collision:** Track D owns the 4→5 bump; Tracks E/I (if ever built) must rebase
  onto it. Track H (ColumnKit) must land *after* schema-5 freezes so `verify`/`--expect-code` keep
  passing.
- **Re-exec environment:** `reproduce` re-runs `current_exe()`; a different binary build legitimately
  may diverge. Mitigated by surfacing the build-identity comparison prominently and by the verdict
  being about output bytes with code-mismatch noted (not silently passed).
- **Over-claim guardrails:** the certificate must render an explicit "re-derivation, not
  correctness" line; the resource line must read "(here: this machine)"; the badge must not imply
  tamper-*proof* (it is tamper-*evident* until signing lands). Honesty is the brand — these are
  acceptance-level, not cosmetic.
