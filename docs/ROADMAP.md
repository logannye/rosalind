# Rosalind — Implementation Roadmap

**Status:** Living engineering-direction document — revised 2026-06-02 (v2). Supersedes the earlier
receipt-first draft with a sharper, honestly-scoped synthesis. Companion to
[`docs/OPEN_PROBLEMS.md`](OPEN_PROBLEMS.md) (the research thesis) and [`CONTRACT.md`](../CONTRACT.md)
(the shipped contract). Audience: builders deciding what to build *on* Rosalind, and contributors
deciding where to help.

This is a *direction*, not a promise of dates. Every claim is scoped to where it is **true**, because
the honesty about where Rosalind does *not* help is the brand.

---

## 1. The organizing principle

Rosalind has one moat, and it is not the caller, the aligner, or a new complexity bound.

> **Memory is a verifiable contract. The receipt is the product; the caller is the demo.**

You declare the RAM you have. Rosalind predicts — *before a byte of compute runs* — whether the job
fits; honors that ceiling (refuse cleanly up front, or fail **loud** the instant a breach is detected,
never a silent OOM-kill); produces a **byte-identical** result; and hands you a **tamper-evident,
content-addressed receipt** you can verify offline, months later, without re-running.

Today that contract is true for **calling, query, and the feature stream**. The one stage where it is
still false is **index construction** (`O(reference)` RAM). **The √t insight is the keystone that makes
the contract whole** — it turns construction into a point on a continuous space/time *curve* a declared
budget selects, so "declare your RAM and it is honored" becomes true *end-to-end*. √t is not a separate
feature to lead with; it is what completes the one thing that is uniquely Rosalind's.

The theory (Williams 2025: any time-*t* computation simulable in `O(√(t log t))` space, buying space by
spending time) is the *justification* for "recompute along a √-shaped curve." It is **not** the runtime
kernel, and **Rosalind makes no new space-complexity claim.** The defensible contribution is the
engineering wrapper: a **declared-budget contract + a measured curve + a verifiable receipt** over a
real memory-bound build.

---

## 2. What is genuinely differentiated — and who feels the pain

Each is tagged **[shipped]** or **[to build]** and tied to the persona who actually needs it.

| Capability | Persona who needs it | Why it's hard for incumbents to match |
|---|---|---|
| **Never a silent OOM — predict, then refuse / fail loud** [shipped] | Field/outbreak operator on a no-swap laptop driving a MinION | GATK `-Xmx` crashes; DeepVariant OOM-kills; pSAscan/Big-BWT/ropebwt3 have no governor. `plan` gives a FITS/REFUSE verdict in ms; the runtime governor (`core/governor.rs`) fails loud at exit 4. |
| **Verify-without-rerun** [shipped] | Air-gapped / regulated / CRO auditor | A non-deterministic caller can't emit a byte-reproducible hash, so it can't anchor an audit at all. `rosalind verify` re-checks inputs/outputs + realized-peak-vs-budget offline in seconds. |
| **Input-size-independent calling peak** [shipped] | Anyone calling on fixed-RAM nodes | Peak ≈ largest-contig reference + a depth-capped active set, *not* BAM size. Emergent-RAM callers can't bound this. |
| **Inherit-the-contract SDK (ColumnKit)** [shipped] | A builder shipping a per-locus analytic | Implement one `ColumnAnalyzer` method → inherit the bounded walk + governor + receipt. No genomics library lets a third party inherit a verifiable memory budget. |
| **The dial for construction — build a genome index on a box that would OOM** [to build, D1a] | T2T/human-index builder on a 16–32 GB workstation, not a fat node | bwa-index/samtools build in *emergent* RAM (you find the ceiling by crashing). External-memory tools exist but at *fixed* points with no declared budget, no feasibility prediction, no degrade-contract, no receipt. |

**Table-stakes — real value, but NOT the moat (say so):**
- **A GIAB F1 number.** Every caller posts one; ours will likely *trail* DeepVariant (SNV-only, no FP
  filters). The contribution is *accuracy paired with a verifiable receipt*, not a leaderboard win.
- **Arrow/Parquet egress.** Universal. The only real claim — a cross-machine-stable hash over the
  shard — is *inherited* from determinism. Sell it as the **on-ramp**, never the moat.
- **Methylation / mod-base / indel-evidence tracks.** modkit/MethylDackel are mature. "Bounded +
  receipted" widens the receipt's *reach* (good for the SDK story); it does not widen the *moat*.

---

## 3. Honest boundaries — stated up front

The decisive test for any workload: **"is the memory bottleneck the computation's working STATE, or its
INPUTS?"** If inputs, the √t dial does nothing.

- **Cohort / joint calling is INPUT-bound** — holding *N* samples is `O(N)` input. Beaten by streaming,
  not √t. Out of scope for the dial.
- **Pangenome-*graph* analysis is INPUT-bound** — the graph *is* the memory. (Indexing a *linearized*
  pangenome reference is in-scope construction; traversing the graph is not.)
- **The `u32` representation ceiling (`MAX_GENOME_LEN ≈ 4.29 Gbp`, `genome_index.rs:20`)** — wheat
  (16 Gbp), conifers, and large pangenomes cannot even be *represented* until a separate `u64`/format-v2
  epic (**D1b**). **Today's credible build market is T2T-human-and-below — not wheat, not pangenome.**
- **Single-threaded** (only the governor spawns a thread) — so a budget-tunable build's time overhead is
  *on top of* a single-threaded baseline; vs *parallel* pSAscan the realistic gap is large. Time is
  telemetry, not a contract — you cannot quote a guaranteed wall-clock up front.
- **Accuracy is 100% simulated today** — the caller is untested on real GIAB. Until that is fixed, the
  regulated/clinical personas (the ones with budget) are unreachable. This is a **precondition, not
  polish** (Phase 2).

---

## 4. Who we build for (ranked by reachability, not by how good the pitch sounds)

1. **Your own stack** — virtual-cell / Tahoe-100M / causal-world-model work. A deterministic, bounded,
   content-addressed feature + eval + provenance backbone is *directly* load-bearing there, with no
   procurement gate. **Build for your own stack first** — it is the soundest adoption thesis, a real
   internal user who needs exactly this.
2. **Reproducibility-conscious researchers** (the fork-and-build-on audience, ~249 stars / 10+ forks).
   The tamper-evident demo + a WASM receipt verifier are top-of-funnel; `reproduce` is the conversion.
3. **Field / edge / instrument** (no-swap, fixed-RAM, never-OOM) — small but genuinely served *today*.
4. **Regulated / biobank — aspirational.** Build the *capability* (signing, chain, cohort) so it is
   ready when there is an *organization* behind the signature; do not build the economics on a buyer a
   solo, pre-1.0 repo cannot transact with.

---

## 5. Done so far (the foundation)

- **Phase A–C + hardening (shipped):** bounded whole-genome germline SNV calling; the memory contract
  (`plan`/`--enforce`/`verify`) with a runtime **governor** (loud exit-4 on breach) + RSS-residual
  telemetry; build-once mmap FM-index; gVCF; region-bounded somatic; deterministic sort; `pack` fleet
  scheduler; ColumnKit SDK; the reproducible feature substrate (TSV); install.sh + the budget Action.
- **Sprint 1 (this session):** the runtime governor (enforced contract); the **tamper-evident
  self-hashing receipt**; forkability hygiene; **MSRV 1.83 + a clippy `-D warnings` CI gate**.
- **Sprint 2.1 (this session):** **genotype-aware GIAB-grade eval** — `eval-germline` now scores
  genotype concordance + decomposes multi-allelic records (synthetic: 40× clean concordance 1.00).

---

## 6. The sequenced roadmap

The through-line: **make the proofs *true* → sell what ships → earn the accuracy gate → extend the
proofs → complete the contract with the √t build.** Effort: **S** = days, **M** = 1–3 weeks, **L** = a
quarter, **XL** = 1–2 quarters.

### Phase 0 — Make the proofs *true* (the soundness foundation; do first)

Everything downstream rests on the proofs being literally correct. Two verified defects + one gap:

- **P0.1 — Sound build-feasibility prediction.** *(S, verified)* `plan --reference` estimates the build
  at a stale **12 B/base** (`genomics/index/report.rs:19`) while D0 *measured* **41 B/base** and the
  code-grounded `BuildMemoryModel` (~45 B/base) sits unused beside it. Today the prediction
  **under-estimates ~3.4×** — a `plan` that says FITS then OOMs, which destroys the one thing the
  contract sells. Wire `BuildMemoryModel` into `plan --reference` so *predicted ≥ realized* is sound
  before any build enforcement ships.
- **P0.2 — Split the receipt into a deterministic CLAIM and a machine-dependent MEASUREMENT.** *(S)*
  `content_hash()` (`provenance/mod.rs:109`) excludes only `manifest_blake3`, so the machine-dependent
  `peak_rss_bytes` / `max_working_set_bytes` are *inside* the hash — the same correct run on two
  machines yields two different receipt hashes. Tamper-evidence still works, but the hash is not a
  cross-machine content-address, which silently breaks chaining/reproduce/cohort. Partition the manifest:
  hash/chain/sign/re-derive the **claim** (inputs, deterministic params, code-identity, output hashes);
  record + budget-check the **measurement** but exclude it from `content_hash()`. Bump
  `MANIFEST_SCHEMA_VERSION` to 2; degrade pre-v2 gracefully (the pre-1.2 self-hash skip is the pattern).
  **DONE** (claim/measurement split shipped; the split is lossless — a second `measurement_blake3` keeps
  the measured cost locally tamper-evident, and a hash-protected `has_measurements` claim marker catches a
  stripped measurement block).
- **P0.2b — Path-normalized / content-only claim.** *(S)* **DONE.** P0.2 removed the measured *cost* from
  the claim hash, but recorded input/output **paths** (`FileHash.path`) were still in the claim, recorded
  verbatim — so two machines with byte-identical data at different paths still hashed differently. The
  schema-3 claim now keys inputs/outputs on their sorted `blake3` digests (paths dropped from the *claim*
  form; the on-disk receipt keeps full paths for humans and for `verify` to re-hash files), making the
  claim hash a true cross-machine content-address. Version-gated (schema <3 reproduces the path-inclusive
  form so pre-P0.2b receipts still self-verify). This is the property chaining/reproduce/cohort need.
- **P0.3 — Build-identity in the receipt.** *(S)* **DONE.** Replaced the `tool_version = "0.1.0"`
  half-measure with `code_git_sha` + `code_dirty` + `rustc_version` + `target_triple` + `deps_lock_blake3`
  (baked at compile time by `build.rs`, each degrading to `"unknown"`), stamped into the *claim* (so they
  are hash-protected and form the reproduction key), plus a `verify --expect-code <sha>` gate (prefix
  match; flags a clean match from a dirty build; rejects a degenerate SHA). Schema 3 → 4. "Exactly this
  code" is no longer a lie. **⇒ Phase 0 complete: the receipt is sound (conservative build estimate),
  cross-machine stable (content-address claim + lossless measurement split), and code-identified.**

### Phase 1 — Sell what ships (harden + market the contract; no √t needed)

The strongest *differentiated* claims are real **today**. Make them legible and viral.

- **P1.1 — Harden + document the contract as the product.** *(S)* Governor fail-loud, offline
  `verify`-without-rerun, `pack` as a k8s/Slurm scheduling input, ColumnKit promoted to a documented SDK
  with 2–3 reference analyzers. Sell what ships.
- **P1.2 — The tamper-evident "caught-you" demo + a WASM receipt verifier.** *(S–M, the front door)*
  An asciinema: `variants` → `verify` OK → a human edits one byte → `verify` exits 5. Then compile
  `provenance/mod.rs` (std + hand-rolled canonical-JSON + blake3, zero htslib) to `wasm32` for a
  drag-a-receipt-in-the-browser verifier — the most viral artifact in the repo and a structural
  capability we're not aware of in incumbent callers. **Honest copy:** catches corruption/casual edits, *not* a motivated
  forger who re-runs `finalize()` — that needs the signature (Phase 3).
- **P1.3 — A scope-boundary benchmark + a README scope table.** *(S, the anti-hype guardrail)* A
  reproducible benchmark showing **both** the in-scope win (construction peak slides down the curve)
  **and** the out-of-scope null (cohort/pangenome memory *unmoved* by the dial). Publish the table:
  *FITS the dial = intermediate-state-bound SA/FM/BWT build ≤ 4.3 Gbp; OUTSIDE = > 4.3 Gbp (until u64),
  pangenome-graph analysis ever, cohort/joint holding.* Turns the honest boundary into a citable
  credibility artifact — a pangenome team adopts Rosalind for the *build* and keeps their graph tooling.

### Phase 2 — Earn the accuracy gate (the precondition for the big markets)

- **P2.1 — One real GIAB HG002 (chr20 + CMRG) number, paired with the receipt.** *(M, the gate)* The
  genotype-aware comparator is ready (Sprint 2.1); run it on real GIAB truth, on an externally-mapped
  BAM, and attach precision/recall/F1 + genotype concordance to the flagship receipt. Frame honestly:
  *"accuracy is now measured, not assumed; the contract is the contribution"* — the realistic first
  result is SNV-competitive, indels trailing. (Blocked in this environment only by missing samtools/
  aligner tooling; runs anywhere they exist.)
- **P2.2 — Fold MAPQ into the germline likelihood.** *(M)* Today MAPQ is a hardcoded placeholder the
  likelihood ignores; weight each read by `P(mismapped) = 10^(-MAPQ/10)`. Behavior-changing → its own
  before/after on synthetic noise; pairs naturally with the real (externally-mapped) GIAB BAM.

### Phase 3 — Extend the proofs (the provenance frontier)

- **P3.1 — Close the provenance CHAIN; `verify --chain`.** *(M, keystone)* `index` writes a `RunManifest`
  sidecar (reuse `provenance/mod.rs`; the header already has `reference_blake3`); downstream receipts pin
  the index by **artifact content hash** (stable across rebuilds), with the upstream receipt hash as a
  secondary attestation. `verify --chain` walks the DAG and checks every node re-hashes + every edge
  resolves — "this VCF came from exactly this index built from exactly this reference," offline.
  **Beats:** nf-core/WDL provenance (path/timestamp text over command lines, not content-addressed bytes
  — and they can't upgrade because the callers underneath are non-deterministic).
- **P3.2 — `rosalind reproduce`: third-party byte re-derivation, CI-fenced.** *(M, the headline proof)*
  Re-run the recorded subcommand against content-located inputs; assert fresh BLAKE3 == recorded output
  hash → REPRODUCED / DIVERGED with a per-field diff. Pair with a frozen golden chain in CI = a
  regression fence on determinism itself. **Why ONLY Rosalind:** GATK/DeepVariant/bwa would report
  DIVERGED on a *correct* run. *Needs P0.2 (cross-machine stable receipt) first.* Add a determinism
  conformance test pinning `rust-htslib` to a deterministic-by-construction mode (htslib has
  multi-threaded bgzf — the I/O layer is C you don't control).
- **P3.3 — `verify-attest`: Ed25519-signed receipts + signed chain root.** *(M, niche)* Sign the
  (now-deterministic) `content_hash()` into a detached `.sig`; `verify-attest --pubkey` checks it.
  Optional, default-features-off. *Build the capability ready; don't bet the economics on the regulated
  buyer.* Key custody is explicitly out of scope.
- **P3.4 — `rosalind cohort`: a verifiable Merkle ledger for biobank-scale re-analysis.** *(L)* Roll
  *N* byte-stable per-sample receipts into a signed root; `cohort verify` re-derives without re-running
  and, on mismatch, reports the **exact changed leaves and which field drifted**; `cohort plan` reuses
  `pack` to show the re-run co-locates within budget. **Beats** GenomicsDB/GLnexus (outputs not byte-stable, no
  receipt, emergent memory). Depends on P0.2 + P3.1 + P3.2.

### Phase 4 — Complete the contract end-to-end (the √t build)

- **P4.1 — D1a: the budget-tunable external-memory blocked SA/BWT build.** *(XL, the second-act
  flagship)* Blocked SA-IS — `sais_u32` each budget-sized block in RAM, spill, deterministic disk-backed
  merge — wired into `index --memory-budget-mb` with refuse/honor/verify + a build receipt. Correctness
  oracle = byte-identical `.idx` vs the full-RAM build on RAM-sized genomes; FM-roundtrip / locate-vs-
  brute-force for genomes too big to build at `b = n`. CI gate: byte-identical `.idx` across ≥3 budget
  points. **Prereqs (gating): P0.1 (sound estimator) + SA-IS adversarial fuzz.**
  > **KILL / REDESIGN CRITERION (hard):** the disk SA merge needs **cross-block suffix comparison**
  > (eSAIS/pSAscan-class) — a *harder* algorithm than the `sort.rs` record-merge; keeping `sais_u32`
  > per-block does **not** de-risk it. If a 4 GB-budget human build exceeds ~5–10× in-RAM wall-clock,
  > **stop and redesign the merge.** Weigh single-threadedness: vs *parallel* pSAscan the gap could be
  > 10–20× ("overnight" → "over a weekend") — **reconsider parallelism as a prerequisite, not a
  > fast-follow.** Run this as a research spike with the kill criterion from day one, *parallel* to the
  > certain-value phases, never on their critical path.
- **P4.2 — The killer demo (once D1a lands).** Build a T2T-CHM13 human FM-index (3.1 Gbp) on a **32 GB
  box where `bwa-index` OOM-kills**; `plan --reference` prints FITS + predicted peak; the blocked build
  spills and completes; the `.idx` is **BLAKE3-byte-identical to the 127 GB full-RAM build** — proven by
  `verify` on a laptop. The byte-identical oracle is the cheap, lethal proof the curve changed *nothing*
  but the RAM.
- **D1b — `u64`/format-v2 for pangenome scale.** *(XL, deferred)* Breaks the 4.29 Gbp ceiling (touches
  the SA element type, sampled SA, on-disk format, every byte-reproducibility golden). Only after D1a
  proves out; inherits all its merge risk.

---

## 7. The flagship demos

1. **Now:** `scripts/flagship_ecoli_demo.sh` — `plan → enforce → verify` on the real *E. coli* genome
   (256 MiB FITS / 22 MiB realized / 8 MiB REFUSED at exit 3). The contract, undeniable, today.
2. **Phase 1:** the tamper-evident "caught-you" + the browser WASM receipt verifier.
3. **Phase 4:** the T2T-on-a-32 GB-box byte-identical build.

Lead live with #1; never invert the order.

---

## 8. What to explicitly NOT do

- **Don't out-call or out-align anyone.** Not DeepVariant on F1, not bwa-mem2/minimap2 on speed. The
  thesis is the substrate, not the leaderboard.
- **Don't lead with the √t build, Parquet egress, or a GIAB number as the headline.** The build is the
  *keystone that completes* the contract, not the front door; Parquet is the on-ramp; a likely-trailing
  GIAB number sold as a win damages the brand more than no number.
- **Don't position methylation / indel / SV / CNV tracks as moat.** They widen reach, not defensibility.
- **Don't chase the input-bound frontiers** (cohort holding, pangenome-graph analysis) *with √t* — they
  need streaming + succinct inputs, orthogonal techniques.
- **Don't add `--threads` as a perf play.** Byte-identical-across-thread-count is a hard *correctness*
  property that forecloses with the speed track; the one exception is evaluating parallelism as a
  *prerequisite* for D1a's build (see the kill criterion), with a byte-identity CI gate as the
  non-negotiable firewall.
- **Don't build the roadmap's economics on the regulated buyer** a solo, pre-1.0 repo cannot yet serve.

---

## 9. The meta-question, answered

**Who actually switches, and the bus factor.** The strongest-*sounding* pitches (clinical chain-of-
custody) are the *least* reachable for a single-maintainer, pre-1.0 tool. The roadmap answers this by
sequencing to the audience that genuinely needs byte-identity and has no procurement gate — **your own
stack first**, then reproducibility researchers, then field/edge — by **fixing the receipt (P0.2) and
the predictor (P0.1) so the proofs are true rather than aspirational**, and by making the lowest-friction
proofs (the tamper demo, the WASM verifier) the front door. The regulated capability is built *ready*,
not bet *upon*.

**The brand is the honesty.** This roadmap names its own me-too risks, scopes claims tightly (positional
not homology leakage; tamper-evident not tamper-proof; a framing not a complexity result), and the
codebase backs the "reuses shipped X" claims under audit. That honesty is rarer than any feature here —
and the one-sentence promise it earns is the whole moat:

> *Rosalind is the only genomics engine where "reproduce this exact result" is a command a stranger can
> run — because determinism, a content-addressed receipt, and a declared memory bound are structural
> properties of every run, not features bolted onto one.*

---

## 10. Build on Rosalind

The fastest path from "I want my own per-locus metric" to "I inherit a machine-checkable memory budget +
a hash-verifiable receipt" is the **ColumnKit SDK** — implement one `ColumnAnalyzer` method and run it
through `run_bounded_whole_genome`. See [`CONTRACT.md`](../CONTRACT.md) and
[`examples/columnkit_coverage.rs`](../examples/columnkit_coverage.rs). The shipped `features` egress is
itself the first `ColumnAnalyzer`, so a third-party analyzer is a first-class citizen.
