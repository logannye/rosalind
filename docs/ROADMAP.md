# Rosalind — Engineering Roadmap

**Status:** Living engineering-direction document — 2026-06-02. Companion to
[`docs/OPEN_PROBLEMS.md`](OPEN_PROBLEMS.md) (the research thesis) and [`CONTRACT.md`](../CONTRACT.md)
(the shipped memory contract). Audience: builders deciding what to build *on* Rosalind, and
contributors deciding where to help.

This document is a *direction*, not a promise of dates. It records **what we will build next and
why**, in priority order, grounded in the actual codebase.

---

## 1. The organizing principle

Rosalind has exactly one moat, and it is not the caller, the aligner, or the `~√t` framing. It is
the **receipt**: a bounded, deterministic run that emits a content-addressed proof of *what it did*
and *what it cost*.

> **The receipt is the product. The caller is the demo.**

Every prioritization decision flows from a single test. A change is worth doing if and only if it
does one of three things to the receipt:

1. **Hardens it** — makes the receipt *unforgeable* and the memory bound it attests *unbreakable*.
2. **Widens it** — makes more lifecycle stages and more artifact types emit a receipt, extending a
   verifiable provenance graph that a non-deterministic tool structurally cannot produce.
3. **Deepens it** — pairs every receipt with a *measured trust number* (accuracy / correctness), so
   an artifact is reproducible **and** right.

Work that does not touch the receipt — raw throughput, structural-variant/CNV tracks, accuracy-parity
races against GATK/DeepVariant, a mismatch-tolerant aligner — is deliberately deprioritized. Those
compete on dimensions Rosalind has declared **non-goals** (see §6), and, more importantly, they do
not *compound*. The receipt does.

---

## 2. Current state (what's real, what's not)

**Shipped and well-tested.** Bounded whole-genome germline *SNV* calling over a sorted BAM (peak
independent of BAM size); a build-once, memory-mapped, byte-reproducible FM-index (`index` / `locate`);
deterministic external-merge `sort` (output provably independent of the chunk boundary); banded gVCF;
region-bounded somatic; the reproducible per-locus feature substrate with BLAKE3 receipts; the
ColumnKit SDK (`ColumnAnalyzer` + `run_bounded_whole_genome`, with `features` as its first impl);
`plan` / `--enforce` / `verify` / `pack`; a self-testing GitHub Action and cross-platform release
tarballs. The historical reference-decode undercount is fixed (`decode_window_arc`), and the
`predicted ≥ realized` working-set inequality is tested.

**Stubbed or simulated — stated plainly.** Accuracy is 100% simulated (no GIAB run yet); the germline
model has none of the standard FP filters (FS/SOR/MQRankSum/ReadPosRankSum). MAPQ is a hardcoded
placeholder the likelihood ignores, so Rosalind's *own* aligner cannot yet be in an accuracy story.
Germline is SNV-only (indels are structurally absent from `Obs`, not merely uncalled). The feature
substrate egresses as TSV through a fully-materializing Python loader — the bounded story currently
breaks at the Python boundary. Published to zero package registries.

**The one structural gap.** Index **build** is `O(reference)` RAM — a single in-RAM SA-IS call peaking
at ~45 B/base. "Declare your RAM" is honored for *calling* and silently false for the *build*, which
is exactly where builders hit the wall first. Closing this is the headline research bet (§5, Sprint 4).

**A quieter, load-bearing gap.** The entire "verifiable upper bound" currently rests on a hardcoded
`PILEUP_IO_RSS_OVERHEAD = 8 MiB` slack constant, with no proof it holds across platforms/allocators,
and **no runtime governor** — so a misprediction could be OOM-killed by the kernel *before* the
post-run exit-4 check fires. That is the exact silent-OOM the thesis forbids. Closing it is the single
highest-ROI item in the repo (Sprint 1).

---

## 3. The asymmetric bets (ranked by ROI, not effort)

| # | Bet | What it does to the receipt | Effort | Why it's asymmetric |
|---|-----|------------------------------|--------|---------------------|
| 1 | Unbreakable bound | Hardens | S–M | Days of work insure the entire thesis against a brand-inversion event. |
| 2 | Unforgeable receipt | Hardens | S | "Reproducible" → "tamper-evident" = the actual clinical/audit bar; schema-version is free now, costly later. |
| 3 | One real GIAB number | Deepens | M | You don't need to *win* GIAB — a number *paired with a memory receipt* is an artifact no one else publishes. |
| 4 | Content-addressed dataset engine | Widens | M | The substrate already emits the bytes; determinism makes "bit-reproducible training data" a category incumbents are locked out of. |
| 5 | Somatic + sort into the envelope | Widens | M | Lifecycle width on existing primitives; closes the "your thesis holds for 2 of 7 stages" gap. |
| 6 | D1a: the build that refuses to OOM | Widens (the last stage) | XL | Category-defining and the literal realization of the √t framing — but high cost, so it runs in parallel, not first. |

---

## 4. The roadmap (sequenced)

Four increments. The **spine** (Sprints 1–3) is certain value that compounds the receipt; the
**headline bet** (D1a) runs in parallel once its cheap prereqs land. Each to-do is shippable as its
own PR via the project's `brainstorm → spec → plan → TDD` flow. Effort: **S** = days, **M** = 1–3
weeks, **L** = a quarter, **XL** = 1–2 quarters.

### Sprint 1 — Harden the receipt (insurance + hygiene)

- **1.1 — Unbreakable bound.** Replace the fixed `PILEUP_IO_RSS_OVERHEAD` with a **measured residual**
  (sample `peak_rss − working_set − baseline`, record p99 in the receipt). Add a `setrlimit(RLIMIT_AS,
  budget)` guard under `--enforce` (Linux) plus a periodic RSS-poll abort fallback (portable), so a
  breach **fails loud at the allocation site (exit 4)** before the kernel OOM-killer can fire. Add an
  upper-bound test at **near-saturation on a chr20-scale fixture**, not the 4 MiB toy.
  *Surfaces:* `core/budget.rs`, `pileup/engine.rs`, `main.rs`, `tests/plan_enforce.rs`. *Effort: S–M.*
- **1.2 — Unforgeable receipt + schema version.** BLAKE3 self-hash over the canonical JSON (digest
  field zeroed); `verify` re-derives it. Remove the verify-budget-from-manifest fallback. Stamp a
  `SCHEMA_VERSION` into `FEATURE_HEADER` + the manifest **now**, before any columnar schema freezes.
  *Surfaces:* `provenance/mod.rs`, `call/features.rs`, `main.rs`. *Effort: S.*
- **1.3 — Forkability hygiene.** `cargo publish`; clippy `-D warnings` + MSRV + `rust-version`; fix the
  known doc-drift defects; run the README quickstart verbatim in CI; add `CONTRIBUTING.md`,
  `ARCHITECTURE.md`, and a deterministic-repro issue template. (Not PyPI/bioconda — the htslib linking
  wall is real work, deferred.) *Effort: S.*

### Sprint 2 — Earn caller trust (∥ de-risk the build bet)

- **2.1 — Genotype-aware GIAB.** Upgrade the comparator: GT/zygosity scoring, multi-allelic
  decomposition (stop skipping comma-containing ALTs), golden left-align cases. Run HG002 chr20
  high-confidence **+ CMRG** through `eval-germline` on an externally-mapped BAM; publish an accuracy
  card *alongside* the memory receipt. Fold MAPQ into the germline likelihood. Frame honestly:
  *"accuracy is now measured, not assumed; the contract is the contribution."*
  *Surfaces:* `genomics/eval/{vcf,compare,normalize}.rs`, `call/germline.rs`. *Effort: M.*
- **2.2 — D1a prereqs (start early; cheap; independent).** Property/adversarial-fuzz `sais_u32` vs a
  naive SA. Swap the coarse 12 B/base build estimate for the code-grounded ~45 B/base
  `BuildMemoryModel` as the `plan --reference` basis, so the build inequality is sound *before*
  enforcement. *Surfaces:* `genomics/suffix_array.rs`, `genomics/index/report.rs`, `call/plan.rs`.
  *Effort: S–M.*

### Sprint 3 — Widen the contract + open the ML surface

- **3.1 — Somatic + sort into the envelope.** `stream_somatic_whole_genome`: two `StreamingBamSource`
  + two depth-capped engines in lockstep, reusing the existing merge-join with **lazy** engines; wire
  `plan`/`--enforce`/`verify`. Give `sort` an enforced (refuse / fail-loud) chunk bound + a build
  receipt. *Surfaces:* `call/{somatic,pipeline}.rs`, `genomics/sort.rs`. *Effort: M.*
- **3.2 — Arrow/Parquet egress + zero-copy bindings.** An `ArrowFeatureAnalyzer` (`ColumnAnalyzer`)
  that flushes bounded 64k-locus Parquet row-groups behind a cargo feature (keep the lean static-musl
  default); add the row-group term to the plan estimator so `predicted ≥ realized` still holds.
  Replace the Python `str.split` loader with `pyarrow → .to_torch()/.to_jax()/__dlpack__`.
  *Surfaces:* `call/{features,columnkit,plan}.rs`, `python/rosalind.py`. *Effort: M.*

### Sprint 4 — The dataset engine + the headline build bet

- **4.1 — Supervised in one command + content-addressed dataset cards.** A `LabeledFeatureAnalyzer`:
  `features --truth truth.vcf --confident-regions conf.bed --emit-label` → bounded `(X, y)`.
  Region-disjoint, content-addressed train/val/test split descriptors; `rosalind dataset-card`;
  `verify` re-checks split-disjointness. The "prove this quarter's model saw exactly last quarter's
  data" artifact. *Surfaces:* `call/features.rs`, `genomics/eval/{normalize,bed}.rs`,
  `provenance/mod.rs`. *Effort: M.*
- **4.2 — D1a: budget-tunable external-memory build (the headline).** Blocked SA-IS (budget-sized
  blocks → spill → deterministic disk merge reusing `sort.rs`'s heap-merge); wire
  `index --memory-budget-mb` into refuse/honor/verify + a build receipt; byte-identical `.idx` oracle
  on RAM-sized genomes, plus an independent FM-roundtrip / `locate`-vs-brute-force check for genomes
  too large to build at `b = n`. Market as *"build a human index on a 2 GB box, or refuse."* Pangenome
  scale (u64 / format-v2, past the current `u32::MAX` ceiling) is **D1b**, a separate deferred epic.
  *Surfaces:* `genomics/{suffix_array,fm_index,genome_index}.rs`, `genomics/index/io.rs`,
  `genomics/sort.rs`. *Effort: XL.*

---

## 5. The single big bet

**D1a — close the build hole, scoped to "human/T2T on a small box," run in parallel with the spine.**

It is the only item that resolves the thesis's defining contradiction: today "declare your RAM" is
honored for calling and silently false for the build — the most memory-hungry stage, and where
builders hit the wall first. No production indexer (bwa, minimap2, samtools) offers a declared-budget,
degrade-to-disk build with a verifiable receipt; they OOM-kill on an emergent peak. Rosalind already
owns every surrounding piece (SA-IS over an integer alphabet, the deterministic external merge in
`sort.rs`, the streamable on-disk block format, and the D0 measurement proving the build is
intermediate-state-bound). It comes with a cheap, powerful oracle — a blocked build must produce a
**byte-identical `.idx`** to the full-RAM build.

Three honest constraints govern it: (1) market it as **human-on-a-laptop, not pangenome** — the
`u32::MAX` ceiling means wheat/conifers/pangenomes cannot even be *represented* until the separate
u64 epic (D1b); (2) the prereqs (SA-IS fuzz, estimator swap) are **gating, not concurrent**; (3) it is
the long pole (XL) and the disk-backed block-boundary SA merge is the genuinely hard part — so it runs
alongside the cheap reach-and-credibility spine, which touches disjoint subsystems.

---

## 6. Non-goals (what we deliberately do not build, and why)

- **No "time/wall-clock curve" as a contract dimension.** A throughput regression fit from past
  receipts has no soundness property, and it is not the √t framing (which is about *space*). Keep
  wall-time as a logged telemetry field only.
- **No SV/CNV evidence track yet.** A read-depth bedgraph + discordant/split-read counts is commodity
  (mosdepth/samtools); "bounded" adds little where depth tracking is already trivially streaming.
- **No full local-assembly haplotype indel calling.** It breaks the bounded contract and is a
  multi-year fight against entrenched, well-validated tools. Bounded-*pileup* indels may ship later as
  a **capability** demonstration, explicitly not accuracy parity, and only after GIAB.
- **No `--threads` as a performance play.** Byte-identical-across-thread-count is a hard *correctness*
  property (SA-merge tie-break, allocator-arena RSS interactions); if pursued, only ever as a
  determinism property with a byte-identity CI gate, and only after the spine and the build bet land.
- **Do not anchor the pitch on `pack` co-location.** It is real but thin — schedulers already pack by
  *requested* memory. It rides on the contract; it is not the headline.

---

## 7. How to build on Rosalind

The fastest path from "I want my own per-locus metric" to "I inherit a machine-checkable memory budget
+ a hash-verifiable receipt" is the **ColumnKit SDK**: implement one `ColumnAnalyzer` method and run it
through `run_bounded_whole_genome`. See [`CONTRACT.md`](../CONTRACT.md) and
[`examples/columnkit_coverage.rs`](../examples/columnkit_coverage.rs). The shipped `features` egress is
itself the first `ColumnAnalyzer`, so a third-party analyzer is a first-class citizen, not a
second-class plugin.
