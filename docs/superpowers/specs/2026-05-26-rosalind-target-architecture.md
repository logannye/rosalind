# Rosalind — Target Architecture & Rebuild Plan

**Status:** Approved direction — 2026-05-26. Living document.

**Context.** Rosalind is a Rust genomics engine (alignment + germline/somatic variant
calling) pitched for *deterministic, bounded-memory, on-prem/edge* genomics, built over a
square-root-space evaluation framework (inspired by Williams 2025 *Simulating Time With
Square-Root Space* and Cook–Mertz 2024 tree evaluation). A structured audit (2026-05) found
genuine, well-built foundations but a large gap between positioning and implementation, plus
several bugs that make outputs **silently wrong on real data**. The repo is gaining forks.
This document defines the architecture we are rebuilding *toward*. Effort is not the
constraint; the quality, uniqueness, and legibility of the engine for the people who build on
it is the sole objective.

---

## 1. North star — the unique bet

**Rosalind is the bounded-memory, deterministic, streaming genomics *kernel* that edge builders
compute on.** Not another monolithic CLI that races GATK on accuracy, but the thing that
exists where the ecosystem has a hole:

- a **library-first** engine whose core primitive — a streaming, CIGAR-aware **pileup column
  stream** — is a public, embeddable substrate;
- where **memory is a contract you declare**, not an outcome you hope for (decisive on edge
  devices with fixed RAM and no swap, where OOM means a dead run);
- where **reproducibility is a verifiable artifact** a third party can check without re-running;
- where **confidence is calibrated and honest** — the engine abstains rather than overcalls
  when evidence is thin (a confident wrong call is the worst outcome in the field, where you
  can't re-run or get a second opinion);
- over which anyone can compute **arbitrary per-locus analytics — including on-device ML
  features — across a whole genome on a laptop, deterministically, within a fixed budget.**

We win by being the trustworthy, embeddable, predictable kernel for genomics at the edge — for
hospital laptops, portable sequencers (MinION in the field), classrooms, and researchers who
want to *build on* a genomics engine rather than shell out to one. Variant calling is the first
consumer of the kernel, not the whole product.

Non-goals (on purpose, so the promises are airtight): out-accuracy-ing GATK/DeepVariant on
benchmarks; cloud/cluster-scale orchestration; supporting every exotic format.

---

## 2. The five promises (each becomes CI-enforced, not asserted)

These are the product. Today most are merely claimed; the rebuild makes each an enforced
invariant.

| Promise | What it means | Enforcement |
|---|---|---|
| **Deterministic** | Byte-identical artifacts for identical inputs+config | Repeat-run byte-equality; shuffled-input → identical calls; (Phase D) `--threads 1` vs `N` equality |
| **Memory is a contract** | You declare a budget; the engine honors it or refuses *up front* with a computed requirement. `rosalind plan` reports the working-set envelope before running | Real-RSS gate while running the CLI on *growing* inputs: assert flat working set + that declared budgets are honored (replaces today's tautological counter). Honest caveat: index *build* is O(reference), documented separately |
| **Standards-compliant** | BAM/VCF the ecosystem accepts | `samtools quickcheck` / `bcftools view` / `bcftools norm` validate emitted artifacts in CI when available |
| **Verifiable reproducibility** | A content-addressed *receipt* (tool version + hashes of inputs, params, outputs) a third party can verify without re-running | JSON manifest emitted by every subcommand; hashed in the determinism test; a `rosalind verify` path re-checks a receipt |
| **Embeddable / composable** | Library-first; the pileup column stream is a public, documented, Python-exposed substrate for arbitrary (ML/QC/coverage) analytics; pipe-native CLI | Public-API test (no htslib types leak); a worked plugin + a Python streaming example; API-stability commitment on `core` types |

---

## 3. Design tenets that make us uniquely valuable

1. **Memory as a contract.** Every streaming stage exposes a provable working-set bound and
   accepts a `MemoryBudget`; the engine adapts blocking to honor it or fails fast with the
   computed requirement. `rosalind plan` answers "will this run in 2 GB?" before you commit.
2. **The pileup is an open substrate.** Variant calling is *a* consumer of the
   `PileupColumn` stream — coverage, QC, methylation, and ML feature extraction are equal
   citizens via the plugin/iterator API, all inheriting the memory + determinism guarantees.
3. **Library-first, embeddable, pipe-native.** Clean typed Rust API + a thin Python binding
   exposing the streaming primitive; the CLI is a thin client; subcommands compose over pipes
   without mandatory big intermediates (field/disk-scarce reality).
4. **Edge-data realism.** Read-length-agnostic and robust to indel-rich **long reads**
   (Nanopore/PacBio — the actual edge-sequencing modality); transparent gzip/bgzf; honest
   handling of real references (N runs, IUPAC, many contigs).
5. **Honest, calibrated confidence over overcalling.** Calls carry true Phred-scaled,
   well-calibrated confidence and an explicit FILTER vocabulary; below-evidence sites are
   flagged/abstained, never emitted as confident calls. (Mirrors the abstention-aware posture:
   don't pretend to predict where evidence is insufficient.)
6. **Verifiable reproducibility.** A portable receipt makes "this artifact came from exactly
   these inputs, params, and tool version" checkable by anyone — the auditability story made
   real and unique.

---

## 4. Target module map

The lingua franca is `core`; every other layer is built over it. Annotated **KEEP** /
**REWRITE** / **NEW** / **SEPARATE**.

```
core/        coordinate model + canonical records — shared by everything
  locus        Contig + Position; per-contig u32, global 64-bit-safe SA offset   NEW
  sequence     2-bit DNA + explicit IUPAC/N handling                             REWRITE (from compressed_dna)
  record       Read / AlignedRead: CIGAR, SAM flags, MAPQ, strand, tags          REWRITE (from types.rs)
  budget       MemoryBudget + working-set model                                  NEW
  error        typed library error taxonomy                                      NEW
io/          the standards boundary (streaming, pipe-native)
  fasta        multi-contig, streaming                                           REWRITE
  fastq        streaming, gzip/bgzf transparent decompression                    REWRITE
  bam          read/write, deterministic sort, .csi index                        KEEP sort + EXTEND
  vcf          header model + record model → one spec-valid emitter              REWRITE
index/       FM-index/SA as a build-once → serialize → mmap-load artifact         KEEP kernels + REWRITE persistence
align/       seed → chain → extend → alignment (real MAPQ, soft-clip, contigs)    KEEP math + REWRITE (later phase)
pileup/      THE single streaming, CIGAR-aware, filtered, bounded engine          REWRITE  ← the kernel; public substrate
call/        pure scoring over columns: germline GL model; somatic LLR; (later) indels  NEW germline + KEEP somatic (re-homed)
eval/        truth-set comparison (precision/recall, normalize, BED)              KEEP (+ fix left-align)
provenance/  content-addressed run manifest + `verify`                           NEW (pulled early)
framework/   bounded map-reduce substrate (space + ledger + tree)                KEEP — for PLUGINS, off the calling path
plugin/      GenomicPlugin extension surface over the column stream              KEEP (elevated to a first-class substrate API)
bindings/py  thin Python API exposing the streaming pileup primitive            NEW (near phase)
theory/      Williams/Cook–Mertz √t Turing-machine simulation demo               SEPARATE (feature `theory`, non-default)
cli (main)   thin orchestration; typed errors → anyhow only here; pipe-native    REWRITE (slim the 1645-line main.rs)
```

**Public-API discipline (hard rule):** no third-party types (e.g. `rust_htslib::bam::Record`)
appear in any public signature. Sources convert to owned `core` types at the boundary. `core`
types carry an API-stability commitment so external/Python/plugin consumers can rely on them.

---

## 5. Core abstractions (get these right; everything follows)

- **`Locus = (Contig, Position)`** — per-contig `u32` positions; the FM-index/SA over the
  *concatenated* genome uses a **64-bit-safe global offset** (standard bwa/minimap2 model).
  Removes the single-contig blocker and the 4.29 Gbp ceiling structurally.
- **`AlignedRead`** (canonical) — contig-aware, CIGAR, SAM flags, MAPQ, strand, tags; SEQ stored
  forward-reference-oriented; `end()` CIGAR-derived. Read-length-agnostic (Illumina or long
  reads).
- **`PileupColumn`** — the public substrate type: one reference position with
  deterministically-ordered observations `{allele, base_qual, strand, mapq}`; size bounded by
  *coverage*. Designed to be consumed by callers, plugins, and (via the Python binding) ML
  pipelines; cheap aggregate views for featurization.
- **`MemoryBudget`** — declared cap; streaming stages report `working_set_bound()` and honor or
  refuse.
- **`Call`** types; **typed library errors**; `anyhow` only at the CLI.

---

## 6. The streaming pileup engine — the kernel

One engine: **contig-aware, CIGAR-aware, read-filtered, strand-aware, deterministic, bounded,
budget-honoring, read-length-agnostic.** Sources: sorted BAM now → `IndexedReader::fetch`
(Phase D) → in-memory sorted reads. Consumers: germline + somatic callers, plugins, Python.
Clean rewrite (not a generalization of htslib-coupled `BamPileupStream`). This is the component
we RSS-gate (promise #2), expose to Python (promise #5/embeddable), and let plugins extend
(tenet #2).

---

## 7. Keep / rewrite / cut

| Preserve (re-home/extend) | Rewrite (structurally wrong) | Separate / cut |
|---|---|---|
| SA-IS suffix array → 64-bit-safe | Pileup engine → one clean kernel | Theory layer → feature-gated `theory/` |
| FM-index math → multi-contig, persistable | Germline caller → calibrated GL model | Unused deps `ff`, `ark-poly` → optional |
| External merge BAM sort (`sort.rs`) | VCF layer → header+record model | Tautological space counter → real RSS gate |
| Eval suite (+ fix left-align) | Core coordinate/record types (`core/`) | Heuristic `ledger.all_merges_complete()` |
| Somatic LLR → `call/somatic` | Index persistence (serialize + mmap) | `unimplemented!()` stubs in `util/` |
| Plugin trait (elevated to substrate API) | Read ingestion (stream, gzip, multi-contig) | |

---

## 8. Locked decisions (2026-05-26)

- **Rebuild approach: evolve toward the target** — rewrite the wrong layers; re-home/extend the
  correct tested kernels; each phase lands green; no throwaway of working code.
- **Theory layer: separate as a feature-gated `theory/` demo**; `ff`/`ark-poly` become optional
  (`theory = ["dep:ff","dep:ark-poly"]`); decoupled from the genomics core.
- **Guiding principle for all open design choices:** maximize *unique* value to edge builders —
  the kernel/substrate, memory-as-contract, calibrated honesty, verifiable reproducibility,
  embeddability (§1, §3).

---

## 9. Sequencing

1. **Phase A — the kernel + calling vertical + receipt.** `core/` (locus, sequence, record,
   budget, `PileupColumn`), the public `pileup/` engine, `call/germline` (calibrated GL model)
   + re-homed `call/somatic`, `io/vcf` (spec-valid), and a **minimal `provenance/` receipt** so
   the first improved outputs are already verifiable. Detailed in
   `2026-05-26-phase-a-unified-pileup-genotype-design.md`.
2. **Phase B — make it real on a genome.** Index persistence (build-once → mmap), multi-contig,
   streaming/gzip/multi-record ingestion, pipe-native composition.
3. **Phase D — fast & honest, memory-as-contract.** rayon (deterministic), `.csi` region
   `fetch`, **real RSS gates + `rosalind plan`** (budget envelope) + `rosalind verify` (receipt
   check), honest memory docs.
4. **Phase C — scoped calling.** Germline indels, richer read QC (mate-overlap dedup,
   strand-bias filter, BQ/BAQ), long-read calling refinements.
5. **Phase P — Python substrate binding.** Expose the `PileupColumn` stream + index/align/call
   to Python with type stubs, pytest, and a wheel matrix (manylinux/macOS/aarch64). Sequenced
   right after the engine API stabilizes (after A/B).
6. **Phase E — interleaved.** Theory separation (§8), CI hardening (clippy `-D warnings`,
   bcftools validation, wheels), README-honesty pass.

---

## 10. Audit findings index (rationale)

Urgent correctness (A): reverse-strand double-complement (`pileup_stream.rs:33-42`,
`somatic/model.rs:325-327`); CIGAR-ignoring/read-dropping pileup (`pileup.rs`,
`pileup_stream.rs:136`, `main.rs parse_cigar` bail); discarded-posterior "caller"
(`statistics.rs:70`); hardcoded MAPQ (`bwt_aligner.rs:206`); non-conformant VCF (`vcf.rs:6`);
empty-position recursion (`pileup_stream.rs:192`). Foundational (B): single-contig
(`main.rs:881`); index rebuilt every run / zero SA samples persisted (`index/io.rs:140`);
O(reference) build RAM. Moat (D): tautological space counter (`space/allocator.rs`),
`rss_budget` caps 50 KB at 2 GB. Trust (E): vestigial theory layer; no clippy; Linux-only CI;
empty bench; Python = RNA-seq demo only.
