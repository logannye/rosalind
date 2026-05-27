# Phase B3c — `rosalind index` + load (design)

**Status:** Spec for review — 2026-05-27. The CLI that turns B3b's persistence machinery into a usable *build-once → query* workflow, under the contract thesis in [`docs/OPEN_PROBLEMS.md`](../../OPEN_PROBLEMS.md). Follows B3b (zero-copy persisted FM-index, merged — PR #16, `990ae07`).

## 1. The capability we are shooting for

> **`rosalind index` builds the reference into a portable artifact once; thereafter `rosalind locate` (and, in B4, `align`/`variants`) memory-maps that artifact and answers queries in milliseconds — it never rebuilds — and the build is byte-identically reproducible.**

B3b delivered this at the library level (`IndexWriter`/`IndexReader` → `FmIndexView`). B3c is the **CLI surface**: a `rosalind index` subcommand that builds + persists a multi-contig index, and a minimal `rosalind locate` that loads + exact-match-queries it without rebuilding. It also lays the **first `MemoryBudget` hook** — a declared budget + a working-set estimate + a plan/receipt line — the visible foundation `rosalind plan`/`verify` (Phase C) builds on. The seed/chain/extend aligner and wiring `align`/`variants` onto the persisted index remain **B4**; budget *enforcement* + `rosalind plan` remain **Phase C**; leaner/sublinear-space construction remains **Phase D**.

## 2. Scope

**In:**
- A **`rosalind index`** subcommand: stream-read all FASTA records → build a multi-contig `GenomeIndex` → persist via `IndexWriter::write_genome_index`. Emits a deterministic **build receipt**.
- A **`rosalind locate`** subcommand: `IndexReader::open` → `GenomeIndexView::locate_exact` → print `contig<TAB>pos` loci. The load+query demo — mmap, no rebuild, exact-match only.
- A **`MemoryBudget` hook**: a `--memory-budget-mb` flag + a coarse build **working-set estimate** + a **plan line** (estimate vs budget, `OK`/`OVER`); realized build peak RSS reported in the receipt. **Record-only — never refuses/aborts.**
- **Pure, unit-testable library helpers**: the working-set estimate and receipt/plan formatting live in the library (return values / `String`s), so the CLI handlers stay thin and the logic is tested without spawning a process.
- **Integration gates** + a short **docs/CLI** note on the build-once → query workflow.

**Out (deferred, by design):**
- Wiring `align`/`variants`/`somatic` onto the persisted index (`--index` instead of `--reference`) → **B4**.
- The seed/chain/extend aligner over the multi-contig index → **B4**.
- `MemoryBudget` *enforcement* (honor-or-refuse, graceful degradation) + `rosalind plan` + `rosalind verify` → **Phase C** (B3c lays the budget *seam* + plan/receipt surface; it does not gate execution on the budget).
- Leaner / sublinear-space construction (the build is still `O(reference)` RAM — the caveat D erases) → **Phase D**.

## 3. `rosalind index` (build + persist)

```
rosalind index --reference <fa> --output <idx> [--memory-budget-mb M]
```

- **Read:** `io::fasta::FastaReader` over the (optionally gzipped, via `io::decompress`) reference, yielding every `FastaRecord { name, sequence }` (sequences already ASCII-uppercased by the reader) — **all contigs**, not just the first. Collected into `Vec<(String, Vec<u8>)>`.
- **Build:** `GenomeIndex::from_named_sequences(&named)` (assembles the `ContigSet` + concatenated reference, validates non-empty / `≤ MAX_GENOME_LEN` / A·C·G·T·N; block size is `GenomeIndex`'s existing internal `√len` heuristic). A configurable block size is **not** exposed in B3c (YAGNI — no consumer needs it yet; trivial to add later).
- **Persist:** `IndexWriter::create(output).write_genome_index(&index)` — the deterministic artifact.
- **Receipt** (stdout, deterministic except the explicitly-marked realized-RSS line): a fixed-format block — index path; per-contig `name\tlength`; contig count + total bp; the reference BLAKE3 (hex); the on-disk index size (bytes); and a realized build peak-RSS line (`util::rss::peak_rss_bytes`, informational). With `--memory-budget-mb`, a **plan line** precedes the build (see §5).

**Determinism:** the **`.idx` artifact** is byte-identical across runs (inherited from B3b's serializer, gated). The receipt's contig/size/BLAKE3 lines are deterministic; the realized-RSS line is per-run informational and is **not** part of any determinism gate.

## 4. `rosalind locate` (load + query)

```
rosalind locate --index <idx> --pattern <ACGT…> [--max-hits N]
```

- `IndexReader::open(idx)` → `ReferenceIndex::genome_view()` → `GenomeIndexView::locate_exact(pattern, max_hits)` → for each `Locus`, print `<contig_name><TAB><pos>` (name resolved via `genome_view().contigs().by_id(locus.contig)`), one per line, in the sorted `(contig, pos)` order `locate_exact` already guarantees. Empty result → a single `no hits` line on stderr (exit 0). `--max-hits` defaults to 1024; the pattern is uppercased like the reference.
- This is **load + exact-match only** — it opens via mmap and never constructs a `BlockedFMIndex` / calls `sais_u32`. It is *not* the aligner (no seeding/chaining/mismatches) — that's B4.

## 5. The `MemoryBudget` hook (record-only)

- `estimate_build_working_set(reference_len: u64) -> WorkingSet` — a **library** function returning a coarse, **documented-as-approximate** estimate of the index *build* peak (dominated by SA-IS over the `u32` text: roughly the text + suffix array + workspace, plus the built structures). It is intentionally a rough upper-ish model; precise accounting is **Phase C**, and the build cost itself is what **Phase D** reduces.
- With `--memory-budget-mb M`: build `MemoryBudget::from_mb(M)`, compute the estimate, and print a **plan line** *before* building: `plan: est. build peak ~<X> / budget <Y>  [OK|OVER]` (via `WorkingSet::fits` / `MemoryBudget::admits`). **The build then proceeds regardless** — B3c never refuses (honor-or-refuse is Phase C). After the build, the receipt's realized-RSS line lets a user compare estimate ↔ realized ↔ budget by eye.
- This threads `MemoryBudget` + `WorkingSet` (which exist in `core/budget.rs` but are currently used nowhere) through a real pipeline for the first time — the seam Phase C formalizes.

## 6. Decomposition (four green sub-stages)

- **B3c.1 — `rosalind index` build + receipt.** Add the `Index` `Commands` variant + a thin `run_index`; multi-record FASTA read → `GenomeIndex::from_named_sequences` → `write_genome_index`; a pure `format_index_receipt(...) -> String` helper. Gate: the artifact is written and **byte-identical across two builds**; the receipt's deterministic fields are correct.
- **B3c.2 — `rosalind locate` load + query.** Add the `Locate` variant + `run_locate`; `open` → `genome_view` → `locate_exact` → print `name\tpos`. Gate: `locate` matches an in-RAM `GenomeIndex::locate_exact` ground truth over a multi-contig + boundary-straddle battery; the load path calls **no** `sais_u32`/`BlockedFMIndex::build`.
- **B3c.3 — `MemoryBudget` hook.** `estimate_build_working_set` (library, documented coarse model) + `--memory-budget-mb` + the plan line + realized peak in the receipt; record-only. Gate: the plan line reads `OK` when the estimate ≤ budget and `OVER` otherwise, and the build **still succeeds when `OVER`** (no enforcement).
- **B3c.4 — integration gates + docs.** `tests/index_cli.rs` (build→locate equivalence, determinism, no-rebuild, bounded residency) + a README/CLI note on the build-once → query workflow.

## 7. Success criteria (gates)

- **Round-trip:** `rosalind index` then `rosalind locate` returns exactly the loci an in-RAM `GenomeIndex::locate_exact` returns for the same patterns (multi-contig + `N`-bearing + boundary-straddle).
- **No rebuild on load:** the `locate` path never calls `sais_u32` / `BlockedFMIndex::build` (structural grep over the load path + behavioral: querying a pre-built index with the in-RAM source absent).
- **Determinism:** two `rosalind index` runs over the same reference produce **byte-identical** `.idx` files.
- **Bounded residency:** `locate` opens + queries via mmap without making the index resident (the B3b view is a borrow; a basic working-set check — the enforced RSS gate is Phase C).
- **Budget seam:** the plan line correctly reports `OK`/`OVER` from `MemoryBudget`/`WorkingSet`, and B3c **never refuses** a build (record-only).

## 8. Testing

- **B3c.1/.3:** pure-helper unit tests (receipt formatting deterministic; `estimate_build_working_set` monotonic in length; the plan line `OK`/`OVER` logic). A CLI/integration test that `index` writes a byte-identical artifact twice.
- **B3c.2:** an integration test that builds a temp index then `locate`s a battery of patterns, comparing to an in-RAM `GenomeIndex::locate_exact` (the ground truth), including a boundary-straddling pattern (expect no hits) and an `N`-bearing fixture.
- **B3c.4:** the gates above wired as `tests/index_cli.rs`, invoking the built binary via `env!("CARGO_BIN_EXE_rosalind")` + `std::process::Command` (no new dependency; the real end-to-end CLI path). The deterministic-artifact and round-trip-vs-in-RAM checks run through the binary; the pure helpers (receipt, estimate) are unit-tested in-library separately.

## 9. Risks & mitigations

- **`main.rs` is already large (~1.7 kLOC).** *Mitigation:* keep `run_index`/`run_locate` thin; put the testable logic (working-set estimate, receipt/plan formatting) in **library** modules with unit tests; do not restructure the existing handlers.
- **The build is `O(reference)` RAM** (SA-IS + the in-RAM `GenomeIndex` before serialization). *Mitigation:* this is the known, documented caveat (the receipt is honest about realized peak; the budget hook reports but does not enforce); **Phase D** is the construction fix. B3c does not claim a bounded *build*.
- **Working-set estimate is approximate.** *Mitigation:* documented as coarse/record-only; the realized-RSS line gives the ground truth; Phase C refines the model. No decision (no enforcement) rides on it in B3c.

## 10. Decisions (this stage)

- **CLI surface = `index` (build) + a minimal `locate` (load+query)** (§3–§4) — the real "build once → query the artifact" demonstration the B3b spec promised, exercising `IndexReader`→`GenomeIndexView` end-to-end, without the B4 aligner.
- **`MemoryBudget` hook is record-only** (§5) — lays the seam + the plan/receipt surface per the roadmap's "lay budget hooks in B," with **no enforcement** (honor-or-refuse is Phase C).
- **Flat subcommands** (`index`, `locate`) matching the existing `Align`/`Variants`/`Sort`/`Somatic` style — not nested `index build`/`index inspect`.
- **Testable logic in the library, thin CLI handlers** — preserves `main.rs` and makes the receipt/estimate unit-testable.
