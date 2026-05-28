# Phase B4 — bounded whole-genome variant calling over the persisted index (design)

**Status:** Spec for review — 2026-05-27. The Phase-B4 caller sub-stage (wire the consumers onto the persisted index), **reprioritized ahead of the aligner** on builder-value grounds (§1), and scoped to deliver **bounded** whole-genome calling (streaming reads + per-contig reference). Follows B4a (`ReferenceView`, merged — PR #18). Under the contract thesis in [`docs/OPEN_PROBLEMS.md`](../../OPEN_PROBLEMS.md).

## 1. The capability we are shooting for

> **`rosalind variants --index <idx> --alignments <sorted.bam>` calls germline variants across *every contig* of a persisted index, in *bounded, predictable memory* — reads stream one record at a time, the reference is read from the index, the pileup is already columnar — emitting a multi-contig VCF plus a reproducibility-and-memory receipt that *shows* the bounded peak.**

This operationalizes Rosalind's core differentiator — **memory as a declared, predictable, verifiable contract** — on the workload builders actually run. Builders align with bwa-mem2/minimap2 and arrive with a sorted BAM; Rosalind's draw is the *caller* (bounded memory, byte-identical reproducibility, calibrated abstention), not its (not-yet-competitive) in-house aligner. `variants` needs no FM-index/aligner — only the reference (self-contained via B4a) plus the streaming pileup + abstention-aware caller already shipped in Phase A. So this is the shortest path to the flagship "call a whole genome from a portable index + a BAM, on a laptop, reproducibly" capability. The aligner (formerly "B4b") is **deferred** until Phase E makes it genuinely good.

## 2. What makes it bounded (the three memory terms)

Whole-genome calling has three memory terms; this stage bounds the two that aren't already:

1. **Reads — bounded NOW (the core add).** A `StreamingBamSource` reads one record at a time from `bam::Reader` and yields it; it **never materializes the BAM**. (Today's `BamSource`/`SliceSource` pre-load + sort the whole file.)
2. **Pileup working set — already bounded.** `PileupEngine::current_working_set()` is "bounded by the active read set (local coverage), independent of total input size" — we surface it (§5).
3. **Reference — per-contig bounded (this stage), O(window) later.** Each contig's reference is decoded from `ReferenceView` (B4a) into an `Arc<[u8]>`; peak = the **largest contig** (~250 MB for human chr1) — predictable and laptop-friendly. On-demand `ref_base` (shrinking this to O(pileup window)) is a documented refinement (§8), not required for "bounded WGS on a laptop."

Net peak ≈ largest contig's reference + the pileup working set + a few buffered reads — **independent of BAM size**.

## 3. Scope

**In:**
- A **`StreamingBamSource`** (`pileup::ReadSource`) that streams a BAM record-by-record (bounded), sharing the existing per-record `Record → AlignedRead` mapping with `read_bam_as_core_reads` (extract it into one helper — DRY).
- **Coordinate-sort enforcement.** Streaming cannot sort in RAM, so the input must be coordinate-sorted. A **monotonicity guard** errors if reads arrive out of `(contig_id, pos)` order (where `contig_id` is resolved against the *index's* `ContigSet`) — this catches both an unsorted BAM and a BAM whose `@SQ` order disagrees with the index. (Optionally also check the `@HD SO:coordinate` header for a fail-fast message.)
- A **`variants --index <idx>`** mode: load the index (`ContigSet` + `ReferenceView`) and call germline variants across **all contigs** via a single sorted streaming pass.
- **Multi-contig VCF** (one `##contig` per contig — `write_germline_vcf` already does this) + the BLAKE3 run manifest, **extended with the realized peak memory** (engine working-set max + process peak RSS).
- A **record-only `--memory-budget-mb`**: prints a plan line and flags if the realized peak exceeds the budget; **never enforces/aborts** (honor-or-refuse is Phase C). Makes the bounded contract *visible + verifiable* on the flagship workload.
- Back-compat: the single-contig **`--reference <fa>`** path is preserved; `--index` and `--reference` are mutually exclusive (exactly one required).

**Out (deferred, by design):**
- The aligner over the persisted index + `align --index` (formerly B4b) and `somatic --index` → later / Phase E.
- **Streaming the SAM path.** BAM (the real WGS input) streams; the legacy SAM path stays materialized for small/targeted inputs — streaming it is a trivial follow-on, not the flagship.
- **On-demand `ref_base`** (§2 term 3 refinement) → Phase C/D.
- **Budget *enforcement* + `rosalind plan`/`verify`** → Phase C (this stage surfaces the receipt, the foundation those build on).
- Indexed random-access (`.csi`) fetch — a single sorted pass suffices; no index file required.

## 4. The bounded multi-contig drive

The BAM is coordinate-sorted ⇒ reads arrive grouped by contig, ascending position. One streaming pass, partitioned per contig, reusing the tested per-contig caller:

```rust
// call_germline_region<S: ReadSource>(source: S, reference: Arc<[u8]>, contig: u32,
//     region: Range<u32>, pileup_params, &germline_params) -> Vec<(Locus, u8, GermlineCall)>
```

- A `StreamingBamSource` yields all reads in `(contig, pos)` order (one at a time), with a **monotonicity guard**.
- A **per-contig adapter** wraps the shared stream: for contig `c` it yields `c`'s reads and returns `None` at the first read of `c+1` (that read is *peeked/buffered*, not dropped). This lets `call_germline_region` consume a per-contig `ReadSource` without materializing.
- The driver iterates the index's contigs in id order; for each `c`: decode `ref_c = ReferenceView::decode_window(c.global_offset, c.global_offset + c.length)`, run `call_germline_region(per_contig(c), ref_c, c.id, 0..c.length, …)`, accumulate `Locus`-tagged rows, track the max engine working set.
- Finally `write_germline_vcf(out, index.contigs(), …, &rows)` → multi-contig VCF; write the manifest (inputs = `.idx` + BAM; realized peak RSS + working-set; budget verdict if `--memory-budget-mb`).

This reuses `call_germline_region`, `write_germline_vcf`, and the manifest writer; the new code is the streaming source, the per-contig adapter + monotonicity guard, the index-load + per-contig loop + reference decode, and the memory receipt.

## 5. Operationalizing the memory contract (the differentiator, made visible)

- The receipt (and a stderr line) report the **realized process peak RSS** (`util::rss::peak_rss_bytes`) and the **maximum engine working set** observed (`PileupEngine::current_working_set`) — a *verifiable* record that the run was bounded, not a claim.
- `--memory-budget-mb M`: print a plan line up front; after the run, flag `within budget` / `EXCEEDED by …` against the realized peak. **Record-only** — the build/call always completes (no refusal; that is Phase C). Reuses the `MemoryBudget`/`WorkingSet` types and the B3c receipt pattern.

## 6. CLI

```
rosalind variants (--index <idx> | --reference <fa>) --alignments <sorted.bam|sam>
    [--mapq-threshold N] [--quality-threshold Q] [--memory-budget-mb M] [-o out.vcf]
```

- `--index` XOR `--reference`, exactly one required (clap group). `--index` → bounded multi-contig calling (reference + contigs from the index); **alignments must be a coordinate-sorted BAM** — SAM under `--index` is an error (the legacy SAM reader is single-contig; use `--reference` for single-contig SAM). `--reference` → today's single-contig path, BAM or SAM (`--chrom`/`--region-start` apply there only).
- `--chrom`/`--region-start` are **an error** under `--index` (they belong to the single-contig path).

## 7. Decomposition (one spec → three green sub-stages)

- **B4-v.1 — streaming sorted source.** `StreamingBamSource` (`ReadSource` over `bam::Reader`, one record at a time) sharing the `Record → AlignedRead` mapping with `read_bam_as_core_reads`; the `(contig_id, pos)` monotonicity guard (errors on unsorted/contig-order-mismatch). Gate: over a sorted fixture it yields exactly the reads the materializing path yields (same set/order); an out-of-order fixture errors; it holds one record at a time (no Vec of all reads).
- **B4-v.2 — bounded multi-contig calling + `variants --index`.** The per-contig adapter + the driver (§4) + index loading + per-contig reference decode + the multi-contig VCF + manifest. Gates: **parity** — on a single-contig index, records match `variants --reference <same fa>`; **multi-contig** — a 2–3-contig sorted BAM yields the union of per-contig calls with correct `##contig` headers + per-record contigs; **self-contained** — calling works with no `--reference` (and after the source FASTA is deleted); **bounded** — calling does not materialize the BAM (streaming source) and never builds an FM-index/SA-IS.
- **B4-v.3 — memory receipt + `--memory-budget-mb`.** Surface realized peak RSS + max working set in the receipt/stderr; the record-only budget plan/flag. Gate: the receipt records a peak; `--memory-budget-mb 0` flags `EXCEEDED` but the call still completes (record-only).

## 8. Success criteria (gates)

- **Bounded reads:** `variants --index` streams the BAM (no full-file `Vec`); peak memory is independent of BAM size (validated structurally — the streaming source — and by the working-set/RSS receipt).
- **Self-contained:** no `--reference` FASTA needed; the reference comes from the `.idx` (calling unaffected by deleting the source FASTA).
- **Parity:** single-contig `variants --index` == `variants --reference <same fa>` (the index's decoded reference is byte-identical to the FASTA, per B4a).
- **Multi-contig / whole-genome:** a multi-contig sorted BAM yields a VCF with per-contig `##contig` headers and correctly-attributed records; per-contig calls equal calling each contig alone.
- **Sort safety:** an unsorted (or index-order-mismatched) BAM is **rejected with a clear error**, never silently mis-piled.
- **Reproducible + verifiable:** identical inputs → byte-identical VCF; the manifest records the `.idx`+BAM BLAKE3s and the realized peak memory.
- **No rebuild:** the `variants --index` path never builds an FM-index / calls SA-IS.

## 9. Testing

- Integration (extend `tests/index_cli.rs` or a new `tests/variants_index.rs`): build an index via `rosalind index`; run `variants --index` on a small **sorted** BAM (or SAM-built-then-sorted fixture); assert parity with `variants --reference` (single contig), multi-contig `##contig` headers + records (2–3 contigs), self-contained (delete the FASTA), and that an out-of-order BAM errors. Assert the receipt records a peak and `--memory-budget-mb 0` flags-but-completes.
- Unit: the `Record → AlignedRead` mapping helper (shared) and the monotonicity guard (out-of-order → Err) in `io/bam.rs`; the per-contig adapter (yields one contig's reads, peeks the next) in the driver module. The caller/pileup are already unit-tested (Phase A).

## 10. Decisions (this stage)

- **Reprioritize the caller ahead of the aligner** (§1) — builder value + unblocked by B4a + sets up Phase C.
- **Deliver bounded WGS now** — fold in the streaming BAM source (term 1) so the "bounded whole-genome on a laptop" capability is real, not aspirational.
- **Coordinate-sorted input required for `--index`** — streaming can't sort; a `(contig_id, pos)` monotonicity guard enforces it with a clear error (standard caller contract; builders `samtools sort`).
- **Per-contig reference (i)** kept; **on-demand `ref_base` (ii)** deferred — streaming reads is the dominant win; ~largest-contig peak is laptop-fine.
- **Surface the memory contract (record-only)** — realized peak RSS + working set in the receipt + `--memory-budget-mb`; **no enforcement** (Phase C). Makes the differentiator visible + verifiable.
- **BAM streams; SAM stays materialized** — BAM is the WGS input; streaming SAM is a trivial deferred follow-on.
- **Germline `variants` only** — `somatic --index` + the aligner deferred.
- **Reuse `call_germline_region` + `write_germline_vcf` + the manifest** unchanged; the new code is the streaming source + per-contig drive + reference decode + memory receipt.
