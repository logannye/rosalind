# Hardening Fast-Follow (design)

**Status:** DESIGN SPEC — 2026-06-02. The cross-cutting hardening items deferred from the keystone
(`2026-06-02-trust-correctness-hardening-keystone-design.md` §7). **Branch:**
`rosalind/hardening-fast-follow` (off `main` `7ffe6b0` = keystone + D0 merged). Source: the 2026-06-02
reflection audit. User chose "finish the hardening fast-follow" as the increment.

## 1. Goal

Close the remaining trust/correctness edges the keystone left, and align the front-door surfaces with
the shipped reality. Four independent fixes, one coherent PR (each small, each tested where it's code).

## 2. The four items

### Item 1 — sort k-way-merge tie-break (finding #5, high)

**File:** `src/genomics/sort.rs`.

`HeapItem::cmp` (sort.rs:188-192) orders only by `SortKey` (tid, pos, is_reverse, qname); `source_idx`
is stored but unused. For records with a fully-equal key (duplicate-marked reads; a primary + an
overlapping mate/split sharing qname+pos+strand), `BinaryHeap::pop` returns an unspecified one,
depending on the chunk partition — i.e. on `--memory-mb`. So the same input sorted at two budgets can
produce byte-different sorted BAMs (violates `docs/determinism.md` Rule 2/3).

**Fix:** add `source_idx` as the final tie-break in `HeapItem::cmp`, reversed to match the max-heap
convention `SortKey::cmp` already uses (so the **lower** `source_idx` pops first). Because records are
read in input order and assigned to chunks sequentially, and each chunk is stable-sorted, "ties broken
by ascending `source_idx`" reproduces **input order** for equal-key records — independent of how they
were partitioned into chunks (i.e. budget-invariant). Make `PartialEq` consistent with the now-total
`Ord` (compare `key` AND `source_idx`); at most one item per chunk is in the heap at once, so
`source_idx` is a strict tie-break and `cmp` never returns `Equal` for distinct heap items.

```rust
impl Ord for HeapItem {
    fn cmp(&self, other: &Self) -> Ordering {
        // Total order: key first (reversed for the min-key-first max-heap), then
        // source_idx (also reversed) so equal-key records pop in input order —
        // making the merge output independent of the chunk partition (--memory-mb).
        self.key
            .cmp(&other.key)
            .then_with(|| other.source_idx.cmp(&self.source_idx))
    }
}
impl PartialEq for HeapItem {
    fn eq(&self, other: &Self) -> bool {
        self.key == other.key && self.source_idx == other.source_idx
    }
}
```

**Test (new `#[cfg(test)] mod tests` in sort.rs):** build two `Record`s with identical (tid, pos,
strand, qname), wrap as `HeapItem::new(0, a)` and `HeapItem::new(3, b)`, push both into a
`BinaryHeap`, pop twice, assert the first popped has `source_idx == 0` (lower pops first = input
order). A third record with a smaller pos pops before both regardless of `source_idx`.

### Item 2 — `@SQ` contig-length cross-check (finding #7, high)

**File:** `src/io/bam.rs`.

`record_to_aligned_read` maps the reference by **name only** (`contigs.by_name`, bam.rs:35) — the BAM
header's `LN` is never compared to the index contig length. A BAM aligned to a different-length
`chr1` (a different assembly/patch) is silently accepted, producing coordinate-shifted/truncated calls
with no error. samtools/bcftools reject `@SQ` length mismatches; Rosalind does not.

**Fix:** in `StreamingBamSource::new` (the bounded `--index` contract path), after opening, validate
every header `@SQ` entry whose name is also in the index `ContigSet`: assert `header.target_len(tid)
== Some(contig.length as u64)`; on mismatch return `CoreError::MalformedRecord` naming the contig and
both lengths. Names in the header but absent from the index keep the existing skip behavior (the
records are skipped downstream). One-time check in `new()` (not per record). The `ContigSet` API is
`by_name(&str) -> Option<&Contig>` with `Contig.length: u32` (core/locus.rs); header lengths via
`HeaderView::target_count()` / `tid2name(tid)` / `target_len(tid) -> Option<u64>`.

```rust
fn validate_contig_lengths(header: &bam::HeaderView, contigs: &ContigSet) -> Result<(), CoreError> {
    for tid in 0..header.target_count() {
        let name = std::str::from_utf8(header.tid2name(tid))
            .map_err(|_| CoreError::MalformedRecord("BAM reference name is not UTF-8".into()))?;
        if let Some(c) = contigs.by_name(name) {
            let hlen = header.target_len(tid);
            if hlen != Some(c.length as u64) {
                return Err(CoreError::MalformedRecord(format!(
                    "BAM @SQ length for contig '{name}' ({}) disagrees with the index ({}); \
                     the alignments were built against a different reference",
                    hlen.map(|l| l.to_string()).unwrap_or_else(|| "missing".into()),
                    c.length
                )));
            }
        }
    }
    Ok(())
}
```
Call it in `StreamingBamSource::new` before `Ok(Self { … })`.

**Test (bam.rs tests):** write a BAM whose `@SQ` `LN` for `chr1` is 2000 against a `ContigSet` where
`chr1` is 1000 → `StreamingBamSource::new` returns `Err`; the matching-length case still constructs Ok
(the existing `streaming_source_*` tests cover the Ok path — they use matching lengths).

### Item 3 — BSD `ru_maxrss` unit fix (finding #6, high)

**File:** `src/util/rss.rs`.

The `#[cfg(not(target_os = "linux"))]` branch returns `ru_maxrss` raw, assuming macOS-style **bytes**.
That is wrong on the BSDs (FreeBSD/NetBSD/OpenBSD/DragonFly), which report **KiB** — a silent 1024×
**under-count** (the dangerous direction: `--enforce` would say "within" / `verify` pass for a job
that blew its budget 1000×). The edge/field/clinical appliances the README courts include BSD-based
NAS/storage boxes.

**Fix (minimal + correct):** invert the cfg. Darwin (macOS/iOS) is the only common platform reporting
bytes; Linux + all BSDs + other unix report KiB. So:

```rust
let raw = usage.ru_maxrss as u64;
// ru_maxrss units differ by platform: bytes on Darwin (macOS/iOS), KiB on Linux
// and the BSDs. Treat Darwin as bytes; everything else as KiB (×1024). For an
// unenumerated target this defaults to KiB — the common case, and the safe
// (never-undercount) direction for the memory contract.
#[cfg(any(target_os = "macos", target_os = "ios"))]
{
    raw
}
#[cfg(not(any(target_os = "macos", target_os = "ios")))]
{
    raw.saturating_mul(1024)
}
```

Update the doc comment accordingly.

**Test (rss.rs tests):** allocate ~64 MiB, touch one byte per 4 KiB page (force residency),
`black_box` it, then assert `peak_rss_bytes() >= 32 * 1024 * 1024`. On the run platform (dev macOS =
bytes; CI Linux = KiB×1024) the value reflects ~64 MiB; a dropped `×1024` on Linux would yield ~64 KiB
< 32 MiB and fail. (Lower-bound assertion only — peak RSS is a process-global high-water mark, so
other tests can only raise it.)

### Item 4 — front-door √t honesty pass (low, but the most jarring inconsistency)

**Files:** `Cargo.toml`, `src/main.rs`, `scale_test_results.txt` (delete), `CONTRACT.md`.

The README/`OPEN_PROBLEMS` honestly demote √t to "not yet load-bearing," but three skimmable surfaces
still sell it as shipped/verified:

- `Cargo.toml:6` description `"Accessible genomics engine with O(√t) space complexity …"` → rewrite
  contract-first, e.g. `"Deterministic, low-memory genomics engine: memory as a verifiable contract
  (declare → predict → honor → verify) for alignment and variant calling"`.
- `src/main.rs:25` CLI about-string `"Genomic analysis engine using O(√t) space"` → e.g.
  `"Deterministic low-memory genomics engine with a verifiable memory contract"`.
- `scale_test_results.txt` (repo root, 278 lines, **unreferenced** by any code/doc/script) — lines
  188-278 are an "O(√t) Space Complexity Verification" with 11× `✓ Space scales as O(√t)`, measured
  on the tautological `SpaceTracker` counter (per the strategy note, the counter is self-incremented,
  so it verifies the counter, not real RSS). A browsing builder reads it as evidence √t is shipped.
  **Delete it** (`git rm`): it is a stray captured-output artifact, not referenced anywhere, and it
  contradicts the README's own honest demotion.
- `CONTRACT.md:107` `"identical inputs produce a byte-identical VCF and a byte-identical manifest"` is
  false — the manifest embeds the machine-dependent realized `peak_rss_bytes` (varies run-to-run).
  Reword to: `"identical inputs produce a byte-identical VCF and a manifest identical except for the
  realized peak_rss_bytes (a machine-dependent measurement)"`.

No test (docs/strings); the build + the existing CLI `--help` smoke remain green.

## 3. Cross-cutting requirements

- Per-item commit; `cargo fmt --check` clean; 0 warnings (debug + release); full `cargo test` green at
  each boundary.
- The sort tie-break must not change output for distinct-key inputs (only equal-key ordering becomes
  deterministic). The `@SQ` check must not break the existing `streaming_source_*` tests (they use
  matching lengths).
- Deleting `scale_test_results.txt`: confirmed unreferenced (grep across `.rs`/`.md`/`.toml`/`.sh`
  found no references) and contradicts the shipped story — safe to remove.

## 4. Out of scope (further fast-follow, not this PR)

- **cgroup-awareness + manifest os/arch provenance** (finding #6 remainder) — read cgroup memory
  limits (where the OOM-killer fires in containers) and stamp os/arch into the receipt. Larger design
  surface (enforcement semantics under a cgroup limit); belongs with a deployment/portability story.
- The bigger strategic bets (adoption on-ramp; ML feature substrate) — separate increments.

## 5. Self-review

- **Coverage:** 4 items, each with file + fix + (where code) a test. ✓
- **Type consistency:** sort `HeapItem` Ord/PartialEq stay consistent; `validate_contig_lengths`
  uses the confirmed `ContigSet`/`HeaderView` APIs; rss cfg branches both return `u64`. ✓
- **No placeholders.** The BSD cfg-inversion judgment call (bytes iff Darwin; unknown → KiB = safe
  direction) is stated explicitly. ✓
- **Scope:** one coherent PR; cgroup/provenance + strategic bets explicitly deferred (§4). ✓
