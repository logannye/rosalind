# ML Feature Substrate — Egress Core (design)

**Status:** DESIGN SPEC — 2026-06-02. **Branch:** `rosalind/ml-feature-substrate` (off `main`
`9506d45`). From the reflection audit (angle 3): the kernel already produces a bounded, deterministic
per-locus feature stream (`PileupColumn`); it is gated from ML builders only by the lack of a feature
**egress** (the `on_row` sink emits germline calls, not features) and a dead-end Python stub.

## 1. Goal

Expose the kernel as a **bounded, deterministic, byte-identical per-locus feature stream** with a
verifiable hash receipt — landing the novel claim no other tool makes: **byte-identical features →
bit-reproducible training inputs.** Reuse the shipped, accuracy-validated, memory-contracted call path.

## 2. Scope (this increment)

The **durable, fully-testable core**: a `rosalind features` CLI + a `call::features` egress module
streaming a per-locus **TSV** under the same memory contract (`plan`/`--enforce`/`verify`), with a
BLAKE3 receipt. **TSV-first** (no new deps, universally readable by pandas/polars/R, trivially
byte-deterministic) — which fully delivers the reproducibility claim; the receipt hashes the feature
file. **Deferred to a follow-up:** Arrow/Parquet egress (efficiency/zero-copy), the real `pyarrow`
Python boundary (replacing the stub), and a reference model demo.

## 3. What exists (reuse)

- `PileupColumn { locus, ref_base, raw_depth, obs }` + `depth()`, `allele_counts() -> [u32;4]`
  (`[A,C,G,T]`), `strand_counts() -> [[u32;2];4]` (`[allele][0=fwd,1=rev]`); `Obs { allele, base_qual,
  mapq, reverse }`. The engine yields one column per **callable** position (obs non-empty), in
  deterministic canonical-obs order (the keystone fix).
- `call_germline_whole_genome` (`src/call/whole_genome.rs`) — the bounded per-contig driver
  (`PerContig` partition + per-contig `decode_window` + `PileupEngine`), returning
  `(WorkingSet, SkipCounts)`. Make `PerContig` `pub(crate)` and reuse it.
- The memory contract: `estimate_variants_working_set` / `--enforce` / `RunManifest` / `peak_rss_bytes`.

## 4. Deliverables

### 4a. `src/call/features.rs` — the egress module

- **`FeatureRow`** (one per callable locus), TSV columns (tab-separated, `\n`-terminated):
  `contig  pos  ref  depth  raw_depth  a  c  g  t  a_fwd  a_rev  c_fwd  c_rev  g_fwd  g_rev  t_fwd
  t_rev  mean_bq  mean_mapq`
  where `pos` is **1-based** (matches VCF POS = `locus.pos.0 + 1`), `contig` is the contig **name**,
  `ref` is the ASCII ref base, `depth` = callable obs, `raw_depth` = covering reads, the 4 `a/c/g/t`
  are `allele_counts`, the 8 strand columns are `strand_counts`, and `mean_bq`/`mean_mapq` are the
  integer means over `obs` formatted `{:.2}` (byte-stable). Header line written once:
  `#contig\tpos\tref\t…`.
- **`feature_row_fields(col: &PileupColumn, contig_name: &str) -> impl Iterator<…>`** (or a
  `write_feature_row<W: Write>(w, contig_name, col)` + `write_feature_header<W: Write>(w)`), mirroring
  `src/io/vcf.rs`'s streaming writer shape. Means: `sum(base_qual as u64)/depth` etc. → `{:.2}`.
- **`stream_features_region<S: ReadSource>(source, reference, contig, region, pileup_params, on_row:
  &mut dyn FnMut(&PileupColumn) -> Result<(), CoreError>) -> Result<(WorkingSet, SkipCounts),
  CoreError>`** — drive `PileupEngine`, call `on_row` per emitted column, track max working set, read
  `skip_counts()` after. (No buffer; bounded.)
- **`stream_features_whole_genome<S>(source, ref_view, contigs, pileup_params, on_row) ->
  Result<(WorkingSet, SkipCounts), CoreError>`** — per-contig loop reusing `PerContig` + per-contig
  `decode_window`, summing `SkipCounts` (mirrors `call_germline_whole_genome`). The sink gets
  `(&PileupColumn, contig_name)` — pass the contig name through (or the sink resolves it).
- **Unit tests:** `feature_row_from_a_known_column` (exact field values for a hand-built column);
  `feature_stream_is_bounded_by_coverage_not_read_count` (working set flat in read count, like the
  caller); `feature_rows_match_per_contig`.

### 4b. `rosalind features` CLI — `src/main.rs`

`Commands::Features { index, alignments, mapq_threshold, max_depth (1000), max_read_len (250),
memory_budget_mb, enforce, output, manifest }` → `run_features`, structurally mirroring
`run_variants_index`:
- `StreamingBamSource` + `stream_features_whole_genome`, streaming each row to a `BufWriter` (TSV
  header once, then one row per column) — no genome-wide buffer.
- The **same** `--enforce` gate (exit 3 refuse / exit 4 breach; the working-set model is identical —
  it is the same pileup engine), the same realized-peak receipt, the same `over_max_depth`/skip
  surfacing. Receipt `subcommand = "features"`, plus a `feature_rows` param (count emitted).
- `plan --index` already predicts this peak (same engine) — no change needed.

### 4c. Tests + docs

- **Integration test** (`tests/features.rs`): run `rosalind features` on the `build_sorted_bam_fixture`
  (or a local fixture); assert the header + ≥1 row with the expected columns; run **twice** and assert
  the two `features.tsv` are **byte-identical** (the reproducibility claim); assert the receipt records
  `feature_rows` and that `verify` passes.
- **README/CONTRACT** short section: "a bounded, deterministic, byte-identical per-locus feature stream
  — bit-reproducible ML training inputs with a verifiable receipt," with the one-liner
  `rosalind features --index ref.idx --alignments sorted.bam -o features.tsv` and a note that pandas
  reads it in one line. Mark Arrow/pyarrow/model as the roadmap follow-up.

## 5. Cross-cutting

- **Determinism is the headline.** `features.tsv` must be byte-identical run-to-run (canonical obs
  order + integer-sum means formatted `{:.2}` + deterministic column stream). The byte-identical
  integration test is the proof.
- **Bounded memory.** Rows stream to disk; one contig's reference resident; active set depth-capped —
  the same guarantee as the caller, inherited by construction.
- Per-item commit; 0 warnings (debug+release); fmt clean; full `cargo test` green.

## 6. Out of scope (follow-up increment)

- **Arrow/Parquet** egress (zero-copy, columnar efficiency) behind a feature flag.
- **`pyarrow` Python boundary** replacing the `python_bindings` stub (`rosalind.features(index, bam,
  region, budget_mb)` yielding RecordBatches) — needs maturin/pyO3 build verification.
- **Reference model** demo (a notebook/script training a tiny SNV classifier on the stream + showing
  the receipts match → bit-reproducible eval).

## 7. Self-review

- **Coverage:** egress module (4a), CLI (4b), tests + docs (4c). ✓
- **Type consistency:** the feature sink is `&mut dyn FnMut(&PileupColumn) -> Result<(), CoreError>`;
  `stream_features_*` return `(WorkingSet, SkipCounts)` like the germline drivers; `PerContig` becomes
  `pub(crate)`. ✓
- **No placeholders:** the TSV schema, means, and 1-based pos are concrete. ✓
- **Reproducibility:** integer means `{:.2}` + canonical obs order + byte-identical test → the claim is
  proven, not asserted. ✓
- **Scope:** TSV core only; Arrow/pyarrow/model explicitly deferred. ✓
