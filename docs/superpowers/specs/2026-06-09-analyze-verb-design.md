# Design: `rosalind analyze <kind>` — a verifiable receipt for any ColumnAnalyzer

**Status:** Approved design — 2026-06-09. Audience: the implementer. Companion to the ColumnKit SDK
(`src/call/columnkit.rs`) and [`CONTRACT.md`](../../../CONTRACT.md).

---

## 1. The goal, in one sentence

Make ColumnKit's promise — *"a verifiable receipt for YOUR metric"* — literally true by wiring
`ColumnAnalyzer::params()` into the receipt through one shared, analyzer-agnostic path.

## 2. Why this, why now

A third-party `ColumnAnalyzer` already inherits the bounded whole-genome walk + the governor +
byte-determinism (via `run_bounded_whole_genome`). But `ColumnAnalyzer::params()` is consumed by
**zero production code** — it is exercised only by a unit test (`columnkit.rs:211`). The shipped
`features` receipt is hand-built for `FeatureAnalyzer` specifically (`run_features` inserts
`feature_rows` directly, `main.rs:~1913`). So a builder who writes their own analyzer inherits the
walk but **not** a receipt for their metric — they'd have to hand-build it. This closes that gap:
`analyze <kind>` runs a registered analyzer through one path that merges `analyzer.params()` into the
content-addressed claim. Combined with the index-receipt/chain work, a stranger's analyzer auto-joins
the verifiable provenance fabric for free.

## 3. Architecture — one shared path (approved)

`run_features`'s ~270 lines of contract machinery (open index → pileup params → predict + refuse →
governor → governed drive → peak/verdict → receipt) become a generic helper, parameterized by the
subcommand label, a receipt-param **prefix**, and the analyzer:

```rust
fn run_bounded_analysis(
    subcommand: &str,        // "features" | "analyze features" | "analyze coverage"
    param_prefix: &str,      // "" for features (byte-identical); "analyzer." for analyze
    analyzer: &mut dyn ColumnAnalyzer,
    index_path: PathBuf,
    alignments_path: PathBuf,
    mapq_threshold: u8,
    memory_budget_mb: Option<u64>,
    max_depth: u32,
    max_read_len: u32,
    enforce: bool,
    output: Option<PathBuf>,
    manifest_out: Option<PathBuf>,
) -> Result<()>
```

It is today's `run_features` body **verbatim**, with exactly two changes:

1. **The analyzer is a parameter** (not a hardcoded `FeatureAnalyzer::default()`).
2. **The hand-inserted `feature_rows` is replaced** by merging `analyzer.params()` into the claim as
   `{param_prefix}{key}`, with a defensive assertion that no `{param_prefix}{key}` lands in
   `MEASUREMENT_KEYS` (so an analyzer can never shadow a measured field out of the claim):

   ```rust
   for (k, v) in analyzer.params() {
       let key = format!("{param_prefix}{k}");
       debug_assert!(
           !rosalind::provenance::MEASUREMENT_KEYS.contains(&key.as_str()),
           "analyzer param {key} collides with a measurement key"
       );
       manifest.params.insert(key, v);
   }
   ```

`CommandCapture::new(subcommand)` records the kind into the replayable `command` (e.g.
`analyze coverage --index @in:… -o @out:…`), so `reproduce` replays the right invocation.
`RunManifest::new(subcommand)` keeps the subcommand field informative.

All other receipt fields (the `CommandCapture` operands/opts, `peak_rss_bytes`,
`predicted_peak_rss_bytes`, `max_working_set_bytes`, `governor`, `baseline_rss_bytes`,
`rss_residual_bytes`, `io_rss_overhead_assumed_bytes`, `over_max_depth`, `reads_skipped_total`,
`contract_verdict`) are inserted exactly as today.

### 3.1 The three call sites

| Caller | `subcommand` | `prefix` | analyzer | receipt param |
|---|---|---|---|---|
| `run_features` (the flagship subcommand, unchanged surface) | `"features"` | `""` | `FeatureAnalyzer` | `feature_rows` — **byte-identical to today** |
| `analyze features` | `"analyze features"` | `"analyzer."` | `FeatureAnalyzer` | `analyzer.feature_rows` |
| `analyze coverage` | `"analyze coverage"` | `"analyzer."` | `CoverageTrack` | `analyzer.analyzer=coverage` |

`run_features` becomes a thin wrapper:
```rust
fn run_features(/* same args */) -> Result<()> {
    let mut analyzer = rosalind::call::FeatureAnalyzer::default();
    run_bounded_analysis("features", "", &mut analyzer, index_path, alignments_path,
        mapq_threshold, memory_budget_mb, max_depth, max_read_len, enforce, output, manifest_out)
}
```

**Byte-identity:** `FeatureAnalyzer::params()` already returns `{"feature_rows": N}`, so `prefix=""`
inserts the exact same `feature_rows` key the hand-built path did. Canonical JSON sorts keys, so the
receipt is byte-identical regardless of insertion order. The output TSV is produced by the analyzer
(untouched), so it is unchanged too.

## 4. A first-party `CoverageTrack` (`src/call/columnkit.rs`)

Add a pub `CoverageTrack` next to `FeatureAnalyzer` — the second registered analyzer and a reference
implementation:

```rust
/// A per-locus coverage track: (contig, 1-based pos, depth). The second first-party
/// ColumnAnalyzer — proof the SDK carries more than the feature egress.
#[derive(Debug, Default)]
pub struct CoverageTrack;

impl ColumnAnalyzer for CoverageTrack {
    fn header(&self) -> Option<String> {
        Some("#contig\tpos\tdepth\n".to_string())
    }
    fn params(&self) -> BTreeMap<String, String> {
        BTreeMap::from([("analyzer".to_string(), "coverage".to_string())])
    }
    fn on_column(&mut self, col: &PileupColumn, contig: &str, out: &mut dyn Write) -> io::Result<()> {
        writeln!(out, "{contig}\t{}\t{}", col.locus.pos.0 + 1, col.depth())
    }
}
```

Re-export it at `src/lib.rs:77` alongside `FeatureAnalyzer`:
`pub use call::{run_bounded_whole_genome, ColumnAnalyzer, CoverageTrack, FeatureAnalyzer};`.

**The example (`examples/columnkit_coverage.rs`) is left untouched** — its inline `CoverageTrack` keeps
its cookbook/teaching value, and not touching the file avoids a merge conflict with the in-flight
PR #79 (which edits the same example). The small duplication (a teaching copy vs the registered copy)
is intentional.

## 5. The CLI

A new `Analyze` subcommand mirroring `features`' flags, with a `kind` positional:

```rust
    /// Run a registered per-locus analyzer over the bounded whole-genome walk, with a
    /// verifiable receipt that records the analyzer's own params.
    Analyze {
        /// Which analyzer to run.
        #[arg(value_enum)]
        kind: AnalyzerKind,
        // --index, --alignments, --mapq-threshold, --memory-budget-mb, --max-depth,
        // --max-read-len, --enforce, -o/--output, --manifest — same as Features.
    },
```
```rust
#[derive(Copy, Clone, Debug, ValueEnum, Eq, PartialEq)]
enum AnalyzerKind { Features, Coverage }
```

Dispatch builds the analyzer and calls the shared helper with `prefix = "analyzer."` and the
kind-aware subcommand label:
```rust
Commands::Analyze { kind, index, alignments, mapq_threshold, memory_budget_mb,
                    max_depth, max_read_len, enforce, output, manifest } => {
    let label = match kind { AnalyzerKind::Features => "analyze features",
                             AnalyzerKind::Coverage => "analyze coverage" };
    match kind {
        AnalyzerKind::Features => {
            let mut a = rosalind::call::FeatureAnalyzer::default();
            run_bounded_analysis(label, "analyzer.", &mut a, index, alignments, mapq_threshold,
                memory_budget_mb, max_depth, max_read_len, enforce, output, manifest)?
        }
        AnalyzerKind::Coverage => {
            let mut a = rosalind::call::CoverageTrack;
            run_bounded_analysis(label, "analyzer.", &mut a, index, alignments, mapq_threshold,
                memory_budget_mb, max_depth, max_read_len, enforce, output, manifest)?
        }
    }
}
```

Adding a future kind (methylation, GC, …) is a one-line enum addition + a match arm — the contributor
on-ramp the SDK is meant to be.

## 6. Testing

- **`features` byte-identical (the load-bearing guard):**
  - The existing `tests/features.rs` suite stays green (output TSV + receipt unchanged).
  - A new assertion: a `features --index … -o out.tsv` receipt's claim still contains the key
    `feature_rows` (and NOT `analyzer.feature_rows`).
- **`analyze coverage` end-to-end (`tests/analyze.rs`, new):** `analyze coverage --index … -o cov.tsv`
  exits 0; the receipt's claim contains `analyzer.analyzer=coverage`; `rosalind verify --manifest …`
  exits 0; the TSV has the `#contig\tpos\tdepth` header.
- **No `MEASUREMENT_KEYS` collision (unit, receipt crate or `columnkit.rs`):** finalize a `RunManifest`
  carrying an `analyzer.`-prefixed param and assert it remains in `params` (the claim) — not relocated
  to `measurements` — proving the prefix protects the claim. (Statically: no `MEASUREMENT_KEY` starts
  with `analyzer.`.)
- **SDK equivalence stays green:** the existing `columnkit.rs` test
  (`feature_analyzer_via_driver_equals_the_direct_features_path`) is untouched and passes.

**Gates:** `cargo test`, `cargo clippy --all-targets -- -D warnings`, `cargo fmt --check`.

## 7. Non-goals / guardrails

- **No dynamic plugin loading** — a compile-time enum registry only (no `.so`/ABI). The SDK is
  first-party-extensible; third parties fork-and-add-a-kind or vendor the crate.
- **Not wired into `variants`/`gvcf`** — those have their own receipt shapes.
- The receipt attests **reproducibility**, never that a metric is biologically correct.
- **Do not modify `examples/columnkit_coverage.rs`** (PR #79 conflict avoidance; the example keeps its
  inline teaching copy).
- `run_features`'s observable behavior (TSV bytes, receipt bytes, exit codes, governor) must not change
  — it is a pure extraction, pinned by `tests/features.rs` + the `feature_rows` key assertion.
