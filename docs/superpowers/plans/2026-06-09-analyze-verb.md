# `rosalind analyze <kind>` Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `rosalind analyze <kind>`, a verb that runs a registered `ColumnAnalyzer` through one shared contract+receipt path that wires `analyzer.params()` into the content-addressed claim — making ColumnKit's "a verifiable receipt for YOUR metric" promise literally true.

**Architecture:** Extract `run_features`'s contract machinery into a generic `run_bounded_analysis(subcommand, param_prefix, analyzer, …)`; `features` becomes a thin wrapper (`prefix=""`, byte-identical), `analyze <kind>` routes through it (`prefix="analyzer."`). Add a first-party `CoverageTrack` and an `Analyze` CLI verb. A compile-time enum registry, no dynamic loading.

**Tech Stack:** Rust 1.83, clap derive, the ColumnKit SDK (`src/call/columnkit.rs`), the `rosalind-receipt` crate. Tests: stdlib `Command` integration + `#[cfg(test)]` units.

**Spec:** `docs/superpowers/specs/2026-06-09-analyze-verb-design.md`.

---

## File Structure

| File | Change | Responsibility |
|---|---|---|
| `src/main.rs` | Rename `run_features` → `run_bounded_analysis` (+ params + 7 body edits); add a thin `run_features` wrapper; add `AnalyzerKind` enum + `Analyze` clap variant + dispatch | The shared path + the CLI |
| `src/call/columnkit.rs` | Add a pub `CoverageTrack` + a unit test | The second first-party analyzer |
| `src/lib.rs` | Add `CoverageTrack` to the `call::{…}` re-export (line 77) | Export it as `rosalind::call::CoverageTrack` |
| `tests/features.rs` | Add a `feature_rows`-not-prefixed assertion | Byte-identity guard |
| `tests/analyze.rs` | Create | `analyze coverage` end-to-end |
| `crates/receipt/src/lib.rs` | Add a `MEASUREMENT_KEYS`/prefix unit test | Claim-protection guard |

**Verified facts:** `tests/features.rs` asserts the receipt has `feature_rows` + `verify: OK` (NOT the stderr summary), so the trailing summary may be generalized. `MEASUREMENT_KEYS` is `pub` in the receipt crate → `rosalind::provenance::MEASUREMENT_KEYS`. `run_features` carries `#[allow(clippy::too_many_arguments)]`.

---

## Task 1: extract `run_bounded_analysis` (the refactor)

**Files:** Modify `src/main.rs` (`run_features`, 1661–1966); Test: `tests/features.rs`

- [ ] **Step 1: Confirm the baseline is green**

Run: `cargo test --test features 2>&1 | tail -5`
Expected: PASS — this is the regression net for the refactor.

- [ ] **Step 2: Apply the extraction (8 surgical edits)**

**Edit 1 — signature** (rename + add 3 params):
```rust
#[allow(clippy::too_many_arguments)]
fn run_bounded_analysis(
    subcommand: &str,
    param_prefix: &str,
    analyzer: &mut dyn rosalind::call::ColumnAnalyzer,
    index_path: PathBuf,
    alignments_path: PathBuf,
    mapq_threshold: u8,
    memory_budget_mb: Option<u64>,
    max_depth: u32,
    max_read_len: u32,
    enforce: bool,
    output: Option<PathBuf>,
    manifest_out: Option<PathBuf>,
) -> Result<()> {
```

**Edit 2 — drop the hardcoded analyzer** (replace the comment + `let mut analyzer = …`):
```rust
    // Drive the analyzer's bounded whole-genome walk straight to the writer (header
    // once, one row per callable locus) — no genome-wide buffer accumulates. This is
    // the ONE analyzer-agnostic path: `features` and `analyze <kind>` both run here.
    // Inline both arms (match arms are exclusive, so moving source/pileup_params in
    // each is fine). Flush on BOTH paths so partial output survives a governed abort;
    // inspect the Result rather than `?`-propagating it.
```
(removing `let mut analyzer = rosalind::call::FeatureAnalyzer::default();`)

**Edit 3 — reborrow in BOTH drive arms** (two occurrences): change `&mut analyzer,` to `&mut *analyzer,` (the param is now `&mut dyn`, so a reborrow is needed to keep `analyzer` usable for `params()` after the drive).

**Edit 4 — delete** `let feature_rows = analyzer.rows();`.

**Edit 5 — generic receipt subcommand:**
```rust
        let mut manifest = RunManifest::new(subcommand);
        let mut cmd = CommandCapture::new(subcommand);
```

**Edit 6 — replace the `feature_rows` insert with the prefixed `params()` merge:**
```rust
        // Merge the analyzer's own params into the CLAIM under the caller's prefix
        // (`""` for features → byte-identical `feature_rows`; `analyzer.` for analyze).
        // The prefix keeps an analyzer from shadowing a measured field out of the claim.
        for (k, v) in analyzer.params() {
            let key = format!("{param_prefix}{k}");
            debug_assert!(
                !rosalind::provenance::MEASUREMENT_KEYS.contains(&key.as_str()),
                "analyzer param {key} collides with a measurement key"
            );
            manifest.params.insert(key, v);
        }
```

**Edit 7 — generalize the trailing summary:**
```rust
    eprintln!(
        "{subcommand}: peak RSS {} MiB; max pileup working set {} KiB",
        peak_rss / (1 << 20),
        max_ws.bytes / 1024
    );
```

**Edit 8 — add the thin `run_features` wrapper** (immediately after `run_bounded_analysis`'s closing brace):
```rust
/// The shipped `features` egress: the FeatureAnalyzer through the one bounded-analysis
/// path. `param_prefix=""` reproduces the historic `feature_rows` claim key exactly, so
/// the receipt stays byte-identical.
#[allow(clippy::too_many_arguments)]
fn run_features(
    index_path: PathBuf,
    alignments_path: PathBuf,
    mapq_threshold: u8,
    memory_budget_mb: Option<u64>,
    max_depth: u32,
    max_read_len: u32,
    enforce: bool,
    output: Option<PathBuf>,
    manifest_out: Option<PathBuf>,
) -> Result<()> {
    let mut analyzer = rosalind::call::FeatureAnalyzer::default();
    run_bounded_analysis(
        "features",
        "",
        &mut analyzer,
        index_path,
        alignments_path,
        mapq_threshold,
        memory_budget_mb,
        max_depth,
        max_read_len,
        enforce,
        output,
        manifest_out,
    )
}
```

- [ ] **Step 3: Build + confirm `features` is unchanged**

Run: `cargo build 2>&1 | tail -5 && cargo test --test features 2>&1 | tail -5`
Expected: PASS — `features` behavior (TSV, receipt, exit codes) is unchanged; `tests/features.rs` stays green.

- [ ] **Step 4: Add the explicit byte-identity guard to `tests/features.rs`**

Find the block that reads the receipt (`let json = std::fs::read_to_string(&manifest).unwrap();`) and add, right after the existing `feature_rows` assertion:
```rust
    assert!(
        !json.contains("analyzer.feature_rows"),
        "features receipt must keep the un-prefixed feature_rows key: {json}"
    );
```

- [ ] **Step 5: Run + commit**

```bash
cargo test --test features 2>&1 | tail -5   # PASS
git add src/main.rs tests/features.rs
git commit -m "refactor(features): extract run_bounded_analysis (one analyzer-agnostic path)

run_features's contract machinery (plan/refuse/governor/drive/receipt) is now a
generic run_bounded_analysis(subcommand, param_prefix, analyzer, ...); features is
a thin wrapper with prefix=\"\" (byte-identical: FeatureAnalyzer.params() already
returns feature_rows). The receipt now merges analyzer.params() into the claim
under the prefix, with a MEASUREMENT_KEYS-collision guard. Sets up the analyze verb.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: a first-party `CoverageTrack`

**Files:** Modify `src/call/columnkit.rs` (after `FeatureAnalyzer`, ≈ line 121) + `src/lib.rs:77`

- [ ] **Step 1: Write the failing unit test**

In `src/call/columnkit.rs`'s `#[cfg(test)] mod tests`, add:
```rust
    #[test]
    fn coverage_track_header_and_params() {
        let t = CoverageTrack;
        assert_eq!(t.header().as_deref(), Some("#contig\tpos\tdepth\n"));
        assert_eq!(
            t.params().get("analyzer").map(String::as_str),
            Some("coverage")
        );
    }
```

- [ ] **Step 2: Run to verify it fails**

Run: `cargo test -p rosalind-bio --lib call::columnkit::tests::coverage_track 2>&1 | tail -8`
Expected: FAIL — compile error, `cannot find type CoverageTrack`.

- [ ] **Step 3: Add `CoverageTrack`** (in `src/call/columnkit.rs`, after the `FeatureAnalyzer` impl, ≈ line 121):
```rust
/// A per-locus coverage track: `(contig, 1-based pos, depth)`. The second first-party
/// [`ColumnAnalyzer`] — proof the SDK carries more than the feature egress. Its
/// `params()` names the analyzer so the receipt records which metric produced the run.
#[derive(Debug, Default)]
pub struct CoverageTrack;

impl ColumnAnalyzer for CoverageTrack {
    fn header(&self) -> Option<String> {
        Some("#contig\tpos\tdepth\n".to_string())
    }

    fn params(&self) -> BTreeMap<String, String> {
        BTreeMap::from([("analyzer".to_string(), "coverage".to_string())])
    }

    fn on_column(
        &mut self,
        col: &PileupColumn,
        contig: &str,
        out: &mut dyn Write,
    ) -> io::Result<()> {
        writeln!(out, "{contig}\t{}\t{}", col.locus.pos.0 + 1, col.depth())
    }
}
```

Re-export it in `src/lib.rs` (line 77):
```rust
pub use call::{run_bounded_whole_genome, ColumnAnalyzer, CoverageTrack, FeatureAnalyzer};
```

- [ ] **Step 4: Run to verify it passes**

Run: `cargo test -p rosalind-bio --lib call::columnkit::tests::coverage_track 2>&1 | tail -6`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/call/columnkit.rs src/lib.rs
git commit -m "feat(columnkit): first-party CoverageTrack analyzer

A pub ColumnAnalyzer emitting (contig, 1-based pos, depth); params() names the
analyzer. Re-exported as rosalind::call::CoverageTrack. The second registered
analyzer for the analyze verb.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 3: the `analyze <kind>` CLI verb

**Files:** Modify `src/main.rs` (`AnalyzerKind` enum + `Analyze` clap variant + dispatch); Test: `tests/analyze.rs` (create)

- [ ] **Step 1: Write the failing integration test** (`tests/analyze.rs`):
```rust
//! `rosalind analyze <kind>`: a registered ColumnAnalyzer with a verifiable receipt
//! whose claim records the analyzer's own params.

use std::env;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{SystemTime, UNIX_EPOCH};

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_rosalind")
}

fn tmpdir() -> PathBuf {
    static COUNTER: AtomicU64 = AtomicU64::new(0);
    let n = COUNTER.fetch_add(1, Ordering::Relaxed);
    let nanos = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_nanos();
    let d = env::temp_dir().join(format!("rosalind-analyze-{nanos}-{n}"));
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn run(args: &[&str]) -> std::process::Output {
    Command::new(bin()).args(args).output().expect("spawn rosalind")
}

fn build_index_and_sorted_bam(dir: &Path) -> (PathBuf, PathBuf) {
    let seq = "ACGTACGTACGTACGTACGTACGTACGTACGT";
    let fa = dir.join("ref.fa");
    std::fs::write(&fa, format!(">chr1\n{seq}\n")).unwrap();
    let fq = dir.join("reads.fq");
    let mut s = String::new();
    for (i, &start) in [0usize, 0, 8].iter().enumerate() {
        let read = &seq[start..start + 16];
        let qual: String = std::iter::repeat_n('I', 16).collect();
        s.push_str(&format!("@r{i}\n{read}\n+\n{qual}\n"));
    }
    std::fs::write(&fq, s).unwrap();
    let idx = dir.join("ref.idx");
    let raw = dir.join("raw.bam");
    let bam = dir.join("sorted.bam");
    assert!(run(&["index", "--reference", fa.to_str().unwrap(), "--output", idx.to_str().unwrap()]).status.success());
    assert!(run(&["align", "--reference", fa.to_str().unwrap(), "--reads", fq.to_str().unwrap(), "--format", "bam", "--output", raw.to_str().unwrap()]).status.success());
    assert!(run(&["sort", "--input", raw.to_str().unwrap(), "--output", bam.to_str().unwrap()]).status.success());
    (idx, bam)
}

#[test]
fn analyze_coverage_writes_a_verifiable_receipt_with_analyzer_params() {
    let d = tmpdir();
    let (idx, bam) = build_index_and_sorted_bam(&d);
    let tsv = d.join("cov.tsv");

    let out = run(&[
        "analyze", "coverage", "--index", idx.to_str().unwrap(),
        "--alignments", bam.to_str().unwrap(), "-o", tsv.to_str().unwrap(),
    ]);
    assert!(out.status.success(), "analyze coverage failed: {}", String::from_utf8_lossy(&out.stderr));

    // The output TSV carries the coverage header.
    let tsv_text = std::fs::read_to_string(&tsv).unwrap();
    assert!(tsv_text.starts_with("#contig\tpos\tdepth"), "tsv header: {tsv_text}");

    // The receipt's CLAIM records the analyzer's params under the analyzer. prefix.
    let manifest = format!("{}.manifest.json", tsv.display());
    let json = std::fs::read_to_string(&manifest).unwrap();
    assert!(
        json.contains("\"analyzer.analyzer\":\"coverage\""),
        "receipt missing analyzer.analyzer=coverage: {json}"
    );

    // And it verifies.
    let v = run(&["verify", "--manifest", &manifest]);
    assert!(
        v.status.success(),
        "verify must pass: {}",
        String::from_utf8_lossy(&v.stderr)
    );

    std::fs::remove_dir_all(&d).ok();
}
```

- [ ] **Step 2: Run to verify it fails**

Run: `cargo test --test analyze 2>&1 | tail -12`
Expected: FAIL — `clap` rejects the unknown `analyze` subcommand (non-zero exit; the assertion on `out.status.success()` fails).

- [ ] **Step 3: Add `AnalyzerKind`** (in `src/main.rs`, next to the `OutputFormat` enum):
```rust
#[derive(Copy, Clone, Debug, ValueEnum, Eq, PartialEq)]
enum AnalyzerKind {
    Features,
    Coverage,
}
```

- [ ] **Step 4: Add the `Analyze` clap variant** (in the `Commands` enum, after the `Features { … }` variant):
```rust
    /// Run a registered per-locus analyzer over the bounded whole-genome walk, with a
    /// verifiable receipt that records the analyzer's own params (under `analyzer.`).
    Analyze {
        /// Which analyzer to run.
        #[arg(value_enum)]
        kind: AnalyzerKind,
        /// Persisted index (`rosalind index`); analyzed over all contigs.
        #[arg(long)]
        index: PathBuf,
        /// Coordinate-sorted alignments (BAM).
        #[arg(long)]
        alignments: PathBuf,
        /// Minimum MAPQ required for a read to be considered.
        #[arg(long, default_value_t = 0)]
        mapq_threshold: u8,
        /// Declared memory budget (MiB). With `--enforce` it is honored (exit 3/4).
        #[arg(long)]
        memory_budget_mb: Option<u64>,
        /// Active-set depth cap (unbiased downsampling). `0` = uncapped.
        #[arg(long, default_value_t = 1000)]
        max_depth: u32,
        /// Max read length assumed by the `--enforce` estimate and enforced at ingest.
        #[arg(long, default_value_t = 250)]
        max_read_len: u32,
        /// Honor the budget: refuse up front (exit 3) / fail after (exit 4).
        #[arg(long, default_value_t = false)]
        enforce: bool,
        /// Output path (stdout if omitted).
        #[arg(short, long)]
        output: Option<PathBuf>,
        /// Where to write the reproducibility receipt (default: `<output>.manifest.json`).
        #[arg(long)]
        manifest: Option<PathBuf>,
    },
```

- [ ] **Step 5: Add the dispatch arm** (in the `match cli.command`, after the `Commands::Features { … } => run_features(…)` arm):
```rust
        Commands::Analyze {
            kind,
            index,
            alignments,
            mapq_threshold,
            memory_budget_mb,
            max_depth,
            max_read_len,
            enforce,
            output,
            manifest,
        } => {
            let label = match kind {
                AnalyzerKind::Features => "analyze features",
                AnalyzerKind::Coverage => "analyze coverage",
            };
            match kind {
                AnalyzerKind::Features => {
                    let mut a = rosalind::call::FeatureAnalyzer::default();
                    run_bounded_analysis(
                        label, "analyzer.", &mut a, index, alignments, mapq_threshold,
                        memory_budget_mb, max_depth, max_read_len, enforce, output, manifest,
                    )?
                }
                AnalyzerKind::Coverage => {
                    let mut a = rosalind::call::CoverageTrack;
                    run_bounded_analysis(
                        label, "analyzer.", &mut a, index, alignments, mapq_threshold,
                        memory_budget_mb, max_depth, max_read_len, enforce, output, manifest,
                    )?
                }
            }
        }
```

- [ ] **Step 6: Run to verify it passes**

Run: `cargo test --test analyze 2>&1 | tail -8`
Expected: PASS — `analyze coverage` writes a coverage TSV + a receipt with `analyzer.analyzer=coverage`, and `verify` exits 0.

- [ ] **Step 7: Commit**

```bash
git add src/main.rs tests/analyze.rs
git commit -m "feat(cli): rosalind analyze <kind> — verifiable receipt for any analyzer

A new analyze verb (kinds: features, coverage) routes a registered ColumnAnalyzer
through run_bounded_analysis with prefix=analyzer., merging the analyzer's own
params() into the content-addressed claim. A stranger's metric now inherits a
chainable, verifiable receipt with zero receipt plumbing.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 4: the `MEASUREMENT_KEYS`/prefix claim-protection guard

**Files:** Modify `crates/receipt/src/lib.rs` (`#[cfg(test)] mod tests`)

- [ ] **Step 1: Write the test**

In `crates/receipt/src/lib.rs`'s `#[cfg(test)] mod tests`, add:
```rust
    #[test]
    fn analyzer_prefixed_params_stay_in_the_claim() {
        // No measurement key starts with the analyzer. prefix, so an analyzer can never
        // shadow a measured field out of the claim.
        assert!(
            MEASUREMENT_KEYS.iter().all(|k| !k.starts_with("analyzer.")),
            "no MEASUREMENT_KEY may start with the analyzer. prefix"
        );
        // And a prefixed param survives finalize in the claim (params), not measurements.
        let mut m = RunManifest::new("analyze coverage");
        m.params
            .insert("analyzer.analyzer".to_string(), "coverage".to_string());
        m.finalize();
        assert_eq!(
            m.params.get("analyzer.analyzer").map(String::as_str),
            Some("coverage"),
            "an analyzer.-prefixed param must remain in the claim"
        );
        assert!(
            !m.measurements.contains_key("analyzer.analyzer"),
            "an analyzer.-prefixed param must NOT be relocated to measurements"
        );
    }
```

- [ ] **Step 2: Run to verify it passes**

Run: `cargo test -p rosalind-receipt analyzer_prefixed 2>&1 | tail -6`
Expected: PASS — the prefix provably protects the claim (it passes immediately; this is a guard codifying the invariant).

- [ ] **Step 3: Commit**

```bash
git add crates/receipt/src/lib.rs
git commit -m "test(receipt): analyzer.-prefixed params stay in the claim

Pins the invariant the analyze verb relies on: no MEASUREMENT_KEY starts with
analyzer., so finalize never relocates an analyzer param out of the claim.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 5: full gate + e2e smoke

**Files:** none (verification only)

- [ ] **Step 1: Full suite + lint + format**

Run:
```bash
cargo test
cargo clippy --all-targets -- -D warnings
cargo fmt --check
```
Expected: PASS on all three. (If `cargo fmt --check` reports diffs, run `cargo fmt` and amend the relevant commit.)

- [ ] **Step 2: End-to-end smoke**

Run:
```bash
BIN=target/debug/rosalind; D=$(mktemp -d)
printf '>chr1\nACGTACGTACGTACGTACGTACGTACGTACGT\n' > "$D/ref.fa"
printf '@r0\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n' > "$D/reads.fq"
$BIN index --reference "$D/ref.fa" --output "$D/ref.idx" >/dev/null 2>&1
$BIN align --reference "$D/ref.fa" --reads "$D/reads.fq" --format bam --output "$D/raw.bam" >/dev/null 2>&1
$BIN sort --input "$D/raw.bam" --output "$D/s.bam" >/dev/null 2>&1
$BIN analyze coverage --index "$D/ref.idx" --alignments "$D/s.bam" -o "$D/cov.tsv" 2>/dev/null
echo "--- tsv head ---"; head -3 "$D/cov.tsv"
echo "--- receipt analyzer param ---"; grep -o '"analyzer.analyzer":"coverage"' "$D/cov.tsv.manifest.json"
echo "--- verify ---"; $BIN verify --manifest "$D/cov.tsv.manifest.json" | grep "verify: OK"
rm -rf "$D"
```
Expected: the TSV starts with `#contig\tpos\tdepth`; the receipt grep prints `"analyzer.analyzer":"coverage"`; `verify` prints `verify: OK`.

- [ ] **Step 3: No commit** (verification only).

---

## Self-Review

**1. Spec coverage:**
- §3 extract `run_bounded_analysis` (verbatim + 2 logical changes: analyzer-param, params-merge; + the trailing-summary generalization) → **Task 1** (Edits 1–8).
- §3.1 three call sites (`features` prefix="", `analyze features`/`analyze coverage` prefix="analyzer.") → **Task 1** (wrapper) + **Task 3** (dispatch).
- §4 `CoverageTrack` + re-export → **Task 2**.
- §5 CLI (`Analyze` variant + `AnalyzerKind` + dispatch) → **Task 3**.
- §6 tests: features byte-identical (Task 1 Step 3/4), analyze coverage e2e (Task 3), no-collision (Task 4), SDK equivalence unchanged (Task 1 leaves `columnkit.rs` test untouched) → covered.
- §7 non-goals: no dynamic loading (enum registry), don't touch the example (untouched), pure extraction (Task 1) → respected.

**2. Placeholder scan:** No TBD/TODO. Every step has literal code or an exact command + expected output.

**3. Type consistency:** `run_bounded_analysis(subcommand: &str, param_prefix: &str, analyzer: &mut dyn rosalind::call::ColumnAnalyzer, …)` — the signature in Task 1 matches every call site (the `run_features` wrapper in Task 1, the two dispatch arms in Task 3). `CoverageTrack` (Task 2) is referenced by `rosalind::call::CoverageTrack` in Task 3's dispatch and re-exported in Task 2. `AnalyzerKind { Features, Coverage }` (Task 3 Step 3) matches the dispatch match (Step 5) and the clap variant (Step 4). `analyzer.params()` returns `BTreeMap<String, String>` (the trait method), iterated in Edit 6. `MEASUREMENT_KEYS` is `rosalind::provenance::MEASUREMENT_KEYS` (Edit 6) and `MEASUREMENT_KEYS` within the receipt crate (Task 4).
