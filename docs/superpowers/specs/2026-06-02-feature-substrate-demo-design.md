# Feature Substrate — Python Boundary + Reproducibility Demo (design)

**Status:** DESIGN SPEC — 2026-06-02. **Branch:** `rosalind/feature-substrate-demo` (off `main`
`e4f43bf`). Follow-up to the feature-substrate egress (PR #28). Goal: **prove the pull** — show an ML
builder can train a model on `rosalind features` output and that the inputs (and therefore the trained
model) are **bit-reproducible**, verifiable by a hash.

## 1. Environment-driven scope (measure-first)

This environment has **numpy** but **no pandas / scikit-learn / pyarrow / maturin**. So a pyO3/pyarrow
in-process binding is **not verifiable here**. This increment delivers the pieces that are both
higher-value-first *and* runnable+verified in-env:

- A **dependency-light Python boundary** (`python/rosalind.py`, stdlib + numpy): a subprocess wrapper
  that runs `rosalind features` and loads the TSV into a numpy array. Meets Python users without pyO3.
- A **numpy-only reference-model demo** that trains a real classifier on the feature matrix and proves
  bit-reproducibility end-to-end.

**Deferred (explicit next step):** the zero-copy pyO3/pyarrow binding (`rosalind.features()` yielding
Arrow RecordBatches) + an Arrow/Parquet egress — built and verified where maturin + pyarrow exist. The
demo proves the pull *before* that investment (the audit's own sequencing).

## 2. Deliverables

### 2a. `python/rosalind.py` — the boundary

`features(index, alignments, *, binary="rosalind", max_depth=1000, max_read_len=250, mapq=0,
memory_budget_mb=None, enforce=False, workdir=None) -> FeatureTable` where `FeatureTable` carries
`columns: list[str]`, `data: np.ndarray` (float64, the numeric columns), `ref: np.ndarray` (the ref
base per row), `contig: np.ndarray`, `pos: np.ndarray`, and `manifest_path`. It runs the binary to a
temp TSV (bounded — the binary streams), parses the header + rows with numpy
(`np.genfromtxt`/manual), and returns the matrix. Pure stdlib + numpy. A `__main__` smoke prints the
shape + column names.

### 2b. `examples/reproducible_features_demo.py` — the demo

Self-contained, numpy-only:
1. Locate the `rosalind` binary (argv[1] or `target/release/rosalind` or PATH); generate the bundled
   `illumina_toy` data if absent; `index` the reference and `sort` the BAM via the CLI.
2. Extract features **twice** to `features_a.tsv` + `features_b.tsv` (with receipts).
3. **Reproducibility proof:** assert the two TSVs are byte-identical; read both manifests and assert
   the recorded **output BLAKE3 hashes match** (the verifiable receipt).
4. **A real supervised task:** label each locus `is_purine = ref in {A,G}`; features = the numeric
   columns (counts, strand counts, mean BQ/MAPQ). Train a **pure-numpy logistic regression**
   (zero-init weights, fixed iterations + learning rate, deterministic standardization from the train
   split, deterministic train/test split by row parity) on extraction A; report test accuracy.
5. **Bit-reproducibility of training:** train the identical model on extraction B; assert the trained
   weight vectors and the test predictions are **bit-identical** (`np.array_equal`) to extraction A's.
6. Print a clear proof block (hashes match; model bit-identical; accuracy). Exit non-zero on any
   assertion failure (so it doubles as a check).

### 2c. Findings doc + README

- `docs/findings/2026-06-02-reproducible-features-demo.md` — the demo's **actual captured output**
  (the matching hashes, the model accuracy, the bit-identical-across-extractions result).
- README: a short pointer under the feature-substrate section to `python/rosalind.py` + the demo.

## 3. Honest notes

- The model task (is-purine-from-counts) is deliberately simple but **genuinely learnable** (read
  counts peak at the ref base) — the headline is **reproducibility**, not biological novelty; stated
  plainly.
- numpy float ops are deterministic for identical inputs + identical operations, so identical feature
  matrices → bit-identical weights. The demo *demonstrates* this rather than assuming it.
- No new Rust code; no new Rust deps. The Rust test suite is unaffected (the demo is a Python script).

## 4. Self-review

- **Coverage:** boundary (2a), demo (2b), findings + README (2c). ✓
- **Verifiable in-env:** numpy-only → the demo runs here and the findings record real output. ✓
- **No placeholders:** the model + reproducibility checks are concrete; pyO3/pyarrow/Arrow explicitly
  deferred with the reason (unverifiable here). ✓
- **Scope:** prove-the-pull + a Python entry point; the efficient binding is the next step. ✓
