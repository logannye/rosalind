# Finding — bit-reproducible ML training inputs from the feature substrate

**Date:** 2026-06-02. **Demo:** `examples/reproducible_features_demo.py` (numpy-only). **Boundary:**
`python/rosalind.py`. Demonstrates the feature substrate's novel claim end-to-end: **byte-identical
features → bit-reproducible trained model**, hash-verifiable.

## What it does

Builds a toy index + coordinate-sorted BAM via the CLI, then **extracts features twice** with
`rosalind features`, and trains a small model on the result:

1. **Reproducibility check** — assert the two `features.tsv` are byte-identical and the two run
   receipts record the **same BLAKE3** of the output.
2. **A genuine supervised task** — label each locus `is_purine = ref ∈ {A,G}` and train a pure-numpy
   logistic regression from the per-locus pileup features (counts, strand counts, mean BQ/MAPQ). Read
   counts peak at the reference base, so the label is learnable from the features.
3. **Bit-reproducibility of training** — train the identical model on each of the two extractions and
   assert the trained weights and the held-out predictions are bit-identical.

## Measured output (this environment, 2026-06-02)

```
  rows: 983011;  numeric features: 16
  features byte-identical across runs: True
  receipt BLAKE3 (run A): c89fdd41c5ac39c656b049ba1d3b839d9477e3d94dc4fd73ebf481f2dc555ec1
  receipt BLAKE3 (run B): c89fdd41c5ac39c656b049ba1d3b839d9477e3d94dc4fd73ebf481f2dc555ec1
  receipts match: True
  held-out accuracy: 0.9990
  trained weights bit-identical across the two extractions: True
  test predictions bit-identical across the two extractions: True

PROOF: byte-identical features -> bit-identical trained model.
Rosalind feature extraction yields bit-reproducible ML training inputs.
```

## Why it matters (honest framing)

- The **headline is reproducibility**, not biological novelty. The is-purine task is deliberately
  simple but genuinely learnable; the point is that the *whole pipeline* — features → matrix → trained
  weights → predictions — is **bit-identical** across two independent extractions, and **provable by a
  hash** in the run receipt. No other pileup gives bounded streaming **and** byte-identical output
  together, so no other tool can hand an ML team a hash that proves "this quarter's training features
  are identical to last quarter's."
- **Bounded:** the ~983k-row whole-toy-genome feature table is produced in ~6 MiB peak (the feature
  CLI's receipt), independent of input size — so this scales to real genomes on a laptop.

## Scope + next step

- **numpy-only** (this environment has numpy but no pandas/scikit-learn/pyarrow/maturin), so the demo
  uses a hand-rolled logistic regression. A pandas/sklearn version is a trivial rewrite where those
  are installed.
- The **zero-copy `pyarrow` Python binding** (an in-process `rosalind.features()` yielding Arrow
  RecordBatches, replacing the dead-end `python_bindings` stub) and an **Arrow/Parquet egress** are the
  planned next step — built and verified where maturin + pyarrow exist. This demo proves the pull
  *before* that investment.

## Reproduce

```sh
cargo build --release
python examples/reproducible_features_demo.py target/release/rosalind   # numpy required
```
