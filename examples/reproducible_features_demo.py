#!/usr/bin/env python3
"""Reproducible-features demo: train a model on `rosalind features` output and
prove the inputs (and therefore the trained model) are bit-reproducible.

The headline is REPRODUCIBILITY, not biological novelty: identical inputs ->
byte-identical features -> bit-identical trained weights, verifiable by the
BLAKE3 receipt. Runs with numpy only (no pandas/sklearn/pyarrow needed).

    python examples/reproducible_features_demo.py [path-to-rosalind-binary]
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "python"))
from rosalind import features  # noqa: E402


def find_binary(argv: list[str]) -> str:
    if len(argv) > 1:
        return argv[1]
    for cand in (os.path.join(REPO, "target", "release", "rosalind"),
                 os.path.join(REPO, "target", "debug", "rosalind")):
        if os.path.exists(cand):
            return cand
    found = shutil.which("rosalind")
    if found:
        return found
    sys.exit("rosalind binary not found — pass it as argv[1] or `cargo build --release`")


def sh(cmd: list[str]) -> None:
    subprocess.run(cmd, check=True)


def output_blake3(manifest_path: str) -> str:
    """Read the recorded BLAKE3 of the feature TSV from the run receipt."""
    with open(manifest_path) as fh:
        m = json.load(fh)
    return m["outputs"][0]["blake3"]


def train_logreg(x: np.ndarray, y: np.ndarray, iters: int = 300, lr: float = 0.2):
    """Deterministic pure-numpy logistic regression (zero init, fixed schedule)."""
    mu = x.mean(0)
    sd = x.std(0)
    sd[sd == 0.0] = 1.0
    xs = (x - mu) / sd
    xb = np.hstack([np.ones((xs.shape[0], 1)), xs])
    w = np.zeros(xb.shape[1])
    n = xb.shape[0]
    for _ in range(iters):
        p = 1.0 / (1.0 + np.exp(-(xb @ w)))
        w -= lr * (xb.T @ (p - y)) / n
    return w, mu, sd


def predict(w, mu, sd, x: np.ndarray) -> np.ndarray:
    sd = sd.copy()
    sd[sd == 0.0] = 1.0
    xs = (x - mu) / sd
    xb = np.hstack([np.ones((xs.shape[0], 1)), xs])
    return (1.0 / (1.0 + np.exp(-(xb @ w)))) >= 0.5


def main() -> int:
    binary = find_binary(sys.argv)
    work = tempfile.mkdtemp(prefix="rosalind-demo-")
    data = os.path.join(work, "toy")

    print("== building a toy dataset + index ==")
    sh([sys.executable, os.path.join(REPO, "scripts", "generate_toy_data.py"), data])
    ref = os.path.join(data, "reference.fa")
    reads = os.path.join(data, "reads_R1.fastq")
    raw = os.path.join(work, "raw.bam")
    bam = os.path.join(work, "sorted.bam")
    idx = os.path.join(work, "ref.idx")
    sh([binary, "align", "--reference", ref, "--reads", reads, "--format", "bam", "--output", raw])
    sh([binary, "sort", "--input", raw, "--output", bam])
    sh([binary, "index", "--reference", ref, "--output", idx])

    print("== extracting features TWICE (independent runs) ==")
    a = features(idx, bam, binary=binary, workdir=os.path.join(work, "a"))
    b = features(idx, bam, binary=binary, workdir=os.path.join(work, "b"))

    # --- Reproducibility proof -------------------------------------------------
    tsv_a = os.path.join(work, "a", "features.tsv")
    tsv_b = os.path.join(work, "b", "features.tsv")
    bytes_identical = open(tsv_a, "rb").read() == open(tsv_b, "rb").read()
    h_a, h_b = output_blake3(a.manifest_path), output_blake3(b.manifest_path)
    local_hash = hashlib.sha256(open(tsv_a, "rb").read()).hexdigest()[:16]

    print(f"  rows: {len(a)};  numeric features: {a.data.shape[1]}")
    print(f"  features byte-identical across runs: {bytes_identical}")
    print(f"  receipt BLAKE3 (run A): {h_a}")
    print(f"  receipt BLAKE3 (run B): {h_b}")
    print(f"  receipts match: {h_a == h_b}")
    assert bytes_identical, "feature TSVs differ across runs"
    assert h_a == h_b, "receipt hashes differ across runs"

    # --- A genuine supervised task: is the reference base a purine (A/G)? -------
    # Read counts peak at the ref base, so this is learnable from the features.
    def label(ft) -> np.ndarray:
        return np.isin(ft.ref, [b"A", b"G"]).astype(np.float64)

    # Deterministic train/test split by position parity.
    def split(ft):
        test = (ft.pos % 2) == 1
        return ~test, test

    ya, yb = label(a), label(b)
    tr_a, te_a = split(a)
    w_a, mu_a, sd_a = train_logreg(a.data[tr_a], ya[tr_a])
    acc_a = float((predict(w_a, mu_a, sd_a, a.data[te_a]) == (ya[te_a] >= 0.5)).mean())

    w_b, mu_b, sd_b = train_logreg(b.data[split(b)[0]], yb[split(b)[0]])
    pred_a = predict(w_a, mu_a, sd_a, a.data[te_a])
    pred_b = predict(w_b, mu_b, sd_b, b.data[split(b)[1]])

    weights_identical = np.array_equal(w_a, w_b)
    preds_identical = np.array_equal(pred_a, pred_b)

    print("== model: predict is-purine(ref) from pileup features (logistic regression) ==")
    print(f"  held-out accuracy: {acc_a:.4f}")
    print(f"  trained weights bit-identical across the two extractions: {weights_identical}")
    print(f"  test predictions bit-identical across the two extractions: {preds_identical}")
    assert weights_identical, "trained weights differ across extractions (not reproducible)"
    assert preds_identical, "predictions differ across extractions"

    print("\nPROOF: byte-identical features (hash %s) -> bit-identical trained model."
          " Rosalind feature extraction yields bit-reproducible ML training inputs." % local_hash)
    shutil.rmtree(work, ignore_errors=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
