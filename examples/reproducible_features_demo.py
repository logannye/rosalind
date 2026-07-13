#!/usr/bin/env python3
"""Reproducible-features demo: train a model on native Arrow feature output and
prove the inputs (and therefore the trained model) are bit-reproducible.

The headline is REPRODUCIBILITY, not biological novelty: identical inputs ->
byte-identical Arrow streams -> bit-identical trained weights, verifiable by the
BLAKE3 receipt. Runs with NumPy and PyArrow (no pandas or sklearn needed).

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
import pyarrow as pa

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

NUMERIC_FEATURES = (
    "depth", "raw_depth", "a", "c", "g", "t",
    "a_fwd", "a_rev", "c_fwd", "c_rev",
    "g_fwd", "g_rev", "t_fwd", "t_rev",
    "mean_bq", "mean_mapq",
)


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
    """Read the recorded BLAKE3 of the Arrow stream from the run receipt."""
    with open(manifest_path) as fh:
        m = json.load(fh)
    return m["outputs"][0]["blake3"]


def extract_features(binary: str, reference: str, bam: str, directory: str):
    """Write one native Arrow artifact and return its table and evidence paths."""
    os.makedirs(directory, exist_ok=True)
    artifact = os.path.join(directory, "features.arrow")
    sh([
        binary, "features", "--reference-pack", reference,
        "--alignments", bam, "--format", "arrow-ipc", "--output", artifact,
    ])
    with pa.memory_map(artifact, "r") as source:
        table = pa.ipc.open_stream(source).read_all()
    return table, artifact, f"{artifact}.manifest.json"


def model_arrays(table: pa.Table):
    """Convert only the explicit training columns after bounded extraction."""
    columns = [
        table[name].to_numpy(zero_copy_only=False).astype(np.float64, copy=False)
        for name in NUMERIC_FEATURES
    ]
    values = np.column_stack(columns)
    positions = table["pos"].to_numpy(zero_copy_only=False)
    references = np.asarray(table["ref"].to_pylist(), dtype="U1")
    return values, positions, references


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

    print("== building a toy dataset + analysis reference pack ==")
    sh([sys.executable, os.path.join(REPO, "scripts", "generate_toy_data.py"), data])
    ref = os.path.join(data, "reference.fa")
    reads = os.path.join(data, "reads_R1.fastq")
    raw = os.path.join(work, "raw.bam")
    bam = os.path.join(work, "sorted.bam")
    reference_pack = os.path.join(work, "ref.rref")
    sh([binary, "align", "--reference", ref, "--reads", reads, "--format", "bam", "--output", raw])
    sh([binary, "sort", "--input", raw, "--output", bam])
    sh([binary, "reference", "build", "--fasta", ref, "--output", reference_pack])

    print("== extracting native Arrow features TWICE (independent runs) ==")
    a, arrow_a, manifest_a = extract_features(
        binary, reference_pack, bam, os.path.join(work, "a")
    )
    b, arrow_b, manifest_b = extract_features(
        binary, reference_pack, bam, os.path.join(work, "b")
    )
    data_a, pos_a, ref_a = model_arrays(a)
    data_b, pos_b, ref_b = model_arrays(b)

    # --- Reproducibility proof -------------------------------------------------
    with open(arrow_a, "rb") as first, open(arrow_b, "rb") as second:
        bytes_a, bytes_b = first.read(), second.read()
    bytes_identical = bytes_a == bytes_b
    h_a, h_b = output_blake3(manifest_a), output_blake3(manifest_b)
    local_hash = hashlib.sha256(bytes_a).hexdigest()[:16]

    print(f"  rows: {a.num_rows};  numeric features: {data_a.shape[1]}")
    print(f"  features byte-identical across runs: {bytes_identical}")
    print(f"  receipt BLAKE3 (run A): {h_a}")
    print(f"  receipt BLAKE3 (run B): {h_b}")
    print(f"  receipts match: {h_a == h_b}")
    assert bytes_identical, "native Arrow streams differ across runs"
    assert h_a == h_b, "receipt hashes differ across runs"

    # --- A genuine supervised task: is the reference base a purine (A/G)? -------
    # Read counts peak at the ref base, so this is learnable from the features.
    def label(references: np.ndarray) -> np.ndarray:
        return np.isin(references, ["A", "G"]).astype(np.float64)

    # Deterministic train/test split by position parity.
    def split(positions: np.ndarray):
        test = (positions % 2) == 1
        return ~test, test

    ya, yb = label(ref_a), label(ref_b)
    tr_a, te_a = split(pos_a)
    w_a, mu_a, sd_a = train_logreg(data_a[tr_a], ya[tr_a])
    acc_a = float((predict(w_a, mu_a, sd_a, data_a[te_a]) == (ya[te_a] >= 0.5)).mean())

    tr_b, te_b = split(pos_b)
    w_b, mu_b, sd_b = train_logreg(data_b[tr_b], yb[tr_b])
    pred_a = predict(w_a, mu_a, sd_a, data_a[te_a])
    pred_b = predict(w_b, mu_b, sd_b, data_b[te_b])

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
