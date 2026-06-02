"""Dependency-light Python boundary for the Rosalind feature substrate.

Runs ``rosalind features`` (the bounded, deterministic per-locus feature stream)
and loads the result into numpy — meeting Python ML pipelines where they live,
with no pyO3 / pyarrow build required (stdlib + numpy only).

    from rosalind import features
    ft = features("ref.idx", "sorted.bam", binary="target/release/rosalind")
    X = ft.data            # (n_loci, n_numeric_features) float64
    y = ft.ref == b"A"     # per-locus labels you derive

The binary streams the table to disk in bounded memory; this loader is the only
part that holds it in RAM. A BLAKE3 receipt is written next to the TSV (``manifest_path``)
so feature extractions are hash-verifiable — identical inputs → byte-identical
features → bit-reproducible training inputs.

The zero-copy in-process binding (pyO3 + pyarrow RecordBatches) is a planned
follow-up; this subprocess boundary is the dependency-light entry point today.
"""

from __future__ import annotations

import os
import subprocess
import tempfile
from dataclasses import dataclass

import numpy as np

# Columns 0..2 are categorical (contig, pos, ref); 3.. are numeric features.
_N_LEADING = 3


@dataclass
class FeatureTable:
    """A per-locus feature matrix loaded from ``rosalind features``."""

    columns: list[str]          # names of the numeric feature columns (len == data.shape[1])
    data: np.ndarray            # (n_loci, n_numeric) float64
    contig: np.ndarray          # (n_loci,) bytes — contig name per locus
    pos: np.ndarray             # (n_loci,) int64 — 1-based position
    ref: np.ndarray             # (n_loci,) bytes — reference base (length-1)
    manifest_path: str          # the BLAKE3 reproducibility receipt for this extraction

    def __len__(self) -> int:
        return self.data.shape[0]


def features(
    index: str,
    alignments: str,
    *,
    binary: str = "rosalind",
    max_depth: int = 1000,
    max_read_len: int = 250,
    mapq: int = 0,
    memory_budget_mb: int | None = None,
    enforce: bool = False,
    workdir: str | None = None,
) -> FeatureTable:
    """Run ``rosalind features`` and load the per-locus table into numpy.

    ``index`` is a ``rosalind index`` artifact; ``alignments`` a coordinate-sorted
    BAM. With ``enforce=True`` and a ``memory_budget_mb`` the run honors the memory
    contract (a breach raises ``CalledProcessError``). Returns a :class:`FeatureTable`.
    """
    tmp = workdir or tempfile.mkdtemp(prefix="rosalind-features-")
    os.makedirs(tmp, exist_ok=True)
    tsv = os.path.join(tmp, "features.tsv")
    manifest = tsv + ".manifest.json"
    cmd = [
        binary, "features",
        "--index", str(index),
        "--alignments", str(alignments),
        "--max-depth", str(max_depth),
        "--max-read-len", str(max_read_len),
        "--mapq-threshold", str(mapq),
        "-o", tsv,
        "--manifest", manifest,
    ]
    if memory_budget_mb is not None:
        cmd += ["--memory-budget-mb", str(memory_budget_mb)]
    if enforce:
        cmd += ["--enforce"]
    subprocess.run(cmd, check=True)
    return load_feature_tsv(tsv, manifest_path=manifest)


def load_feature_tsv(path: str, *, manifest_path: str = "") -> FeatureTable:
    """Parse a feature TSV emitted by ``rosalind features`` into a FeatureTable."""
    with open(path, "r", encoding="ascii") as fh:
        header = fh.readline().lstrip("#").rstrip("\n").split("\t")
        rows = [line.rstrip("\n").split("\t") for line in fh if line.strip()]
    if not rows:
        empty = np.empty((0, len(header) - _N_LEADING), dtype=np.float64)
        return FeatureTable(header[_N_LEADING:], empty, np.array([]), np.array([], dtype=np.int64), np.array([]), manifest_path)
    contig = np.array([r[0] for r in rows], dtype="S32")
    pos = np.array([int(r[1]) for r in rows], dtype=np.int64)
    ref = np.array([r[2] for r in rows], dtype="S1")
    data = np.array([r[_N_LEADING:] for r in rows], dtype=np.float64)
    return FeatureTable(header[_N_LEADING:], data, contig, pos, ref, manifest_path)


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 4:
        print("usage: python rosalind.py <binary> <index> <sorted.bam>", file=sys.stderr)
        raise SystemExit(2)
    ft = features(sys.argv[2], sys.argv[3], binary=sys.argv[1])
    print(f"loaded {len(ft)} loci x {ft.data.shape[1]} numeric features")
    print("columns:", ft.columns)
    print("receipt:", ft.manifest_path)
