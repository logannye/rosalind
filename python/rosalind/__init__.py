"""Typed Python API for the Rosalind binary analysis contract."""

from importlib.metadata import PackageNotFoundError, version

from .features import FeatureRun, RunResult, collect_features, features, iter_features
from .receipts import Receipt, inspect_receipt
from .evidence import EvidenceProcessError, EvidenceResult, EvidenceRun, iter_evidence, materialize_evidence, panel_qc

try:
    __version__ = version("rosalind-bio")
except PackageNotFoundError:
    # Source-tree fallback; release wheels derive this value from Cargo metadata.
    __version__ = "0.5.0"

__all__ = [
    "EvidenceProcessError",
    "EvidenceResult",
    "EvidenceRun",
    "iter_evidence",
    "materialize_evidence",
    "panel_qc",
    "FeatureRun",
    "Receipt",
    "RunResult",
    "collect_features",
    "features",
    "inspect_receipt",
    "iter_features",
]
