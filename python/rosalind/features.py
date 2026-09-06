"""Bounded native Arrow feature iteration through the bundled Rosalind CLI."""

from __future__ import annotations

import shutil
import re
import subprocess
import sysconfig
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Optional, Union

try:
    import pyarrow as pa
except ImportError:  # Allows receipt-only use from a source checkout.
    pa = None  # type: ignore[assignment]


def _require_pyarrow():
    if pa is None:
        raise RuntimeError("feature iteration requires pyarrow>=15")
    return pa


@dataclass(frozen=True)
class RunResult:
    """Completed process and receipt information."""

    manifest_path: Path
    returncode: int


def _package_version() -> str:
    from . import __version__

    return __version__


def _normalized_version(version: str) -> str:
    """Map Cargo's numbered RC spelling to the corresponding wheel version."""
    return re.sub(r"^(\d+\.\d+\.\d+)-rc\.(\d+)$", r"\1rc\2", version)


def _resolve_binary(binary: Optional[Union[str, Path]], allow_version_mismatch: bool) -> Path:
    if binary is None:
        candidate = Path(sysconfig.get_path("scripts")) / "rosalind"
        if not candidate.is_file():
            located = shutil.which("rosalind")
            if located is None:
                raise FileNotFoundError("the bundled rosalind executable is not available")
            candidate = Path(located)
    else:
        candidate = Path(binary)
        if not candidate.is_file():
            located = shutil.which(str(binary))
            if located is None:
                raise FileNotFoundError(f"rosalind executable not found: {binary}")
            candidate = Path(located)
    version = subprocess.run(
        [str(candidate), "--version"], capture_output=True, text=True, check=True
    ).stdout.strip().split()[-1]
    if not allow_version_mismatch and _normalized_version(version) != _normalized_version(_package_version()):
        raise RuntimeError(
            f"Python package {_package_version()} does not match Rosalind binary {version}"
        )
    return candidate


class FeatureRun:
    """One bounded Arrow IPC subprocess and its batch iterator."""

    def __init__(self, process: subprocess.Popen[bytes], manifest_path: Path):
        self._process = process
        self._reader: Optional[pa.RecordBatchReader] = None
        self.manifest_path = manifest_path
        self.result: Optional[RunResult] = None
        self._closed = False

    def batches(self) -> Iterator[pa.RecordBatch]:
        """Yield native fixed-size record batches without collecting the stream."""
        arrow = _require_pyarrow()
        if self._reader is not None or self._closed:
            raise RuntimeError("FeatureRun batches can only be consumed once")
        assert self._process.stdout is not None
        try:
            self._reader = arrow.ipc.open_stream(self._process.stdout)
            for batch in self._reader:
                yield batch
            returncode = self._process.wait()
            self.result = RunResult(self.manifest_path, returncode)
            if returncode:
                raise subprocess.CalledProcessError(returncode, self._process.args)
        finally:
            self.close()

    def __iter__(self) -> Iterator[pa.RecordBatch]:
        return self.batches()

    def close(self) -> None:
        """Terminate an incompletely consumed feature process."""
        if self._closed:
            return
        if self._reader is not None:
            self._reader.close()
        if self._process.stdout is not None:
            self._process.stdout.close()
        if self._process.poll() is None:
            self._process.terminate()
        self._process.wait()
        self._closed = True

    def __enter__(self) -> "FeatureRun":
        return self

    def __exit__(self, *_: object) -> None:
        self.close()


def iter_features(
    reference: Union[str, Path],
    alignments: Union[str, Path],
    *,
    binary: Optional[Union[str, Path]] = None,
    allow_version_mismatch: bool = False,
    max_depth: int = 1000,
    max_read_len: int = 250,
    mapq: int = 0,
    memory_budget_mb: Optional[int] = None,
    enforce: bool = False,
    region: Optional[str] = None,
    regions: Optional[Union[str, Path]] = None,
    shard_count: Optional[int] = None,
    shard_index: Optional[int] = None,
    workdir: Optional[Union[str, Path]] = None,
) -> FeatureRun:
    """Start native Arrow feature extraction and return a lazy batch stream."""
    executable = _resolve_binary(binary, allow_version_mismatch)
    directory = Path(workdir) if workdir is not None else Path(tempfile.mkdtemp(prefix="rosalind-features-"))
    directory.mkdir(parents=True, exist_ok=True)
    manifest = directory / "features.manifest.json"
    reference = Path(reference)
    reference_flag = "--reference-pack" if reference.suffix == ".rref" else "--index"
    command = [
        str(executable),
        "features",
        reference_flag,
        str(reference),
        "--alignments",
        str(alignments),
        "--format",
        "arrow-ipc",
        "--max-depth",
        str(max_depth),
        "--max-read-len",
        str(max_read_len),
        "--mapq-threshold",
        str(mapq),
        "--manifest",
        str(manifest),
    ]
    selections = sum(
        value is not None
        for value in (region, regions, shard_count)
    )
    if selections > 1 or (shard_count is None) != (shard_index is None):
        raise ValueError("region, regions, and a complete shard pair are mutually exclusive")
    if region is not None:
        command += ["--region", region]
    if regions is not None:
        command += ["--regions", str(regions)]
    if shard_count is not None and shard_index is not None:
        command += ["--shard-count", str(shard_count), "--shard-index", str(shard_index)]
    if memory_budget_mb is not None:
        command += ["--memory-budget-mb", str(memory_budget_mb)]
    if enforce:
        command.append("--enforce")
    process: subprocess.Popen[bytes] = subprocess.Popen(command, stdout=subprocess.PIPE)
    return FeatureRun(process, manifest)


def collect_features(*args: object, **kwargs: object) -> pa.Table:
    """Explicitly materialize a complete feature stream as a PyArrow table."""
    arrow = _require_pyarrow()
    run = iter_features(*args, **kwargs)
    batches = list(run.batches())
    if batches:
        return arrow.Table.from_batches(batches)
    assert run._reader is not None
    return arrow.Table.from_batches([], schema=run._reader.schema)


def features(*args: object, **kwargs: object):
    """Compatibility alias for explicit collection; prefer ``collect_features``."""
    return collect_features(*args, **kwargs)
