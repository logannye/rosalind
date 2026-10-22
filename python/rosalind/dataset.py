"""Portable evidence queries and Parquet export through the verified native reader.

The dataset owns scientific filter and sample settings. A query can project fields
or select covered loci; it cannot infer evidence discarded during extraction.
Native budgets exclude Arrow batches and tables retained by Python consumers.
"""
from __future__ import annotations

import json
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence

from .evidence import EvidenceProcessError, EvidenceResult, EvidenceRun, PathLike
from .features import _resolve_binary


@dataclass(frozen=True)
class ParquetExport:
    """A committed directory of uint64-preserving Parquet parts and its receipt."""

    directory: Path
    manifest_path: Path


@dataclass(frozen=True)
class EvidenceDataset:
    """A lazy handle; every native invocation verifies metadata and consumed parts.

    Original alignments and reference files need not be available. Their recorded
    identities remain provenance; offline queries do not claim to rehash them.
    Keep the dataset directory immutable throughout each invocation.
    """

    manifest_path: Path
    binary: Optional[PathLike] = None
    allow_version_mismatch: bool = False
    max_dataset_metadata_bytes: int = 33_554_432

    def _command(self, action: str, *, memory_budget_mb: Optional[int] = None) -> list[str]:
        command = [str(_resolve_binary(self.binary, self.allow_version_mismatch)), "dataset", action, "--dataset", str(self.manifest_path)]
        command += ["--max-dataset-metadata-bytes", str(self.max_dataset_metadata_bytes)]
        if memory_budget_mb is not None:
            if isinstance(memory_budget_mb, bool) or memory_budget_mb <= 0:
                raise ValueError("memory_budget_mb must be a positive integer")
            command += ["--memory-budget-mb", str(memory_budget_mb)]
        return command

    def inspect(self, *, memory_budget_mb: Optional[int] = None) -> dict:
        """Verify and return metadata; does not scan evidence partition bodies."""
        return _json(self._command("inspect", memory_budget_mb=memory_budget_mb))

    def verify(self, *, memory_budget_mb: Optional[int] = None) -> dict:
        """Verify every stored evidence row and partition, including zero-depth rows."""
        return _json(self._command("verify", memory_budget_mb=memory_budget_mb))

    def coverage(self, *, regions: Optional[PathLike] = None, sites: Optional[PathLike] = None,
                 fields: Optional[Sequence[str]] = None, memory_budget_mb: Optional[int] = None,
                 max_output_bytes: int = 4_194_304) -> dict:
        """Split normalized requested loci into covered and missing metadata.

        The ``rosalind-dataset-coverage-v1`` response contains the canonical query,
        its native digest, parent/descriptor byte identities, coverage counts and
        selections. Coordinates are zero-based; SNV bases are ASCII byte values.
        This verifies metadata only: partition bodies and requested REF alleles
        remain unchecked until consumed. It never reads original source files.

        Responses above ``max_output_bytes`` are refused before JSON is returned;
        callers must submit smaller selection windows. The maximum is 32 MiB.
        This is an output bound, not a maximum process RSS. Strict saved queries
        still refuse missing loci; partial consumers must query the intersection.
        """
        if (isinstance(max_output_bytes, bool) or not isinstance(max_output_bytes, int)
                or not 1 <= max_output_bytes <= 33_554_432):
            raise ValueError("max_output_bytes must be an integer in 1..=33554432")
        command = self._command("coverage", memory_budget_mb=memory_budget_mb)
        command += _query(regions, sites, fields) + ["--max-output-bytes", str(max_output_bytes)]
        result = _json(command, max_bytes=max_output_bytes)
        if not isinstance(result, dict) or result.get("schema") != "rosalind-dataset-coverage-v1":
            raise ValueError("unsupported native dataset coverage schema")
        return result

    def batches(self, *, workdir: Optional[PathLike] = None, regions: Optional[PathLike] = None,
                sites: Optional[PathLike] = None, fields: Optional[Sequence[str]] = None,
                memory_budget_mb: Optional[int] = None) -> EvidenceRun:
        """Return a lazy, single-use iterator over native batches of at most 1,024 rows.

        Use as a context manager for early termination. ``run.result`` is set
        only after successful stream exhaustion; its receipt has no persisted
        artifact. Use materialize() when byte verification and replay are needed.
        """
        command = self._command("extract", memory_budget_mb=memory_budget_mb)
        command += _query(regions, sites, fields)
        directory = Path(workdir) if workdir else Path(tempfile.mkdtemp(prefix="rosalind-dataset-"))
        directory.mkdir(parents=True, exist_ok=True)
        manifest = directory / "query.manifest.json"
        return EvidenceRun(command + ["--format", "arrow-ipc", "--manifest", str(manifest)], manifest)

    def plan(self, *, regions: Optional[PathLike] = None, sites: Optional[PathLike] = None,
             fields: Optional[Sequence[str]] = None, memory_budget_mb: Optional[int] = None) -> dict:
        """Plan an Arrow query, including source decoder and projected output costs."""
        return _json(self._command("extract", memory_budget_mb=memory_budget_mb)
                     + _query(regions, sites, fields) + ["--plan"])

    def materialize(self, output: PathLike, *, format: str = "arrow-ipc", force: bool = False,
                    regions: Optional[PathLike] = None, sites: Optional[PathLike] = None,
                    fields: Optional[Sequence[str]] = None,
                    memory_budget_mb: Optional[int] = None) -> EvidenceResult:
        """Persist a covered query and a receipt suitable for byte verification/replay."""
        if format not in ("arrow-ipc", "tsv"):
            raise ValueError("format must be arrow-ipc or tsv")
        command = self._command("extract", memory_budget_mb=memory_budget_mb)
        command += _query(regions, sites, fields) + ["--format", format]
        return _persist(command, output, force)

    def panel_qc(self, regions: PathLike, output: PathLike, *, min_callable_depth: Optional[int] = None,
                 force: bool = False, memory_budget_mb: Optional[int] = None) -> EvidenceResult:
        """Compute exact target summaries; missing stored positions cause refusal."""
        command = self._command("panel-qc", memory_budget_mb=memory_budget_mb)
        command += ["--regions", str(regions)]
        if min_callable_depth is not None:
            if isinstance(min_callable_depth, bool) or min_callable_depth < 0:
                raise ValueError("min_callable_depth must be non-negative")
            command += ["--min-callable-depth", str(min_callable_depth)]
        return _persist(command, output, force)

    def export_parquet(self, directory: PathLike, *, regions: Optional[PathLike] = None,
                       sites: Optional[PathLike] = None, fields: Optional[Sequence[str]] = None,
                       memory_budget_mb: Optional[int] = None) -> ParquetExport:
        """Publish a new directory of bounded Parquet parts; never replace a directory.

        Parquet preserves unsigned 64-bit counts/sums. The receipt verifies physical
        files and lineage. Use Arrow/TSV materialization for native byte replay.
        """
        directory = Path(directory)
        command = self._command("export", memory_budget_mb=memory_budget_mb)
        command += _query(regions, sites, fields) + ["-o", str(directory)]
        _run(command)
        return ParquetExport(directory, directory / "parquet-export.manifest.json")


def open_dataset(manifest: PathLike, *, binary: Optional[PathLike] = None,
                 allow_version_mismatch: bool = False,
                 max_dataset_metadata_bytes: int = 33_554_432) -> EvidenceDataset:
    """Return a lazy portable dataset handle; native reads verify actual file bytes."""
    if isinstance(max_dataset_metadata_bytes, bool) or not isinstance(max_dataset_metadata_bytes, int) or max_dataset_metadata_bytes <= 0:
        raise ValueError("max_dataset_metadata_bytes must be a positive integer")
    return EvidenceDataset(Path(manifest), binary, allow_version_mismatch, max_dataset_metadata_bytes)


def _query(regions: Optional[PathLike], sites: Optional[PathLike], fields: Optional[Sequence[str]]) -> list[str]:
    if regions is not None and sites is not None:
        raise ValueError("supply at most one of regions or sites")
    command = []
    if regions is not None:
        command += ["--regions", str(regions)]
    if sites is not None:
        command += ["--sites", str(sites)]
    if fields is not None:
        if isinstance(fields, str):
            raise TypeError("fields must be a sequence of field names, not a string")
        command += ["--fields", ",".join(fields) if fields else "none"]
    return command


def _run(command: list[str], stdout=None) -> None:
    with tempfile.TemporaryFile() as diagnostics:
        code = subprocess.run(command, stdout=stdout, stderr=diagnostics, check=False).returncode
        if code:
            diagnostics.seek(0, 2)
            diagnostics.seek(max(0, diagnostics.tell() - 65536))
            raise EvidenceProcessError(code, diagnostics.read().decode("utf-8", errors="replace"))


def _json(command: list[str], *, max_bytes: Optional[int] = None) -> dict:
    with tempfile.TemporaryFile() as output:
        _run(command, stdout=output)
        if max_bytes is not None:
            output.seek(0, 2)
            if output.tell() > max_bytes:
                raise ValueError("native JSON response exceeds max_output_bytes")
        output.seek(0)
        return json.load(output)


def _persist(command: list[str], output: PathLike, force: bool) -> EvidenceResult:
    output = Path(output)
    manifest = Path(str(output) + ".manifest.json")
    command += ["-o", str(output), "--manifest", str(manifest)]
    if force:
        command.append("--force")
    _run(command)
    return EvidenceResult(manifest, output)
