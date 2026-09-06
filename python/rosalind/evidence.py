"""Exact short-read evidence in bounded native Arrow batches.

The native memory budget excludes memory retained by Python consumers. Use a
context manager when stopping iteration early; successful stream exhaustion
finalizes ``run.result``. Materialization produces a byte-verifiable artifact.
"""

from __future__ import annotations

import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Optional, Sequence, Union

from .features import _require_pyarrow, _resolve_binary

PathLike = Union[str, Path]


@dataclass(frozen=True)
class EvidenceResult:
    """Finalized native run; streamed runs have no persisted artifact."""

    manifest_path: Path
    artifact_path: Optional[Path]
    returncode: int = 0
    annotated_variants_path: Optional[Path] = None


class EvidenceProcessError(RuntimeError):
    """Native refusal, resource failure, invalid input, or analyzer failure."""

    def __init__(self, returncode: int, diagnostics: str):
        self.returncode = returncode
        super().__init__(f"Rosalind exited {returncode}: {diagnostics.strip()}")


class EvidenceRun:
    """Single-use Arrow iterator with explicit child-process ownership."""

    def __init__(self, command: list[str], manifest: Path):
        _require_pyarrow()
        self.manifest_path = manifest
        self.result: Optional[EvidenceResult] = None
        self._command = command
        self._process = None
        self._diagnostics = None
        self._consumed = False
        self._closed = False

    def batches(self) -> Iterator:
        """Start on first consumption and yield at most 1,024 rows per batch."""
        if self._consumed or self._closed:
            raise RuntimeError("EvidenceRun can only be consumed once")
        self._consumed = True
        arrow = _require_pyarrow()
        # File-backed stderr prevents a pipe deadlock and unbounded Python memory.
        self._diagnostics = tempfile.TemporaryFile()
        self._process = subprocess.Popen(
            self._command, stdout=subprocess.PIPE, stderr=self._diagnostics
        )
        try:
            try:
                with arrow.ipc.open_stream(self._process.stdout) as reader:
                    for batch in reader:
                        yield batch
            except (arrow.ArrowInvalid, arrow.ArrowIOError, OSError) as error:
                # A malformed producer may still be writing to a full stdout
                # pipe. Never wait indefinitely after the consumer rejects it.
                try:
                    code = self._process.wait(timeout=0.5)
                except subprocess.TimeoutExpired:
                    self.close()
                    raise error
                if code:
                    raise EvidenceProcessError(code, self._read_diagnostics()) from error
                raise
            code = self._process.wait()
            if code:
                raise EvidenceProcessError(code, self._read_diagnostics())
            self.result = EvidenceResult(self.manifest_path, None)
        finally:
            self.close()

    def _read_diagnostics(self) -> str:
        if self._diagnostics is None:
            return ""
        self._diagnostics.seek(0, 2)
        self._diagnostics.seek(max(0, self._diagnostics.tell() - 65536))
        return self._diagnostics.read().decode("utf-8", errors="replace")

    def __iter__(self) -> Iterator:
        return self.batches()

    def close(self) -> None:
        """Cancel early consumption and release pipes; no successful result is set."""
        if self._closed:
            return
        self._closed = True
        if self._process is not None:
            if self._process.stdout is not None:
                self._process.stdout.close()
            if self._process.poll() is None:
                self._process.terminate()
                try:
                    self._process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    self._process.kill()
            self._process.wait()
        if self._diagnostics is not None:
            self._diagnostics.close()

    def __enter__(self) -> EvidenceRun:
        return self

    def __exit__(self, *_: object) -> None:
        self.close()


def _command(
    reference: Optional[PathLike], alignments: PathLike, *,
    sites: Optional[PathLike] = None, regions: Optional[PathLike] = None,
    binary: Optional[PathLike] = None, allow_version_mismatch: bool = False,
    mapq: int = 20, base_quality: int = 20, memory_budget_mb: Optional[int] = None,
    max_read_len: int = 250, tile_bases: int = 16384, workers: int = 1,
    cache_dir: Optional[PathLike] = None, resume: bool = False,
    cram_reference: Optional[PathLike] = None, panel: bool = False,
    min_callable_depth: Optional[int] = None,
    sample: Optional[str] = None, pool_samples: bool = False,
    fields: Optional[Union[str, Sequence[str]]] = None,
) -> list[str]:
    if (sites is None) == (regions is None):
        raise ValueError("supply exactly one of sites (VCF) or regions (BED)")
    if panel and sites is not None:
        raise ValueError("panel QC requires BED regions")
    if not panel and reference is None:
        raise ValueError("SNV evidence requires an explicit reference")
    if sample is not None and pool_samples:
        raise ValueError("sample and pool_samples are mutually exclusive")
    if sample is not None and not sample:
        raise ValueError("sample must be a nonempty declared sample name")
    if not panel and min_callable_depth is not None:
        raise ValueError("min_callable_depth is a panel QC option")
    field_names = None
    if fields is not None:
        field_names = fields.split(",") if isinstance(fields, str) else list(fields)
        groups = {"depths", "alleles", "strands", "quality-sums", "quality-histograms", "read-position", "allele-quality"}
        if not field_names:
            field_names = ["none"]
        if field_names not in (["all"], ["all-supported"], ["none"]) and any(name not in groups for name in field_names):
            raise ValueError("fields must contain known evidence groups, or all/none alone")
        if panel and field_names not in (["all"], ["all-supported"]) and not {"depths", "quality-sums"}.issubset(field_names):
            raise ValueError("panel QC fields must include depths and quality-sums")
    executable = _resolve_binary(binary, allow_version_mismatch)
    command = [str(executable), "analyze", "panel-qc" if panel else "evidence",
               "--alignments", str(alignments), "--mapq-threshold", str(mapq),
               "--base-quality-threshold", str(base_quality),
               "--max-read-len", str(max_read_len), "--tile-bases", str(tile_bases),
               "--workers", str(workers)]
    command += ["--sites", str(sites)] if sites is not None else ["--regions", str(regions)]
    if reference is not None:
        command += ["--reference", str(reference)]
    if cram_reference is not None:
        command += ["--cram-reference", str(cram_reference)]
    if memory_budget_mb is not None:
        command += ["--memory-budget-mb", str(memory_budget_mb)]
    if cache_dir is not None:
        command += ["--cache-dir", str(cache_dir)]
    if resume:
        command.append("--resume")
    if sample is not None:
        command += ["--sample", sample]
    if pool_samples:
        command.append("--pool-samples")
    if field_names is not None:
        command += ["--fields", ",".join(field_names)]
    # The native executable owns the shared default; Python only passes overrides.
    if panel and min_callable_depth is not None:
        command += ["--min-callable-depth", str(min_callable_depth)]
    return command


def iter_evidence(
    reference: PathLike, alignments: PathLike, *, workdir: Optional[PathLike] = None,
    **options: object,
) -> EvidenceRun:
    """Return a lazy exact evidence run; requires exactly one sites/regions option.

    Histograms, counts and sums use the versioned shortread-dna-readcount-v1
    profile. Default MAPQ and base quality thresholds are both 20. Overlapping
    mates count as separate reads. A supplied budget governs the native process.
    Use ``sample="name"`` for a declared SM sample or ``pool_samples=True`` for
    deliberate pooling. Ambiguous sample scope is refused by default.
    Pass ``fields=["depths", "alleles"]`` to allocate and encode only those
    groups. Omitted columns are absent, never filled with zeros.
    """
    command = _command(reference, alignments, **options)
    directory = Path(workdir) if workdir else Path(tempfile.mkdtemp(prefix="rosalind-evidence-"))
    directory.mkdir(parents=True, exist_ok=True)
    manifest = directory / "evidence.manifest.json"
    return EvidenceRun(command + ["--format", "arrow-ipc", "--manifest", str(manifest)], manifest)


def materialize_evidence(
    reference: PathLike, alignments: PathLike, output: PathLike, *,
    format: str = "arrow-ipc", force: bool = False,
    annotated_variants: Optional[PathLike] = None, **options: object,
) -> EvidenceResult:
    """Persist evidence and its receipt; optionally annotate original SNV records.

    ``annotated_variants`` requires ``sites`` (VCF, VCF.gz, or BCF) and a
    .vcf/.vcf.gz/.bcf destination. Existing annotations, allele/record order and
    genotypes are preserved. INFO fields describe the chosen alignment sample
    scope. Request ``allele-quality`` for exact per-allele quality/position sums.
    """
    if format not in ("arrow-ipc", "tsv"):
        raise ValueError("format must be arrow-ipc or tsv")
    command = _command(reference, alignments, **options)
    if annotated_variants is not None:
        if options.get("sites") is None:
            raise ValueError("annotated_variants requires sites (VCF, VCF.gz, or BCF)")
        command += ["--annotated-variants", str(annotated_variants)]
    result = _materialize(command, output, format, force)
    return EvidenceResult(result.manifest_path, result.artifact_path, result.returncode,
                          Path(annotated_variants) if annotated_variants is not None else None)


def panel_qc(
    alignments: PathLike, regions: PathLike, output: PathLike, *,
    reference: Optional[PathLike] = None, force: bool = False, **options: object,
) -> EvidenceResult:
    """Persist per-target TSV summaries with explicit coverage denominators.

    The matching native executable supplies the shared 10x callability default;
    pass ``min_callable_depth`` to override it explicitly.
    """
    command = _command(reference, alignments, regions=regions, panel=True, **options)
    return _materialize(command, output, "tsv", force)


def _materialize(command: list[str], output: PathLike, format: str, force: bool) -> EvidenceResult:
    output = Path(output)
    manifest = Path(str(output) + ".manifest.json")
    command += ["--format", format, "-o", str(output), "--manifest", str(manifest)]
    if force:
        command.append("--force")
    with tempfile.TemporaryFile() as diagnostics:
        code = subprocess.run(command, stderr=diagnostics, check=False).returncode
        if code:
            diagnostics.seek(0, 2)
            diagnostics.seek(max(0, diagnostics.tell() - 65536))
            raise EvidenceProcessError(code, diagnostics.read().decode("utf-8", errors="replace"))
    return EvidenceResult(manifest, output)
