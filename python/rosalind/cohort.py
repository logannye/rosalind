"""Preview local cohort queries through Rosalind's native verified evidence engine.

Handles are lazy and snapshot identities are immutable. Saved-only queries never
open original alignments. Extension requires explicit content-matched sources.
Native budgets exclude batches/tables retained by the Python consumer.
"""
from __future__ import annotations

import json
import subprocess
import tempfile
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Iterator, Optional, Sequence

from .dataset import _json
from .evidence import EvidenceProcessError, EvidenceResult, EvidenceRun, PathLike
from .features import _require_pyarrow, _resolve_binary


class CohortRun(EvidenceRun):
    """Lazy, single-use disk-backed Arrow batches with child-process ownership.

    The native query materializes and verifies a complete artifact before the
    first batch is yielded. Iteration retains at most one native canonical batch
    (1,024 rows) unless the consumer saves batches. Use a context manager when
    stopping early. The artifact and receipt remain in ``workdir`` (a new
    temporary directory by default); ``result`` is set on full exhaustion only.
    """

    def __init__(self, command: list[str], manifest: Path, artifact: Path):
        super().__init__(command, manifest)
        self.artifact_path = artifact
        self._reader = None
        self._transport = None

    def batches(self) -> Iterator:
        if self._consumed or self._closed:
            raise RuntimeError("CohortRun can only be consumed once")
        self._consumed = True
        arrow = _require_pyarrow()
        self._diagnostics = tempfile.TemporaryFile()
        try:
            self._transport = _CommandTransport(self._command)
            self._process = subprocess.Popen(
                self._transport.prepare(), stdout=subprocess.DEVNULL, stderr=self._diagnostics
            )
            code = self._process.wait()
            if code:
                raise EvidenceProcessError(code, self._read_diagnostics())
            if self._closed:
                return
            with arrow.ipc.open_stream(self.artifact_path) as reader:
                self._reader = reader
                for batch in reader:
                    yield batch
            self.result = EvidenceResult(self.manifest_path, self.artifact_path)
        finally:
            self.close()

    def close(self) -> None:
        if self._reader is not None:
            self._reader.close()
            self._reader = None
        super().close()
        if self._transport is not None:
            self._transport.close()


@dataclass(frozen=True)
class CohortExtension:
    """New immutable handle plus native extension measurements; parent is unchanged."""

    cohort: EvidenceCohort
    report: dict


@dataclass(frozen=True)
class EvidenceCohort:
    """A lazy handle to one local immutable snapshot (next-minor preview).

    Member IDs/metadata are user assertions, not verified biological identities.
    Inspection and queries verify native metadata and consumed saved evidence;
    ``verify`` scans every stored row. Keep the cohort immutable during each run.
    """

    directory: Path
    snapshot_id: str
    binary: Optional[PathLike] = None
    allow_version_mismatch: bool = False
    max_snapshot_bytes: int = 8_388_608
    max_dataset_metadata_bytes: int = 33_554_432

    def _command(self, action: str, *, memory_budget_mb: Optional[int] = None,
                 require_os_limit: bool = False) -> list[str]:
        command = [str(_resolve_binary(self.binary, self.allow_version_mismatch)),
                   "cohort", action, "--cohort", str(self.directory),
                   "--snapshot", self.snapshot_id,
                   "--max-snapshot-bytes", str(self.max_snapshot_bytes),
                   "--max-dataset-metadata-bytes", str(self.max_dataset_metadata_bytes)]
        if memory_budget_mb is not None:
            _positive(memory_budget_mb, "memory_budget_mb")
            command += ["--memory-budget-mb", str(memory_budget_mb), "--enforce"]
        if require_os_limit:
            if memory_budget_mb is None:
                raise ValueError("require_os_limit requires memory_budget_mb")
            command.append("--require-os-limit")
        return command

    def inspect(self, *, memory_budget_mb: Optional[int] = None,
                require_os_limit: bool = False) -> dict:
        """Return verified snapshot/leaf metadata without scanning Arrow rows."""
        return _json(self._command("inspect", memory_budget_mb=memory_budget_mb,
                                   require_os_limit=require_os_limit))

    def verify(self, *, memory_budget_mb: Optional[int] = None,
               require_os_limit: bool = False) -> dict:
        """Verify every stored row and declared inventory, without original sources."""
        return _json(self._command("verify", memory_budget_mb=memory_budget_mb,
                                   require_os_limit=require_os_limit))

    def _query(self, action: str, sites: PathLike, *, members: Optional[Sequence[str]] = None,
               fields: Optional[Sequence[str]] = None, missing: str = "strict",
               min_callable_depth: int = 10, memory_budget_mb: Optional[int] = None,
               require_os_limit: bool = False, tile_bases: int = 16_384,
               max_candidate_sites: int = 1_000_000,
               max_variant_header_bytes: int = 8_388_608,
               max_variant_record_bytes: int = 1_048_576,
               max_receipt_bytes: int = 33_554_432) -> list[str]:
        if missing not in ("strict", "partial"):
            raise ValueError("missing must be strict or partial")
        for name, value in (("min_callable_depth", min_callable_depth),
                            ("tile_bases", tile_bases),
                            ("max_candidate_sites", max_candidate_sites),
                            ("max_variant_header_bytes", max_variant_header_bytes),
                            ("max_variant_record_bytes", max_variant_record_bytes),
                            ("max_receipt_bytes", max_receipt_bytes)):
            _positive(value, name)
        command = self._command(action, memory_budget_mb=memory_budget_mb,
                                require_os_limit=require_os_limit)
        command += ["--sites", str(sites), "--missing", missing,
                    "--min-callable-depth", str(min_callable_depth),
                    "--tile-bases", str(tile_bases),
                    "--max-candidate-sites", str(max_candidate_sites),
                    "--max-variant-header-bytes", str(max_variant_header_bytes),
                    "--max-variant-record-bytes", str(max_variant_record_bytes),
                    "--max-receipt-bytes", str(max_receipt_bytes)]
        if members is not None:
            if isinstance(members, str):
                raise TypeError("members must be a sequence of IDs, not a string")
            ids = list(members)
            if any(not isinstance(member, str) or not member for member in ids):
                raise ValueError("member IDs must be nonempty strings")
            if len(ids) != len(set(ids)):
                raise ValueError("member IDs must be unique")
            # This native operand preserves the distinction between None and [].
            command += ["--cohort-members", json.dumps(ids, separators=(",", ":"))]
        if fields is not None:
            if isinstance(fields, str):
                raise TypeError("fields must be a sequence of names, not a string")
            names = list(fields)
            if not names or any(not isinstance(name, str) or not name for name in names):
                raise ValueError("fields must contain nonempty native field names")
            command += ["--fields", ",".join(names)]
        return command

    def plan(self, *, sites: PathLike, operation: str = "extract",
             sources: Optional[PathLike] = None, workdir: Optional[PathLike] = None,
             pairs: Optional[PathLike] = None, max_pairs: int = 65_536,
             max_pair_table_bytes: int = 8_388_608, **options) -> dict:
        """Plan extract/summarize/extend/compare-pairs without publishing or opening raw sources.

        Query options shared by all methods: members (None=all, []=none), fields,
        missing (strict/partial), min_callable_depth (default technical screen 10),
        memory_budget_mb (enforced cooperatively), require_os_limit, tile_bases,
        max_candidate_sites, max_variant_header_bytes, max_variant_record_bytes,
        and max_receipt_bytes. Scientific incompatibilities remain native errors
        or structured blocking plan issues; partial only permits unmeasured loci.
        """
        if operation not in ("extract", "summarize", "extend", "compare-pairs"):
            raise ValueError("operation must be extract, summarize, extend or compare-pairs")
        command = self._query(operation, sites, **options)
        if operation == "compare-pairs":
            command += _pairs(pairs, max_pairs, max_pair_table_bytes)
        elif pairs is not None:
            raise ValueError("pairs applies only to compare-pairs")
        if operation == "extend":
            command += _extension(sources, workdir)
        elif sources is not None or workdir is not None:
            raise ValueError("sources/workdir apply only to extension plans")
        return _native_json(command + ["--plan"])

    def materialize(self, output: PathLike, *, sites: PathLike,
                    format: str = "arrow-ipc", force: bool = False,
                    **options) -> EvidenceResult:
        """Persist sample-by-candidate evidence and a verifiable/replayable receipt.

        See ``plan`` for query options. Strict coverage is the default. Partial
        coverage represents unmeasured cells as nulls; measured zero stays zero.
        """
        _format(format)
        return _materialize(self._query("extract", sites, **options) + ["--format", format], output, force)

    def summarize(self, output: PathLike, *, sites: PathLike,
                  format: str = "arrow-ipc", force: bool = False,
                  **options) -> EvidenceResult:
        """Persist native candidate summaries with explicit sample denominators.

        ``min_callable_depth`` defaults to 10, a technical screen rather than
        confidence. ALT-support proportions are neither genotypes nor population
        allele frequencies. See ``plan`` for shared query options.
        """
        _format(format)
        return _materialize(self._query("summarize", sites, **options) + ["--format", format], output, force)

    def compare_pairs(self, output: PathLike, *, pairs: PathLike, sites: PathLike,
                      format: str = "arrow-ipc", force: bool = False,
                      max_pairs: int = 65_536, max_pair_table_bytes: int = 8_388_608,
                      **options) -> EvidenceResult:
        """Persist explicit pair comparisons, in table order then candidate order.

        ``pairs`` is a TSV with exactly id, left, right columns. Direction is
        right-minus-left observed ALT fraction. Both sides retain exact uint64
        counts, missingness and technical eligibility. Difference magnitude and
        denominator are decimal strings (uint128); zero-depth or unmeasured
        sides yield null differences. Names/metadata never imply pairing.
        ``fields`` must be depths and alleles; see ``plan`` for shared options.
        """
        _format(format)
        command = self._query("compare-pairs", sites, **options)
        command += _pairs(pairs, max_pairs, max_pair_table_bytes)
        return _materialize(command + ["--format", format], output, force)

    def batches(self, *, sites: PathLike, workdir: Optional[PathLike] = None,
                **options) -> CohortRun:
        """Lazily materialize once, then read bounded Arrow batches; see CohortRun."""
        command = self._query("extract", sites, **options)
        directory = Path(workdir) if workdir is not None else Path(tempfile.mkdtemp(prefix="rosalind-cohort-"))
        directory.mkdir(parents=True, exist_ok=True)
        artifact = directory / "cohort.arrow"
        manifest = directory / "cohort.arrow.manifest.json"
        return CohortRun(command + ["--format", "arrow-ipc", "-o", str(artifact),
                                    "--manifest", str(manifest)], manifest, artifact)

    def extend(self, sources: PathLike, *, sites: PathLike,
               workdir: Optional[PathLike] = None, **options) -> CohortExtension:
        """Publish missing loci into a child snapshot; retain this parent unchanged.

        ``sources`` is a TSV with id, role, path columns (paths relative to that
        table). Supply exactly the affected members and every recorded scientific
        source role. Native code verifies original hashes, scope and filters.
        Existing loci cannot silently gain fields. ``workdir`` must exist outside
        the cohort (defaults to the current directory). A no-op returns this ID.
        """
        report = _native_json(self._query("extend", sites, **options) + _extension(sources, workdir))
        return CohortExtension(replace(self, snapshot_id=report["snapshot_id"]), report)


def open_cohort(directory: PathLike, snapshot: str, *, binary: Optional[PathLike] = None,
                allow_version_mismatch: bool = False, max_snapshot_bytes: int = 8_388_608,
                max_dataset_metadata_bytes: int = 33_554_432) -> EvidenceCohort:
    """Return a lazy handle; no binary, cohort or original-source files are opened."""
    if not isinstance(snapshot, str) or len(snapshot) != 64 or any(c not in "0123456789abcdef" for c in snapshot):
        raise ValueError("snapshot must be a lowercase 64-character content identity")
    _positive(max_snapshot_bytes, "max_snapshot_bytes")
    _positive(max_dataset_metadata_bytes, "max_dataset_metadata_bytes")
    return EvidenceCohort(Path(directory), snapshot, binary, allow_version_mismatch,
                          max_snapshot_bytes, max_dataset_metadata_bytes)


def _pairs(path: Optional[PathLike], max_pairs: int, max_table_bytes: int) -> list[str]:
    if path is None:
        raise ValueError("compare-pairs requires an explicit pair table")
    _positive(max_pairs, "max_pairs")
    _positive(max_table_bytes, "max_pair_table_bytes")
    return ["--pairs", str(path), "--max-pairs", str(max_pairs),
            "--max-pair-table-bytes", str(max_table_bytes)]


def _extension(sources: Optional[PathLike], workdir: Optional[PathLike]) -> list[str]:
    if sources is None:
        raise ValueError("extension requires an explicit source mapping table")
    command = ["--sources", str(sources)]
    if workdir is not None:
        command += ["--work-dir", str(workdir)]
    return command


def _materialize(command: list[str], output: PathLike, force: bool) -> EvidenceResult:
    output = Path(output)
    manifest = Path(str(output) + ".manifest.json")
    command += ["-o", str(output), "--manifest", str(manifest)]
    if force:
        command.append("--force")
    # Consume the native status JSON instead of printing it from a library call.
    _native_json(command)
    return EvidenceResult(manifest, output)


def _positive(value: int, name: str) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(f"{name} must be a positive integer")


def _format(value: str) -> None:
    if value not in ("arrow-ipc", "tsv"):
        raise ValueError("format must be arrow-ipc or tsv")


class _CommandTransport:
    """Own a private bounded argument file when member selection exceeds argv."""

    argv_bytes = 65_536
    request_bytes = 33_554_432

    def __init__(self, command: list[str]):
        self.command = command
        self.path = None

    def prepare(self) -> list[str]:
        if sum(len(token.encode("utf-8")) + 1 for token in self.command) <= self.argv_bytes:
            return self.command
        if self.command[1:3] not in (["cohort", "extract"], ["cohort", "summarize"], ["cohort", "compare-pairs"]):
            raise ValueError("large explicit member selections are supported for extract/summarize/compare-pairs only; extension can select all members and use its bounded source table")
        payload = json.dumps(self.command[1:], ensure_ascii=False, separators=(",", ":")).encode("utf-8")
        if len(payload) > self.request_bytes:
            raise ValueError("cohort argument request exceeds the native 32 MiB envelope")
        with tempfile.NamedTemporaryFile(prefix="rosalind-cohort-request-", suffix=".json", delete=False) as request:
            self.path = Path(request.name)
            request.write(payload)
        command = [self.command[0], "cohort", "replay", "--request", str(self.path)]
        if "--memory-budget-mb" in self.command:
            budget = int(self.command[self.command.index("--memory-budget-mb") + 1]) << 20
            command += ["--memory-budget-bytes", str(budget)]
        for flag in ("--enforce", "--require-os-limit"):
            if flag in self.command:
                command.append(flag)
        return command

    def close(self) -> None:
        if self.path is not None:
            self.path.unlink(missing_ok=True)
            self.path = None


def _native_json(command: list[str]) -> dict:
    transport = _CommandTransport(command)
    try:
        return _json(transport.prepare())
    finally:
        transport.close()
