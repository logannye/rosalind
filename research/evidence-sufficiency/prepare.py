"""Bounded P0 acquisition and template-preserving preparation.

``prepare_inputs(resources_path, protocol, out)`` requires a new output directory.
Only its explicit call imports/uses pysam; import and unit tests perform no native
execution. Run the entry point on the hosted study worker. Required protocol keys
are ``intervals``, ``downsample.{seed,numerator,denominator}``, and
``limits.{download_bytes,scratch_bytes}``. Optional ``sample`` defaults to HG002;
``limits.analysis_seconds`` defaults to 7200 and ``network_timeout_seconds`` to 30.
``require_nonempty_benchmark_coverage`` optionally refuses any selected block with
zero benchmark-BED coverage, before reference expansion or native execution.
Deadlines and scratch limits are cooperative checks; the hosted job must impose a
hard outer timeout. Failed preparation retains bounded partial evidence but never
publishes preparation.json. No shared cache or existing output is overwritten.
"""

from __future__ import annotations

import bisect
import gzip
import hashlib
import importlib
import json
import os
from pathlib import Path
import re
import shutil
import stat
import time
import urllib.error
import urllib.parse
import urllib.request


CHUNK_BYTES = 1 << 20
REQUIRED = ("reads_bam", "reads_bai", "reference", "truth_vcf", "truth_bed")


def _positive_int(value, name):
    if type(value) is not int or value <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return value


def _regular(path):
    path = Path(path)
    if not stat.S_ISREG(path.lstat().st_mode):
        raise ValueError(f"input must be a regular file, not a symlink: {path}")
    return path


def _read_regular(path):
    path = _regular(path)
    fd = os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
    if not stat.S_ISREG(os.fstat(fd).st_mode):
        os.close(fd)
        raise ValueError(f"input changed to a nonregular file: {path}")
    return os.fdopen(fd, "rb")


def _digest(path, budget=None):
    digest = hashlib.sha256()
    with _read_regular(path) as source:
        while True:
            if budget:
                budget.check_time()
            chunk = source.read(CHUNK_BYTES)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def _source_lock(path):
    with _read_regular(path) as source:
        raw = source.read((1 << 20) + 1)
    if len(raw) > 1 << 20:
        raise ValueError("resource lock exceeds 1 MiB")
    rows, names = {}, set()
    for line in raw.decode("utf-8").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) != 4:
            raise ValueError("resource lock requires four tab-separated fields")
        key, name, sha, url = parts
        if not re.fullmatch(r"[A-Za-z0-9_-]+", key) or key in rows:
            raise ValueError("invalid or duplicate resource identifier")
        if (not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", name)
                or len(name.encode()) > 255 or name.casefold() in names):
            raise ValueError("resource filename must be a unique safe basename")
        if not re.fullmatch(r"[0-9a-f]{64}", sha):
            raise ValueError("resource requires lowercase SHA256")
        parsed = urllib.parse.urlsplit(url)
        if (parsed.scheme != "https" or not parsed.hostname or parsed.username is not None
                or parsed.password is not None or parsed.fragment
                or any(c.isspace() or ord(c) < 32 for c in url)):
            raise ValueError("resource requires an HTTPS URL without credentials or fragment")
        names.add(name.casefold())
        rows[key] = {"id": key, "filename": name, "sha256": sha, "url": url}
    if any(key not in rows for key in REQUIRED):
        raise ValueError(f"resource lock must contain {REQUIRED}")
    return [rows[key] for key in REQUIRED], raw


def _validate_protocol(protocol):
    if type(protocol.get("require_nonempty_benchmark_coverage", False)) is not bool:
        raise ValueError("require_nonempty_benchmark_coverage must be boolean")
    limits = protocol["limits"]
    download = _positive_int(limits["download_bytes"], "download_bytes")
    scratch = _positive_int(limits["scratch_bytes"], "scratch_bytes")
    if download > 4 * (1 << 30) or scratch > 30 * (1 << 30):
        raise ValueError("P0 ceilings are 4 GiB acquisition and 30 GiB scratch")
    _positive_int(limits.get("analysis_seconds", 7200), "analysis_seconds")
    _positive_int(limits.get("network_timeout_seconds", 30), "network_timeout_seconds")
    downsample = protocol["downsample"]
    numerator = _positive_int(downsample["numerator"], "numerator")
    denominator = _positive_int(downsample["denominator"], "denominator")
    if numerator >= denominator:
        raise ValueError("downsample fraction must be strictly between zero and one")
    seed = downsample["seed"]
    if type(seed) not in (int, str) or isinstance(seed, str) and not seed:
        raise ValueError("seed must be an integer or nonempty string")
    sample = protocol.get("sample", "HG002")
    if not isinstance(sample, str) or not sample or any(c.isspace() for c in sample):
        raise ValueError("sample must be a nonempty name without whitespace")
    intervals = protocol["intervals"]
    if not isinstance(intervals, list) or not 1 <= len(intervals) <= 1024:
        raise ValueError("intervals must contain between 1 and 1024 entries")
    result, block_ids = [], set()
    for item in intervals:
        if item["contig"] != "chr20":
            raise ValueError("P0 preparation accepts chr20 intervals only")
        start, end = item["start"], item["end"]
        if type(start) is not int or type(end) is not int or not 0 <= start < end < 1 << 63:
            raise ValueError("intervals require integer zero-based half-open coordinates")
        for key in ("partition", "block_id"):
            value = item[key]
            if not isinstance(value, str) or not re.fullmatch(r"[A-Za-z0-9_.-]+", value):
                raise ValueError(f"invalid {key}")
        if item["block_id"] in block_ids:
            raise ValueError("duplicate block_id")
        block_ids.add(item["block_id"])
        result.append(dict(item))
    result.sort(key=lambda x: (x["start"], x["end"]))
    if any(left["end"] > right["start"] for left, right in zip(result, result[1:])):
        raise ValueError("intervals must not overlap across or within partitions")
    return result


class _Budget:
    def __init__(self, root, limits):
        self.root = Path(root)
        self.download_limit = limits["download_bytes"]
        self.scratch_limit = limits["scratch_bytes"]
        self.started = time.monotonic()
        self.deadline = self.started + limits.get("analysis_seconds", 7200)
        self.timeout = limits.get("network_timeout_seconds", 30)
        self.download_bytes = 0
        self.scratch_high_water = 0
        self.attempts = 0

    def check_time(self):
        if time.monotonic() >= self.deadline:
            raise TimeoutError("preparation deadline exceeded")

    def check_scratch(self, incoming=0):
        self.check_time()
        total = 0
        for directory, dirs, files in os.walk(self.root, followlinks=False):
            for name in dirs:
                if Path(directory, name).is_symlink():
                    raise ValueError("symlink directory in preparation output")
            for name in files:
                path = _regular(Path(directory, name))
                total += path.stat().st_size
        self.scratch_high_water = max(self.scratch_high_water, total)
        if total + incoming > self.scratch_limit:
            raise ValueError("scratch byte budget exceeded")
        if incoming and shutil.disk_usage(self.root).free < incoming:
            raise ValueError("insufficient free disk space")
        return total


class _HttpsRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        if urllib.parse.urlsplit(newurl).scheme != "https":
            raise ValueError("download redirect must remain HTTPS")
        return super().redirect_request(req, fp, code, msg, headers, newurl)


def _download(row, root, budget, opener=None):
    destination = root / row["filename"]
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(destination)
    opener = opener or urllib.request.build_opener(_HttpsRedirect()).open
    partial = destination.with_name(destination.name + ".partial")
    started = time.monotonic()
    for attempt in range(3):
        budget.check_time()
        budget.attempts += 1
        try:
            request = urllib.request.Request(row["url"], headers={"Accept-Encoding": "identity"})
            timeout = max(0.001, min(budget.timeout, budget.deadline - time.monotonic()))
            with opener(request, timeout=timeout) as response, partial.open("xb") as output:
                length = response.headers.get("Content-Length")
                if length is not None:
                    length = int(length)
                    if length < 0 or length > budget.download_limit - budget.download_bytes:
                        raise ValueError("download byte budget exceeded by Content-Length")
                    budget.check_scratch(length)
                digest, received = hashlib.sha256(), 0
                while True:
                    budget.check_time()
                    remaining = budget.download_limit - budget.download_bytes
                    if remaining == 0:
                        if length is not None and received == length:
                            break
                        raise ValueError("download byte budget exhausted")
                    chunk = response.read(min(CHUNK_BYTES, remaining))
                    if not chunk:
                        break
                    budget.download_bytes += len(chunk)
                    received += len(chunk)
                    budget.check_scratch(len(chunk))
                    output.write(chunk)
                    digest.update(chunk)
                output.flush()
                os.fsync(output.fileno())
                if length is not None and received != length:
                    raise OSError("incomplete download")
                if digest.hexdigest() != row["sha256"]:
                    raise ValueError(f"download SHA256 mismatch: {row['id']}")
            os.link(partial, destination)
            partial.unlink()
            budget.check_scratch()
            return {**row, "path": str(destination), "bytes": received,
                    "download_seconds": time.monotonic() - started, "attempts": attempt + 1}
        except (OSError, urllib.error.URLError):
            partial.unlink(missing_ok=True)
            budget.check_time()
            if attempt == 2:
                raise
        except BaseException:
            partial.unlink(missing_ok=True)
            raise
    raise AssertionError("unreachable")


def _keep_template(seed, source_sha256, read_group, name, numerator, denominator):
    if not isinstance(name, str) or not name:
        raise ValueError("template-preserving downsampling requires a QNAME")
    key = json.dumps([seed, source_sha256, read_group, name],
                     ensure_ascii=True, separators=(",", ":")).encode()
    value = int.from_bytes(hashlib.sha256(key).digest(), "big")
    return value * denominator < numerator * (1 << 256)


def _overlaps(start, end, starts, ends):
    index = bisect.bisect_left(starts, end) - 1
    return index >= 0 and ends[index] > start


def benchmark_coverage(truth_bed: Path, intervals: list, budget=None) -> dict:
    """Intersect a merged benchmark BED with selected zero-based half-open blocks.

    Reads coordinates only, never the truth VCF or candidate observations. Plain
    and gzip BEDs are supported; overlapping/duplicate/unsorted mask rows are
    merged before counting. Only mask rows overlapping selected blocks are kept.
    ``intervals`` must already have passed protocol validation.
    """
    selected = {}
    for item in intervals:
        selected.setdefault(item["contig"], []).append(item)
    bounds = {}
    masks = {contig: [] for contig in selected}
    for contig, blocks in selected.items():
        blocks.sort(key=lambda item: (item["start"], item["end"]))
        bounds[contig] = ([item["start"] for item in blocks],
                          [item["end"] for item in blocks])
    with _read_regular(truth_bed) as source:
        compressed = source.read(2) == b"\x1f\x8b"
        source.seek(0)
        with gzip.GzipFile(fileobj=source, mode="rb") if compressed else source as bed:
            line_number = 0
            while True:
                if budget:
                    budget.check_time()
                line = bed.readline(CHUNK_BYTES + 1)
                if not line:
                    break
                line_number += 1
                if len(line) > CHUNK_BYTES:
                    raise ValueError(f"benchmark BED line {line_number} exceeds 1 MiB")
                fields = line.split()
                if not fields or fields[0].startswith(b"#") or fields[0] in (b"track", b"browser"):
                    continue
                if (len(fields) < 3 or not re.fullmatch(rb"[0-9]+", fields[1])
                        or not re.fullmatch(rb"[0-9]+", fields[2])):
                    raise ValueError(f"invalid benchmark BED coordinates at line {line_number}")
                contig = fields[0].decode("utf-8")
                start, end = int(fields[1]), int(fields[2])
                if not 0 <= start <= end < 1 << 63:
                    raise ValueError(f"invalid benchmark BED interval at line {line_number}")
                if contig in selected and start < end:
                    starts, ends = bounds[contig]
                    if _overlaps(start, end, starts, ends):
                        masks[contig].append((max(start, starts[0]), min(end, ends[-1])))
    merged = {}
    for contig, rows in masks.items():
        union = []
        for start, end in sorted(rows):
            if union and start <= union[-1][1]:
                union[-1] = (union[-1][0], max(end, union[-1][1]))
            else:
                union.append((start, end))
        merged[contig] = (union, [end for _, end in union])
    blocks = {}
    for item in intervals:
        union, ends = merged[item["contig"]]
        confident = 0
        index = bisect.bisect_right(ends, item["start"])
        while index < len(union) and union[index][0] < item["end"]:
            start, end = union[index]
            confident += min(item["end"], end) - max(item["start"], start)
            index += 1
        selected_bases = item["end"] - item["start"]
        blocks[item["block_id"]] = {key: item[key] for key in
                                    ("contig", "start", "end", "partition")}
        blocks[item["block_id"]].update({"selected_bases": selected_bases,
                                        "confident_bases": confident,
                                        "confident_fraction": confident / selected_bases})
    return {"by_block": blocks,
            "selected_bases": sum(block["selected_bases"] for block in blocks.values()),
            "confident_bases": sum(block["confident_bases"] for block in blocks.values()),
            "empty_block_ids": [key for key, block in blocks.items() if block["confident_bases"] == 0]}


def _audit_inputs(source, fasta, sample, intervals):
    if tuple(zip(source.references, source.lengths)) != tuple(zip(fasta.references, fasta.lengths)):
        raise ValueError("full ordered BAM/reference dictionaries differ")
    groups = source.header.to_dict().get("RG", [])
    if (not groups or any(not group.get("ID") or group.get("SM") != sample for group in groups)
            or len({group["ID"] for group in groups}) != len(groups)):
        raise ValueError("BAM must have distinct read groups with the sole named selected SM")
    dictionary = dict(zip(source.references, source.lengths))
    if any(item["end"] > dictionary.get(item["contig"], -1) for item in intervals):
        raise ValueError("interval outside reference dictionary")
    if not source.has_index():
        raise ValueError("BAM index is unavailable")
    return {group["ID"] for group in groups}


def _downsample(source, destination, pysam, intervals, sample_groups, source_sha, spec, budget):
    starts = [item["start"] for item in intervals]
    ends = [item["end"] for item in intervals]
    metrics = {"traversed_records": 0, "selected_records": 0, "retained_records": 0,
               "selected_query_bases": 0, "retained_query_bases": 0}
    with pysam.AlignmentFile(str(destination), "wb", header=source.header) as output:
        # A single indexed iterator yields each physical record once. Never use
        # QNAME/content sets: legitimate duplicate records must remain present.
        for record in source.fetch("chr20"):
            metrics["traversed_records"] += 1
            if metrics["traversed_records"] % 4096 == 0:
                budget.check_scratch()
            start = record.reference_start
            end = record.reference_end
            if start < 0 or end is None or end <= start:
                continue
            if not _overlaps(start, end, starts, ends):
                continue
            if not record.has_tag("RG") or record.get_tag("RG") not in sample_groups:
                raise ValueError("selected record has absent or unknown RG")
            metrics["selected_records"] += 1
            metrics["selected_query_bases"] += record.query_length or 0
            if _keep_template(spec["seed"], source_sha, record.get_tag("RG"), record.query_name,
                              spec["numerator"], spec["denominator"]):
                output.write(record)
                metrics["retained_records"] += 1
                metrics["retained_query_bases"] += record.query_length or 0
    budget.check_scratch()
    return metrics


def verify_prepared_inputs(state: dict) -> dict:
    """Recheck every prepared input against its original byte identity.

    Call after scientific consumption and before publication. The same pathname
    may occur as both a source and a derived-input alias; hash it once only after
    checking that all aliases declare identical sizes and digests. This does not
    reauthenticate the biological material or establish a new identity from the
    files' current bytes. Regular files are required, including for empty inputs.
    """
    identities = state["identities"]
    sources = identities["sources"]
    if not isinstance(sources, dict) or any(key not in sources for key in REQUIRED):
        raise ValueError("prepared identities must retain all five locked sources")
    entries = [(f"source:{key}", value) for key, value in sources.items()]
    for key in ("original_bam", "reference", "reference_fai", "downsampled_bam",
                "downsampled_bai", "regions"):
        if key not in identities:
            raise ValueError(f"missing prepared identity: {key}")
        entries.append((key, identities[key]))
    expected = {}
    for label, record in entries:
        if not isinstance(record, dict):
            raise ValueError(f"invalid prepared identity: {label}")
        path, sha, size = record.get("path"), record.get("sha256"), record.get("bytes")
        if (not isinstance(path, str) or not path or not Path(path).is_absolute()
                or not isinstance(sha, str) or not re.fullmatch(r"[0-9a-f]{64}", sha)
                or type(size) is not int or size < 0):
            raise ValueError(f"incomplete prepared byte identity: {label}")
        normalized = Path(os.path.abspath(path))
        value = (sha, size)
        if normalized in expected and expected[normalized] != value:
            raise ValueError(f"conflicting prepared identities for {normalized}")
        expected[normalized] = value
    checked_bytes = 0
    for path, (sha, size) in expected.items():
        before = _regular(path).stat()
        if before.st_size != size:
            raise ValueError(f"prepared input size changed: {path}")
        observed = _digest(path)
        after = _regular(path).stat()
        stable = lambda value: (value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns)
        if observed != sha or stable(before) != stable(after):
            raise ValueError(f"prepared input SHA256 changed or file mutated during verification: {path}")
        checked_bytes += size
    return {"files_checked": len(expected), "bytes_checked": checked_bytes}


def prepare_inputs(resources_path: Path, protocol: dict, out: Path) -> dict:
    """Acquire five locked resources and create one reduced-depth interval BAM.

    Source BAM bytes remain unchanged. Metrics count physical records, not unique
    templates. The reference's compressed indices are intentionally not acquired:
    a new full, uncompressed FASTA and its own faidx are produced and hashed.
    """
    intervals = _validate_protocol(protocol)
    rows, lock_bytes = _source_lock(resources_path)
    out = Path(out).absolute()
    out.mkdir(parents=True, exist_ok=False)
    budget = _Budget(out, protocol["limits"])
    sample = protocol.get("sample", "HG002")
    (out / "resources.tsv").write_bytes(lock_bytes)
    sources = {row["id"]: _download(row, out, budget) for row in rows}
    coverage = benchmark_coverage(Path(sources["truth_bed"]["path"]), intervals, budget)
    if protocol.get("require_nonempty_benchmark_coverage", False) and coverage["empty_block_ids"]:
        raise ValueError("benchmark coverage preflight refused: selected blocks have no confident bases: "
                         + ", ".join(coverage["empty_block_ids"]))
    pysam = importlib.import_module("pysam")
    reference = out / "reference.fa"
    with _read_regular(sources["reference"]["path"]) as compressed:
        with gzip.GzipFile(fileobj=compressed, mode="rb") as source, reference.open("xb") as output:
            while True:
                budget.check_time()
                chunk = source.read(CHUNK_BYTES)
                if not chunk:
                    break
                budget.check_scratch(len(chunk))
                output.write(chunk)
    budget.check_scratch()
    pysam.faidx(str(reference))
    budget.check_scratch()
    original_bam = Path(sources["reads_bam"]["path"])
    reduced_bam = out / "downsampled.bam"
    with pysam.AlignmentFile(str(original_bam), "rb", index_filename=sources["reads_bai"]["path"]) as source:
        with pysam.FastaFile(str(reference)) as fasta:
            groups = _audit_inputs(source, fasta, sample, intervals)
        dictionary_count = len(source.references)
        metrics = _downsample(source, reduced_bam, pysam, intervals, groups,
                              sources["reads_bam"]["sha256"], protocol["downsample"], budget)
    budget.check_scratch()
    pysam.index(str(reduced_bam))
    budget.check_scratch()
    regions = out / "regions.bed"
    with regions.open("x") as output:
        for item in intervals:
            output.write("\t".join(str(item[key]) for key in
                         ("contig", "start", "end", "partition", "block_id")) + "\n")
    identities = {"resource_lock_sha256": hashlib.sha256(lock_bytes).hexdigest(),
                  "sources": sources, "sample": sample, "dictionary_contigs": dictionary_count,
                  "downsample": {**protocol["downsample"], "algorithm": "sha256-json-template-v1"},
                  "benchmark_coverage": {"truth_bed_sha256": sources["truth_bed"]["sha256"],
                                         "required_nonempty": protocol.get(
                                             "require_nonempty_benchmark_coverage", False),
                                         **coverage},
                  "pysam_version": pysam.__version__}
    for key, path in (("original_bam", original_bam), ("reference", reference),
                      ("reference_fai", Path(str(reference) + ".fai")),
                      ("downsampled_bam", reduced_bam),
                      ("downsampled_bai", Path(str(reduced_bam) + ".bai")), ("regions", regions)):
        identities[key] = {"path": str(path), "sha256": _digest(path, budget), "bytes": path.stat().st_size}
    for value in sources.values():
        if _digest(value["path"], budget) != value["sha256"]:
            raise ValueError(f"source changed during preparation: {value['id']}")
    budget.check_scratch()
    metrics.update({"download_bytes": budget.download_bytes, "download_attempts": budget.attempts,
                    "scratch_high_water_bytes": budget.scratch_high_water,
                    "wall_seconds": time.monotonic() - budget.started,
                    "benchmark_coverage": coverage,
                    "interval_bases": sum(x["end"] - x["start"] for x in intervals)})
    result = {"reference": str(reference), "original_bam": str(original_bam),
              "downsampled_bam": str(reduced_bam), "sample": sample,
              "truth_vcf": sources["truth_vcf"]["path"], "truth_bed": sources["truth_bed"]["path"],
              "regions": str(regions), "identities": identities, "metrics": metrics}
    encoded = (json.dumps(result, indent=2, sort_keys=True) + "\n").encode()
    budget.check_scratch(len(encoded))
    partial = out / ".preparation.json.partial"
    with partial.open("xb") as output:
        output.write(encoded)
        output.flush()
        os.fsync(output.fileno())
    budget.check_time()
    os.link(partial, out / "preparation.json")
    partial.unlink()
    return result
