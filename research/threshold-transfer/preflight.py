"""Hosted, outcome-blind metadata admission for the exact GIAB no-RG source.

Only the pinned BAI and bounded BAM header/indexed-read ranges are acquired.
No reference sequence, truth labels or model are loaded. Read bases and qualities
are not inspected, used, or retained in the metadata receipt.
The finite probe establishes input compatibility, not global identity uniqueness.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import threading
import time

HEADER_SHA256 = "2a507e6a549c0e6b62bdccb14da7c4b10dec3c7882fb04e357782672ceabc856"
DICTIONARY_SHA256 = "64416f59da1f276495e898e8e1d5a6f071f463bd081aa0f561068f6c62714da7"
EXPECTED_BAM = {
    "url": "https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/"
           "NIST_NA12878_HG001_HiSeq_300x/NHGRI_Illumina300X_novoalign_bams/"
           "HG001.GRCh38_full_plus_hs38d1_analysis_set_minus_alts.300x.bam",
    "bytes": 589702020566, "md5": "f48bf321c34f32b7c08feb17b9e5b220",
    "last_modified": "Tue, 08 Mar 2016 18:53:34 GMT",
}
EXPECTED_BAI = {
    "url": EXPECTED_BAM["url"] + ".bai", "bytes": 12006264,
    "md5": "00fd8ac75a92daf341e79008e7e9d900",
}
MAX_DOWNLOAD_BYTES = 64 << 20
MAX_WALL_SECONDS = 120
MAX_PROBE_RECORDS = 64
PROBE_BP = 1000
QNAME_RE = re.compile(r"[A-Za-z0-9_.-]+:[0-9]+:[A-Za-z0-9_.-]+:[1-9][0-9]*:[0-9]+:-?[0-9]+:-?[0-9]+")


def _adapter():
    # Deferred import avoids a cycle when the acquisition adapter imports us.
    import prepare_fresh
    return prepare_fresh


def _canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def source_identity_sha256(source):
    if source != EXPECTED_BAM:
        raise ValueError("no-RG admission requires the exact reviewed source object")
    return hashlib.sha256(_canonical(source).encode()).hexdigest()


def technical_read_group_id(source):
    """A technical source namespace, never an inferred original library ID."""
    return "rosalind_source_" + source_identity_sha256(source)


def validate_header_identity(header, source):
    identity = source_identity_sha256(source)
    raw = str(header).encode()
    sha = hashlib.sha256(raw).hexdigest()
    dictionary = header.to_dict()
    ordered = [[row["SN"], int(row["LN"])] for row in dictionary.get("SQ", [])]
    ordered_sha = hashlib.sha256(_canonical(ordered).encode()).hexdigest()
    if sha != HEADER_SHA256:
        raise ValueError("source SAM header differs from the exact previously observed header")
    if dictionary.get("RG"):
        raise ValueError("reviewed no-RG source unexpectedly has header read groups")
    if ordered_sha != DICTIONARY_SHA256 or ["chr20", 64444167] not in ordered:
        raise ValueError("source ordered dictionary or chr20 metadata differs")
    return {"source_header_sha256": sha, "source_header_bytes": len(raw),
            "source_identity_sha256": identity, "dictionary_sha256": ordered_sha,
            "contigs": len(ordered), "chr20_length": 64444167,
            "header_read_groups": [], "declared_sort_order": dictionary.get("HD", {}).get("SO"),
            "technical_read_group": technical_read_group_id(source),
            "technical_group_scope": "source-level selection namespace; original library/read-group identity unavailable",
            "reference_sequence_verified": False}


def validate_cluster_qname(name):
    """Validate syntax only; return the entire unchanged physical-cluster name."""
    if not isinstance(name, str) or len(name) > 254 or QNAME_RE.fullmatch(name) is None:
        raise ValueError("record QNAME lacks the complete seven-field Illumina cluster identity")
    return name


def validate_no_rg_record_identity(record):
    """Return the unmodified complete physical-cluster QNAME or reject."""
    if record.has_tag("RG"):
        raise ValueError("record carries RG despite the admitted no-RG source header")
    name = validate_cluster_qname(record.query_name)
    if not record.is_secondary and not record.is_supplementary:
        if not record.is_paired or bool(record.is_read1) == bool(record.is_read2):
            raise ValueError("primary record must be paired with exactly one read1/read2 flag")
    return name


class _MetadataBudget:
    """Reserve from both the 64 MiB metadata cap and the enclosing study cap."""
    def __init__(self, root, parent=None):
        adapter = _adapter()
        self.local = adapter._Budget(root, {"download_bytes": MAX_DOWNLOAD_BYTES,
            "scratch_bytes": 128 << 20, "analysis_seconds": MAX_WALL_SECONDS - 10,
            "network_timeout_seconds": 5})
        self.parent = parent
        self.lock = threading.Lock()
        if parent is not None:
            self.local.deadline = min(self.local.deadline, parent.deadline - 10)

    def __getattr__(self, name):
        return getattr(self.local, name)

    def read(self, response, size):
        with self.lock:
            self.local.check_time()
            remaining = self.local.download_limit - self.local.download_bytes
            if remaining <= 0:
                raise ValueError("metadata preflight acquisition byte budget exceeded")
            if self.parent is None:
                return self.local.read(response, min(size, remaining))
            chunk = self.parent.read(response, min(size, remaining))
            self.local.download_bytes += len(chunk)
            return chunk

    def check_scratch(self):
        self.local.check_scratch()
        if self.parent is not None:
            self.parent.check_scratch()


def _probe_worker(request_path):
    import pysam
    request = json.loads(Path(request_path).read_text())
    out = Path(request["out"])
    report = {"accepted": False, "records": [], "outcome_data_loaded": False,
              "scope": "capped indexed metadata probe; not global template-identity uniqueness or reference verification"}
    try:
        with pysam.AlignmentFile(request["url"], "rb", index_filename=request["bai"]) as source:
            (out / "source-header.sam").write_text(str(source.header))
            report["observed_header_sha256"] = hashlib.sha256(str(source.header).encode()).hexdigest()
            report.update(validate_header_identity(source.header, request["source"]))
            if not source.has_index():
                raise ValueError("pinned BAI cannot index the source BAM")
            region = request["region"]
            report["probe_region"] = region
            previous = -1
            primary_ends = {}
            errors = []
            iterator = source.fetch(region["contig"], region["start"], region["end"])
            for _ in range(MAX_PROBE_RECORDS):
                try:
                    record = next(iterator)
                except StopIteration:
                    break
                metadata = {"qname": record.query_name, "flag": record.flag,
                    "has_rg": record.has_tag("RG"),
                    "rg_value": record.get_tag("RG") if record.has_tag("RG") else None,
                    "reference_start": record.reference_start,
                    "is_secondary": record.is_secondary, "is_supplementary": record.is_supplementary,
                    "is_read1": record.is_read1, "is_read2": record.is_read2}
                report["records"].append(metadata)
                try:
                    validate_no_rg_record_identity(record)
                    if record.reference_start < previous:
                        raise ValueError("indexed probe coordinates are not monotone")
                except ValueError as exc:
                    errors.append({"record_index": len(report["records"]) - 1, "reason": str(exc)})
                previous = record.reference_start
                if not record.is_secondary and not record.is_supplementary:
                    ends = primary_ends.setdefault(record.query_name, [])
                    ends.append(1 if record.is_read1 else 2 if record.is_read2 else 0)
            report["probe_records"] = len(report["records"])
            report["read_identity_errors"] = errors
            report["primary_templates_both_ends_observed"] = sum(set(ends) == {1, 2} for ends in primary_ends.values())
            report["primary_templates_missing_mate_in_capped_probe"] = sum(set(ends) != {1, 2} for ends in primary_ends.values())
            report["repeated_primary_end_observations"] = sum(len(ends) - len(set(ends)) for ends in primary_ends.values())
            report["physical_records_deduplicated"] = False
            if not report["records"]:
                raise ValueError("fixed indexed metadata probe returned no records; no adaptive relocation")
            if not primary_ends:
                raise ValueError("fixed metadata probe has no primary paired-record identity support")
            if errors:
                raise ValueError("indexed metadata probe has incompatible record identities")
            report["accepted"] = True
    except BaseException as exc:
        report["failure"] = {"type": type(exc).__name__, "reason": str(exc)}
        raise
    finally:
        _write(out / "read-metadata.json", report)


def run_preflight(protocol, out, budget=None):
    """Return a metadata receipt with a reusable verified BAI; failures persist."""
    if os.environ.get("GITHUB_ACTIONS") != "true":
        raise RuntimeError("metadata acquisition and native inspection require hosted execution")
    adapter = _adapter()
    intervals, _ = adapter._validate(protocol)
    sources = protocol["fresh"].get("sources", adapter.DEFAULT_SOURCES)
    if sources["reads_bam"] != EXPECTED_BAM or sources["reads_bai"] != EXPECTED_BAI:
        raise ValueError("preflight requires the pinned BAM and matching BAI")
    out = Path(out).absolute()
    out.mkdir(parents=True, exist_ok=False)
    metadata_budget = _MetadataBudget(out, budget)
    started = time.monotonic()
    ranges, errors = [], []
    report = {"accepted": False, "source": sources["reads_bam"],
        "limits": {"download_bytes": MAX_DOWNLOAD_BYTES, "wall_seconds": MAX_WALL_SECONDS,
                   "probe_records": MAX_PROBE_RECORDS, "probe_bp": PROBE_BP},
        "reference_sequence_verified": False, "outcome_data_loaded": False,
        "remote_bam_full_md5_verified": False,
        "scope": "metadata compatibility only; full extraction must revalidate every record"}
    try:
        report["bai_download"] = adapter._download(sources["reads_bai"], out / "source.bam.bai",
            metadata_budget, out / "download-reads_bai.json")
        first = intervals[0]
        region = {"contig": first["contig"], "start": first["start"],
                  "end": min(first["end"], first["start"] + PROBE_BP)}
        with adapter._range_relay(sources["reads_bam"], metadata_budget) as (url, ranges, errors):
            request = {"url": url, "source": sources["reads_bam"],
                       "bai": report["bai_download"]["path"], "out": str(out), "region": region}
            request_path = out / "probe-request.json"
            _write(request_path, request)
            subprocess.run([sys.executable, str(Path(__file__).resolve()), "--probe-worker", str(request_path)],
                           check=True, timeout=max(0.001, metadata_budget.deadline - time.monotonic()))
        if errors:
            raise ValueError("metadata range relay reported acquisition errors")
        worker = json.loads((out / "read-metadata.json").read_text())
        if not worker.get("accepted"):
            raise ValueError("metadata worker did not admit the source")
        if time.monotonic() - started > MAX_WALL_SECONDS:
            raise TimeoutError("metadata preflight exceeded wall-clock limit")
        report.update({key: worker[key] for key in ("source_header_sha256", "source_identity_sha256",
            "technical_read_group", "dictionary_sha256", "chr20_length", "probe_records")})
        report["read_metadata_path"] = str(out / "read-metadata.json")
        report["accepted"] = True
    except BaseException as exc:
        report["failure"] = {"type": type(exc).__name__, "reason": str(exc)}
        raise
    finally:
        report["wall_seconds"] = time.monotonic() - started
        report["download_bytes"] = metadata_budget.download_bytes
        report["received_bam_bytes"] = sum(item["received_bytes"] for item in ranges)
        report["receipt_path"] = str(out / "preflight.json")
        _write(out / "bam-byte-ranges.json", ranges)
        _write(out / "range-errors.json", errors)
        _write(out / "preflight.json", report)
    files = {path.name: _file_identity(path) for path in sorted(out.iterdir()) if path.is_file()}
    _write(out / "completion.json", {"schema": "rosalind-metadata-preflight-v1", "files": files})
    return report


def _file_identity(path):
    path = Path(path)
    if path.is_symlink() or not path.is_file():
        raise ValueError("preflight receipt must be a regular file")
    sha, md5, size = hashlib.sha256(), hashlib.md5(), 0
    with path.open("rb") as handle:
        while chunk := handle.read(1 << 20):
            size += len(chunk)
            sha.update(chunk)
            md5.update(chunk)
    return {"bytes": size, "sha256": sha.hexdigest(), "md5": md5.hexdigest()}


def validate_preflight_receipt(directory, protocol):
    """Revalidate exact metadata artifacts before the expensive acquisition."""
    directory = Path(directory).absolute()
    manifest = json.loads((directory / "completion.json").read_text())
    if manifest.get("schema") != "rosalind-metadata-preflight-v1":
        raise ValueError("preflight completion manifest missing or unsupported")
    required = {"source.bam.bai", "source-header.sam", "read-metadata.json", "download-reads_bai.json",
                "probe-request.json", "bam-byte-ranges.json", "range-errors.json", "preflight.json"}
    files = manifest.get("files", {})
    if set(files) != required:
        raise ValueError("preflight manifest does not bind every required receipt")
    for name, expected in files.items():
        if _file_identity(directory / name) != expected:
            raise ValueError(f"preflight receipt identity changed: {name}")
    adapter = _adapter()
    intervals, _ = adapter._validate(protocol)
    sources = protocol["fresh"].get("sources", adapter.DEFAULT_SOURCES)
    report = json.loads((directory / "preflight.json").read_text())
    worker = json.loads((directory / "read-metadata.json").read_text())
    download = json.loads((directory / "download-reads_bai.json").read_text())
    request = json.loads((directory / "probe-request.json").read_text())
    ranges = json.loads((directory / "bam-byte-ranges.json").read_text())
    errors = json.loads((directory / "range-errors.json").read_text())
    first = intervals[0]
    region = {"contig": first["contig"], "start": first["start"],
              "end": min(first["end"], first["start"] + PROBE_BP)}
    if (sources["reads_bam"] != EXPECTED_BAM or sources["reads_bai"] != EXPECTED_BAI
            or report.get("source") != EXPECTED_BAM or request.get("source") != EXPECTED_BAM
            or request.get("region") != region or worker.get("probe_region") != region):
        raise ValueError("preflight receipt source or probe interval mismatch")
    expected = {"source_header_sha256": HEADER_SHA256, "dictionary_sha256": DICTIONARY_SHA256,
                "chr20_length": 64444167, "source_identity_sha256": source_identity_sha256(EXPECTED_BAM),
                "technical_read_group": technical_read_group_id(EXPECTED_BAM)}
    if any(report.get(key) != value or worker.get(key) != value for key, value in expected.items()):
        raise ValueError("preflight receipt source identity does not match exact admitted metadata")
    if files["source-header.sam"]["sha256"] != HEADER_SHA256:
        raise ValueError("retained source header differs from previously observed header")
    if (report.get("accepted") is not True or worker.get("accepted") is not True
            or report.get("failure") or worker.get("failure") or errors or worker.get("read_identity_errors")
            or report.get("outcome_data_loaded") is not False or worker.get("outcome_data_loaded") is not False
            or report.get("reference_sequence_verified") is not False
            or report.get("remote_bam_full_md5_verified") is not False):
        raise ValueError("preflight does not contain an accepted metadata-only probe")
    records = worker.get("records", [])
    if not 0 < len(records) <= MAX_PROBE_RECORDS or report.get("probe_records") != len(records):
        raise ValueError("preflight record cap or count mismatch")
    previous, primary = -1, 0
    for record in records:
        validate_cluster_qname(record.get("qname"))
        if record.get("has_rg") is not False or record.get("rg_value") is not None:
            raise ValueError("preflight metadata includes unexplained record RG")
        if record["reference_start"] < previous:
            raise ValueError("preflight indexed metadata coordinates are not monotone")
        previous = record["reference_start"]
        flag = record["flag"]
        if not flag & (256 | 2048):
            primary += 1
            if not flag & 1 or bool(flag & 64) == bool(flag & 128):
                raise ValueError("preflight primary paired-record flags are invalid")
    if not primary:
        raise ValueError("preflight has no primary paired-record identity support")
    bai = files["source.bam.bai"]
    if (bai["bytes"] != EXPECTED_BAI["bytes"] or bai["md5"] != EXPECTED_BAI["md5"]
            or report["bai_download"].get("sha256") != bai["sha256"]
            or report["bai_download"].get("bytes") != bai["bytes"]
            or report["bai_download"].get("published_md5_verified") is not True
            or download.get("source") != EXPECTED_BAI
            or download.get("published_md5_verified") is not True
            or download.get("received_sha256") != bai["sha256"]
            or download.get("received_md5") != bai["md5"]
            or download.get("received_bytes") != bai["bytes"] or download.get("failure")):
        raise ValueError("preflight BAI identity or published checksum mismatch")
    received = sum(row["received_bytes"] for row in ranges)
    if (report.get("received_bam_bytes") != received
            or report.get("download_bytes") != bai["bytes"] + received
            or not 0 < report["download_bytes"] <= MAX_DOWNLOAD_BYTES
            or not 0 <= report["wall_seconds"] <= MAX_WALL_SECONDS):
        raise ValueError("preflight exceeded its acquisition boundary or accounting differs")
    return dict(report, bai_path=str(directory / "source.bam.bai"), verified_preflight_files=files)


def main():
    if len(sys.argv) == 3 and sys.argv[1] == "--probe-worker":
        if os.environ.get("GITHUB_ACTIONS") != "true":
            raise RuntimeError("native metadata worker requires hosted execution")
        _probe_worker(sys.argv[2])
        return
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    run_preflight(json.loads(args.protocol.read_text()), args.output)


if __name__ == "__main__":
    main()
