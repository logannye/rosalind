"""Hosted-only, bounded acquisition and de novo evidence for independent HG001.

``prepare_fresh(protocol, out, binary)`` returns records, candidate_recall, metrics,
identities, and state. It accepts ``protocol['fresh']`` with sources, test-only
intervals, conditions (name/numerator/denominator), and seed. Shared candidate,
extraction, rich_fields and limits keys follow the P0 protocol. Imports are inert;
native reads, extraction, caller execution and scientific calculations are hosted.
"""
from __future__ import annotations

from contextlib import contextmanager
import gzip
import hashlib
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import threading
import time
import urllib.request

HERE = Path(__file__).resolve().parent
P0 = HERE.parent / "evidence-sufficiency"
sys.path.insert(0, str(P0))
import prepare
import run_pilot as pilot

READ_BASE = ("https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/"
             "NIST_NA12878_HG001_HiSeq_300x/NHGRI_Illumina300X_novoalign_bams/")
BAM_NAME = "HG001.GRCh38_full_plus_hs38d1_analysis_set_minus_alts.300x.bam"
TRUTH_BASE = ("https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/"
              "NA12878_HG001/NISTv4.2.1/GRCh38/")
REFERENCE_URL = ("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/"
                 "GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/"
                 "GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz")
DEFAULT_SOURCES = {
    "reads_bam": {"url": READ_BASE + BAM_NAME, "bytes": 589702020566,
                  "md5": "f48bf321c34f32b7c08feb17b9e5b220",
                  "last_modified": "Tue, 08 Mar 2016 18:53:34 GMT"},
    "reads_bai": {"url": READ_BASE + BAM_NAME + ".bai", "bytes": 12006264,
                  "md5": "00fd8ac75a92daf341e79008e7e9d900"},
    "reference": {"url": REFERENCE_URL, "bytes": 874637925,
                  "md5": "a056c57649f3c9964c68aead3849bbf8"},
    "truth_vcf": {"url": TRUTH_BASE + "HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz",
                  "bytes": 125932193, "md5": "ec61d76ad725d808e3a9e1c8acec402a"},
    "truth_bed": {"url": TRUTH_BASE + "HG001_GRCh38_1_22_v4.2.1_benchmark.bed",
                  "bytes": 15479939, "md5": "e1ea90c830f981d29a46043eb618e7bf"},
}


def _identity(path):
    path = Path(path).absolute()
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": pilot.digest(path)}


def verify_input_snapshot(snapshot):
    """Recheck identities captured before scientific consumption, never rebase."""
    total = 0
    for entry in snapshot.values():
        path = prepare._regular(entry["path"])
        if path.stat().st_size != entry["bytes"] or pilot.digest(path) != entry["sha256"]:
            raise ValueError(f"fresh consumed input changed: {path.name}")
        total += entry["bytes"]
    return {"passed": True, "files_verified": len(snapshot), "bytes_verified": total}


def _validate(protocol):
    fresh = protocol["fresh"]
    intervals = fresh["intervals"]
    if not intervals or any(r["partition"] != "test" or r["contig"] != "chr20" for r in intervals):
        raise ValueError("fresh HG001 intervals must be chr20 test-only")
    ordered = sorted(intervals, key=lambda r: r["start"])
    if any(type(r["start"]) is not int or type(r["end"]) is not int
           or not 0 <= r["start"] < r["end"] <= 64444167 for r in ordered):
        raise ValueError("fresh interval outside GRCh38 chr20")
    if len({r["block_id"] for r in ordered}) != len(ordered):
        raise ValueError("fresh interval block identifiers must be distinct")
    if any(a["end"] > b["start"] for a, b in zip(ordered, ordered[1:])):
        raise ValueError("fresh intervals overlap")
    conditions = fresh["conditions"]
    names = [c["name"] for c in conditions]
    if len(names) != len(set(names)) or not conditions:
        raise ValueError("conditions must have distinct names")
    for condition in conditions:
        if not re.fullmatch(r"[a-z][a-z0-9_]*", condition["name"]):
            raise ValueError("unsafe condition name")
        if (type(condition["numerator"]) is not int or type(condition["denominator"]) is not int
                or not 0 < condition["numerator"] < condition["denominator"]):
            raise ValueError("condition retention must be strictly between zero and one")
    if any(c["denominator"] != 300 for c in conditions):
        raise ValueError("source retention denominator must be the prespecified nominal 300x")
    if fresh.get("sources", DEFAULT_SOURCES) != DEFAULT_SOURCES:
        raise ValueError("fresh sources differ from reviewed authoritative metadata")
    for key in ("download_bytes", "scratch_bytes", "analysis_seconds"):
        prepare._positive_int(protocol["limits"][key], key)
    if protocol["limits"]["download_bytes"] > 4 * (1 << 30):
        raise ValueError("fresh combined acquisition ceiling is 4 GiB")
    if protocol["limits"]["scratch_bytes"] > 30 * (1 << 30):
        raise ValueError("fresh scratch ceiling is 30 GiB")
    return ordered, conditions


class _Budget(prepare._Budget):
    def __init__(self, root, limits):
        super().__init__(root, limits)
        self.lock = threading.Lock()

    def read(self, response, size):
        """Reserve/read under one lock, so concurrent ranges cannot overshoot."""
        with self.lock:
            self.check_time()
            remaining = self.download_limit - self.download_bytes
            if remaining <= 0:
                raise ValueError("combined fresh acquisition byte budget exceeded")
            chunk = response.read(min(size, remaining))
            self.download_bytes += len(chunk)
            return chunk


def _request(url, **kwargs):
    if not url.startswith("https://"):
        raise ValueError("source transport must be HTTPS")
    return urllib.request.build_opener(prepare._HttpsRedirect()).open(
        urllib.request.Request(url, **kwargs), timeout=30)


def _download(source, path, budget):
    started = time.monotonic()
    md5, sha = hashlib.md5(), hashlib.sha256()
    with _request(source["url"], headers={"Accept-Encoding": "identity"}) as response:
        if response.status != 200 or int(response.headers.get("Content-Length", -1)) != source["bytes"]:
            raise ValueError("official resource length or status changed")
        with path.open("xb") as handle:
            while chunk := budget.read(response, 1 << 20):
                handle.write(chunk)
                md5.update(chunk)
                sha.update(chunk)
                if handle.tell() % (128 << 20) < 1 << 20:
                    budget.check_scratch()
    if path.stat().st_size != source["bytes"] or md5.hexdigest() != source["md5"]:
        raise ValueError("official published resource MD5 or size mismatch")
    budget.check_scratch()
    return {**source, "path": str(path), "sha256": sha.hexdigest(),
            "published_md5_verified": True, "wall_seconds": time.monotonic() - started}


def _check_bam_headers(headers, source, start=None, end=None):
    if headers.get("Last-Modified") != source["last_modified"]:
        raise ValueError("remote BAM Last-Modified changed")
    if start is None:
        if int(headers.get("Content-Length", -1)) != source["bytes"]:
            raise ValueError("remote BAM size changed")
    elif headers.get("Content-Range") != f"bytes {start}-{end}/{source['bytes']}":
        raise ValueError("remote BAM range identity mismatch")


@contextmanager
def _range_relay(source, budget):
    """Account every received BAM byte while htslib performs indexed retrieval.

    Native htslib reads a loopback URL. Each upstream request is a closed HTTPS
    range of at most 8 MiB and must return 206 with exact Content-Range. Full-body
    fallback is rejected. Byte limits include bytes read before client disconnect.
    The published whole-BAM MD5 is retained as metadata, never claimed verified.
    """
    with _request(source["url"], method="HEAD", headers={"Accept-Encoding": "identity"}) as response:
        if response.status != 200 or response.headers.get("Accept-Ranges") != "bytes":
            raise ValueError("remote BAM does not advertise byte ranges")
        _check_bam_headers(response.headers, source)
    logs, errors, stopping = [], [], threading.Event()
    log_lock = threading.Lock()

    class Handler(BaseHTTPRequestHandler):
        protocol_version = "HTTP/1.0"

        def log_message(self, *args):
            pass

        def do_HEAD(self):
            if self.path != "/source.bam":
                self.send_error(404)
                return
            self.send_response(200)
            self.send_header("Content-Length", source["bytes"])
            self.send_header("Accept-Ranges", "bytes")
            self.send_header("Last-Modified", source["last_modified"])
            self.end_headers()

        def do_GET(self):
            if self.path != "/source.bam":
                self.send_error(404)
                return
            incoming = self.headers.get("Range")
            match = re.fullmatch(r"bytes=([0-9]+)-([0-9]*)", incoming or "bytes=0-")
            if not match:
                self.send_error(416)
                return
            start = int(match[1])
            end = int(match[2]) if match[2] else source["bytes"] - 1
            if not 0 <= start <= end < source["bytes"]:
                self.send_error(416)
                return
            self.send_response(206 if incoming else 200)
            self.send_header("Content-Length", end - start + 1)
            self.send_header("Accept-Ranges", "bytes")
            if incoming:
                self.send_header("Content-Range", f"bytes {start}-{end}/{source['bytes']}")
            self.end_headers()
            try:
                while start <= end and not stopping.is_set():
                    stop = min(end, start + (8 << 20) - 1)
                    sha, received = hashlib.sha256(), 0
                    request_headers = {"Range": f"bytes={start}-{stop}",
                                       "If-Unmodified-Since": source["last_modified"],
                                       "Accept-Encoding": "identity"}
                    try:
                        with _request(source["url"], headers=request_headers) as response:
                            if response.status != 206:
                                raise ValueError("refusing non-range BAM response; whole genome transfer forbidden")
                            _check_bam_headers(response.headers, source, start, stop)
                            while not stopping.is_set():
                                remaining = stop - start + 1 - received
                                if remaining == 0:
                                    break
                                chunk = budget.read(response, min(65536, remaining))
                                if not chunk:
                                    break
                                received += len(chunk)
                                sha.update(chunk)
                                self.wfile.write(chunk)
                            if not stopping.is_set() and received != stop - start + 1:
                                raise OSError("short BAM byte range")
                    finally:
                        with log_lock:
                            logs.append({"start": start, "requested_end_inclusive": stop,
                                         "received_bytes": received, "received_sha256": sha.hexdigest(),
                                         "complete": received == stop - start + 1})
                    start = stop + 1
            except (BrokenPipeError, ConnectionResetError):
                pass  # htslib ends a range stream when it seeks or finishes.
            except BaseException as exc:
                with log_lock:
                    errors.append({"type": type(exc).__name__, "reason": str(exc)})

    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    server.daemon_threads = False
    server.block_on_close = True
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        yield f"http://127.0.0.1:{server.server_port}/source.bam", logs, errors
    finally:
        stopping.set()
        server.shutdown()
        server.server_close()
        thread.join(timeout=2)


def _extract_bams_worker(request_path):
    """Separate process keeps htslib network waits independent of relay threads."""
    import pysam
    request = json.loads(Path(request_path).read_text())
    intervals, conditions = request["intervals"], request["conditions"]
    paths = {c["name"]: str(Path(request["out"]) / (c["name"] + ".bam")) for c in conditions}
    metrics = {name: {"records": 0, "query_bases": 0} for name in paths}
    traversed, selected = 0, 0
    source_identity = hashlib.sha256(pilot.canonical(request["source"]).encode()).hexdigest()
    with pysam.AlignmentFile(request["url"], "rb", index_filename=request["bai"]) as source:
        groups = source.header.to_dict().get("RG", [])
        samples = {group.get("SM") for group in groups}
        if len(samples) != 1 or not samples <= {"HG001", "NA12878"}:
            raise ValueError("fresh BAM sample is not the single declared HG001/NA12878 specimen")
        sample = next(iter(samples))
        with pysam.FastaFile(request["reference"]) as fasta:
            group_ids = prepare._audit_inputs(source, fasta, sample, intervals)
            chr20 = fasta.fetch("chr20").upper().encode()
            chr20_md5 = hashlib.md5(chr20).hexdigest()
            header_m5 = next(sq for sq in source.header.to_dict()["SQ"] if sq["SN"] == "chr20").get("M5")
            if header_m5 and header_m5.lower() != chr20_md5:
                raise ValueError("chr20 source header M5 differs from acquired reference sequence")
        writers = {name: pysam.AlignmentFile(path, "wb", header=source.header) for name, path in paths.items()}
        try:
            for region_index, region in enumerate(intervals):
                for record in source.fetch(region["contig"], region["start"], region["end"]):
                    traversed += 1
                    start, end = record.reference_start, record.reference_end
                    if start < 0 or end is None or end <= start:
                        continue
                    # A read spanning more than one selected interval belongs to
                    # its first overlapping interval. Do not deduplicate by QNAME.
                    if any(start < previous["end"] and end > previous["start"] for previous in intervals[:region_index]):
                        continue
                    if not record.has_tag("RG") or record.get_tag("RG") not in group_ids:
                        raise ValueError("fresh selected record has missing or unexpected RG")
                    selected += 1
                    for condition in conditions:
                        if prepare._keep_template(request["seed"], source_identity, record.get_tag("RG"),
                                                  record.query_name, condition["numerator"], condition["denominator"]):
                            name = condition["name"]
                            writers[name].write(record)
                            metrics[name]["records"] += 1
                            metrics[name]["query_bases"] += record.query_length or 0
        finally:
            for writer in writers.values():
                writer.close()
    for path in paths.values():
        pysam.index(path)
    result = {"bams": paths, "sample": sample, "conditions": metrics,
              "traversed_records": traversed, "selected_records": selected,
              "chr20_md5": chr20_md5, "bam_chr20_header_m5": header_m5,
              "full_ordered_bam_reference_dictionary_match": True,
              "template_identity_sha256": source_identity,
              "downsampling_algorithm": "sha256-json-template-v1; same hash nested at all fractions"}
    pilot.write_json(Path(request["out"]) / "bam-preparation.json", result)


def counts_features(row, alt):
    """Compute exactly F0 from measured depths+alleles, without rich-field defaults."""
    from features import BASES, _base, _required, _nullable, _ratio
    if row.get("observed") is False:
        raise ValueError("counts features require observed data")
    ref, alt = _base(row.get("ref"), "ref"), _base(alt, "alt")
    if ref == alt:
        raise ValueError("ALT must differ from REF")
    depth = _required(row, "callable_depth")
    counts = [_required(row, base.lower()) for base in BASES]
    if sum(counts) != depth:
        raise ValueError("counts do not reconcile with callable depth")
    ref_count, alt_count = counts[BASES.index(ref)], counts[BASES.index(alt)]
    other_count = depth - ref_count - alt_count
    features = {}
    for name, value in (("callable_depth", depth), ("ref_count", ref_count),
                        ("alt_count", alt_count), ("other_count", other_count)):
        features[name], features[name + "_log1p"] = float(value), math.log1p(value)
    for name, count in (("ref", ref_count), ("alt", alt_count), ("other", other_count)):
        _nullable(features, name + "_fraction", _ratio(count, depth))
    for base in BASES:
        features["ref_" + base], features["alt_" + base] = float(ref == base), float(alt == base)
    return features


def prepare_fresh(protocol, out, binary):
    if os.environ.get("GITHUB_ACTIONS") != "true":
        raise RuntimeError("fresh scientific acquisition and execution require hosted workflow")
    import pysam
    from labels import TruthIndex
    intervals, conditions = _validate(protocol)
    fresh = protocol["fresh"]
    out, binary = Path(out).absolute(), Path(binary).resolve(strict=True)
    out.mkdir(parents=True, exist_ok=False)
    inputs = out / "inputs"
    inputs.mkdir()
    budget = _Budget(out, protocol["limits"])
    sources = fresh.get("sources", DEFAULT_SOURCES)
    acquired = {}
    paths = {"reads_bai": "source.bam.bai", "reference": "reference.fa.gz",
             "truth_vcf": "truth.vcf.gz", "truth_bed": "truth.bed"}
    for name, filename in paths.items():
        acquired[name] = _download(sources[name], inputs / filename, budget)
    coverage = prepare.benchmark_coverage(inputs / "truth.bed", intervals, budget)
    if coverage["empty_block_ids"]:
        raise ValueError("fresh selected interval has no benchmark confident coverage; no outcome-driven relocation")
    reference = inputs / "reference.fa"
    with gzip.open(inputs / "reference.fa.gz", "rb") as source, reference.open("xb") as dest:
        while chunk := source.read(1 << 20):
            dest.write(chunk)
            if dest.tell() % (128 << 20) < 1 << 20:
                budget.check_scratch()
    pysam.faidx(str(reference))
    regions = inputs / "regions.bed"
    with regions.open("x") as handle:
        for region in intervals:
            handle.write("\t".join(str(region[k]) for k in ("contig", "start", "end")) + "\n")
    with _range_relay(sources["reads_bam"], budget) as (url, ranges, errors):
        request = {"url": url, "source": sources["reads_bam"], "bai": str(inputs / "source.bam.bai"),
                   "reference": str(reference), "out": str(inputs), "intervals": intervals,
                   "conditions": conditions, "seed": fresh.get("seed", 42)}
        request_path = inputs / "extraction-request.json"
        pilot.write_json(request_path, request)
        subprocess.run([sys.executable, str(Path(__file__).resolve()), "--extract-bams", str(request_path)],
                       check=True, timeout=max(1, budget.deadline - time.monotonic()))
    if errors:
        pilot.write_json(out / "range-errors.json", errors)
        raise ValueError("bounded BAM range relay reported acquisition errors")
    pilot.write_json(out / "bam-byte-ranges.json", ranges)
    bam_preparation = json.loads((inputs / "bam-preparation.json").read_text())
    expected_m5 = fresh.get("expected_chr20_md5")
    if expected_m5 and expected_m5 != bam_preparation["chr20_md5"]:
        raise ValueError("fresh chr20 differs from pinned HG002 reference sequence")
    input_snapshot = {str(p.relative_to(inputs)): _identity(p) for p in inputs.rglob("*") if p.is_file()}
    for name, entry in acquired.items():
        before = input_snapshot[Path(entry["path"]).name]
        if before["sha256"] != entry["sha256"] or before["bytes"] != entry["bytes"]:
            raise ValueError("acquired source changed during preparation")
    pilot.write_json(out / "input-snapshot-before-analysis.json", input_snapshot)
    base_state = {"reference": str(reference), "sample": bam_preparation["sample"],
                  "regions": str(regions), "truth_vcf": str(inputs / "truth.vcf.gz"),
                  "truth_bed": str(inputs / "truth.bed")}
    local_protocol = dict(protocol, intervals=intervals)
    truth = TruthIndex(base_state["truth_vcf"], base_state["truth_bed"], intervals,
                       reference_path=base_state["reference"])
    records, candidate_recall, costs, leaves = [], {}, {}, {}
    try:
        for condition_index, condition in enumerate(conditions):
            name = condition["name"]
            state = dict(base_state, original_bam=bam_preparation["bams"][name],
                         downsampled_bam=bam_preparation["bams"][name])
            folder = out / name
            folder.mkdir()
            started = time.monotonic()
            loci = pilot.scan_candidates(state, local_protocol, binary, folder)
            scan_seconds = time.monotonic() - started
            if not loci:
                raise ValueError(f"no de novo candidates at {name}; no replacement regions")
            sites = folder / "candidates.vcf"
            pilot.write_sites(loci, state, sites)
            observations, extraction = {}, {}
            # Counterbalance shared OS-cache ordering across the prespecified depths.
            order = ("f0", "f1") if condition_index % 2 else ("f1", "f0")
            for group in order:
                group_protocol = dict(local_protocol, rich_fields=["depths", "alleles"] if group == "f0"
                                      else protocol["rich_fields"])
                started = time.monotonic()
                leaf = pilot.extract_leaf(state, group_protocol, binary, sites, group, folder)
                observations[group] = pilot.read_native(leaf, binary, group_protocol, folder / ("read-" + group))
                leaves[name, group] = (leaf, group_protocol,
                                      hashlib.sha256(pilot.canonical(observations[group]).encode()).hexdigest())
                extraction[group] = {"wall_seconds": time.monotonic() - started,
                                     "retained_bytes": pilot.directory_size(folder / group),
                                     "scope": "materialization, verification and first saved read at identical candidate sites"}
            rich_by_locus = {(r["contig"], r["pos"]): r for r in observations["f1"]}
            count_by_locus = {(r["contig"], r["pos"]): r for r in observations["f0"]}
            if rich_by_locus.keys() != count_by_locus.keys():
                raise ValueError("standalone counts/rich candidate universes differ")
            for key, row in count_by_locus.items():
                if any(value != rich_by_locus[key].get(field) for field, value in row.items()):
                    raise ValueError("standalone counts observations differ from rich projection")
            pilot.audit_counts(observations["f1"], state, name, folder)
            parquet = pilot.parquet_roundtrip(observations["f1"], folder / "parquet")
            feature_started = time.monotonic()
            native_records = pilot.feature_records(observations["f1"], loci, name, truth)
            rich_feature_seconds = time.monotonic() - feature_started
            if native_records != pilot.feature_records(parquet, loci, name, truth):
                raise ValueError("fresh native/Parquet features differ")
            feature_started = time.monotonic()
            for row in native_records:
                projected = counts_features(count_by_locus[row["contig"], row["pos0"] + 1], row["alt"])
                if projected != row["f0"]:
                    raise ValueError("standalone counts F0 differs from frozen feature schema")
                row["specimen"] = "HG001"
                row["sample_id"] = "HG001"
            count_feature_seconds = time.monotonic() - feature_started
            proposed = {(r["contig"], r["pos0"], r["ref"], r["alt"]) for r in native_records}
            generation = {}
            for region in intervals:
                keys = {key for key in truth.positives if key[0] == region["contig"]
                        and region["start"] <= key[1] < region["end"]}
                generation[region["block_id"]] = {"truth_snvs": len(keys), "proposed_truth_snvs": len(keys & proposed),
                    "candidate_generation_recall": len(keys & proposed) / len(keys) if keys else None}
            candidate_recall[name] = generation
            started = time.monotonic()
            caller = pilot.caller_scores(state, sites, folder)
            caller_seconds = time.monotonic() - started
            caller_labels = [truth.label(*key)[0] for key in caller]
            for row in native_records:
                row["caller_score"] = caller.get((row["contig"], row["pos0"], row["ref"], row["alt"]))
            records.extend(native_records)
            pilot.write_jsonl(folder / "records.jsonl", native_records)
            costs[name] = {"candidate_scan_seconds": scan_seconds, "extraction_order": order,
                           "matched_candidate_extraction": extraction, "caller_seconds": caller_seconds,
                           "feature_derivation_seconds": {"f0": count_feature_seconds, "f1": rich_feature_seconds},
                           "feature_derivation_scope": "F1 includes candidate assembly/label lookup; F0 includes exact parity check; not a matched microbenchmark",
                           "caller_unthresholded_snv_calls": {"true": caller_labels.count(1), "false": caller_labels.count(0),
                                                               "unknown": caller_labels.count(None)},
                           "candidate_rows": len(native_records), "unknown_labels": sum(r["label"] is None for r in native_records),
                           "input_bam_bytes": Path(state["original_bam"]).stat().st_size}
            budget.check_scratch()
        pilot.write_json(out / "label-domain.json", truth.domain_report())
    finally:
        truth.close()
    # Verify saved replay on every separately extracted condition/feature path.
    hidden = out / "inputs.hidden-for-replay"
    inputs.rename(hidden)
    try:
        for (name, group), (leaf, group_protocol, original_sha) in leaves.items():
            replay = pilot.read_native(leaf, binary, group_protocol, out / name / ("replay-" + group))
            if hashlib.sha256(pilot.canonical(replay).encode()).hexdigest() != original_sha:
                raise ValueError("fresh saved-only replay differs")
    finally:
        hidden.rename(inputs)
    retained = sum(p.stat().st_size for p in out.rglob("*") if p.is_file() and inputs not in p.parents)
    if retained > protocol["limits"]["retained_evidence_bytes"]:
        raise ValueError("fresh retained evidence byte budget exceeded")
    recheck = verify_input_snapshot(input_snapshot)
    pilot.write_json(out / "input-recheck.json", recheck)
    identities = {"sources": acquired, "remote_bam": sources["reads_bam"],
                  "remote_bam_full_md5_verified": False,
                  "remote_bam_identity_scope": "official URL/size/Last-Modified plus received byte-range hashes; full BAM not acquired",
                  "inputs": input_snapshot, "post_analysis_input_recheck": recheck,
                  "bam_preparation": bam_preparation,
                  "cross_specimen_chr20_sequence_match": True if expected_m5 else "not separately verified"}
    metrics = {"wall_seconds": time.monotonic() - budget.started, "download_bytes": budget.download_bytes,
               "retained_bytes": retained, "scratch_high_water_bytes": budget.scratch_high_water,
               "condition_costs": costs, "benchmark_coverage": coverage,
               "correctness": {"independent_counts": True, "native_parquet_features": True,
                               "standalone_counts_f0_equal": True, "saved_only_replay": True},
               "cost_limitations": "Single counterbalanced extraction per arm/depth, shared acquisition and scan; no selective fallback or user-review workflow measured"}
    result = {"records": records, "candidate_recall": candidate_recall, "metrics": metrics,
              "identities": identities, "state": base_state | {"condition_bams": bam_preparation["bams"]}}
    pilot.write_jsonl(out / "records.jsonl", records)
    pilot.write_json(out / "candidate-recall.json", candidate_recall)
    pilot.write_json(out / "fresh-preparation.json", {key: value for key, value in result.items() if key != "records"})
    return result


if __name__ == "__main__":
    if len(sys.argv) == 3 and sys.argv[1] == "--extract-bams":
        if os.environ.get("GITHUB_ACTIONS") != "true":
            raise RuntimeError("native extraction worker requires hosted workflow")
        _extract_bams_worker(sys.argv[2])
    else:
        raise SystemExit("Import prepare_fresh; standalone only supports the hosted extraction worker")
