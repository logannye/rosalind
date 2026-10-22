#!/usr/bin/env python3
"""Bounded, diagnostic B/C/D pilot. Run only with an explicitly built binary.

No downloads, implicit builds, real participants, cold-cache claims or timing
claims are embedded. Each measured worker is a separate process, including setup.
D independently counts the same BAM reads; it never consumes analyzer output.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import statistics
import subprocess
import sys
import threading
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "python"))


def sha256(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            value.update(block)
    return value.hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, sort_keys=True, indent=2) + "\n")


def size(path):
    return sum(p.stat().st_size for p in Path(path).rglob("*") if p.is_file())


def write_sites(path, candidates, length):
    Path(path).write_text("##fileformat=VCFv4.2\n"
        f"##contig=<ID=synthetic,length={length}>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
        "".join(f"synthetic\t{p}\t.\tA\t{alt}\t.\t.\t.\n" for p, alt in candidates))


def prepare(root, config):
    import pysam
    if pysam.__version__ != "0.23.3":
        raise ValueError("pilot requires pinned pysam==0.23.3")
    limits = config["limits"]
    if not 1 <= config["members"] <= limits["max_members"] <= 32:
        raise ValueError("member bound")
    if not 1 <= config["loci"] <= limits["max_loci"] <= 8192:
        raise ValueError("candidate bound")
    if not 1 <= len(config["requests"]) <= limits["max_requests"] <= 20:
        raise ValueError("request bound")
    positions = [1 + index * config["spacing"] for index in range(config["loci"])]
    if positions[-1] > config["reference_length"]:
        raise ValueError("reference bound")
    raw = root / "raw"
    raw.mkdir()
    reference = raw / "reference.fa"
    reference.write_text(">synthetic\n" + "A" * config["reference_length"] + "\n")
    pysam.faidx(str(reference))
    members = []
    for member_index in range(config["members"]):
        specimen = f"synthetic-{member_index:03d}"
        stored = [p for index, p in enumerate(positions)
                  if member_index == 0 or (index + member_index) % config["missing_modulus"]]
        bed = raw / (specimen + ".bed")
        bed.write_text("".join(f"synthetic\t{p-1}\t{p}\n" for p in stored))
        bam = raw / (specimen + ".bam")
        header = {"HD": {"VN": "1.6", "SO": "coordinate"},
                  "SQ": [{"SN": "synthetic", "LN": config["reference_length"]}],
                  "RG": [{"ID": specimen, "SM": specimen}]}
        with pysam.AlignmentFile(str(bam), "wb", header=header) as stream:
            for index, pos in enumerate(positions):
                # Distinct deterministic source bytes/scopes; zero-depth observations
                # remain stored. Extra excluded reads exercise quality/flag semantics.
                depth = (index * 3 + member_index * 7) % 24
                for ordinal in range(depth + 3):
                    read = pysam.AlignedSegment(stream.header)
                    read.query_name = f"{specimen}-{index}-{ordinal}"
                    read.query_sequence = "ACGT"[(ordinal + index + member_index) % 4]
                    read.query_qualities = [10 if ordinal == depth + 1 else 30]
                    read.flag = 1024 if ordinal == depth else 0
                    read.reference_id = 0
                    read.reference_start = pos - 1
                    read.mapping_quality = 10 if ordinal == depth + 2 else 60
                    read.cigarstring = "1M"
                    read.set_tag("RG", specimen)
                    stream.write(read)
        pysam.index(str(bam))
        members.append({"specimen_id": specimen, "bam": str(bam), "regions": str(bed),
                        "stored": stored})
    requests = []
    for request in config["requests"]:
        indices = sorted({index for start, end in request["indices"] for index in range(start, end)})
        if not indices or indices[-1] >= len(positions) or indices[0] < 0:
            raise ValueError("invalid frozen request indices")
        candidates = [(positions[index], request["alt"]) for index in indices]
        path = raw / (request["id"] + ".vcf")
        write_sites(path, candidates, config["reference_length"])
        requests.append({**request, "sites": str(path), "candidates": candidates})
    state = {"version": 1, "root": str(root), "reference": str(reference), "members": members,
             "requests": requests, "config": config, "synthetic": True}
    write_json(root / "state.json", state)
    return state


def checked(command, log):
    with Path(log).open("ab") as output:
        subprocess.run(list(map(str, command)), check=True, stdout=output, stderr=output)


def source_identities(state):
    paths = [Path(state["reference"]), Path(state["reference"] + ".fai")]
    for member in state["members"]:
        paths.extend([Path(member["bam"]), Path(member["bam"] + ".bai"), Path(member["regions"])])
    return {str(path): sha256(path) for path in paths}


def setup_native(state, binary, output):
    from rosalind import open_dataset
    output.mkdir()
    members = []
    parts = {}
    for member in state["members"]:
        artifact = output / (member["specimen_id"] + ".arrow")
        checked([binary, "analyze", "evidence", "--reference", state["reference"],
                 "--alignments", member["bam"], "--sample", member["specimen_id"],
                 "--regions", member["regions"], "--fields", "depths,alleles", "--format", "arrow-ipc",
                 "--memory-budget-mb", state["config"]["native_memory_budget_mb"],
                 "--cache-dir", output / (member["specimen_id"] + "-cache"), "--output", artifact],
                output / "commands.log")
        receipt = json.loads(Path(str(artifact) + ".manifest.json").read_text())
        dataset = Path(receipt["measurements"]["execution.evidence_dataset_manifest"])
        handle = open_dataset(dataset, binary=binary)
        handle.verify(memory_budget_mb=state["config"]["native_memory_budget_mb"])
        count = len(list(dataset.parent.glob("*/evidence.arrow")))
        if count < 2:
            raise AssertionError("pilot must exercise more than one physical partition per specimen")
        parts[member["specimen_id"]] = count
        members.append({"specimen_id": member["specimen_id"], "dataset": str(dataset)})
    write_json(output / "members.json", {"version": 1, "members": members})
    write_json(output / "setup.json", {"physical_partitions": parts, "sources": source_identities(state),
                                     "scope": "native extraction, materialization and initial full dataset verification"})


def independent_rows(member):
    """Independent sequential 1M read-count oracle; no native rows or pileup defaults.

    Each source is generated above and coordinate sorted. Only one locus is kept
    at a time. The explicit saved mask is applied, including stored zero rows.
    """
    import pysam
    stored = iter(member["stored"])
    wanted = next(stored, None)
    counts = [0, 0, 0, 0]
    with pysam.AlignmentFile(member["bam"], "rb") as alignment:
        for read in alignment:
            if read.cigarstring != "1M" or len(read.query_sequence) != 1:
                raise ValueError("synthetic independent oracle only supports 1M")
            position = read.reference_start + 1
            while wanted is not None and wanted < position:
                yield {"pos": wanted, "callable_depth": sum(counts),
                       **dict(zip("acgt", counts))}
                counts = [0, 0, 0, 0]
                wanted = next(stored, None)
            if wanted != position:
                continue
            if read.flag & (0x100 | 0x800 | 0x200 | 0x400) or read.mapping_quality == 255 or read.mapping_quality < 20:
                continue
            if read.query_qualities[0] < 20 or read.query_sequence not in "ACGT":
                continue
            counts["ACGT".index(read.query_sequence)] += 1
    while wanted is not None:
        yield {"pos": wanted, "callable_depth": sum(counts), **dict(zip("acgt", counts))}
        counts = [0, 0, 0, 0]
        wanted = next(stored, None)


def setup_columnar(state, output):
    import pyarrow as pa
    import pyarrow.parquet as pq
    output.mkdir()
    sources = source_identities(state)
    schema = pa.schema([(name, pa.uint64()) for name in ["pos", "callable_depth", "a", "c", "g", "t"]])
    inventory = []
    for member in state["members"]:
        partition, buffer = None, []

        def flush():
            if not buffer:
                return
            path = output / f"member={member['specimen_id']}" / f"ownership={partition}" / "counts.parquet"
            path.parent.mkdir(parents=True, exist_ok=True)
            pq.write_table(pa.Table.from_pylist(buffer, schema=schema), path,
                           row_group_size=1024, compression="zstd", write_statistics=True)
            # Full independent readback validates materialized values and uint64 type.
            actual = pq.ParquetFile(path).read()
            if actual.to_pylist() != buffer or actual.schema != schema:
                raise AssertionError("columnar setup readback mismatch")
            inventory.append({"path": str(path), "sha256": sha256(path), "member": member["specimen_id"],
                              "ownership": partition, "min_pos": buffer[0]["pos"], "max_pos": buffer[-1]["pos"]})
            buffer.clear()

        for row in independent_rows(member):
            owner = (row["pos"] - 1) // 16384
            if partition != owner:
                flush()
                partition = owner
            buffer.append(row)
        flush()
    if sources != source_identities(state):
        raise AssertionError("source changed during independent extraction")
    write_json(output / "inventory.json", {"version": 1, "sources": sources, "files": inventory,
        "scope": "independent 1M BAM counting, materialization, uint64 readback and SHA256; no native receipts"})


def counter_rows(request, members):
    return {(p, alt, threshold): [0, 0, 0, 0, 0, 0, 0]
            for p, alt in request["candidates"] for threshold in request["thresholds"]}


def add(counters, request, row):
    position = row["pos"]
    alt = request["alt"]
    depth, allele = row["callable_depth"], row[alt.lower()]
    if (type(position) is not int or position <= 0
            or type(depth) is not int or type(allele) is not int
            or not 0 <= allele <= depth <= (1 << 64) - 1):
        raise ValueError("invalid exact count")
    for threshold in request["thresholds"]:
        key = (position, alt, threshold)
        if key not in counters:
            raise ValueError("reader returned an unrequested row")
        value = counters[key]
        value[0] += 1
        value[3] += depth
        value[4] += allele
        if depth >= threshold:
            value[1] += 1
            value[2] += int(allele > 0)
            value[5] += depth
            value[6] += allele
        if max(value) > (1 << 64) - 1:
            raise OverflowError("summary outside uint64")


def finish(counters, member_count):
    result = []
    for (pos, alt, threshold), (observed, eligible, supported, depth, allele, ed, ea) in sorted(counters.items()):
        result.append({"contig": "synthetic", "pos": pos, "ref": "A", "alt": alt,
            "min_callable_depth": threshold, "n_requested": member_count, "n_observed": observed,
            "n_depth_eligible": eligible, "n_alt_supported": supported,
            "total_callable_depth": depth, "total_alt_count": allele,
            "eligible_callable_depth": ed, "eligible_alt_count": ea,
            "alt_fraction_numerator": allele if depth else None, "alt_fraction_denominator": depth or None,
            "eligible_alt_fraction_numerator": ea if ed else None, "eligible_alt_fraction_denominator": ed or None,
            "reference_validation": "checked_on_observed_rows" if observed else "unavailable"})
    return result


def saved_request(state, binary, request, output):
    from rosalind import open_dataset
    members = json.loads((Path(state["root"]) / "native/members.json").read_text())["members"]
    counters = counter_rows(request, members)
    launches = 0
    sources = {source["specimen_id"]: source for source in state["members"]}
    if (len(sources) != len(state["members"]) or len(members) != len(sources)
            or {member["specimen_id"] for member in members} != set(sources)):
        raise ValueError("saved baseline member identities differ from the original masks")
    for member in members:
        source = sources[member["specimen_id"]]
        covered = [(pos, alt) for pos, alt in request["candidates"] if pos in source["stored"]]
        if not covered:
            continue
        selection = output / (member["specimen_id"] + ".vcf")
        write_sites(selection, covered, state["config"]["reference_length"])
        dataset = open_dataset(member["dataset"], binary=binary)
        with dataset.batches(sites=selection, fields=["depths", "alleles"],
                             memory_budget_mb=state["config"]["native_memory_budget_mb"],
                             workdir=output / member["specimen_id"]) as stream:
            launches += 1
            for batch in stream:
                for row in batch.to_pylist():
                    add(counters, request, row)
            if stream.result is None:
                raise AssertionError("saved reader did not complete")
    return finish(counters, len(members)), {"native_query_invocations": launches, "native_version_probes": launches,
        "native_process_launches": 2 * launches,
        "integrity_scope": "native verified consumed inputs per invocation; no outer union recheck"}


def columnar_request(state, request, output, integrity):
    import pyarrow.dataset as ds
    root = Path(state["root"]) / "columnar"
    inventory_path = root / "inventory.json"
    before_inventory = sha256(inventory_path) if integrity else None
    inventory = json.loads(inventory_path.read_text())
    counters = counter_rows(request, state["members"])
    selected = sorted({pos for pos, _ in request["candidates"]})
    consumed = []
    rows_read = 0
    for member in state["members"]:
        relevant = [part for part in inventory["files"] if part["member"] == member["specimen_id"]
                    and any(part["min_pos"] <= p <= part["max_pos"] for p in selected)]
        if not relevant:
            continue
        if integrity:
            for part in relevant:
                if sha256(part["path"]) != part["sha256"]:
                    raise ValueError("columnar consumed content changed")
        consumed.extend(relevant)
        # Explicit file pruning plus Parquet row-group statistics and projection.
        dataset = ds.dataset([part["path"] for part in relevant], format="parquet")
        scanner = dataset.scanner(columns=["pos", "callable_depth", "a", "c", "g", "t"],
            filter=ds.field("pos").isin(selected), batch_size=1024,
            batch_readahead=1, fragment_readahead=1, use_threads=False)
        for batch in scanner.to_batches():
            for row in batch.to_pylist():
                rows_read += 1
                add(counters, request, row)
    if integrity:
        for part in consumed:
            if sha256(part["path"]) != part["sha256"]:
                raise ValueError("columnar consumed content changed before publication")
        if sha256(inventory_path) != before_inventory:
            raise ValueError("columnar inventory changed")
    return finish(counters, len(state["members"])), {"native_process_launches": 0, "rows_read": rows_read,
        "consumed_files": len(consumed), "bytes_hashed": 2 * (sum(Path(p["path"]).stat().st_size for p in consumed) + inventory_path.stat().st_size) if integrity else 0,
        "integrity_scope": "consumed SHA256 before/after and inventory identity" if integrity else "performance-only; no content hashing"}


def worker(args):
    state = json.loads(args.state.read_text())
    if args.worker == "setup-native":
        setup_native(state, args.binary, args.output)
        return
    if args.worker == "setup-columnar":
        setup_columnar(state, args.output)
        return
    args.output.mkdir()
    answers, metrics = [], []
    # The harness knows the frozen family; individual reducers receive only the
    # current request. They cannot coalesce or prefetch a later query.
    visible_state = {key: value for key, value in state.items() if key != "requests"}
    for index, request in enumerate(state["requests"][:args.prefix]):
        destination = args.output / f"request-{index:02d}"
        if args.worker == "analyzer":
            command = [sys.executable, str(HERE / "analyze.py"), "--members", str(Path(state["root"]) / "native/members.json"),
                "--sites", request["sites"], "--missingness", "partial", "--binary", str(args.binary),
                "--native-memory-budget-mb", str(state["config"]["native_memory_budget_mb"]),
                "--parent-memory-budget-mb", str(state["config"]["parent_memory_budget_mb"]),
                "--window-size", str(args.window), "--output", str(destination)]
            for threshold in request["thresholds"]:
                command += ["--min-callable-depth", str(threshold)]
            # Every path retains one Python interpreter for this K-prefix.
            # Native subprocesses remain fresh, and only this request is passed.
            import analyze
            with (args.output / "commands.log").open("a") as log:
                log.write(json.dumps({"equivalent_cli": command, "execution": "in-process run_analysis"}) + "\n")
            analyze.run_analysis(analyze.build_parser().parse_args(command[2:]))
            answer = [json.loads(line) for line in (destination / "summary.jsonl").read_text().splitlines()]
            metric = json.loads((destination / "metrics.json").read_text()) if (destination / "metrics.json").exists() else {}
        else:
            destination.mkdir()
            if args.worker == "saved":
                answer, metric = saved_request(visible_state, args.binary, request, destination)
            else:
                answer, metric = columnar_request(visible_state, request, destination, args.integrity)
        answers.append(answer)
        metrics.append(metric)
    write_json(args.output / "answers.json", answers)
    write_json(args.output / "metrics.json", metrics)


def measure(command, destination):
    import psutil
    destination.parent.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    peak, samples = 0, 0
    stopped = threading.Event()
    with destination.with_suffix(".stdout.log").open("wb") as stdout, destination.with_suffix(".stderr.log").open("wb") as stderr:
        process = subprocess.Popen(command, stdout=stdout, stderr=stderr)

        def sample():
            nonlocal peak, samples
            try:
                parent = psutil.Process(process.pid)
            except psutil.Error:
                return
            while not stopped.is_set():
                try:
                    processes = [parent, *parent.children(recursive=True)]
                    rss = 0
                    for item in processes:
                        try:
                            rss += item.memory_info().rss
                        except psutil.Error:
                            pass
                    peak = max(peak, rss)
                    samples += 1
                except psutil.Error:
                    pass
                stopped.wait(0.01)

        monitor = threading.Thread(target=sample, daemon=True)
        monitor.start()
        try:
            code = process.wait()
            # Capture completion immediately, not at the next 10-ms RSS poll.
            elapsed = time.perf_counter() - started
        finally:
            stopped.set()
            monitor.join()
    return {"argv": command, "exit_code": code, "outer_wall_seconds": elapsed,
            "sampled_process_tree_peak_rss_bytes": peak, "rss_samples": samples,
            "rss_sample_interval_ms": 10, "rss_scope": "simultaneous worker and current descendants; sampling may miss short peaks",
            "memory_budget_semantics": "cooperative native/parent admission; no OS RSS cap",
            "output_bytes": size(destination) if destination.exists() else 0}


def run(args):
    import importlib.metadata
    config = json.loads(args.config.read_text())
    if config["repetitions"] < 3:
        raise ValueError("at least three counterbalanced repetitions required")
    args.output.mkdir(parents=True, exist_ok=False)
    root = args.output.resolve()
    report = {"version": 1, "status": "running", "synthetic": True,
        "purpose": "diagnostic G2 pilot; not final application performance evidence",
        "config": config, "config_sha256": sha256(args.config),
        "binary": {"path": str(args.binary), "sha256": sha256(args.binary),
                   "version": subprocess.check_output([str(args.binary), "--version"], text=True).strip()},
        "source_revision": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "source_status": subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT, text=True),
        "environment": {"platform": platform.platform(), "machine": platform.machine(),
                        "processor": platform.processor(), "cpu_count": os.cpu_count(), "python": sys.version,
                        "dependencies": {name: importlib.metadata.version(name) for name in ("pyarrow", "pysam", "psutil", "blake3")}},
        "setup": {}, "trials": [], "failures": [],
        "boundaries": ["Synthetic generation is outside both setup timings; both paths read the identical generated BAM observations.",
                       "Setup measured once per representation; this is a diagnostic estimate, not a repeated setup distribution.",
                       "D content checks are separately labeled and do not implement native receipt semantics.",
                       "No OS cache eviction, no cold-cache or p99 claim. All trial order and failures retained.",
                       "No raw/fresh A comparison or real-read conclusion in this pilot.",
                       "All query paths retain one Python worker per K-prefix; each analyzer request is an in-process API call. Native children remain fresh."]}
    try:
        state = prepare(root, config)
        base = [sys.executable, str(Path(__file__).resolve()), "--state", str(root / "state.json"), "--binary", str(args.binary)]
        for worker_name, name in [("setup-native", "native"), ("setup-columnar", "columnar")]:
            command = [*base, "--worker", worker_name, "--output", str(root / name)]
            result = measure(command, root / name)
            report["setup"][name] = result
            if result["exit_code"]:
                raise RuntimeError(f"{name} setup failed; see retained stderr")
        paths = [("B-saved", "saved", None, False)]
        paths += [(f"C-window-{window}", "analyzer", window, False) for window in config["windows"]]
        paths += [("D-performance-only", "columnar", None, False), ("D-content-checked", "columnar", None, True)]
        reference = {}
        for repetition in range(config["repetitions"]):
            order = paths[repetition % len(paths):] + paths[:repetition % len(paths)]
            for prefix in config["prefixes"]:
                if not 1 <= prefix <= len(state["requests"]):
                    raise ValueError("prefix outside request family")
                for name, worker_name, window, integrity in order:
                    destination = root / "trials" / f"r{repetition}-k{prefix}-{name}"
                    command = [*base, "--worker", worker_name, "--prefix", str(prefix), "--output", str(destination)]
                    if window is not None:
                        command += ["--window", str(window)]
                    if integrity:
                        command.append("--integrity")
                    result = measure(command, destination)
                    result.update({"path": name, "repetition": repetition, "prefix": prefix})
                    report["trials"].append(result)
                    if result["exit_code"]:
                        report["failures"].append({"path": name, "prefix": prefix, "repetition": repetition,
                                                   "reason": "nonzero worker exit; retained stderr"})
                    else:
                        answers = json.loads((destination / "answers.json").read_text())
                        reference.setdefault(prefix, answers)
                        result["exact_agreement"] = answers == reference[prefix]
                        result["answers_sha256"] = sha256(destination / "answers.json")
                        if not result["exact_agreement"]:
                            report["failures"].append({"path": name, "prefix": prefix, "repetition": repetition,
                                                       "reason": "scientific mismatch"})
                    write_json(root / "results.json", report)
                    if size(root) > config["limits"]["max_output_bytes"]:
                        raise RuntimeError("declared retained-output byte bound exceeded")
        medians = []
        for name, worker_name, _window, _integrity in paths:
            for prefix in config["prefixes"]:
                trials = [r for r in report["trials"] if r["path"] == name and r["prefix"] == prefix and r["exit_code"] == 0]
                values = [r["outer_wall_seconds"] for r in trials]
                if not values:
                    continue
                setup = report["setup"]["columnar" if worker_name == "columnar" else "native"]["outer_wall_seconds"]
                medians.append({"path": name, "prefix": prefix, "trials": len(trials),
                    "query_wall_median_seconds": statistics.median(values), "query_wall_range_seconds": [min(values), max(values)],
                    "single_setup_plus_query_median_seconds": setup + statistics.median(values),
                    "sampled_tree_rss_median_bytes": statistics.median(r["sampled_process_tree_peak_rss_bytes"] for r in trials),
                    "sampled_tree_rss_range_bytes": [min(r["sampled_process_tree_peak_rss_bytes"] for r in trials),
                                                      max(r["sampled_process_tree_peak_rss_bytes"] for r in trials)]})
        report["medians"] = medians
        report["storage"] = {"raw_bytes": size(root / "raw"), "native_retained_bytes": size(root / "native"),
                             "columnar_retained_bytes": size(root / "columnar"), "all_trial_bytes": size(root / "trials"),
                             "scope": "logical file bytes; includes retained disposable query outputs, not sampled peak temporary storage"}
        report["status"] = "failed" if report["failures"] else "passed"
        report["decision"] = "Inspect measured totals and memory before choosing any optimization. No optimization selected automatically."
    except BaseException as exc:
        report["status"] = "failed"
        report["failures"].append({"reason": repr(exc)})
        write_json(root / "results.json", report)
        raise
    write_json(root / "results.json", report)
    if report["status"] != "passed":
        raise SystemExit("pilot failed; all trials and failures retained in results.json")
    print(root / "results.json")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", required=True, type=lambda value: Path(value).resolve(strict=True))
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--config", type=Path, default=HERE / "workloads/pilot-v1.json")
    parser.add_argument("--state", type=Path)
    parser.add_argument("--worker", choices=["setup-native", "setup-columnar", "saved", "analyzer", "columnar"])
    parser.add_argument("--prefix", type=int, default=1)
    parser.add_argument("--window", type=int, default=256)
    parser.add_argument("--integrity", action="store_true")
    args = parser.parse_args()
    if args.worker:
        if args.state is None:
            parser.error("internal worker requires --state")
        worker(args)
    else:
        run(args)


if __name__ == "__main__":
    main()
