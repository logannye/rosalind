#!/usr/bin/env python3
"""Bounded research-only candidate summaries over one saved leaf per specimen."""
from __future__ import annotations

import argparse
from itertools import groupby
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))

from contract import normalized_candidates
from planner import (FIELDS, MAX_CHILDREN, MAX_FILES, MAX_RETAINED_BYTES, MIB,
                     covered_keys, prepare, windows, write_sites)
from reduce import WindowReducer
from run_manifest import (SCHEMA, ProcessTreeSampler, canonical, dependency_versions,
                          digest_json, file_digest, load_json, output_inventory,
                          publish_new, source_identity, verify_run, write_json)


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--members", required=True, type=Path)
    parser.add_argument("--sites", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--binary", type=Path, default=os.environ.get("ROSALIND_BINARY"))
    parser.add_argument("--min-callable-depth", action="append", type=int)
    parser.add_argument("--missingness", choices=("strict", "partial"), default="strict")
    parser.add_argument("--native-memory-budget-mb", type=int, default=512)
    parser.add_argument("--parent-memory-budget-mb", type=int, default=128)
    parser.add_argument("--window-size", type=int, default=256, help="maximum distinct loci per window")
    parser.add_argument("--plan", action="store_true", help="publish metadata-only plan; no successful run manifest")
    return parser


def _receipt_dependencies(prepared, index, sites, receipt):
    """Accept only actual native dependencies matching the pinned descriptor."""
    member = prepared.members[index]
    entry = prepared.plan["members"][index]
    descriptor = load_json(prepared.directory / "metadata" / f"{index:04d}.json")
    allowed = {str(Path(entry[name]["path"]).resolve()): entry[name]["blake3"]
               for name in ("parent", "descriptor")}
    partition_paths = set()
    for part in descriptor["partitions"]:
        for role in ("receipt", "arrow"):
            identity = part[role]
            path = str((member.dataset.parent / identity["path"]).resolve())
            allowed[path] = identity["blake3"]
            if role == "arrow":
                partition_paths.add(path)
    allowed[str(sites.resolve())] = file_digest(sites, "blake3")[0]
    inputs = receipt.get("inputs")
    if not isinstance(inputs, list) or len(inputs) > MAX_FILES:
        raise ValueError("invalid native child dependency inventory")
    seen, consumed = set(), []
    for identity in inputs:
        path = str(Path(identity["path"]).resolve(strict=True))
        if path not in allowed or identity["blake3"] != allowed[path]:
            raise ValueError("child dependency does not match pinned metadata/query identity")
        prepared.ledger.pin(path, identity["blake3"], "blake3")
        seen.add(path)
        if path in partition_paths:
            consumed.append(path)
    required = {str(Path(entry[name]["path"]).resolve()) for name in ("parent", "descriptor")}
    if not required <= seen or str(sites.resolve()) not in seen:
        raise ValueError("successful child receipt omits expected metadata/selection dependencies")
    return consumed


def _stream(prepared, index, window, child_directory, metrics):
    import pyarrow as pa

    coverage = load_json(prepared.coverage_paths[index], prepared.plan["max_metadata_file_bytes"])
    keys = covered_keys(coverage)
    covered = [candidate for candidate in window if (candidate.contig_id, candidate.position0) in keys]
    if not covered:
        return None  # an omitted sites flag would erroneously query the full leaf
    child_directory.mkdir()
    sites = child_directory / "sites.vcf"
    write_sites(sites, covered, prepared.plan["contigs"])
    expected_loci = len({(c.contig, c.pos) for c in covered})
    dataset = prepared.datasets[index]
    run = dataset.batches(sites=sites, fields=FIELDS, workdir=child_directory,
                          memory_budget_mb=prepared.plan["native_memory_budget_mb"])
    metrics["native_query_invocations"] += 1
    columns = ("contig", "pos", "ref", "callable_depth", "a", "c", "g", "t")
    def consume(reducer):
        with run:
            for batch in run:
                if batch.num_rows > 1024:
                    raise ValueError("native batch exceeds the declared 1024-row envelope")
                for name in columns[3:]:
                    if batch.schema.field(name).type != pa.uint64():
                        raise ValueError("native counts must be physical uint64 columns")
                arrays = [batch.column(batch.schema.get_field_index(name)) for name in columns]
                for row_index in range(batch.num_rows):
                    row = {name: array[row_index].as_py() for name, array in zip(columns, arrays)}
                    if (row["contig"], row["pos"]) not in expected_keys:
                        raise ValueError("native row outside covered intersection")
                    reducer.observe(row)
        if run.result is None or run.result.returncode != 0:
            raise ValueError("native stream did not complete successfully")
        receipt = load_json(run.result.manifest_path, 2 * MIB)
        prepared.ledger.pin(run.result.manifest_path)
        if receipt.get("params", {}).get("run_status") != "completed" or receipt.get("outputs"):
            raise ValueError("invalid completion/scope for streamed native receipt")
        if int(receipt.get("measurements", {}).get("execution.emitted_loci", -1)) != expected_loci:
            raise ValueError("native completed row count differs from covered selection")
        native_query = json.loads(receipt["params"]["dataset.query"])
        normalized = normalized_candidates(native_query["selection"], prepared.plan["contigs"])
        if (native_query.get("version") != 1 or native_query.get("fields") != 3
                or [c.as_dict() for c in normalized] != [c.as_dict() for c in covered]):
            raise ValueError("native receipt query differs from the consumed window")
        consumed = _receipt_dependencies(prepared, index, sites, receipt)
        for path in consumed:
            metrics["partition_decode_visits"] += 1
            metrics["partition_visit_counts"][path] = metrics["partition_visit_counts"].get(path, 0) + 1
        build = {key: receipt["params"].get(key) for key in
                 ("code_git_sha", "code_dirty", "rustc_version", "target_triple", "producer.version")}
        if metrics["native_build_identity"] is not None and metrics["native_build_identity"] != build:
            raise ValueError("native build identity changed during run")
        metrics["native_build_identity"] = build
        peak = int(receipt.get("measurements", {}).get("peak_rss_bytes", 0))
        metrics["largest_native_reported_peak_rss_bytes"] = max(
            peak, metrics["largest_native_reported_peak_rss_bytes"])
        return expected_loci, run.result.manifest_path
    expected_keys = {(c.contig, c.pos) for c in covered}
    return consume


def run_analysis(args, *, after_child=None):
    from rosalind.features import _resolve_binary

    started = time.perf_counter()
    if args.native_memory_budget_mb <= 0 or args.native_memory_budget_mb > (1 << 40):
        raise ValueError("native memory budget must be a positive bounded MiB integer")
    output = Path(args.output).absolute()
    if output.exists() or output.is_symlink():
        raise FileExistsError(f"create-new output already exists: {output}")
    if not output.parent.is_dir():
        raise ValueError("output parent directory must already exist")
    binary = _resolve_binary(args.binary, False).resolve(strict=True)
    binary_hash, binary_size = file_digest(binary)
    binary_version = subprocess.check_output([str(binary), "--version"], text=True).strip()
    source = source_identity(ROOT)
    staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.staging-", dir=output.parent))
    metrics = {"native_query_invocations": 0, "partition_decode_visits": 0,
               "partition_visit_counts": {}, "largest_native_reported_peak_rss_bytes": 0,
               "native_build_identity": None, "summary_rows": 0}
    try:
        with ProcessTreeSampler() as sampler:
            prepared = prepare(args, staging, binary)
            prepared.ledger.pin(binary, binary_hash)
            planning_finished = time.perf_counter()
            if args.plan:
                # Plan identity and output are explicitly metadata-only, and
                # cannot be mistaken for a completed aggregation by verify.
                prepared.ledger.recheck()
                publish_new(staging, output)
                return {"status": "planned", "output": str(output)}
            receipts = staging / "children"
            receipts.mkdir()
            scientific = {"schema": "rosalind-cohort-scientific-v1",
                          "comparison": prepared.plan["comparison"], "query": prepared.plan["query"],
                          "thresholds": prepared.plan["thresholds"], "missingness": args.missingness,
                          "members": [{"specimen_id": member.specimen_id,
                                       "metadata": member.metadata,
                                       "parent_blake3": entry["parent"]["blake3"],
                                       "descriptor_blake3": entry["descriptor"]["blake3"]}
                                      for member, entry in zip(prepared.members, prepared.plan["members"])],
                          "fraction_policy": "exact numerator/denominator; both null at zero denominator",
                          "duplicate_policy": "native union of identical REF loci and distinct ALT alleles",
                          "read_unit": "read bases; overlapping mates count separately"}
            write_json(staging / "scientific.json", scientific)
            write_sites(staging / "candidates.vcf", prepared.candidates, prepared.plan["contigs"])
            retained_estimate = sum(p.stat().st_size for p in staging.rglob("*") if p.is_file())
            with (staging / "summary.jsonl").open("x") as summary:
                for window_index, window in enumerate(windows(prepared.candidates, args.window_size)):
                    reducer = WindowReducer(window, prepared.plan["thresholds"], len(prepared.members))
                    for index, member in enumerate(prepared.members):
                        reducer.begin_member(member.specimen_id)
                        child = receipts / f"{window_index:06d}-{index:04d}"
                        consume = _stream(prepared, index, window, child, metrics)
                        expected = 0
                        if consume is not None:
                            expected, receipt_path = consume(reducer)
                            retained_estimate += sum(p.stat().st_size for p in child.iterdir() if p.is_file())
                            if retained_estimate > MAX_RETAINED_BYTES or metrics["native_query_invocations"] > MAX_CHILDREN:
                                raise ValueError("retained child evidence envelope exceeded")
                            if after_child is not None:
                                after_child(prepared, index, receipt_path)
                        if reducer.end_member() != expected:
                            raise ValueError("observed rows differ from covered preflight denominator")
                    for row in reducer.finish(missingness=args.missingness):
                        line = canonical(row) + "\n"
                        summary.write(line)
                        retained_estimate += len(line.encode())
                        metrics["summary_rows"] += 1
                        if retained_estimate > MAX_RETAINED_BYTES:
                            raise ValueError("retained scientific output envelope exceeded")
            if metrics["summary_rows"] != prepared.plan["limits"]["output_rows"]:
                raise ValueError("canonical output cardinality mismatch")
            prepared.ledger.recheck()
            current_source = source_identity(ROOT)
            if any(current_source[key] != source[key] for key in ("git_revision", "python_sources_sha256")):
                raise ValueError("analyzer source changed during execution")
            write_json(staging / "input-identities.json", prepared.ledger.as_list())
            metrics["unique_consumed_partitions"] = len(metrics["partition_visit_counts"])
            metrics["repeated_partition_decodes"] = (metrics["partition_decode_visits"]
                                                       - metrics["unique_consumed_partitions"])
            metrics["parent_bytes_hashed_before_finalization"] = prepared.ledger.bytes_hashed
            metrics["final_union_recheck_bytes"] = sum(e["bytes"] for e in prepared.ledger.entries.values())
            metrics["native_metadata_invocations"] = prepared.plan["native_metadata_invocations"]
            metrics["native_version_probes"] = metrics["native_metadata_invocations"] + metrics["native_query_invocations"] + 2
            metrics["native_process_launches"] = (metrics["native_version_probes"] +
                                                   metrics["native_metadata_invocations"] + metrics["native_query_invocations"])
            metrics["planning_wall_seconds"] = planning_finished - started
            metrics.update(sampler.report())
            metrics["rss_scope"] = "simultaneous parent and descendants through reduction and first final rehash; excludes final publication phase; benchmark outer sampler covers full workflow"
            # Outer benchmark owns full wall time including final publication.
            metrics["wall_seconds_before_finalization"] = time.perf_counter() - started
            metrics["timing_scope"] = "phase observations may overlap; outer caller measures publication and process startup"
            write_json(staging / "metrics.json", metrics)
            outputs, retained = output_inventory(staging, MAX_FILES + MAX_CHILDREN * 2, MAX_RETAINED_BYTES)
            manifest = {"schema": SCHEMA, "status": "completed", "analyzer": source,
                        "native_binary": {"path": str(binary), "sha256": binary_hash,
                                          "bytes": binary_size, "version": binary_version,
                                          "build": metrics["native_build_identity"]},
                        "dependencies": dependency_versions(), "scientific_sha256": digest_json(scientific),
                        "outputs": outputs, "retained_bytes_excluding_manifest": retained,
                        "native_receipt_scope": "individual completed stdout queries, no persisted intermediate artifact",
                        "trust_scope": "unsigned content integrity and recorded provenance; not biological correctness",
                        "mutation_scope": "final consumed-file union rehash; no snapshot isolation or adversarial race guarantee"}
            manifest["manifest_sha256"] = digest_json(manifest)
            write_json(staging / "run-manifest.json", manifest)
            # Every earlier child is covered; do not silently scan untouched parts.
            prepared.ledger.recheck()
            publish_new(staging, output)
        return {"status": "completed", "output": str(output), "summary_rows": metrics["summary_rows"],
                "wall_seconds_including_publication": time.perf_counter() - started}
    except BaseException as error:
        if staging.exists():
            (staging / "run-manifest.json").unlink(missing_ok=True)
            failure = {"schema": "rosalind-cohort-failure-v1", "status": "failed",
                       "error_type": type(error).__name__, "message": str(error)[:8192],
                       "elapsed_seconds": time.perf_counter() - started,
                       "successful_output_published": False}
            try:
                write_json(staging / "failure.json", failure)
            except OSError:
                pass
        raise


def replay(argv):
    parser = argparse.ArgumentParser(description="Replay saved research summary and compare canonical science")
    parser.add_argument("run", type=Path)
    parser.add_argument("--members", required=True, type=Path, help="explicit relocated or original member catalog")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--binary", type=Path, default=os.environ.get("ROSALIND_BINARY"))
    parser.add_argument("--window-size", type=int)
    args = parser.parse_args(argv)
    original = verify_run(args.run)
    plan = load_json(args.run / "plan.json")
    options = argparse.Namespace(members=args.members, sites=args.run / "candidates.vcf", output=args.output,
                                 binary=args.binary, min_callable_depth=plan["thresholds"],
                                 missingness=plan["missingness"],
                                 native_memory_budget_mb=plan["native_memory_budget_mb"],
                                 parent_memory_budget_mb=plan["limits"]["parent_memory_budget_bytes"] // MIB,
                                 window_size=args.window_size or plan["limits"]["window_loci"], plan=False)
    result = run_analysis(options)
    repeated = verify_run(args.output)
    equal = (file_digest(args.run / "summary.jsonl")[0] == file_digest(args.output / "summary.jsonl")[0]
             and original["scientific_sha256"] == repeated["scientific_sha256"])
    report = {"schema": "rosalind-cohort-replay-v1", "verdict": "reproduced" if equal else "diverged",
              "original_manifest_sha256": original["manifest_sha256"],
              "replayed_manifest_sha256": repeated["manifest_sha256"],
              "scope": "canonical summary and scientific metadata; timings and execution settings excluded"}
    # Certificate is a separate create-new sibling: completed outputs stay intact.
    write_json(Path(str(args.output) + ".replay.json"), report)
    if not equal:
        raise ValueError("replay diverged; the independently completed rerun is retained")
    return {**result, "replay": report}


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    def cancelled(signum, frame):
        raise KeyboardInterrupt(f"cancelled by signal {signum}")
    signal.signal(signal.SIGTERM, cancelled)
    try:
        if argv and argv[0] == "verify":
            parser = argparse.ArgumentParser(description="Verify retained research-run content integrity")
            parser.add_argument("run", type=Path)
            args = parser.parse_args(argv[1:])
            manifest = verify_run(args.run)
            result = {"status": "verified", "manifest_sha256": manifest["manifest_sha256"],
                      "scope": "research manifest and retained output bytes; not native outer-aggregation verification"}
        elif argv and argv[0] == "replay":
            result = replay(argv[1:])
        else:
            result = run_analysis(build_parser().parse_args(argv))
        print(canonical(result))
        return 0
    except KeyboardInterrupt as error:
        print(f"cancelled: {error}", file=sys.stderr)
        return 130
    except (ValueError, OSError, RuntimeError, OverflowError) as error:
        print(f"refused/failed: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
