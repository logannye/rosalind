"""Same-candidate saved-feature cost measurements with a columnar baseline.

No manual-review-time or net end-user-cost inference is made from these timings.
Extraction/caller measurements belong to the independently timed fresh adapter.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import statistics
import time


def _summarize(values):
    return {"rounds": values, "median_wall_seconds": statistics.median(v["wall_seconds"] for v in values),
            "median_cpu_seconds": statistics.median(v["cpu_seconds"] for v in values)}


def measure_costs(records, bundle, manifest, output, protocol):
    import pyarrow as pa
    import pyarrow.parquet as pq
    import thresholds
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    rows = sorted([r for r in records if r["condition"] == protocol["primary_condition"]],
                  key=lambda r: r["candidate_id"])
    if not rows:
        raise ValueError("Cannot benchmark an empty primary candidate pool")
    metadata = ("candidate_id", "contig", "pos0", "ref", "alt", "condition", "partition", "block_id")
    rounds = protocol["timing"]["rounds"]
    if rounds < 3:
        raise ValueError("Counterbalanced timing requires at least three rounds")
    paths, storage, serialization = {}, {}, {}
    for label, group in (("counts", "f0"), ("rich", "f1")):
        start, cpu = time.perf_counter(), time.process_time()
        names = bundle["feature_names"][group]
        columns = {key: [r[key] for r in rows] for key in metadata}
        table = pa.table({**{key: pa.array(vals) for key, vals in columns.items()},
                          **{name: pa.array([r[group][name] for r in rows], type=pa.float64()) for name in names}})
        path = output / f"features-{group}.parquet"
        pq.write_table(table, path, compression="zstd", row_group_size=1024)
        serialization[label] = {"wall_seconds": time.perf_counter()-start, "cpu_seconds": time.process_time()-cpu}
        paths[label] = path
        storage[label] = {"bytes": path.stat().st_size, "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                          "feature_columns": len(names), "candidate_rows_including_unknown": len(rows)}

    expected = {}
    structural = {}
    for label in ("counts", "rich"):
        expected[label] = thresholds.predict_policy(rows, bundle, manifest, label)
        structural[label] = thresholds._structural(manifest["policies"][manifest["selected"][label]])

    def measure(label, parse):
        start, cpu = time.perf_counter(), time.process_time()
        group = "f0" if label == "counts" else "f1"
        source = rows
        if parse:
            source = []
            for batch in pq.ParquetFile(paths[label]).iter_batches(batch_size=1024, use_threads=False):
                for record in batch.to_pylist():
                    features = {name: record.pop(name) for name in bundle["feature_names"][group]}
                    record[group] = features
                    if group == "f1":
                        record["f0"] = {name: features[name] for name in bundle["feature_names"]["f0"]}
                    source.append(record)
            if [r["candidate_id"] for r in source] != [r["candidate_id"] for r in rows]:
                raise ValueError("Columnar timing candidate identity/order changed")
        decisions = thresholds.predict_policy(source, bundle, manifest, label)
        wall, cpu_seconds = time.perf_counter()-start, time.process_time()-cpu
        if decisions != expected[label]:
            raise ValueError("Saved-columnar decisions differ from identical in-memory features")
        return {"wall_seconds": wall, "cpu_seconds": cpu_seconds}

    measurements = {mode: {label: [] for label in ("counts", "rich")} for mode in ("warm_inference", "parquet_parse_and_inference")}
    for mode in measurements:
        parse = mode == "parquet_parse_and_inference"
        for label in ("counts", "rich"):
            measure(label, parse)
        for index in range(rounds):
            for label in (("counts", "rich") if index % 2 == 0 else ("rich", "counts")):
                measurements[mode][label].append(measure(label, parse))
    summarized = {mode: {label: _summarize(values) for label, values in arms.items()}
                  for mode, arms in measurements.items()}
    report = {
        "schema": "rosalind-threshold-transfer-costs-v1", "primary_condition": protocol["primary_condition"],
        "candidate_rows_including_unknown": len(rows), "rounds": rounds,
        "order": "counts/rich alternating with rich/counts, after one warmup per arm",
        "columnar_format": "Parquet zstd, row_group_size1024, one thread, identical candidate order",
        "storage": storage, "serialization": serialization, "measurements": summarized,
        "selected_structural_retain_all": structural,
        "structural_policy_cost_scope": "A structural retain-all policy needs neither model scoring nor stored features. Its serialized-table costs are optional representation overhead, not necessary retention-policy cost.",
        "repeated_saved_queries": {str(k): {label: k*summarized["parquet_parse_and_inference"][label]["median_wall_seconds"]
                                              for label in ("counts", "rich")}
                                   for k in protocol["timing"]["repeated_queries"]},
        "repeated_saved_queries_scope": "Arithmetic projection of measured median repeated-query costs, not observed end-to-end crossover",
        "native_parquet_prediction_parity": True,
        "manual_review_seconds": None, "net_end_user_benefit_established": False,
        "scope": "Standalone feature storage, serialization and repeated same-candidate matrix/model/policy execution. Acquisition/extraction/caller costs reported separately; no native-vs-Parquet speed claim.",
    }
    with (output / "report.json").open("x") as handle:
        json.dump(report, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
    return report
