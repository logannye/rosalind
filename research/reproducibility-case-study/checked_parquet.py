"""Bounded, checked Parquet representation of the fixed legacy evidence stream.

This is a representation comparator, not an independent measurement engine.
The caller verifies the native source receipt and independent scientific oracle.
Coordinates are the native Arrow contract: one-based uint32 positions.  The
inventory digest must be anchored by the outer run/bundle identity for trust;
self-consistent replacement of an inventory and payload is not authentication.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import time


MAX_ROWS = 1000
BATCH_ROWS = 1024
MAX_DATA_BYTES = 16 * 1024 * 1024
MAX_INVENTORY_BYTES = 2 * 1024 * 1024
FORMAT = "rosalind-case-study-checked-parquet-v1"
SCALARS = (
    "prefilter_depth", "aligned_depth", "callable_depth", "a", "c", "g", "t",
    "a_fwd", "a_rev", "c_fwd", "c_rev", "g_fwd", "g_rev", "t_fwd", "t_rev",
    "base_quality_sum", "mapping_quality_sum", "read_position_sum", "read_length_sum",
    "filtered_secondary", "filtered_supplementary", "filtered_qc_fail", "filtered_duplicate",
    "filtered_unavailable_mapq", "filtered_low_mapq", "filtered_unavailable_base_quality",
    "filtered_low_base_quality", "filtered_ambiguous_base",
)
HISTOGRAMS = {"base_quality_histogram": 94, "mapping_quality_histogram": 255}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def legacy_schema():
    """Exact mask-63 schema from src/evidence/encoding.rs (including uint64)."""
    import pyarrow as pa
    fields = [pa.field("contig", pa.string(), nullable=False),
              pa.field("pos", pa.uint32(), nullable=False),
              pa.field("ref", pa.string(), nullable=False),
              pa.field("requested_alts", pa.string(), nullable=False)]
    fields += [pa.field(name, pa.uint64(), nullable=False) for name in SCALARS]
    fields += [pa.field(name, pa.list_(pa.field("item", pa.uint64(), nullable=False), size),
                        nullable=False) for name, size in HISTOGRAMS.items()]
    return pa.schema(fields, metadata={b"rosalind.evidence.schema": b"1"})


def _schema_check(schema):
    if not schema.equals(legacy_schema(), check_metadata=True):
        raise ValueError("expected complete native legacy mask-63 Arrow schema")


def _integer(value, maximum=(1 << 64) - 1):
    if type(value) is not int or not 0 <= value <= maximum:
        raise ValueError("expected exact unsigned integer")
    return value


def _contig_ranks(contig_order):
    if (not isinstance(contig_order, list) or not 1 <= len(contig_order) <= 4096
            or any(not isinstance(name, str) or not 1 <= len(name) <= 4096 for name in contig_order)
            or len(set(contig_order)) != len(contig_order)):
        raise ValueError("invalid reference contig order")
    return {name: rank for rank, name in enumerate(contig_order)}


def _locus(row, ranks):
    reference = row.get("ref")
    if (row.get("contig") not in ranks or not isinstance(reference, str)
            or len(reference) != 1 or reference not in "ACGT"):
        raise ValueError("invalid contig or REF")
    pos = _integer(row.get("pos"), (1 << 32) - 1)
    if pos == 0:
        raise ValueError("native Arrow positions are one-based")
    return ranks[row["contig"]], pos


def _canonical_bytes(value):
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n").encode()


def _full_rows(batches, ranks, expected_rows):
    """Validate/hash bounded batches, retaining only bounded coverage metadata."""
    digest = hashlib.sha256()
    coverage = []
    previous = None
    for batch in batches:
        _schema_check(batch.schema)
        if not 0 <= batch.num_rows <= BATCH_ROWS:
            raise ValueError("Arrow batch exceeds 1024-row bound")
        if len(coverage) + batch.num_rows > expected_rows:
            raise ValueError("evidence exceeds declared row bound")
        for row in batch.to_pylist():
            key = _locus(row, ranks)
            if previous is not None and key <= previous:
                raise ValueError("evidence loci must be sorted and unique")
            previous = key
            alts = row["requested_alts"]
            if (not isinstance(alts, str) or len(alts) > 3
                    or alts != "".join(sorted(set(alts)))
                    or any(base not in "ACGT" or base == row["ref"] for base in alts)):
                raise ValueError("invalid native requested ALTs")
            for name in SCALARS:
                _integer(row[name])
            for name, size in HISTOGRAMS.items():
                values = row[name]
                if not isinstance(values, list) or len(values) != size:
                    raise ValueError("invalid fixed-size histogram")
                for value in values:
                    _integer(value)
            digest.update(_canonical_bytes(row))
            coverage.append({name: row[name] for name in ("contig", "pos", "ref")})
    if len(coverage) != expected_rows:
        raise ValueError("missing evidence rows")
    return digest.hexdigest(), coverage


def _relative_file(root: Path, name: str) -> Path:
    """Reject lexical escapes, noncanonical names and all symlink components."""
    if (not isinstance(name, str) or not name or "\\" in name or ":" in name or "\0" in name
            or PurePosixPath(name).is_absolute()
            or any(part in ("", ".", "..") for part in name.split("/"))):
        raise ValueError("inventory paths must be canonical bundle-relative paths")
    root = Path(root)
    if root.is_symlink():
        raise ValueError("inventory root must not be a symlink")
    path = root
    for part in name.split("/"):
        path = path / part
        if path.is_symlink():
            raise ValueError("inventory paths must not traverse symlinks")
    if not path.is_file():
        raise ValueError("missing inventory payload")
    path.resolve(strict=True).relative_to(root.resolve(strict=True))
    return path


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError("duplicate JSON inventory key")
        result[key] = value
    return result


def _digest(value):
    if not isinstance(value, str) or re.fullmatch(r"[0-9a-f]{64}", value) is None:
        raise ValueError("invalid SHA-256 identity")
    return value


def _load_inventory(directory, expected_inventory_sha256=None):
    root = Path(directory)
    inventory_path = _relative_file(root, "inventory.json")
    size = inventory_path.stat().st_size
    if not 0 < size <= MAX_INVENTORY_BYTES:
        raise ValueError("inventory exceeds bounded metadata envelope")
    raw = inventory_path.read_bytes()
    if len(raw) != size or len(raw) > MAX_INVENTORY_BYTES:
        raise ValueError("inventory changed while reading bounded metadata")
    before = hashlib.sha256(raw).hexdigest()
    if expected_inventory_sha256 is not None and before != _digest(expected_inventory_sha256):
        raise ValueError("inventory differs from anchored digest")
    inventory = json.loads(raw, object_pairs_hook=_unique_object)
    if inventory.get("format") != FORMAT or inventory.get("coordinates") != "one-based":
        raise ValueError("unsupported checked Parquet inventory")
    rows = _integer(inventory.get("rows"), MAX_ROWS)
    if rows == 0:
        raise ValueError("empty case-study evidence")
    ranks = _contig_ranks(inventory.get("contig_order"))
    _digest(inventory.get("full_rows_sha256"))
    files = inventory.get("files")
    if not isinstance(files, list) or len(files) != 1:
        raise ValueError("case study requires exactly one sorted Parquet file")
    part = files[0]
    path = _relative_file(root, part.get("path"))
    _digest(part.get("sha256"))
    if not 0 < _integer(part.get("bytes"), MAX_DATA_BYTES):
        raise ValueError("invalid payload size")
    coverage = inventory.get("coverage")
    if not isinstance(coverage, list) or len(coverage) != rows:
        raise ValueError("stored coverage differs from row count")
    keys = [_locus(row, ranks) for row in coverage]
    if keys != sorted(set(keys)):
        raise ValueError("stored coverage must be sorted and unique")
    if not isinstance(inventory.get("semantics"), dict) or not inventory["semantics"]:
        raise ValueError("missing sample/selection/filter semantics")
    return inventory, path, inventory_path, before, size


def _check_payload(path, part):
    if path.stat().st_size != part["bytes"] or sha256(path) != part["sha256"]:
        raise ValueError("consumed Parquet payload changed")


def _check_after(path, part, inventory_path, inventory_digest):
    _check_payload(path, part)
    if sha256(inventory_path) != inventory_digest:
        raise ValueError("Parquet inventory changed before publication")


def _readback(path, ranks, expected_rows):
    import pyarrow.parquet as pq
    parquet = pq.ParquetFile(path)
    try:
        _schema_check(parquet.schema_arrow)
        for index in range(parquet.metadata.num_row_groups):
            group = parquet.metadata.row_group(index)
            if group.num_rows > BATCH_ROWS:
                raise ValueError("Parquet row group exceeds 1024-row bound")
            if any(group.column(column).compression != "ZSTD" for column in range(group.num_columns)):
                raise ValueError("Parquet payload is not fully zstd compressed")
        return _full_rows(parquet.iter_batches(batch_size=BATCH_ROWS, use_threads=False), ranks, expected_rows)
    finally:
        parquet.close()


def convert_arrow(arrow_path, destination, *, semantics, contig_order, expected_rows=1000):
    """Publish one checked zstd file and inventory after full-schema readback.

    Input must be a verified, sorted native Arrow IPC stream with <=1,000 rows.
    ``semantics`` retains the caller's sample, selection and immutable filters.
    Failure preserves the new directory but never publishes an inventory.
    """
    import pyarrow as pa
    import pyarrow.parquet as pq
    if not 1 <= _integer(expected_rows, MAX_ROWS):
        raise ValueError("case-study row count must be in 1..1000")
    if not isinstance(semantics, dict) or not semantics:
        raise ValueError("sample/selection/filter semantics are required")
    ranks = _contig_ranks(contig_order)
    source = Path(arrow_path)
    if not 0 < source.stat().st_size <= MAX_DATA_BYTES:
        raise ValueError("Arrow input exceeds bounded file envelope")
    source_bytes, source_digest = source.stat().st_size, sha256(source)
    root = Path(destination)
    root.mkdir(parents=True, exist_ok=False)
    payload = root / "evidence.parquet"
    with source.open("rb") as handle, pa.ipc.open_stream(handle) as reader:
        _schema_check(reader.schema)
        with pq.ParquetWriter(payload, reader.schema, compression="zstd", use_dictionary=False,
                              write_statistics=True, use_compliant_nested_type=False) as writer:
            def written_batches():
                for batch in reader:
                    # Admission and row validation occur before this generator
                    # resumes to write and then advances the bounded input stream.
                    yield batch
                    writer.write_batch(batch, row_group_size=BATCH_ROWS)
            expected_digest, coverage = _full_rows(written_batches(), ranks, expected_rows)
    payload_digest = sha256(payload)
    actual_digest, actual_coverage = _readback(payload, ranks, expected_rows)
    if actual_digest != expected_digest or actual_coverage != coverage:
        raise ValueError("full-schema Parquet readback differs from native Arrow")
    if source.stat().st_size != source_bytes or sha256(source) != source_digest:
        raise ValueError("native Arrow source changed during conversion")
    if sha256(payload) != payload_digest:
        raise ValueError("Parquet payload changed during full-schema readback")
    inventory = {"format": FORMAT, "coordinates": "one-based", "rows": expected_rows,
                 "contig_order": contig_order, "coverage": coverage, "semantics": semantics,
                 "full_rows_sha256": expected_digest,
                 "source_arrow": {"sha256": source_digest, "bytes": source_bytes},
                 "files": [{"path": "evidence.parquet", "sha256": payload_digest,
                            "bytes": payload.stat().st_size}],
                 "compression": "zstd", "row_group_rows_max": BATCH_ROWS,
                 "integrity_scope": "complete native schema and values; independent oracle is a separate gate"}
    encoded = _canonical_bytes(inventory)
    if len(encoded) > MAX_INVENTORY_BYTES or payload.stat().st_size > MAX_DATA_BYTES:
        raise ValueError("converted cache exceeds bounded envelope")
    with (root / "inventory.json").open("xb") as handle:
        handle.write(encoded)
    return inventory


def verify(directory, *, expected_inventory_sha256=None):
    """Recheck full schema/values/coverage plus before/after payload identities."""
    inventory, path, inventory_path, before, size = _load_inventory(directory, expected_inventory_sha256)
    part = inventory["files"][0]
    _check_payload(path, part)
    digest, coverage = _readback(path, _contig_ranks(inventory["contig_order"]), inventory["rows"])
    if digest != inventory["full_rows_sha256"] or coverage != inventory["coverage"]:
        raise ValueError("full-schema Parquet content differs from inventory")
    _check_after(path, part, inventory_path, before)
    return {"rows": inventory["rows"], "full_schema_verified": True,
            "inventory_sha256": before, "bytes_hashed": 2 * (part["bytes"] + size)}


def query(directory, requests, *, thresholds=(), expected_inventory_sha256=None):
    """Return mixed-ALT answers only after all integrity and missingness checks.

    Each request is {contig, pos, ref, alt}; output preserves request order.
    A caller may publish returned values only after this function succeeds.
    """
    import pyarrow.dataset as ds
    if not isinstance(requests, list) or not 1 <= len(requests) <= MAX_ROWS:
        raise ValueError("request must contain 1..1000 loci")
    if tuple(thresholds) not in ((), (5, 10, 20)) or any(type(v) is not int for v in thresholds):
        raise ValueError("only the frozen depth screens 5/10/20 are supported")
    check_started = time.perf_counter()
    inventory, path, inventory_path, before, size = _load_inventory(directory, expected_inventory_sha256)
    verification_seconds = time.perf_counter() - check_started
    part = inventory["files"][0]
    ranks = _contig_ranks(inventory["contig_order"])
    selected = {}
    grouped = {}
    alts = set()
    for request in requests:
        _locus(request, ranks)
        alt = request.get("alt")
        if not isinstance(alt, str) or len(alt) != 1 or alt not in "ACGT" or alt == request["ref"]:
            raise ValueError("invalid requested ALT")
        key = request["contig"], request["pos"]
        if key in selected:
            raise ValueError("duplicate requested locus")
        selected[key] = request
        grouped.setdefault(key[0], []).append(key[1])
        alts.add(alt.lower())
    check_started = time.perf_counter()
    _check_payload(path, part)
    dataset = ds.dataset(str(path), format="parquet")
    _schema_check(dataset.schema)
    verification_seconds += time.perf_counter() - check_started
    columns = ["contig", "pos", "ref", "callable_depth"]
    for base in sorted(alts):
        columns.extend([base, f"{base}_fwd", f"{base}_rev"])
    predicate = None
    for contig, positions in grouped.items():
        clause = (ds.field("contig") == contig) & ds.field("pos").isin(positions)
        predicate = clause if predicate is None else predicate | clause
    scanner = dataset.scanner(columns=columns, filter=predicate, batch_size=BATCH_ROWS,
                              batch_readahead=1, fragment_readahead=1, use_threads=False)
    answers = {}
    for batch in scanner.to_batches():
        if batch.num_rows > BATCH_ROWS:
            raise ValueError("scanner exceeded bounded batch size")
        for row in batch.to_pylist():
            key = row["contig"], row["pos"]
            if key not in selected or key in answers:
                raise ValueError("unexpected or duplicate consumed locus")
            request = selected[key]
            if row["ref"] != request["ref"]:
                raise ValueError("requested REF mismatch")
            base = request["alt"].lower()
            answer = {name: request[name] for name in ("contig", "pos", "ref", "alt")}
            answer.update(callable_depth=_integer(row["callable_depth"]),
                          alt_count=_integer(row[base]), alt_fwd=_integer(row[f"{base}_fwd"]),
                          alt_rev=_integer(row[f"{base}_rev"]))
            answer.update({f"depth_ge_{threshold}": answer["callable_depth"] >= threshold
                           for threshold in thresholds})
            answers[key] = answer
    if len(answers) != len(selected):
        raise ValueError("missing requested locus; unavailable evidence is not zero")
    check_started = time.perf_counter()
    _check_after(path, part, inventory_path, before)
    verification_seconds += time.perf_counter() - check_started
    return [answers[(request["contig"], request["pos"])] for request in requests], {
        "native_process_launches": 0, "rows_read": len(answers), "consumed_files": 1,
        "bytes_hashed": 2 * (part["bytes"] + size), "inventory_sha256": before,
        "verification_seconds": verification_seconds,
        "verification_scope": "inventory loading/anchor and metadata validation; pre-scan payload hash and native schema check; post-scan payload and inventory hashes",
        "projected_columns": columns, "batch_size": BATCH_ROWS, "batch_readahead": 1,
        "fragment_readahead": 1, "integrity_scope": "consumed SHA256 before/after and inventory identity"}
