"""Metadata-only admission and native-owned coordinate normalization."""
from __future__ import annotations

from dataclasses import dataclass
from itertools import groupby
from pathlib import Path

from contract import load_members, normalized_candidates, validate_comparison
from run_manifest import InputLedger, canonical, digest_json, load_json, write_json

FIELDS = ["depths", "alleles"]
MIB = 1 << 20
MAX_MEMBERS = 128
MAX_CANDIDATES = 100000
MAX_CHILDREN = 4096
MAX_FILES = 16384
MAX_RETAINED_BYTES = 256 * MIB
MAX_OUTPUT_ROWS = 1000000


@dataclass
class Prepared:
    plan: dict
    members: list
    candidates: list
    ledger: InputLedger
    datasets: list
    coverage_paths: list
    directory: Path


def windows(candidates, size):
    """Window by unique locus; never divide the ALT rows for one locus."""
    window = []
    count = 0
    for _, group in groupby(candidates, key=lambda c: (c.contig_id, c.position0)):
        if count == size:
            yield window
            window, count = [], 0
        window.extend(group)
        count += 1
    if window:
        yield window


def write_sites(path, candidates, contigs):
    """Serialize already native-normalized sites; this is not an input parser."""
    with Path(path).open("x", encoding="ascii") as stream:
        stream.write("##fileformat=VCFv4.2\n")
        for contig in contigs:
            stream.write(f'##contig=<ID={contig["name"]},length={contig["length"]}>\n')
        stream.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for _, group in groupby(candidates, key=lambda c: (c.contig_id, c.position0)):
            rows = list(group)  # at most three normalized SNV ALTs
            first = rows[0]
            stream.write(f"{first.contig}\t{first.position0 + 1}\t.\t{first.ref}\t"
                         f"{','.join(c.alt for c in rows)}\t.\t.\t.\n")


def covered_keys(coverage):
    return {(site["contig"], site["position"]) for site in coverage["covered"]["sites"]}


def admit(parent_mb, metadata_bytes, candidate_rows, thresholds, member_count,
          window_size, unique_loci, inventory_files):
    if not isinstance(parent_mb, int) or parent_mb <= 0:
        raise ValueError("parent memory budget must be positive")
    if not isinstance(window_size, int) or window_size <= 0 or window_size > MAX_CANDIDATES:
        raise ValueError("window size must be in 1..100000")
    if len(thresholds) > 8 or not thresholds:
        raise ValueError("one to eight depth thresholds are supported")
    if candidate_rows > MAX_CANDIDATES or member_count > MAX_MEMBERS or inventory_files > MAX_FILES:
        raise ValueError("candidate/member/dependency admission limit exceeded")
    child_bound = member_count * ((unique_loci + window_size - 1) // window_size)
    if child_bound > MAX_CHILDREN:
        raise ValueError("window would exceed the 4096 retained child-receipt limit; use a larger window")
    output_rows = candidate_rows * len(thresholds)
    if output_rows > MAX_OUTPUT_ROWS:
        raise ValueError("summary output cardinality limit exceeded")
    # Deliberately conservative admitted-state model, not an OS allocation cap.
    # Covers JSON parse/encoding coexistence, Python object amplification, one
    # coverage record, one projected batch, bounded ledger and output inventory.
    parent_bytes = (64 * MIB + metadata_bytes * 32 + candidate_rows * 768
                    + inventory_files * 2048 + child_bound * 2048
                    + min(unique_loci, window_size) * 3 * len(thresholds) * 2048)
    if parent_bytes > parent_mb * MIB:
        raise ValueError(f"parent state resource refusal: modeled {parent_bytes} bytes > {parent_mb} MiB")
    return {"model": "cohort-python-admission-v1", "modeled_parent_bytes": parent_bytes,
            "parent_memory_budget_bytes": parent_mb * MIB,
            "native_child_invocations_upper_bound": child_bound, "output_rows": output_rows,
            "max_members": MAX_MEMBERS, "max_candidates": MAX_CANDIDATES,
            "max_child_receipts": MAX_CHILDREN, "max_consumed_files": MAX_FILES,
            "max_retained_bytes": MAX_RETAINED_BYTES, "window_loci": window_size,
            "assurance": "bounded admitted state model; cooperative native budget; not an OS RSS cap"}


def prepare(args, directory, binary):
    from rosalind import open_dataset

    directory = Path(directory)
    metadata_dir = directory / "metadata"
    metadata_dir.mkdir()
    coverage_dir = directory / "coverage"
    coverage_dir.mkdir()
    ledger = InputLedger(MAX_FILES)
    members_path, sites_path = Path(args.members).absolute(), Path(args.sites).absolute()
    # Refuse before parsing input/catalog or launching a native selection parser.
    max_input_bytes = min(4 * MIB, max(0, (args.parent_memory_budget_mb * MIB - 64 * MIB) // 32))
    if max_input_bytes <= 0:
        raise ValueError("parent budget is below the 64 MiB process/stream reserve")
    for path in (members_path, sites_path):
        if path.stat().st_size > max_input_bytes:
            raise ValueError(f"input exceeds admitted {max_input_bytes} byte limit: {path}")
        ledger.pin(path)
    # A compressed file's byte size cannot bound the native parser's decoded
    # selection. Keep this prototype to byte-bounded plain VCF, while letting
    # the native parser exclusively validate its coordinates and alleles.
    with sites_path.open("rb") as stream:
        if not stream.read(17).startswith(b"##fileformat=VCFv"):
            raise ValueError("research selection requires uncompressed VCF; compressed VCF/BCF is outside this envelope")
    members = load_members(members_path, max_manifest_bytes=max_input_bytes, max_members=MAX_MEMBERS)
    datasets, descriptors, coverage_paths, entries = [], [], [], []
    metadata_bytes = members_path.stat().st_size + sites_path.stat().st_size
    total_inventory = 2
    candidates, query, contigs = None, None, None
    thresholds = sorted(set(args.min_callable_depth or [10]))
    if any(type(t) is not int or t <= 0 or t > (1 << 64) - 1 for t in thresholds):
        raise ValueError("depth thresholds must be positive uint64 integers")
    if len(thresholds) > 8:
        raise ValueError("at most eight distinct thresholds are supported")
    with (directory / "coverage.jsonl").open("x") as coverage_report:
        for index, member in enumerate(members):
            descriptor_path = member.dataset.parent / "dataset.descriptor.json"
            for path in (member.dataset, descriptor_path):
                size = path.stat().st_size
                if size > max_input_bytes:
                    raise ValueError("dataset metadata exceeds admitted per-file envelope")
                metadata_bytes += size
                ledger.pin(path)
            # Check cumulative metadata before launching a parser or retaining it.
            admit(args.parent_memory_budget_mb, metadata_bytes, 0, thresholds, len(members),
                  args.window_size, 0, total_inventory)
            dataset = open_dataset(member.dataset, binary=binary,
                                   max_dataset_metadata_bytes=max_input_bytes)
            descriptor = dataset.inspect(memory_budget_mb=args.native_memory_budget_mb)
            descriptors.append((member, descriptor))
            datasets.append(dataset)
            total_inventory += 2 + 2 * len(descriptor["partitions"])
            coverage = dataset.coverage(sites=sites_path, fields=FIELDS,
                                        memory_budget_mb=args.native_memory_budget_mb,
                                        max_output_bytes=max_input_bytes)
            if coverage.get("schema") != "rosalind-dataset-coverage-v1":
                raise ValueError("unsupported native coverage adapter; build the changed revision")
            for identity in (coverage["parent"], coverage["descriptor"]):
                ledger.pin(identity["path"], identity["blake3"], "blake3")
            if candidates is None:
                query, contigs = coverage["query"], descriptor["contigs"]
                candidates = normalized_candidates(query["selection"], contigs, max_candidates=MAX_CANDIDATES)
            elif query != coverage["query"]:
                raise ValueError("native normalized requests differ between members")
            path = coverage_dir / f"{index:04d}.json"
            write_json(path, coverage)
            # Execution rereads this split on every member/window pass. Bind it
            # now so final output hashing cannot silently accept a changed mask.
            ledger.pin(path)
            total_inventory += 1
            coverage_paths.append(path)
            # Serialized split selections are disk-backed, with bounded parsing
            # one member at a time during execution, not a retained M x L mask.
            metadata_bytes += path.stat().st_size
            descriptor_copy = metadata_dir / f"{index:04d}.json"
            write_json(descriptor_copy, descriptor)
            ledger.pin(descriptor_copy)
            total_inventory += 1
            entry = {"specimen_id": member.specimen_id, "dataset": str(member.dataset),
                     "metadata": dict(member.metadata),
                     "parent": coverage["parent"], "descriptor": coverage["descriptor"],
                     "sample_scope": descriptor["sample_scope"], "sources": descriptor["sources"],
                     "physical_fields": descriptor["fields"], "required_fields": 3,
                     "requested_loci": coverage["requested_loci"],
                     "covered_loci": coverage["covered_loci"], "missing_loci": coverage["missing_loci"]}
            entries.append(entry)
            coverage_report.write(canonical({"specimen_id": member.specimen_id,
                                            "requested_loci": coverage["requested_loci"],
                                            "covered_loci": coverage["covered_loci"],
                                            "missing_loci": coverage["missing_loci"],
                                            "missing": coverage["missing"]}) + "\n")
    comparison = validate_comparison(descriptors)
    del descriptors
    unique_loci = sum(1 for _ in groupby(candidates, key=lambda c: (c.contig_id, c.position0)))
    bounds = admit(args.parent_memory_budget_mb, metadata_bytes, len(candidates), thresholds,
                   len(members), args.window_size, unique_loci, total_inventory)
    if args.missingness == "strict" and any(entry["missing_loci"] for entry in entries):
        raise ValueError("strict coverage refusal: one or more member/candidate loci are absent")
    native_plans = []
    for index, (entry, dataset, path) in enumerate(zip(entries, datasets, coverage_paths)):
        coverage = load_json(path, max_input_bytes)
        if not entry["covered_loci"]:
            native_plans.append(None)  # never turn empty coverage into all positions
            continue
        covered = normalized_candidates(coverage["covered"], contigs, max_candidates=MAX_CANDIDATES)
        sites = metadata_dir / f"{index:04d}-covered.vcf"
        write_sites(sites, covered, contigs)
        native_plans.append(dataset.plan(sites=sites, fields=FIELDS,
                                        memory_budget_mb=args.native_memory_budget_mb))
    plan = {"schema": "rosalind-cohort-research-plan-v1", "verification_level": "metadata_only",
            "partition_bodies_verified": False, "reference_validation": "deferred_to_observed_rows",
            "missingness": args.missingness, "thresholds": thresholds, "members": entries,
            "comparison": comparison, "query": query, "query_sha256": digest_json(query),
            "candidate_alt_rows": len(candidates), "candidate_loci": unique_loci,
            "contigs": contigs, "native_plans": native_plans, "limits": bounds,
            "metadata_bytes_charged": metadata_bytes, "max_metadata_file_bytes": max_input_bytes,
            "native_memory_budget_mb": args.native_memory_budget_mb,
            "native_metadata_invocations": len(members) * 2 + sum(p is not None for p in native_plans)}
    ledger.recheck()
    write_json(directory / "plan.json", plan)
    return Prepared(plan, members, candidates, ledger, datasets, coverage_paths, directory)
