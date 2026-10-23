#!/usr/bin/env python3
"""Bounded hosted P0 execution. Never builds or installs native code.

The protocol is committed before acquisition. Completed scientific results require
native stream completion, independent count checks, saved-data replay, and exact
Parquet feature equivalence. The final decision is separate from job completion.
"""
from __future__ import annotations

import argparse
from collections import Counter
from contextlib import contextmanager
import hashlib
import importlib.metadata
import importlib.util
import json
import os
from pathlib import Path
import platform
import signal
import subprocess
import sys
import threading
import time
import traceback

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "python"))


def digest(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def write_jsonl(path, records):
    with Path(path).open("w") as out:
        for row in records:
            out.write(canonical(row) + "\n")


def directory_size(path):
    return sum(p.stat().st_size for p in Path(path).rglob("*") if p.is_file())


def locate(intervals, contig, pos0):
    matches = [r for r in intervals if r["contig"] == contig and r["start"] <= pos0 < r["end"]]
    if len(matches) != 1:
        raise ValueError(f"coordinate has {len(matches)} assigned roles: {contig}:{pos0}")
    return matches[0]


def validate_protocol(protocol):
    intervals = protocol["intervals"]
    if {r["partition"] for r in intervals} != {"train", "calibration", "threshold", "test"}:
        raise ValueError("protocol must define all four evaluation roles")
    if len({r["block_id"] for r in intervals}) != len(intervals):
        raise ValueError("duplicate genomic block ID")
    ordered = sorted(intervals, key=lambda r: (r["contig"], r["start"]))
    for i, region in enumerate(ordered):
        if not 0 <= region["start"] < region["end"]:
            raise ValueError("invalid genomic interval")
        if i and region["contig"] == ordered[i-1]["contig"] and region["start"] < ordered[i-1]["end"]:
            raise ValueError("overlapping evaluation blocks")
    if len([r for r in intervals if r["partition"] == "train"]) != 2:
        raise ValueError("P0 requires exactly two training blocks")


def candidate_alts(row, rule):
    ref = row["ref"]
    if ref not in ("A", "C", "G", "T"):
        return []
    depth = row["callable_depth"]
    return [base for base in "ACGT" if base != ref
            and row[base.lower()] >= rule["minimum_alt_count"]
            and row[base.lower()] * rule["fraction_denominator"] >= depth * rule["fraction_numerator"]]


def advancement_allowed(gate, correctness):
    """Unknown statuses and incomplete external checks never unlock P1."""
    return (gate.get("status") == "proceed_to_p1"
            and gate.get("proceed_to_p1") is True
            and gate.get("test_outcomes_used") is False
            and all(correctness.get(key) is True for key in (
                "independent_counts", "native_parquet_features", "saved_only_replay", "prediction_equivalence")))


def truth_recovery(study, generation, intervals):
    denominators = Counter()
    for region in intervals:
        denominators[region["partition"]] += generation[region["block_id"]]["truth_snvs"]
    return {role:{condition:{arm:{"accepted_true_snvs":metrics["true_positives"],
                 "assessed_truth_snvs":denominators[role],
                 "truth_snv_recovery":metrics["true_positives"]/denominators[role] if denominators[role] else None}
                 for arm,metrics in arms.items()} for condition,arms in conditions.items()}
                 for role,conditions in study["metrics"].items()}


def freeze_files(directories):
    return {str(path):digest(path) for directory in directories
            for path in sorted(Path(directory).rglob("*")) if path.is_file()}


class Telemetry:
    def __init__(self, output, protocol):
        self.output, self.protocol = output, protocol
        self.start = time.monotonic()
        self.phases, self.peak, self.scratch_peak = {}, 0, 0
        self.stop = threading.Event()
        self.resource_error = None

    def __enter__(self):
        import psutil
        root = psutil.Process()
        def sample():
            last_disk = 0
            while not self.stop.wait(0.1):
                try:
                    self.peak = max(self.peak, sum(p.memory_info().rss for p in [root, *root.children(recursive=True)] if p.is_running()))
                    if self.peak > self.protocol["limits"]["process_tree_rss_bytes"]:
                        self.resource_error = "process tree RSS budget exceeded"
                        os.kill(os.getpid(), signal.SIGTERM)
                        return
                    if time.monotonic() - last_disk > 5:
                        self.scratch_peak = max(self.scratch_peak, directory_size(self.output))
                        last_disk = time.monotonic()
                        if self.scratch_peak > self.protocol["limits"]["scratch_bytes"]:
                            self.resource_error = "scratch budget exceeded"
                            os.kill(os.getpid(), signal.SIGTERM)
                except (psutil.Error, FileNotFoundError):
                    pass
        self.thread = threading.Thread(target=sample, daemon=True)
        self.thread.start()
        return self

    @contextmanager
    def phase(self, name):
        print(f"START {name}", flush=True)
        started = time.monotonic()
        try:
            yield
        finally:
            self.phases[name] = time.monotonic() - started
            self.scratch_peak = max(self.scratch_peak, directory_size(self.output))
            write_json(self.output / "telemetry.json", self.report())
            print(f"END {name}: {self.phases[name]:.2f}s", flush=True)

    def report(self):
        return {"elapsed_seconds":time.monotonic()-self.start, "phases_seconds":self.phases,
                "sampled_process_tree_peak_rss_bytes":self.peak, "rss_sample_seconds":0.1,
                "sampled_scratch_peak_bytes":self.scratch_peak, "scratch_sample_seconds":5,
                "resource_error":self.resource_error,
                "limitations":"sampled peaks may miss brief maxima; cooperative admission plus job timeout"}

    def __exit__(self, *args):
        self.stop.set()
        self.thread.join()
        write_json(self.output / "telemetry.json", self.report())


def source_identity(binary, protocol):
    return {"git_revision":subprocess.check_output(["git","rev-parse","HEAD"], cwd=ROOT, text=True).strip(),
            "native_sha256":digest(binary), "protocol_sha256":digest(protocol),
            "native_version":subprocess.check_output([str(binary),"--version"], text=True).strip(),
            "platform":platform.platform(), "python":sys.version, "logical_cpus":os.cpu_count(),
            "github_run_id":os.environ.get("GITHUB_RUN_ID"),
            "dependencies":{p:importlib.metadata.version(p) for p in
                ("numpy","pyarrow","pysam","psutil","blake3","scipy","scikit-learn","joblib","threadpoolctl")},
            "source_sha256":{str(p.relative_to(ROOT)):digest(p) for p in
                sorted(HERE.glob("*.py")) + sorted((ROOT/"python/rosalind").glob("*.py"))}}


def native_options(protocol, state, binary):
    return {**protocol["extraction"], "sample":state["sample"], "binary":binary}


def scan_candidates(state, protocol, binary, out):
    from rosalind.evidence import iter_evidence
    loci, totals, seen = [], Counter(), 0
    opts = native_options(protocol, state, binary)
    with iter_evidence(state["reference"], state["original_bam"], regions=state["regions"],
                       fields=["depths","alleles"], workdir=out/"candidate-scan", **opts) as stream:
        for batch in stream:
            for row in batch.to_pylist():
                pos0 = int(row["pos"])-1
                region = locate(protocol["intervals"], row["contig"], pos0)
                totals[region["block_id"]+".positions"] += 1
                totals[region["block_id"]+".callable_depth_sum"] += int(row["callable_depth"])
                alts = candidate_alts(row, protocol["candidate_rule"])
                if alts:
                    loci.append({"contig":row["contig"], "pos0":pos0, "ref":row["ref"], "alts":alts,
                                 "partition":region["partition"], "block_id":region["block_id"]})
                    if len(loci) > protocol["limits"]["candidate_loci"]:
                        raise ValueError("candidate-locus budget exceeded; no outcome-based truncation")
                seen += 1
            if seen % 200704 == 0:
                print(f"Candidate scan: {seen} positions, {len(loci)} candidate loci", flush=True)
        if stream.result is None:
            raise ValueError("candidate scan did not complete successfully")
    if seen != sum(r["end"]-r["start"] for r in protocol["intervals"]):
        raise ValueError("candidate scan omitted measured positions")
    if sum(len(x["alts"]) for x in loci) > protocol["limits"]["candidate_rows"]:
        raise ValueError("candidate-row budget exceeded")
    write_jsonl(out/"candidates.jsonl", loci)
    write_json(out/"candidate-scan.json", {"rows":seen,"candidate_loci":len(loci),"block_totals":dict(totals),
                                          "truth_labels_used":False,"benchmark_mask_used_for_candidate_selection":False,
                                          "receipt":str(stream.result.manifest_path)})
    return loci


def write_sites(loci, state, path):
    import pysam
    with pysam.FastaFile(state["reference"]) as fasta, path.open("w") as out:
        out.write("##fileformat=VCFv4.2\n")
        for contig in dict.fromkeys(x["contig"] for x in loci):
            out.write(f"##contig=<ID={contig},length={fasta.get_reference_length(contig)}>\n")
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for row in loci:
            out.write(f"{row['contig']}\t{row['pos0']+1}\t.\t{row['ref']}\t{','.join(row['alts'])}\t.\tPASS\t.\n")


def extract_leaf(state, protocol, binary, sites, condition, out):
    from rosalind.evidence import materialize_evidence
    from rosalind import open_dataset
    folder = out/condition
    folder.mkdir()
    bam = state["original_bam"] if condition == "original" else state["downsampled_bam"]
    result = materialize_evidence(state["reference"], bam, folder/"evidence.arrow", sites=sites,
                fields=protocol["rich_fields"], cache_dir=folder/"cache", **native_options(protocol,state,binary))
    receipt = json.loads(result.manifest_path.read_text())
    leaf = Path(receipt["measurements"]["execution.evidence_dataset_manifest"])
    open_dataset(leaf, binary=binary).verify(memory_budget_mb=protocol["extraction"]["memory_budget_mb"])
    return leaf


def read_native(leaf, binary, protocol, workdir):
    from rosalind import open_dataset
    rows = []
    with open_dataset(leaf,binary=binary).batches(workdir=workdir,
            fields=protocol["rich_fields"], memory_budget_mb=protocol["extraction"]["memory_budget_mb"]) as stream:
        for batch in stream:
            rows.extend(batch.to_pylist())
            if len(rows) > protocol["limits"]["candidate_loci"]:
                raise ValueError("saved reader exceeded admitted rows")
        if stream.result is None:
            raise ValueError("saved reader did not complete")
    return rows


def parquet_roundtrip(rows, directory):
    import pyarrow as pa
    import pyarrow.parquet as pq
    directory.mkdir()
    # All source fields keep their native exact integer/vector types.
    scalar_names = set(rows[0]) - {"contig","ref","requested_alts","pos"} if rows else set()
    arrays = {}
    for name in rows[0] if rows else []:
        vals = [r[name] for r in rows]
        typ = pa.list_(pa.uint64(),4) if name.startswith("allele_") else (pa.uint64() if name in scalar_names else (pa.uint32() if name == "pos" else pa.string()))
        arrays[name] = pa.array(vals,type=typ)
    table = pa.table(arrays)
    path = directory/"evidence.parquet"
    pq.write_table(table,path,compression="zstd",row_group_size=1024)
    before = digest(path)
    result = []
    for batch in pq.ParquetFile(path).iter_batches(batch_size=1024,use_threads=False):
        result.extend(batch.to_pylist())
    if before != digest(path) or rows != result:
        raise ValueError("Parquet raw observation identity mismatch")
    write_json(directory/"manifest.json", {"sha256":before,"rows":len(rows),"bytes":path.stat().st_size,
               "scope":"content checked columnar projection of identical native observations; not independent extraction"})
    return result


def feature_records(rows, loci, condition, truth):
    from features import feature_vectors
    indexed = {(r["contig"],int(r["pos"])-1):r for r in rows}
    if len(indexed) != len(rows) or set(indexed) != {(r["contig"],r["pos0"]) for r in loci}:
        raise ValueError("candidate observations missing/duplicated; absent rows cannot become zeros")
    result = []
    for locus in loci:
        row = indexed[locus["contig"],locus["pos0"]]
        if row["ref"] != locus["ref"]:
            raise ValueError("candidate/evidence reference mismatch")
        for alt in locus["alts"]:
            label, reason = truth.label(locus["contig"],locus["pos0"],locus["ref"],alt)
            result.append({k:locus[k] for k in ("contig","pos0","ref","partition","block_id")} |
                          {"candidate_id":f"{locus['contig']}:{locus['pos0']}:{locus['ref']}:{alt}",
                           "alt":alt,"condition":condition,"label":label,"label_reason":reason,
                           **feature_vectors(row,alt)})
    return result


def audit_counts(rows, state, condition, out):
    import pysam
    spec=importlib.util.spec_from_file_location("independent_cigar_oracle", ROOT/"benchmarks/evidence/pysam_evidence.py")
    oracle=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(oracle)
    chosen=sorted(rows,key=lambda r:hashlib.sha256(f"{r['contig']}:{r['pos']}".encode()).hexdigest())[:64]
    bam_path=state["original_bam"] if condition=="original" else state["downsampled_bam"]
    with pysam.AlignmentFile(bam_path,"rb") as bam, pysam.FastaFile(state["reference"]) as fasta:
        for row in chosen:
            observed=oracle.row_for(bam,fasta,row["contig"],int(row["pos"])-1,row["ref"],row["requested_alts"])
            for name,value in zip(oracle.HEADER[4:-2],observed[4:-2]):
                if row[name] != value:
                    raise ValueError(f"independent count mismatch at {row['contig']}:{row['pos']} {name}")
    report={"passed":True,"audited_loci":len(chosen),"fields":oracle.HEADER[4:-2],
            "oracle_sha256":digest(ROOT/"benchmarks/evidence/pysam_evidence.py"),
            "scope":"independent CIGAR scalar/strand/filter count check; per-allele sums have independent fixture tests"}
    write_json(out/f"count-audit-{condition}.json",report)
    return report


def caller_scores(state, sites, out):
    import pysam
    import pysam.bcftools
    raw=out/"caller-pileup.bcf"
    calls=out/"caller.vcf"
    arguments=["-f",state["reference"],"-R",state["regions"],"-q","20","-Q","20","-d","100000","-Ou","-o",str(raw),state["original_bam"]]
    pysam.bcftools.mpileup(*arguments,catch_stdout=False)
    pysam.bcftools.call("-m","-v","-Ov","-o",str(calls),str(raw),catch_stdout=False)
    scores={}
    with pysam.VariantFile(str(calls)) as vcf:
        for record in vcf:
            gt=next(iter(record.samples.values())).get("GT",())
            for i,alt in enumerate(record.alts or (),1):
                if i in gt and len(record.ref)==1 and len(alt)==1 and alt in "ACGT":
                    scores[(record.contig,record.start,record.ref,alt)]=float(record.qual or 0)
    raw.unlink()  # Identified regenerable pileup; retain called VCF and invocation.
    write_json(out/"caller-identity.json",{"pysam":pysam.__version__,"samtools":pysam.__samtools_version__,
          "mpileup_argv":arguments,"call_argv":["-m","-v","-Ov"],
          "semantics":"bcftools default BAQ/overlap treatment, MAPQ/BQ20; separate caller baseline, not exact-count parity"})
    return scores


def run(args):
    from prepare import prepare_inputs, verify_prepared_inputs
    from labels import TruthIndex
    from study import run_study, predict_records
    out=args.output.absolute()
    out.mkdir(parents=True,exist_ok=False)
    protocol=json.loads(args.protocol.read_text())
    validate_protocol(protocol)
    binary=args.binary.resolve(strict=True)
    write_json(out/"protocol.json",protocol)
    identity=source_identity(binary,args.protocol)
    identity["resource_registry_sha256"]=digest(args.resources)
    write_json(out/"identity.json",identity)
    def interrupted(signum,frame):
        raise TimeoutError(f"analysis interrupted by signal {signum}; no completed output")
    signal.signal(signal.SIGALRM,interrupted)
    signal.signal(signal.SIGTERM,interrupted)
    signal.alarm(protocol["limits"]["analysis_seconds"])
    report={"schema":"rosalind-evidence-sufficiency-result-v1","status":"running","p1_authorized_by_gate":False}
    with Telemetry(out,protocol) as telemetry:
        try:
            with telemetry.phase("prepare"):
                state=prepare_inputs(args.resources,protocol,out/"inputs")
                write_json(out/"prepared.json",state)
            with telemetry.phase("candidate_generation"):
                loci=scan_candidates(state,protocol,binary,out)
                if not loci:
                    raise ValueError("no candidates; study is data-limited")
                sites=out/"candidates.vcf"
                write_sites(loci,state,sites)
                sites_digest=digest(sites)
            leaves, original_rows, parquet_rows, audits = {}, {}, {}, {}
            for condition in ("original","downsampled"):
                with telemetry.phase(f"extract_{condition}"):
                    leaves[condition]=extract_leaf(state,protocol,binary,sites,condition,out)
                    original_rows[condition]=read_native(leaves[condition],binary,protocol,out/f"read-{condition}")
                with telemetry.phase(f"audit_{condition}"):
                    audits[condition]=audit_counts(original_rows[condition],state,condition,out)
                with telemetry.phase(f"parquet_{condition}"):
                    parquet_rows[condition]=parquet_roundtrip(original_rows[condition],out/f"parquet-{condition}")
            retained=sum(directory_size(out/c) for c in ("original","downsampled","parquet-original","parquet-downsampled"))
            consumed_files=freeze_files(out/c for c in ("original","downsampled","parquet-original","parquet-downsampled"))
            if retained>protocol["limits"]["retained_evidence_bytes"]:
                raise ValueError("retained evidence budget exceeded")
            with telemetry.phase("labels_and_features"):
                truth=TruthIndex(state["truth_vcf"],state["truth_bed"],protocol["intervals"],reference_path=state["reference"])
                records=[]
                parquet_records=[]
                for condition in ("original","downsampled"):
                    records.extend(feature_records(original_rows[condition],loci,condition,truth))
                    parquet_records.extend(feature_records(parquet_rows[condition],loci,condition,truth))
                if records != parquet_records:
                    raise ValueError("native/Parquet feature equivalence failed")
                proposed={(r["contig"],r["pos0"],r["ref"],alt) for r in loci for alt in r["alts"]}
                generation={}
                for region in protocol["intervals"]:
                    keys={k for k in truth.positives if k[0]==region["contig"] and region["start"]<=k[1]<region["end"]}
                    generation[region["block_id"]]={"truth_snvs":len(keys),"proposed_truth_snvs":len(keys&proposed),
                                                    "candidate_generation_recall":len(keys&proposed)/len(keys) if keys else None}
                write_json(out/"candidate-recall.json",generation)
                write_json(out/"label-domain.json",truth.domain_report())
                truth.close()
            with telemetry.phase("caller_baseline"):
                scores=caller_scores(state,sites,out)
                for collection in (records,parquet_records):
                    for r in collection:
                        if r["condition"]=="original":
                            r["caller_score"]=scores.get((r["contig"],r["pos0"],r["ref"],r["alt"]))
                write_jsonl(out/"records.jsonl",records)
            with telemetry.phase("saved_only_replay"):
                hidden=out/"inputs.raw-unavailable"
                (out/"inputs").rename(hidden)
                try:
                    replay={}
                    for condition in ("original","downsampled"):
                        replay[condition]=read_native(leaves[condition],binary,protocol,out/f"replay-{condition}")
                        if replay[condition]!=original_rows[condition]:
                            raise ValueError("saved-only replay differs from original evidence")
                finally:
                    hidden.rename(out/"inputs")
                write_json(out/"replay.json",{"passed":True,"original_raw_paths_unavailable":True,"conditions":list(replay)})
            with telemetry.phase("model_study"):
                del original_rows, parquet_rows, replay
                study=run_study(records,out/"study",protocol)
                write_json(out/"truth-recovery.json",truth_recovery(study,generation,protocol["intervals"]))
            model_path=out/"study/frozen-models.pkl"
            if model_path.exists():
                with telemetry.phase("model_backend_equivalence"):
                    native_predictions=predict_records(model_path,records)
                    native_hash=hashlib.sha256(canonical(native_predictions).encode()).hexdigest()
                    del native_predictions
                    columnar_predictions=predict_records(model_path,parquet_records)
                    columnar_hash=hashlib.sha256(canonical(columnar_predictions).encode()).hexdigest()
                    del columnar_predictions
                    if native_hash!=columnar_hash:
                        raise ValueError("model predictions differ by storage path")
                    write_json(out/"prediction-equivalence.json",{"passed":True,"sha256":native_hash})
            if digest(binary)!=identity["native_sha256"] or digest(args.protocol)!=identity["protocol_sha256"]:
                raise ValueError("source executable or protocol changed during analysis")
            if digest(args.resources)!=identity["resource_registry_sha256"]:
                raise ValueError("resource registry changed during analysis")
            for relative,expected in identity["source_sha256"].items():
                if digest(ROOT/relative)!=expected:
                    raise ValueError("analysis source changed during execution")
            with telemetry.phase("final_identity_check"):
                verified=verify_prepared_inputs(state)
                if digest(sites)!=sites_digest:
                    raise ValueError("candidate sites changed during execution")
                for path,expected in consumed_files.items():
                    if digest(path)!=expected:
                        raise ValueError("consumed saved evidence changed during execution")
                write_json(out/"input-recheck.json",verified)
            gate=json.loads((out/"study/gate.json").read_text())
            report.update(status="completed",study=study,gate=gate,correctness={"independent_counts":True,
                         "native_parquet_features":True,"saved_only_replay":True,"prediction_equivalence":model_path.exists()},
                         retained_evidence_bytes=retained,scope=protocol["primary_scope"])
            # Explicit recognized outcomes only: never promote an unknown gate shape.
            report["p1_authorized_by_gate"]=advancement_allowed(gate,report["correctness"])
            write_json(out/"report.json",report)
            inventory={str(p.relative_to(out)):{"sha256":digest(p),"bytes":p.stat().st_size}
                       for p in out.rglob("*") if p.is_file() and "inputs" not in p.relative_to(out).parts and p.name not in ("run-manifest.json","telemetry.json")}
            write_json(out/"run-manifest.json",{"status":"completed","identity":identity,"files":inventory,
                       "scope":"retained generated results; original input identities in prepared.json; no biological authentication"})
            print(canonical({"status":report["status"],"gate":gate,"p1_authorized_by_gate":report["p1_authorized_by_gate"]}),flush=True)
        except BaseException as error:
            (out/"report.json").unlink(missing_ok=True)
            (out/"run-manifest.json").unlink(missing_ok=True)
            report.update(status="invalid",error=repr(error),traceback=traceback.format_exc())
            write_json(out/"failure.json",report)
            raise
        finally:
            signal.alarm(0)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary",type=Path,required=True)
    parser.add_argument("--protocol",type=Path,default=HERE/"protocol-p0.json")
    parser.add_argument("--resources",type=Path,default=ROOT/"benchmarks/giab/resources.tsv")
    parser.add_argument("--output",type=Path,required=True)
    run(parser.parse_args())


if __name__=="__main__":
    main()
