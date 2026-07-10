#!/usr/bin/env python3
"""Rosalind claims harness — re-run it, don't trust it.

Each check below is a property a skeptic can re-derive on their own machine. Any failure
makes the whole run exit non-zero, so it doubles as a CI regression gate and can't be
cherry-picked. Every row names the exact `rosalind` subcommand it runs.

This harness is about the **verifiable memory contract** and **byte-reproducibility** —
explicitly NOT speed (the engine is single-threaded) and NOT real-world accuracy
(calling here is on bundled, simulated, SNV-only toy data). The index *build* is
O(reference) and out of scope; only the `variants --index` / `features --index`
streaming paths are bounded. Machine-dependent numbers (peak RSS) are *recorded* but the
PASS/FAIL is on the qualitative property (predicted >= realized, exit codes, byte
equality, verdicts), never a fixed MiB value.
"""

import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

BIN = os.environ.get("ROSALIND_BIN", "target/release/rosalind")
ROOT = Path.cwd().resolve()
DATA = Path("examples/data/illumina_toy")
RESULTS = Path("benchmarks/results.json")
MiB = 1 << 20

rows = []


def add(claim, check, expected, observed, ok):
    rows.append(
        {
            "claim": claim,
            "check": check,
            "expected": expected,
            "observed": observed,
            "pass": bool(ok),
        }
    )


def run(args, cwd=None):
    return subprocess.run([BIN, *args], cwd=cwd, capture_output=True, text=True)


def manifest_of(vcf):
    return json.loads(Path(str(vcf) + ".manifest.json").read_text())


def main():
    global BIN
    if Path(BIN).exists() is False and shutil.which(BIN) is None:
        sys.exit(f"binary not found: {BIN} (build it: cargo build --release)")
    BIN = str(Path(BIN).resolve()) if Path(BIN).exists() else shutil.which(BIN)
    if not (DATA / "reference.fa").exists():
        sys.exit(f"bundled data missing: {DATA} (run from the repo root)")

    work = Path(tempfile.mkdtemp(prefix="rosalind-bench-"))
    try:
        idx = work / "ref.idx"
        sbam = work / "sorted.bam"
        if run(["index", "--reference", str(DATA / "reference.fa"), "--output", str(idx)]).returncode != 0:
            sys.exit("setup failed: rosalind index")
        if run(["sort", "--input", str(DATA / "alignments.bam"), "--output", str(sbam)]).returncode != 0:
            sys.exit("setup failed: rosalind sort")

        # 1) predicted peak (from the index header, no BAM) >= realized peak.
        plan = json.loads(run(["plan", "--index", str(idx), "--budget-mb", "4096", "--json"]).stdout)
        predicted = int(plan["predicted_peak_rss_bytes"])
        out1 = work / "c1.vcf"
        r1 = run(["variants", "--index", str(idx), "--alignments", str(sbam),
                  "--memory-budget-mb", "4096", "--enforce", "-o", str(out1)])
        realized = int(manifest_of(out1)["measurements"]["peak_rss_bytes"])
        add(
            "predicted peak is a conservative upper bound on realized peak",
            "rosalind plan --json predicted_peak_rss_bytes >= variants receipt peak_rss_bytes (same --max-depth, predicted before the BAM is read)",
            "predicted >= realized",
            f"predicted={predicted/MiB:.1f} MiB, realized={realized/MiB:.1f} MiB, margin={(predicted-realized)/MiB:.1f} MiB",
            r1.returncode == 0 and predicted >= realized,
        )

        # 2) honor-or-refuse before any work — never a silent OOM.
        v_fit = json.loads(run(["plan", "--index", str(idx), "--budget-mb", "512", "--json"]).stdout)["verdict"]
        v_ref = json.loads(run(["plan", "--index", str(idx), "--budget-mb", "1", "--json"]).stdout)["verdict"]
        refused = work / "refused.vcf"
        r2 = run(["variants", "--index", str(idx), "--alignments", str(sbam),
                  "--memory-budget-mb", "1", "--enforce", "-o", str(refused)])
        add(
            "the declared budget is honored or refused up front (never a silent OOM)",
            "plan verdict flips fits@512MiB -> refuse@1MiB; variants --enforce@1MiB exits 3 and writes NO output file",
            "fits / refuse / exit 3 + no output",
            f"plan@512={v_fit}, plan@1={v_ref}, enforce@1 exit={r2.returncode}, output_written={refused.exists()}",
            v_fit == "fits" and v_ref == "refuse" and r2.returncode == 3 and not refused.exists(),
        )

        # 3) identical inputs -> byte-identical TEXT output (VCF, TSV).
        vcfs, rec_hashes = [], []
        for i in range(3):
            d = work / f"run{i}"
            d.mkdir()
            o = d / "calls.vcf"
            run(["variants", "--index", str(idx), "--alignments", str(sbam), "-o", str(o)])
            vcfs.append(o.read_bytes())
            rec_hashes.append(manifest_of(o)["outputs"][0]["blake3"])
        tsvs = []
        for i in range(2):
            d = work / f"feat{i}"
            d.mkdir()
            o = d / "f.tsv"
            run(["features", "--index", str(idx), "--alignments", str(sbam), "-o", str(o)])
            tsvs.append(o.read_bytes())
        vcf_same = vcfs[0] == vcfs[1] == vcfs[2]
        add(
            "identical inputs produce byte-identical text outputs",
            "3 variants runs -> identical VCF bytes + one recorded blake3; 2 features runs -> identical TSV bytes",
            "all identical",
            f"vcf_bytes_identical={vcf_same}, recorded_blake3={rec_hashes[0][:10]}… (distinct={len(set(rec_hashes))}), tsv_bytes_identical={tsvs[0]==tsvs[1]}",
            vcf_same and len(set(rec_hashes)) == 1 and tsvs[0] == tsvs[1],
        )

        # 4) re-derive byte-for-byte from the receipt; tamper-evident.
        rep_out = work / "rep.vcf"
        run(["variants", "--index", str(idx), "--alignments", str(sbam), "-o", str(rep_out)])
        man_path = Path(str(rep_out) + ".manifest.json")
        rep = run(["reproduce", "--manifest", str(man_path), "--inputs", str(work)])
        reproduced = rep.returncode == 0 and "REPRODUCED" in rep.stdout
        man = json.loads(man_path.read_text())
        mb = man["params"]["manifest_blake3"]
        flipped = ("1" if mb[0] == "0" else "0") + mb[1:]
        tampered = work / "tampered.manifest.json"
        tampered.write_text(man_path.read_text().replace(mb, flipped, 1))
        tv = run(["verify", "--manifest", str(tampered)])
        add(
            "a recorded result re-derives byte-for-byte, and the receipt is tamper-evident",
            "reproduce -> exit 0 REPRODUCED; flip one byte of manifest_blake3 -> verify exits 5 (TAMPERED)",
            "REPRODUCED (0) / TAMPERED (5)",
            f"reproduce_exit={rep.returncode} (REPRODUCED={reproduced}), tampered_verify_exit={tv.returncode}",
            reproduced and tv.returncode == 5,
        )

        # 5) pack: additive predicted peaks -> a placement decision, run-free.
        jobs = work / "jobs.tsv"
        jobs.write_text(f"{idx}\t1000\t250\n{idx}\t1000\t250\n")
        pk = json.loads(run(["pack", "--jobs", str(jobs), "--node-mb", "64000", "--json"]).stdout)
        cap = int(pk["node_mb"]) * MiB
        within = all(0 < int(n["used_bytes"]) <= cap for n in pk["nodes"])
        toosmall = run(["pack", "--jobs", str(jobs), "--node-mb", "1", "--nodes", "1"])
        max_used = max(int(n["used_bytes"]) for n in pk["nodes"])
        add(
            "pack shows a co-location fits within capacity by additive predicted peaks, before launching a byte",
            "every node's summed predicted peak <= capacity; an impossible packing refuses (exit 3)",
            "within capacity / exit 3",
            f"nodes={len(pk['nodes'])}, max_node={max_used/MiB:.1f} MiB <= {pk['node_mb']} MiB ({within}), impossible_pack_exit={toosmall.returncode}",
            within and toosmall.returncode == 3,
        )

        # 6) A downstream binary inherits the complete contract without copying the CLI.
        project = work / "external-analyzer"
        scaffold = run(["new", "analyzer", "external-analyzer", "--output", str(project)])
        cargo_toml = project / "Cargo.toml"
        if scaffold.returncode == 0:
            cargo_lines = []
            for line in cargo_toml.read_text().splitlines():
                if line.startswith("rosalind-bio ="):
                    line = f'rosalind-bio = {{ path = "{ROOT}", features = ["contract-testkit"] }}'
                elif line.startswith("rosalind-build-info ="):
                    line = f'rosalind-build-info = {{ path = "{ROOT / "crates/build-info"}" }}'
                cargo_lines.append(line)
            cargo_toml.write_text("\n".join(cargo_lines) + "\n")
        built = subprocess.run(
            ["cargo", "build", "--quiet", "--offline"],
            cwd=project,
            capture_output=True,
            text=True,
        ) if scaffold.returncode == 0 else scaffold
        external_bin = project / "target" / "debug" / "external-analyzer"
        external_outputs = []
        for label, scale in (("a", 1), ("b", 1), ("scaled", 2)):
            output = work / f"external-{label}.tsv"
            external_outputs.append(output)
            if built.returncode == 0:
                subprocess.run(
                    [str(external_bin), "run", "--index", str(idx), "--alignments", str(sbam),
                     "--scale", str(scale), "--output", str(output)],
                    capture_output=True,
                    text=True,
                )
        ext_manifest = Path(str(external_outputs[0]) + ".manifest.json")
        verified = run(["verify", "--manifest", str(ext_manifest), "--json"]) if ext_manifest.exists() else None
        reproduced_ext = run([
            "reproduce", "--manifest", str(ext_manifest), "--inputs", str(work),
            "--binary", str(external_bin), "--json",
        ]) if ext_manifest.exists() else None
        scaled_manifest = Path(str(external_outputs[2]) + ".manifest.json")
        causal = run(["diff", str(ext_manifest), str(scaled_manifest)]) \
            if ext_manifest.exists() and scaled_manifest.exists() else None
        inherited = (
            scaffold.returncode == 0
            and built.returncode == 0
            and verified is not None and verified.returncode == 0
            and reproduced_ext is not None and reproduced_ext.returncode == 0
            and causal is not None and causal.returncode == 1
            and external_outputs[0].read_bytes() == external_outputs[1].read_bytes()
            and "analyzer.scale" in causal.stdout
        )
        add(
            "an external analyzer inherits the bounded, receipted, replayable contract",
            "scaffold -> offline build -> run -> verify -> reproduce --binary -> causal parameter diff",
            "all stages succeed; diff exits 1 with a parameter cause",
            f"scaffold={scaffold.returncode}, build={built.returncode}, verify={getattr(verified, 'returncode', None)}, reproduce={getattr(reproduced_ext, 'returncode', None)}, diff={getattr(causal, 'returncode', None)}",
            inherited,
        )

        # 7) The front-door demo is packaged, offline, and completes end to end.
        demo_dir = work / "demo"
        demo = run(["demo", "--output-dir", str(demo_dir), "--json"])
        try:
            demo_json = json.loads(demo.stdout)
        except json.JSONDecodeError:
            demo_json = {}
        demo_complete = (
            demo.returncode == 0
            and demo_json.get("ok") is True
            and (demo_dir / "calls.vcf").exists()
            and (demo_dir / "calls.vcf.manifest.json").exists()
            and (demo_dir / "calls.vcf.manifest.json.repro.json").exists()
        )
        add(
            "the packaged offline demo completes the full trust journey",
            "demo --json runs embedded index -> align -> sort -> plan -> enforce -> verify -> reproduce -> chain",
            "exit 0, ok=true, output + receipt + certificate exist",
            f"exit={demo.returncode}, ok={demo_json.get('ok')}, receipt={bool(demo_dir / 'calls.vcf.manifest.json')}",
            demo_complete,
        )
    finally:
        shutil.rmtree(work, ignore_errors=True)

    passed = sum(r["pass"] for r in rows)
    total = len(rows)
    RESULTS.parent.mkdir(exist_ok=True)
    RESULTS.write_text(
        json.dumps(
            {
                "schema": 1,
                "scope": "bundled simulated toy data; SNV-only; single-threaded; single-process; verifiability not speed/accuracy",
                "passed": passed,
                "total": total,
                "checks": rows,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )

    print("\nRosalind claims harness — verifiability, not speed or accuracy")
    print("scope: bundled simulated toy data · SNV-only · single-threaded · single-process\n")
    for r in rows:
        print(f"  [{'PASS' if r['pass'] else 'FAIL'}] {r['claim']}")
        print(f"         {r['observed']}")
    print(f"\n{passed}/{total} checks passed  →  wrote {RESULTS}")
    sys.exit(0 if passed == total else 1)


if __name__ == "__main__":
    main()
