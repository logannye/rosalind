#!/usr/bin/env python3
"""Run and record a reproducible, public-data evidence-reuse demonstration."""
import argparse
import csv
import hashlib
import io
import json
import platform
from pathlib import Path
import shlex
import shutil
import subprocess
import sys


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(65536), b""):
            digest.update(block)
    return digest.hexdigest()


def table(path):
    return list(csv.DictReader(io.StringIO(path.read_text().lstrip("#")), delimiter="\t"))


def normalize(text, replacements):
    # Both spellings matter on macOS, where /tmp commonly resolves to /private/tmp.
    for value, label in sorted(replacements.items(), key=lambda item: len(item[0]), reverse=True):
        text = text.replace(value, label)
    return text


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True, help="exact native executable; never auto-selected")
    parser.add_argument("--expected-version", required=True, help="exact expected `rosalind --version` version")
    parser.add_argument("--output", type=Path, required=True, help="new directory; retained for inspection")
    parser.add_argument("--label", choices=("source-preview", "release-candidate", "stable"), default="source-preview",
                        help="operator-supplied distribution label, not proof of publication")
    args = parser.parse_args(argv)
    binary = args.binary.resolve(strict=True)
    actual_version = subprocess.check_output([str(binary), "--version"], text=True).strip().split()[-1]
    if actual_version != args.expected_version:
        parser.error(f"binary version {actual_version!r} does not match expected {args.expected_version!r}")
    import pysam
    if pysam.__version__ != "0.23.3":
        parser.error("preparation requires pysam==0.23.3 in the selected Python environment")
    repository = Path(__file__).resolve().parents[1]
    output = args.output.absolute()
    output.mkdir(parents=True, exist_ok=False)
    output = output.resolve()
    replacements = {str(binary): "$ROSALIND", str(args.binary.absolute()): "$ROSALIND",
                    str(output): "$DEMO", str(args.output.absolute()): "$DEMO",
                    str(repository): "$REPOSITORY", str(Path(sys.executable).absolute()): "$PYTHON"}
    transcript = []
    commands = []

    def emit(message):
        message = normalize(str(message), replacements)
        print(message, flush=True)
        transcript.append(message)
        (output / "transcript.txt").write_text("\n".join(transcript) + "\n")

    def run(*arguments, executable=binary):
        command = [str(executable), *map(str, arguments)]
        emit("$ " + shlex.join(command))
        result = subprocess.run(command, cwd=output, text=True, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, check=False)
        if result.stderr.strip():
            emit(result.stderr.rstrip())
        if result.stdout.strip():
            emit(result.stdout.rstrip())
        commands.append({"argv": [normalize(value, replacements) for value in command],
                         "returncode": result.returncode})
        if result.returncode:
            raise subprocess.CalledProcessError(result.returncode, command)
        return result.stdout

    report = {"schema": 1, "distribution_label": args.label,
              "label_is_operator_assertion": True, "binary_version": actual_version,
              "binary_sha256": sha256(binary), "platform": platform.system(),
              "architecture": platform.machine(), "python": platform.python_version(),
              "pysam": pysam.__version__, "commands": commands}
    emit(f"Rosalind evidence reuse demonstration ({args.label})")
    emit(f"Version: {actual_version}; executable SHA256: {report['binary_sha256']}")
    emit("Public NA18507 tutorial slices; research evidence, not a calling-accuracy or speed benchmark.")
    try:
        emit("\n1. Prepare pinned public inputs and inspect four supplied candidate SNVs.")
        inputs = output / "inputs"
        run(repository / "examples/research-filter/prepare.py", inputs, executable=sys.executable)
        preparation = json.loads((inputs / "preparation.json").read_text())
        report["preparation"] = json.loads(normalize(json.dumps(preparation), replacements))
        common = ("--reference", inputs / "ex1.fa", "--alignments", inputs / "sample.bam",
                  "--memory-budget-mb", "256")
        run("analyze", "evidence", *common, "--sites", inputs / "candidates.vcf",
            "--output", output / "candidate-evidence.tsv")
        run("verify", "--manifest", output / "candidate-evidence.tsv.manifest.json")
        candidates = table(output / "candidate-evidence.tsv")
        expected = [("seq1", "548", 36, 17), ("seq1", "1294", 34, 16),
                    ("seq2", "505", 44, 22), ("seq2", "1344", 27, 12)]
        observed = [(r["contig"], r["pos"], int(r["callable_depth"]),
                     sum(int(r[base.lower()]) for base in r["requested_alts"].split(",")))
                    for r in candidates]
        if observed != expected:
            raise ValueError(f"pinned candidate evidence changed: {observed!r}")
        emit("contig\tposition\tcallable_depth\tALT_read_observations")
        for row in observed:
            emit("\t".join(map(str, row)))
        report["candidate_evidence"] = [list(row) for row in observed]

        emit("\n2. Save every target position for later candidate queries and panel QC.")
        run("analyze", "evidence", *common, "--regions", inputs / "targets.bed",
            "--fields", "depths,alleles,quality-sums", "--format", "arrow-ipc",
            "--cache-dir", output / "cache", "--output", output / "target-evidence.arrow")
        receipt = json.loads((output / "target-evidence.arrow.manifest.json").read_text())
        manifest = Path(receipt["measurements"]["execution.evidence_dataset_manifest"])
        relocated = output / "relocated"
        relocated.mkdir()
        selection = relocated / "second-candidates.vcf"
        lines = (inputs / "candidates.vcf").read_text().splitlines(keepends=True)
        headers = [line for line in lines if line.startswith("#")]
        records = [line for line in lines if not line.startswith("#")]
        if len(records) != 4:
            raise ValueError("expected four supplied candidate records")
        selection.write_text("".join(headers + records[:2]))
        shutil.copy2(inputs / "targets.bed", relocated / "targets.bed")

        emit("\n3. Keep fresh answers to check the changed two-candidate shortlist and panel summary.")
        run("analyze", "evidence", *common, "--sites", selection, "--fields", "depths,alleles",
            "--output", output / "fresh-candidates.tsv")
        run("analyze", "panel-qc", *common, "--regions", inputs / "targets.bed",
            "--min-callable-depth", "10", "--output", output / "fresh-panel.tsv")
        # This is an ordinary move of the entire portable dataset, not a manifest rewrite.
        shutil.move(str(manifest.parent), str(relocated / "portable"))
        # Delete only the public inputs created by this invocation under a new output directory.
        shutil.rmtree(inputs)
        if inputs.exists():
            raise ValueError("generated input removal failed")
        emit("Moved the complete portable dataset to $DEMO/relocated/portable.")
        emit("Removed only this run's generated $DEMO/inputs (including BAM, FASTA, and indexes).")

        emit("\n4. Answer both questions from the relocated evidence alone.")
        dataset = relocated / "portable/evidence-dataset.manifest.json"
        verified = json.loads(run("dataset", "verify", "--dataset", dataset))
        if verified["verified_loci"] != 3159:
            raise ValueError("pinned target position count changed")
        run("dataset", "extract", "--dataset", dataset, "--sites", selection,
            "--fields", "depths,alleles", "--memory-budget-mb", "256",
            "--format", "tsv", "--output", output / "reused-candidates.tsv")
        run("dataset", "panel-qc", "--dataset", dataset, "--regions", relocated / "targets.bed",
            "--min-callable-depth", "10", "--memory-budget-mb", "256",
            "--output", output / "reused-panel.tsv")
        report["verified_loci"] = verified["verified_loci"]
        report["outputs"] = {}
        for fresh, reused in (("fresh-candidates.tsv", "reused-candidates.tsv"),
                              ("fresh-panel.tsv", "reused-panel.tsv")):
            if (output / fresh).read_bytes() != (output / reused).read_bytes():
                raise ValueError(f"fresh/reused evidence differs: {reused}")
            run("verify", "--manifest", output / (reused + ".manifest.json"))
            report["outputs"][reused] = {"sha256": sha256(output / reused), "equals_fresh_bytes": True}
            emit(f"PASS: {reused} equals {fresh} byte-for-byte.")
        panel = table(output / "reused-panel.tsv")
        if [(int(r["length"]), int(r["callable_positions"])) for r in panel] != [(1575, 1493), (1584, 1517)]:
            raise ValueError("pinned panel denominators or depth-eligible counts changed")
        report["panel_summary"] = [{k: r[k] for k in ("target", "length", "callable_positions", "callable_threshold")}
                                   for r in panel]

        emit("\n5. Replay the saved-evidence candidate result with the explicitly selected binary.")
        replay = json.loads(run("reproduce", "--manifest", output / "reused-candidates.tsv.manifest.json",
                                "--inputs", relocated, "--binary", binary, "--json", "--no-attest"))
        if replay.get("verdict") != "REPRODUCED" or replay.get("code_identity_matches") is not True:
            raise ValueError("saved-evidence replay did not reproduce with matching code identity")
        if sha256(binary) != report["binary_sha256"]:
            raise ValueError("selected executable changed during the demonstration")
        report["selected_executable_unchanged"] = True
        report["replay"] = replay
        report["original_inputs_absent"] = not inputs.exists()
        report["passed"] = True
        emit("PASS: four candidate records; 3,159 saved positions; two-candidate reuse and panel QC match fresh extraction.")
        emit("Depth 10 is an illustrative technical screen. Stored zero depth differs from an unmeasured locus.")
        emit("Retained outputs, receipts, sanitized transcript, and demo-report.json are in $DEMO.")
    except Exception as error:
        report["passed"] = False
        report["error"] = normalize(str(error), replacements)
        emit(f"FAILED: {error}")
        raise
    finally:
        (output / "demo-report.json").write_text(normalize(json.dumps(report, indent=2), replacements) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
