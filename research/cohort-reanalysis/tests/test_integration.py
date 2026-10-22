"""Native integration gates; require an explicitly selected current-source binary.

These execute only when ROSALIND_BINARY is set. Hosted CI builds that binary from
its checked-out revision. This file never builds or installs a native package.
The authored fixture and its independent validator remain unchanged.
"""
from __future__ import annotations

import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
FIXTURE = ROOT / "examples" / "cohort-reanalysis"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "python"))


def rows(path):
    return [json.loads(line) for line in path.read_text().splitlines() if line]


def sites(path, candidates):
    path.write_text("##fileformat=VCFv4.2\n##contig=<ID=synthetic,length=64>\n"
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
                    "".join(f"synthetic\t{p}\t.\t{ref}\t{alt}\t.\t.\t.\n"
                            for p, ref, alt in candidates))
    return path


@unittest.skipUnless(os.environ.get("ROSALIND_BINARY"), "select a built binary with ROSALIND_BINARY; no implicit build")
class CohortIntegration(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory(prefix="rosalind-research-integration-")
        cls.root = Path(cls.temp.name)
        cls.binary = Path(os.environ["ROSALIND_BINARY"]).resolve(strict=True)
        cls.prepared = cls.root / "fixture"
        subprocess.run([sys.executable, str(FIXTURE / "prepare.py"), str(cls.prepared)], check=True)
        cls.expected = json.loads((FIXTURE / "expected.json").read_text())
        cls.members = []
        for specimen in cls.expected["members"]:
            output = cls.root / (specimen + ".arrow")
            subprocess.run([str(cls.binary), "analyze", "evidence", "--reference",
                            str(cls.prepared / "raw/reference.fa"), "--alignments",
                            str(cls.prepared / "raw" / (specimen + ".bam")), "--regions",
                            str(FIXTURE / (specimen + ".bed")), "--fields", "depths,alleles",
                            "--cache-dir", str(cls.root / (specimen + "-cache")),
                            "--memory-budget-mb", "256", "--format", "arrow-ipc",
                            "--output", str(output)], check=True, capture_output=True)
            receipt = json.loads(Path(str(output) + ".manifest.json").read_text())
            leaf = Path(receipt["measurements"]["execution.evidence_dataset_manifest"])
            cls.members.append({"specimen_id": specimen, "dataset": str(leaf)})
        cls.manifest = cls.root / "members.json"
        cls.manifest.write_text(json.dumps({"version": 1, "members": cls.members}))

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def setUp(self):
        self.directory = Path(tempfile.mkdtemp(dir=self.root))

    def command(self, output, *, members=None, selection=None, window=2,
                thresholds=(10,), policy="partial", native=256, parent=128):
        command = [sys.executable, str(HERE / "analyze.py"), "--binary", str(self.binary),
                   "--members", str(members or self.manifest), "--sites",
                   str(selection or FIXTURE / "candidates.vcf"), "--output", str(output),
                   "--missingness", policy, "--window-size", str(window),
                   "--native-memory-budget-mb", str(native), "--parent-memory-budget-mb", str(parent)]
        for threshold in thresholds:
            command += ["--min-callable-depth", str(threshold)]
        return command

    def run_cli(self, output, succeeds=True, **kwargs):
        result = subprocess.run(self.command(output, **kwargs), text=True, capture_output=True)
        if succeeds:
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertTrue((output / "run-manifest.json").is_file())
        else:
            self.assertNotEqual(result.returncode, 0)
            self.assertFalse((output / "run-manifest.json").exists())
        return result

    def assert_oracle(self, output, thresholds=(10,), candidates=None):
        candidates = candidates or [(int(p), "A", a) for p, a in self.expected["candidates"].items()]
        actual = rows(output / "summary.jsonl")
        expected = []
        for pos, ref, alt in sorted(candidates, key=lambda r: (r[0], r[2])):
            measured = [r for m in self.expected["members"].values()
                        for r in m["full_rows"] if r[0] == pos and pos in m["stored_positions"]]
            for threshold in sorted(thresholds):
                eligible = [r for r in measured if r[3] >= threshold]
                depth = sum(r[3] for r in measured)
                ac = sum(r[4 + "ACGT".index(alt)] for r in measured)
                ed = sum(r[3] for r in eligible)
                ea = sum(r[4 + "ACGT".index(alt)] for r in eligible)
                expected.append({"contig": "synthetic", "pos": pos, "ref": ref, "alt": alt,
                    "min_callable_depth": threshold, "n_requested": 3, "n_observed": len(measured),
                    "n_depth_eligible": len(eligible),
                    "n_alt_supported": sum(r[4 + "ACGT".index(alt)] > 0 for r in eligible),
                    "total_callable_depth": depth, "total_alt_count": ac,
                    "eligible_callable_depth": ed, "eligible_alt_count": ea,
                    "alt_fraction_numerator": ac if depth else None,
                    "alt_fraction_denominator": depth or None,
                    "eligible_alt_fraction_numerator": ea if ed else None,
                    "eligible_alt_fraction_denominator": ed or None,
                    "reference_validation": "checked_on_observed_rows" if measured else "unavailable"})
        self.assertEqual(actual, expected)

    def test_partial_fixture_and_strict_refusal(self):
        output = self.directory / "partial"
        self.run_cli(output)
        self.assert_oracle(output)
        self.run_cli(self.directory / "strict", policy="strict", succeeds=False)
        selection = sites(self.directory / "complete.vcf", [(10, "A", "C"), (20, "A", "G")])
        self.run_cli(self.directory / "complete", selection=selection, policy="strict")
        self.assert_oracle(self.directory / "complete", candidates=[(10, "A", "C"), (20, "A", "G")])

    def test_windows_budgets_thresholds_order_and_duplicate_invariance(self):
        reference = None
        reversed_members = self.directory / "reversed-members.json"
        reversed_members.write_text(json.dumps({"version": 1, "members": list(reversed(self.members))}))
        selection = sites(self.directory / "reversed.vcf", [(40, "A", "C"), (30, "A", "T"),
                          (20, "A", "G"), (10, "A", "C"), (10, "A", "C")])
        for window, native, parent in [(1, 128, 96), (2, 256, 128), (4, 512, 128)]:
            output = self.directory / f"run-{window}"
            self.run_cli(output, members=reversed_members, selection=selection, window=window,
                         native=native, parent=parent, thresholds=(20, 5, 10))
            self.assert_oracle(output, thresholds=(5, 10, 20))
            content = (output / "summary.jsonl").read_bytes()
            if reference is not None:
                self.assertEqual(content, reference)
            reference = content

    def test_multialt_and_completely_unobserved_reference_annotation(self):
        candidates = [(10, "A", "G"), (10, "A", "C"), (50, "G", "T")]
        selection = sites(self.directory / "new-alt.vcf", candidates)
        output = self.directory / "new-alt"
        self.run_cli(output, selection=selection, window=1, thresholds=(5, 10, 20))
        self.assert_oracle(output, candidates=candidates, thresholds=(5, 10, 20))

    def test_duplicate_source_and_observed_ref_mismatch_refuse(self):
        members = self.directory / "aliases.json"
        aliases = self.members + [{"specimen_id": "another-name", "dataset": self.members[0]["dataset"]}]
        members.write_text(json.dumps({"version": 1, "members": aliases}))
        self.run_cli(self.directory / "duplicate", members=members, succeeds=False)
        selection = sites(self.directory / "wrong-ref.vcf", [(10, "G", "C")])
        self.run_cli(self.directory / "wrong-ref", selection=selection, succeeds=False)

    def test_offline_relocation_replay_and_tampered_summary_verification(self):
        original = self.directory / "original"
        self.run_cli(original)
        relocated = []
        for member in self.members:
            source = Path(member["dataset"])
            destination = self.directory / "portable" / member["specimen_id"]
            shutil.copytree(source.parent, destination)
            relocated.append({"specimen_id": member["specimen_id"], "dataset": str(destination / source.name)})
        manifest = self.directory / "relocated.json"
        manifest.write_text(json.dumps({"version": 1, "members": relocated}))
        raw = self.prepared / "raw"
        unavailable = self.prepared / "generated-raw-unavailable"
        raw.rename(unavailable)
        try:
            replay = self.directory / "replayed"
            result = subprocess.run([sys.executable, str(HERE / "analyze.py"), "replay", str(original),
                                     "--members", str(manifest), "--output", str(replay),
                                     "--binary", str(self.binary), "--window-size", "1"],
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertEqual((original / "summary.jsonl").read_bytes(), (replay / "summary.jsonl").read_bytes())
        finally:
            unavailable.rename(raw)
        good = subprocess.run([sys.executable, str(HERE / "analyze.py"), "verify", str(original)],
                              capture_output=True, text=True)
        self.assertEqual(good.returncode, 0, good.stderr)
        with (original / "summary.jsonl").open("a") as stream:
            stream.write("{}\n")
        bad = subprocess.run([sys.executable, str(HERE / "analyze.py"), "verify", str(original)], capture_output=True)
        self.assertNotEqual(bad.returncode, 0)

    def load_analysis(self):
        import analyze
        return analyze

    def test_prior_child_partition_mutation_prevents_publication(self):
        analyze = self.load_analysis()
        output = self.directory / "mutated"
        args = analyze.build_parser().parse_args(self.command(output, window=4)[2:])
        changed = []

        def mutate(*_args, **_kwargs):
            if not changed:
                # The first successful child has consumed specimen-a's only part.
                part = next(Path(self.members[0]["dataset"]).parent.glob("*/evidence.arrow"))
                before = part.read_bytes()
                changed.append((part, before))
                part.write_bytes(before + b"changed-after-successful-child")

        try:
            with self.assertRaisesRegex(ValueError, "input changed before publication"):
                analyze.run_analysis(args, after_child=mutate)
            self.assertTrue(changed, "must mutate after a successful child, not fail before execution")
            self.assertFalse((output / "run-manifest.json").exists())
        finally:
            for part, before in changed:
                part.write_bytes(before)

    def test_future_member_coverage_mutation_prevents_publication(self):
        analyze = self.load_analysis()
        output = self.directory / "coverage-mutated"
        args = analyze.build_parser().parse_args(self.command(output, window=4)[2:])
        changed = []

        def mutate(prepared, member_index, _receipt_path):
            if not changed:
                self.assertEqual(member_index, 0)
                # A later member's saved mask is an execution input. Removing a
                # site could otherwise turn real stored evidence into missingness.
                coverage_path = prepared.coverage_paths[1]
                coverage = json.loads(coverage_path.read_text())
                self.assertTrue(coverage["covered"]["sites"])
                coverage["covered"]["sites"].pop(0)
                coverage_path.write_text(json.dumps(coverage))
                changed.append(coverage_path)

        with self.assertRaisesRegex(ValueError, "input changed before publication"):
            analyze.run_analysis(args, after_child=mutate)
        self.assertTrue(changed, "must mutate after the first completed child")
        self.assertFalse(output.exists())
        self.assertFalse((output / "run-manifest.json").exists())

    def test_consumer_failure_and_cancellation_do_not_publish(self):
        analyze = self.load_analysis()
        for failure in (RuntimeError("simulated downstream consumer failure"), KeyboardInterrupt()):
            output = self.directory / type(failure).__name__
            args = analyze.build_parser().parse_args(self.command(output)[2:])
            called = []

            def fail(*_args, **_kwargs):
                called.append(True)
                raise failure

            with self.assertRaises(type(failure)):
                analyze.run_analysis(args, after_child=fail)
            self.assertTrue(called)
            self.assertFalse((output / "run-manifest.json").exists())

    def test_active_child_failure_and_cancellation_reap_native_processes(self):
        analyze = self.load_analysis()
        original_popen = subprocess.Popen
        for failure in (RuntimeError("consumer failed while native stream active"), KeyboardInterrupt()):
            output = self.directory / ("active-" + type(failure).__name__)
            args = analyze.build_parser().parse_args(self.command(output)[2:])
            processes = []

            def capture(*arguments, **kwargs):
                process = original_popen(*arguments, **kwargs)
                processes.append(process)
                return process

            with patch.object(subprocess, "Popen", side_effect=capture):
                with patch.object(analyze.WindowReducer, "observe", side_effect=failure) as observe:
                    with self.assertRaises(type(failure)):
                        analyze.run_analysis(args)
                    self.assertTrue(observe.called, "failure must occur while consuming a native batch")
            self.assertTrue(processes)
            self.assertTrue(all(process.poll() is not None for process in processes),
                            "active native child must be reaped after consumer failure/cancellation")
            self.assertFalse((output / "run-manifest.json").exists())

    def test_corrupted_partition_fails_without_success(self):
        source = Path(self.members[0]["dataset"])
        destination = self.directory / "corrupt-leaf"
        shutil.copytree(source.parent, destination)
        part = next(destination.glob("*/evidence.arrow"))
        with part.open("ab") as stream:
            stream.write(b"corruption")
        changed = [dict(m) for m in self.members]
        changed[0]["dataset"] = str(destination / source.name)
        manifest = self.directory / "corrupt-members.json"
        manifest.write_text(json.dumps({"version": 1, "members": changed}))
        self.run_cli(self.directory / "corrupt", members=manifest, succeeds=False)


if __name__ == "__main__":
    unittest.main()
