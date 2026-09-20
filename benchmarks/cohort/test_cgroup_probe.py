import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

_spec = importlib.util.spec_from_file_location("cohort_cgroup_probe", Path(__file__).with_name("cgroup_probe.py"))
probe = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(probe)


class CohortProbeTests(unittest.TestCase):
    def observed(self, expected="completed", **changes):
        case = {"expected": expected, "hard_limit_bytes": 512 << 20,
                "require_os_limit": expected == "completed"}
        item = {"case": case, "analysis_exit_code": 0, "container_exit_code": 0,
                "docker_oom_killed": False, "events_before": {"oom_kill": 0},
                "events_after": {"oom_kill": 0}, "timed_out": False,
                "memory_max_bytes": 512 << 20, "swap_max_bytes": 0, "memory_peak_bytes": 64 << 20,
                "success_artifact": True, "success_receipt": True, "partial_artifact": False,
                "partial_receipt": False, "receipt_status": "completed", "verified": True,
                "assurance": "cgroup-v2", "files": {"result.tsv": {"sha256": "expected"}},
                "measurements": {"execution.original_alignment_records_decoded": "0",
                                 "resource.os_limit_bytes": str(512 << 20)}}
        item.update(changes)
        item["errors"] = probe.evidence.check_case(case, item)
        return item

    def test_completes_only_with_verified_equal_bytes_and_observed_hard_limits(self):
        self.assertEqual(probe.validate_case(self.observed(), "expected"), [])
        for changes in [dict(memory_max_bytes=256 << 20), dict(swap_max_bytes=1),
                        dict(memory_peak_bytes=None), dict(events_after=None), dict(events_before={}),
                        dict(verified=False), dict(assurance="observed-only"),
                        dict(measurements={"execution.original_alignment_records_decoded": "1"}),
                        dict(container_exit_code=137), dict(docker_oom_killed=True)]:
            self.assertTrue(probe.validate_case(self.observed(**changes), "expected"), changes)
        self.assertTrue(probe.validate_case(self.observed(), "different"))

    def test_tiny_budget_must_refuse_without_any_success_or_partial(self):
        refused = dict(analysis_exit_code=3, container_exit_code=3,
                       success_artifact=False, success_receipt=False, files={})
        self.assertEqual(probe.validate_case(self.observed("refused", **refused), "unused"), [])
        for changes in [dict(success_artifact=True), dict(success_receipt=True),
                        dict(partial_artifact=True), dict(partial_receipt=True),
                        dict(analysis_exit_code=4, container_exit_code=4), dict(docker_oom_killed=True)]:
            self.assertTrue(probe.validate_case(self.observed("refused", **dict(refused, **changes)), "unused"))

    def test_container_has_only_saved_inputs_network_disabled_and_swap_zero(self):
        argv = probe.container_arguments("unique", probe.DEFAULT_IMAGE, "/binary", "/saved-cohort",
                                         "/candidates", "/result", ["/inputs/rosalind", "--version"])
        self.assertEqual(argv[argv.index("--network") + 1], "none")
        self.assertEqual(argv[argv.index("--memory") + 1], "512m")
        self.assertEqual(argv[argv.index("--memory-swap") + 1], "512m")
        self.assertIn("--read-only", argv)
        mounts = [argv[index + 1] for index, value in enumerate(argv) if value == "--mount"]
        self.assertEqual(len(mounts), 4)
        self.assertTrue(all("readonly" in value for value in mounts[:3]))
        self.assertFalse(any("alignments" in value or "reference" in value for value in mounts))
        self.assertRegex(probe.DEFAULT_IMAGE, r"@sha256:[0-9a-f]{64}$")

    def test_actual_docker_independent_failure_is_retained_not_claimed_as_pass(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            binary = root / "non-linux"
            binary.write_bytes(b"not an ELF binary")
            output = root / "report"
            code = probe.main(["--binary", str(binary), "--demo-report", str(root / "missing.json"),
                               "--output", str(output)])
            self.assertEqual(code, 1)
            report = json.loads((output / "report.json").read_text())
            self.assertEqual(report["status"], "failed")
            self.assertIn("Linux/amd64 ELF", report["failure"])
            self.assertEqual(report["cases"], [])
            self.assertTrue((output / "harness.py").is_file())
            self.assertTrue((output / "evidence_cgroup_helpers.py").is_file())

    def test_native_argv_requires_os_assurance_only_for_admitted_case(self):
        admitted = probe.command("binary", "summarize", "cohort", "id", "sites", "result", os_limit=True)
        refused = probe.command("binary", "extract", "cohort", "id", "sites", "result", budget=1)
        self.assertIn("--require-os-limit", admitted)
        self.assertNotIn("--require-os-limit", refused)
        self.assertEqual(refused[refused.index("--memory-budget-mb") + 1], "1")
        self.assertIn("--enforce", refused)

    def test_pairs_use_only_an_explicit_readonly_table_and_preserve_direction(self):
        argv = probe.command("binary", "compare-pairs", "cohort", "id", "sites", "result",
                             os_limit=True, pairs="/inputs/pairs.tsv")
        self.assertEqual(argv[argv.index("--pairs") + 1], "/inputs/pairs.tsv")
        self.assertIn("--require-os-limit", argv)
        container = probe.container_arguments("unique", probe.DEFAULT_IMAGE, "/binary", "/saved-cohort",
            "/candidates", "/result", argv, pairs="/ordered-pairs.tsv")
        mounts = [container[index + 1] for index, value in enumerate(container) if value == "--mount"]
        self.assertIn("type=bind,src=/ordered-pairs.tsv,dst=/inputs/pairs.tsv,readonly", mounts)
        self.assertEqual(len(mounts), 5)
        for operation, pairs in [("compare-pairs", None), ("extract", "/table")]:
            with self.assertRaisesRegex(ValueError, "explicit pair table"):
                probe.command("binary", operation, "cohort", "id", "sites", "result", pairs=pairs)


if __name__ == "__main__":
    unittest.main()
