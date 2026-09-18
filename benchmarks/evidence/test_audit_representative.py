#!/usr/bin/env python3
"""Fast post-run audit regressions; no extraction, network, or representative data."""
import json
from pathlib import Path
import tempfile
import unittest

import audit_representative as auditor


class AuditTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.root = Path(self.directory.name)
        (self.root / "artifacts").mkdir()
        (self.root / "raw").mkdir()
        (self.root / "cache").mkdir()
        self.source = self.root / "reference.fa"
        self.binary = self.root / "rosalind"
        self.driver = self.root / "representative.py"
        for path, data in ((self.source, ">chr1\nA\n"), (self.binary, "binary"), (self.driver, "driver")):
            path.write_text(data)
        self.manifest = self.root / "workloads.json"
        manifest = dict(workloads=[dict(id="w", expected_rows=1, gates=dict(min_effective_tiles=1))],
                        cases=[dict(id="c", budget_mib=1)])
        self.manifest.write_text(json.dumps(manifest))
        artifact = self.root / "artifacts" / "w.tsv"
        artifact.write_text("#contig\tpos\nchr1\t1\n")
        verified = self.root / "raw" / "w.verify.json"
        verified.write_text('{"verified": true}\n')
        self.receipt_path = Path(str(artifact) + ".manifest.json")
        self.values = {"execution.emitted_loci": "1", "execution.microtile_bases": "128",
                       "execution.microtiles": "1", "execution.record_visits": "3",
                       "peak_rss_bytes": "500", "predicted_peak_rss_bytes": "900", "contract_verdict": "within"}
        self.receipt = dict(params=dict(memory_budget_mb="1", run_status="completed"), measurements=self.values)
        self.write_receipt()
        self.row = dict(label="w", kind="rosalind", case="c", workload="w", phase="fresh", repeat=1,
                        status="completed", exit_code=0, valid=True,
                        artifact=str(artifact.relative_to(self.root)), output_sha256=auditor.fingerprint(artifact)["sha256"],
                        output_bytes=artifact.stat().st_size, declared_budget_mib=1, peak_rss_bytes=600,
                        argv=[str(self.binary), "analyze", "evidence", "--manifest", str(self.receipt_path),
                              "--memory-budget-mb", "1"], receipt_measurements=self.values,
                        verification=dict(exit_code=0, stdout_sha256=auditor.fingerprint(verified)["sha256"], wall_seconds=0.2))
        self.report = dict(status="passed", measurements=[self.row], job_inventory=[dict(label="w")], manifest=manifest,
                           inputs=[dict(role="reference", path=str(self.source), sha256=auditor.fingerprint(self.source)["sha256"])],
                           environment=dict(binary_sha256=auditor.fingerprint(self.binary)["sha256"],
                               manifest_sha256=auditor.fingerprint(self.manifest)["sha256"],
                               harness_sha256={self.driver.name: auditor.fingerprint(self.driver)["sha256"]}))
        self.report_path = self.root / "report.json"

    def tearDown(self):
        self.directory.cleanup()

    def write_receipt(self):
        self.receipt_path.write_text(json.dumps(self.receipt))

    def audit(self):
        self.report_path.write_text(json.dumps(self.report))
        before = self.report_path.read_bytes()
        result = auditor.audit(self.report_path, self.manifest, harness_dir=self.root)
        self.assertEqual(self.report_path.read_bytes(), before)
        return result

    def test_identity_resource_work_and_storage_audit(self):
        (self.root / "cache" / "evidence.arrow").write_bytes(b"arrow")
        result = self.audit()
        self.assertEqual(result["status"], "passed", result["issues"])
        self.assertEqual(len(result["identities"]), 4)
        self.assertEqual(result["execution"][0]["observations"][0]["alignment_record_visits"], 3)
        self.assertEqual(result["retained_storage"][0]["logical_bytes"], 5)
        self.assertEqual(result["verification_summaries"][0]["wall_seconds"]["median"], 0.2)

    def test_changed_startup_files_fail_individually(self):
        for path in (self.source, self.binary, self.driver, self.manifest):
            with self.subTest(path=path.name):
                previous = path.read_bytes()
                path.write_bytes(previous + b"changed")
                result = self.audit()
                self.assertEqual(result["status"], "failed")
                self.assertTrue(any("startup identity mismatch" in issue for issue in result["issues"]))
                path.write_bytes(previous)

    def test_prediction_underestimate_is_not_a_budget_breach(self):
        self.values["predicted_peak_rss_bytes"] = "100"
        self.write_receipt()
        result = self.audit()
        self.assertEqual(result["status"], "passed", result["issues"])
        self.assertEqual(len(result["prediction_underestimates"]), 1)
        self.assertFalse(result["prediction_underestimates"][0]["budget_breach"])
        self.row["peak_rss_bytes"] = 1048577
        result = self.audit()
        self.assertEqual(result["status"], "failed")
        self.assertTrue(result["measurements"][0]["budget_breach"])

    def test_changed_receipt_budget_and_observed_denominator_fail(self):
        self.receipt["params"]["memory_budget_mb"] = "2"
        self.write_receipt()
        self.report["manifest"]["workloads"][0]["expected_rows"] = 2
        result = self.audit()
        self.assertEqual(result["status"], "failed")
        self.assertIn("denominator", str(result["issues"]))
        self.report["manifest"]["workloads"][0]["expected_rows"] = 1
        result = self.audit()
        self.assertIn("receipt budget differs", str(result["issues"]))

    def test_failed_and_malformed_receipt_rows_are_retained(self):
        self.row["status"], self.row["exit_code"], self.row["valid"] = "refused", 3, False
        result = self.audit()
        self.assertEqual(result["status"], "failed")
        self.assertEqual(result["measurements"][0]["status"], "refused")
        self.row["status"], self.row["exit_code"], self.row["valid"] = "completed", 0, True
        self.receipt_path.write_text("bad json")
        result = self.audit()
        self.assertEqual(result["status"], "failed")
        self.assertEqual(len(result["measurements"]), 1)
        self.assertTrue(result["measurements"][0]["issues"])

    def test_running_and_missing_invocations_fail(self):
        self.report["status"] = "running"
        self.assertIn("still running", str(self.audit()["issues"]))
        self.report["status"] = "passed"
        self.report["job_inventory"].append(dict(label="unrecorded"))
        self.assertEqual(self.audit()["inventory"]["missing"], ["unrecorded"])

    def test_plan_widths_alone_do_not_prove_execution_changed(self):
        self.report["manifest"]["workloads"][0]["gates"]["min_effective_tiles"] = 2
        result = self.audit()
        self.assertEqual(result["status"], "failed")
        self.assertIn("no observed work difference", str(result["issues"]))

    def test_bounded_metadata_and_no_storage_symlink_traversal(self):
        self.source.write_text("x" * 100)
        with self.assertRaisesRegex(ValueError, "exceeds"):
            auditor.read_json(self.source, 50)
        (self.root / "cache" / "outside").symlink_to(self.source)
        with self.assertRaisesRegex(ValueError, "symlink"):
            auditor.storage_inventory(self.root / "cache")


if __name__ == "__main__":
    unittest.main()
