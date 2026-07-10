#!/usr/bin/env python3

import copy
import json
from pathlib import Path
import tempfile
import unittest

from compare import compare


def evidence(score: float = 1.0) -> dict:
    return {
        "calls_filter_all": {"precision": score},
        "calls_filter_pass": {"precision": score},
        "external_happy_vcfeval": {"container": "image@sha256:" + "a" * 64, "f1": score},
        "memory": {"predicted_peak_rss_bytes": 10, "peak_rss_bytes": 9},
        "receipt_claim": "b" * 64,
        "producer_identity": {"producer.version": "0.4.0"},
        "command_argv": ["rosalind", "variants"],
        "data_manifest": {"schema": 1},
    }


class ComparisonTests(unittest.TestCase):
    def run_compare(self, baseline: dict, latest: dict):
        temporary = tempfile.TemporaryDirectory()
        root = Path(temporary.name)
        baseline_path, latest_path = root / "baseline.json", root / "latest.json"
        baseline_path.write_text(json.dumps(baseline))
        latest_path.write_text(json.dumps(latest))
        before = baseline_path.read_bytes()
        status = compare(baseline_path, latest_path, root)
        self.assertEqual(before, baseline_path.read_bytes())
        return temporary, root, status

    def test_missing_baseline_emits_review_candidate_and_succeeds(self):
        temporary, root, status = self.run_compare(
            {"status": "not-yet-established"}, evidence()
        )
        self.addCleanup(temporary.cleanup)
        self.assertEqual(status, "candidate-pending-review")
        self.assertTrue((root / "baseline-candidate.json").is_file())

    def test_exact_evidence_reproduces_without_candidate(self):
        accepted = evidence()
        temporary, root, status = self.run_compare(accepted, copy.deepcopy(accepted))
        self.addCleanup(temporary.cleanup)
        self.assertEqual(status, "reproduced")
        self.assertFalse((root / "baseline-candidate.json").exists())

    def test_poor_or_changed_metrics_are_recorded_honestly(self):
        temporary, root, status = self.run_compare(evidence(), evidence(0.01))
        self.addCleanup(temporary.cleanup)
        self.assertEqual(status, "diverged")
        self.assertEqual(
            json.loads((root / "baseline-candidate.json").read_text())["calls_filter_all"]["precision"],
            0.01,
        )
        self.assertIn("diverged", (root / "credibility.md").read_text())

    def test_invalid_candidate_is_rejected(self):
        latest = evidence()
        del latest["receipt_claim"]
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            baseline, candidate = root / "baseline.json", root / "latest.json"
            baseline.write_text(json.dumps({"status": "not-yet-established"}))
            candidate.write_text(json.dumps(latest))
            with self.assertRaises(ValueError):
                compare(baseline, candidate, root)


if __name__ == "__main__":
    unittest.main()
