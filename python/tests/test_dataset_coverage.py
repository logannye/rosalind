"""Coverage wrapper contracts without a native build or optional Arrow runtime."""
import json
import unittest
from pathlib import Path
from unittest.mock import patch

from rosalind import open_dataset
from rosalind.dataset import _json


class DatasetCoverageWrapperTest(unittest.TestCase):
    def test_native_query_and_output_envelope_are_explicit(self):
        response = {"schema": "rosalind-dataset-coverage-v1", "covered_loci": 2**53 + 1}
        with patch("rosalind.dataset._resolve_binary", return_value=Path("/native/rosalind")), \
                patch("rosalind.dataset._json", return_value=response) as read:
            result = open_dataset("leaf/manifest.json").coverage(
                sites="sites.vcf", fields=["depths", "alleles"],
                memory_budget_mb=512, max_output_bytes=8192)
        self.assertEqual(result["covered_loci"], 2**53 + 1)
        read.assert_called_once_with([
            "/native/rosalind", "dataset", "coverage", "--dataset", "leaf/manifest.json",
            "--max-dataset-metadata-bytes", "33554432", "--memory-budget-mb", "512",
            "--sites", "sites.vcf", "--fields", "depths,alleles", "--max-output-bytes", "8192",
        ], max_bytes=8192)

    def test_output_limit_rejects_invalid_values_before_resolving_binary(self):
        with patch("rosalind.dataset._resolve_binary") as resolve:
            for limit in [True, False, 0, -1, 1.5, "1024", 33_554_433, None]:
                with self.subTest(limit=limit), self.assertRaises(ValueError):
                    open_dataset("manifest.json").coverage(max_output_bytes=limit)
        resolve.assert_not_called()

    def test_conflicting_selections_and_field_strings_are_rejected(self):
        with patch("rosalind.dataset._resolve_binary", return_value=Path("/native/rosalind")):
            with self.assertRaises(ValueError):
                open_dataset("manifest.json").coverage(regions="a.bed", sites="b.vcf")
            with self.assertRaises(TypeError):
                open_dataset("manifest.json").coverage(fields="depths")

    def test_unknown_schema_is_not_silently_accepted(self):
        with patch("rosalind.dataset._resolve_binary", return_value=Path("/native/rosalind")):
            for response in [{}, {"schema": "rosalind-dataset-coverage-v2"}, []]:
                with self.subTest(response=response), \
                        patch("rosalind.dataset._json", return_value=response), \
                        self.assertRaisesRegex(ValueError, "unsupported native"):
                    open_dataset("manifest.json").coverage()

    def test_json_output_is_bounded_before_parsing(self):
        def oversized(_command, stdout):
            stdout.write(b"not valid JSON, and larger than the admitted limit")
        with patch("rosalind.dataset._run", side_effect=oversized), \
                self.assertRaisesRegex(ValueError, "exceeds max_output_bytes"):
            _json(["unused"], max_bytes=4)

    def test_json_at_the_bound_preserves_integer_values(self):
        encoded = json.dumps({"n": 2**64 - 1}).encode()
        def completed(_command, stdout):
            stdout.write(encoded)
        with patch("rosalind.dataset._run", side_effect=completed):
            self.assertEqual(_json(["unused"], max_bytes=len(encoded)), {"n": 2**64 - 1})


if __name__ == "__main__":
    unittest.main()
