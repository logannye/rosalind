"""Offline structural/link checks for the blank adoption kit (no participants)."""
import importlib.util
import json
from pathlib import Path
import re
import subprocess
import unittest

ROOT = Path(__file__).resolve().parents[1]
KIT = ROOT / "examples/adoption"
SPEC = importlib.util.spec_from_file_location("onboarding", ROOT / "scripts/onboarding.py")
ONBOARDING = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ONBOARDING)


class AdoptionKitTests(unittest.TestCase):
    def test_template_has_no_participant_results_or_return_claim(self):
        record = json.loads((KIT / "session.template.json").read_text())
        self.assertEqual(record["schema"], 1)
        self.assertEqual(record["record_kind"], "template")
        self.assertIsNone(record["task"])
        self.assertFalse(record["participant"]["consent_to_publish"])
        self.assertTrue(all(value is None for key, value in record["participant"].items()
                            if key != "consent_to_publish"))
        self.assertEqual(record["attempts"], [])
        self.assertEqual(record["assessment"]["state"], "not-started")
        self.assertIsNone(record["assessment"]["real_task_completed"])
        self.assertEqual(record["followup_30_days"]["status"], "not-scheduled")
        self.assertIsNone(record["followup_30_days"]["returned_to_real_work"])
        self.assertTrue(all(value is None for value in record["effort"].values()))

    def test_form_and_schema_match_without_extending_release_partner_fields(self):
        record = json.loads((KIT / "session.template.json").read_text())
        schema = json.loads((KIT / "session.schema.json").read_text())
        partner = json.loads((ROOT / "release/schemas/design-partner-v1.schema.json").read_text())

        def check_keys(value, shape):
            if isinstance(value, dict):
                self.assertFalse(shape["additionalProperties"])
                self.assertEqual(set(value), set(shape["required"]))
                self.assertEqual(set(value), set(shape["properties"]))
                for key, child in value.items():
                    check_keys(child, shape["properties"][key])
        check_keys(record, schema)
        self.assertNotIn("task", partner["properties"])
        self.assertNotIn("effort", partner["properties"])
        self.assertFalse(partner["additionalProperties"])
        self.assertEqual(set(schema["properties"]["task"]["enum"]), {
            None, "researcher-candidate-evidence", "builder-external-analyzer", "workflow-integration"})

    def test_kit_links_and_shell_blocks_are_valid(self):
        pending = [ROOT / "docs/adoption-validation.md", KIT / "README.md"]
        visited = set()
        while pending:
            document = pending.pop()
            if document in visited:
                continue
            visited.add(document)
            for target in ONBOARDING.local_links(document, ROOT):
                if target.is_dir() and (target / "README.md").is_file():
                    pending.append(target / "README.md")
                elif target.suffix.lower() == ".md":
                    pending.append(target)
        document = (ROOT / "docs/adoption-validation.md").read_text()
        for block in re.findall(r"```sh\n(.*?)\n```", document, re.S):
            result = subprocess.run(["bash", "-n"], input=block, text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
        for task in ("researcher", "builder", "workflow"):
            self.assertEqual(document.count(f"<!-- adoption:{task} -->"), 1)
        # The researcher joins text columns; dataset extract defaults to Arrow.
        self.assertRegex(document, r"dataset extract[\s\S]*?--format tsv")
        for relative in ("examples/research-filter/prepare.py", "examples/research-filter/join.py",
                         "integrations/nextflow/examples/evidence/main.nf",
                         "integrations/nextflow/examples/evidence/nextflow.test.config"):
            self.assertTrue((ROOT / relative).is_file(), relative)


if __name__ == "__main__":
    unittest.main()
