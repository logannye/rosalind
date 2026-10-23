"""Hosted-only integration for source-only technical read-group derivation.

The tiny header is authorized by patching its identity validator only. Real
source pinning is independently tested by test_preflight; no real reads are
downloaded and no scientific model is fitted by this fixture.
"""
from __future__ import annotations

import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import preflight
import prepare_fresh as fresh


QNAME = "HISEQ:22:FC1:1:1101:100:200"


def _core(record):
    """All SAM core fields plus unchanged original optional tags, excluding RG."""
    return {
        "qname": record.query_name, "flag": record.flag,
        "reference_id": record.reference_id, "reference_start": record.reference_start,
        "reference_end": record.reference_end, "mapq": record.mapping_quality,
        "cigar": record.cigartuples, "mate_reference_id": record.next_reference_id,
        "mate_start": record.next_reference_start, "tlen": record.template_length,
        "sequence": record.query_sequence,
        "qualities": list(record.query_qualities) if record.query_qualities is not None else None,
        "original_optional_tags": [(key, value, kind) for key, value, kind in record.get_tags(with_value_type=True)
                                   if key != "RG"],
    }


@unittest.skipUnless(importlib.util.find_spec("pysam"), "Run native pysam integration on hosted validation")
class SourceOnlyWorkerTests(unittest.TestCase):
    def fixture(self, root, *, conflicting_repeat=False, reciprocal=True):
        import pysam
        inputs = root / "inputs"
        inputs.mkdir()
        reference = root / "reference.fa"
        reference.write_text(">chr20\nACGTACGTACGT\n")
        pysam.faidx(str(reference))
        header = {"HD": {"VN": "1.6", "SO": "unknown"},
                  "SQ": [{"SN": "chr20", "LN": 12}],
                  "CO": ["Original no-RG fixture metadata must remain present"]}
        source = root / "fixture.bam"
        with pysam.AlignmentFile(str(source), "wb", header=header) as writer:
            for index, first in enumerate((True, True, False)):
                record = pysam.AlignedSegment()
                record.query_name = QNAME
                record.query_sequence = "TCGT" if conflicting_repeat and index == 1 else "ACGT"
                record.flag = 99 if first else 147
                record.reference_id = 0
                record.reference_start = 0 if first else 4
                record.mapping_quality = 60
                record.cigartuples = [(0, 4)]
                record.next_reference_id = 0
                record.next_reference_start = 4 if first else (0 if reciprocal else 1)
                record.template_length = 8 if first else -8
                record.query_qualities = pysam.qualitystring_to_array("IIII")
                record.set_tag("NM", 0)
                writer.write(record)
        pysam.index(str(source))
        with pysam.AlignmentFile(str(source), "rb") as handle:
            original_header = handle.header.to_dict()
            original_header_text = str(handle.header)
            before_records = [_core(record) for record in handle.fetch("chr20")]
        source_spec = copy.deepcopy(fresh.DEFAULT_SOURCES["reads_bam"])
        bound = {"source_header_sha256": hashlib.sha256(original_header_text.encode()).hexdigest(),
                 "source_identity_sha256": preflight.source_identity_sha256(source_spec),
                 "dictionary_sha256": hashlib.sha256(b'[["chr20",12]]').hexdigest(),
                 "technical_read_group": preflight.technical_read_group_id(source_spec),
                 "header_read_groups": [], "contigs": 1, "chr20_length": 12}
        receipt = dict(bound, accepted=True)
        request = {"url": str(source), "source": source_spec, "bai": str(source) + ".bai",
                   "reference": str(reference), "out": str(inputs), "audit_out": str(root),
                   "seed": 42, "preflight_receipt": receipt,
                   "intervals": [{"contig": "chr20", "start": 0, "end": 12}],
                   "conditions": [{"name": "original", "numerator": 35, "denominator": 300},
                                  {"name": "downsampled", "numerator": 10, "denominator": 300},
                                  {"name": "novel_low", "numerator": 5, "denominator": 300}]}
        request_path = root / "request.json"
        request_path.write_text(json.dumps(request))
        return {"inputs": inputs, "source": source, "request_path": request_path,
                "request": request, "bound": bound, "header": original_header,
                "records": before_records, "sha256": fresh.pilot.digest(source)}

    def test_source_only_worker_preserves_every_physical_record_and_hash_input(self):
        import pysam
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fixture = self.fixture(root)
            request, bound = fixture["request"], fixture["bound"]
            with patch.object(preflight, "validate_header_identity", return_value=bound) as header_check, \
                    patch.object(fresh.prepare, "_keep_template", return_value=True) as keep:
                fresh._extract_bams_worker(fixture["request_path"])
            self.assertEqual(header_check.call_count, 1)
            self.assertEqual(header_check.call_args.args[1], fresh.DEFAULT_SOURCES["reads_bam"])
            self.assertEqual(fresh.pilot.digest(fixture["source"]), fixture["sha256"])
            identity = hashlib.sha256(fresh.pilot.canonical(request["source"]).encode()).hexdigest()
            expected_hash_arguments = [(42, identity, bound["technical_read_group"], QNAME,
                                       condition["numerator"], condition["denominator"])
                                      for _ in fixture["records"] for condition in request["conditions"]]
            self.assertEqual([call.args for call in keep.call_args_list], expected_hash_arguments)
            expected_header = copy.deepcopy(fixture["header"])
            expected_header["HD"]["SO"] = "coordinate"
            expected_header["RG"] = [{"ID": bound["technical_read_group"], "SM": "HG001",
                "DS": "Synthetic source-object group; original library and read-group identity unavailable"}]
            for condition in request["conditions"]:
                path = fixture["inputs"] / (condition["name"] + ".bam")
                self.assertTrue(Path(str(path) + ".bai").is_file())
                with pysam.AlignmentFile(str(path), "rb") as derived:
                    # SAM section order may be normalized by htslib (RG before
                    # CO). Compare all section contents without treating the
                    # top-level OrderedDict insertion order as biological data.
                    self.assertEqual(dict(derived.header.to_dict()), dict(expected_header))
                    rows = list(derived.fetch("chr20"))
                    self.assertEqual([_core(record) for record in rows], fixture["records"])
                    self.assertEqual([r.get_tag("RG") for r in rows], [bound["technical_read_group"]] * 3)
                    self.assertEqual([r.flag for r in rows], [99, 99, 147])
            with pysam.AlignmentFile(str(fixture["source"]), "rb") as original:
                self.assertEqual(original.header.to_dict(), fixture["header"])
                self.assertTrue(all(not record.has_tag("RG") for record in original.fetch("chr20")))
            result = json.loads((fixture["inputs"] / "bam-preparation.json").read_text())
            audit = json.loads((root / "template-identity-audit.json").read_text())
            admission = json.loads((root / "source-sample-admission.json").read_text())
            self.assertEqual(result["template_identity_audit"], audit)
            self.assertTrue(result["source_only_technical_group"])
            self.assertTrue(result["coordinate_order_verified"])
            self.assertTrue(result["full_ordered_bam_reference_dictionary_match"])
            self.assertEqual(result["selected_records"], 3)
            self.assertEqual(result["traversed_records"], 3)
            self.assertEqual(audit["physical_records_observed"], 3)
            self.assertEqual(audit["repeated_primary_end_signatures"], 1)
            self.assertEqual(audit["mapped_primary_pairs_reciprocity_checked"], 1)
            self.assertEqual(audit["physical_records_deduplicated"], 0)
            self.assertEqual(audit["qnames_rewritten"], 0)
            self.assertEqual(audit["technical_synthetic_group"], bound["technical_read_group"])
            self.assertTrue(admission["synthetic_source_group"])
            self.assertEqual(admission["original_read_groups"], [])
            for condition in request["conditions"]:
                self.assertEqual(result["conditions"][condition["name"]]["records"], 3)

    def test_conflicting_primary_repeat_rejects_before_preparation_receipt(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fixture = self.fixture(root, conflicting_repeat=True)
            with patch.object(preflight, "validate_header_identity", return_value=fixture["bound"]), \
                    patch.object(fresh.prepare, "_keep_template", return_value=True):
                with self.assertRaisesRegex(ValueError, "conflicting primary-end identity"):
                    fresh._extract_bams_worker(fixture["request_path"])
            self.assertFalse((fixture["inputs"] / "bam-preparation.json").exists())
            self.assertFalse((root / "template-identity-audit.json").exists())
            self.assertFalse(any(fixture["inputs"].glob("*.bai")))
            self.assertEqual(fresh.pilot.digest(fixture["source"]), fixture["sha256"])

    def test_nonreciprocal_primary_mates_reject_before_preparation_receipt(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fixture = self.fixture(root, reciprocal=False)
            with patch.object(preflight, "validate_header_identity", return_value=fixture["bound"]), \
                    patch.object(fresh.prepare, "_keep_template", return_value=True):
                with self.assertRaisesRegex(ValueError, "nonreciprocal"):
                    fresh._extract_bams_worker(fixture["request_path"])
            self.assertFalse((fixture["inputs"] / "bam-preparation.json").exists())
            self.assertFalse((root / "template-identity-audit.json").exists())
            self.assertEqual(fresh.pilot.digest(fixture["source"]), fixture["sha256"])

    def test_preflight_receipt_cannot_admit_a_different_header_identity(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fixture = self.fixture(root)
            request = fixture["request"]
            request["preflight_receipt"]["dictionary_sha256"] = "0" * 64
            fixture["request_path"].write_text(json.dumps(request))
            with patch.object(preflight, "validate_header_identity", return_value=fixture["bound"]):
                with self.assertRaisesRegex(ValueError, "accepted metadata preflight"):
                    fresh._extract_bams_worker(fixture["request_path"])
            self.assertFalse((fixture["inputs"] / "bam-preparation.json").exists())
            self.assertFalse(any(fixture["inputs"].glob("*.bam")))
            self.assertEqual(fresh.pilot.digest(fixture["source"]), fixture["sha256"])


if __name__ == "__main__":
    unittest.main()
