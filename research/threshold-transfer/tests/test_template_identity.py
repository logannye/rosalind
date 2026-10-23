"""Source-local template audit boundaries; execute on hosted validation only."""
from __future__ import annotations

import copy
from pathlib import Path
import sys
import unittest

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
from template_identity import TemplateIdentityAudit


class Record:
    def __init__(self, *, name="HISEQ:22:FC1:1:1101:100:200", first=True,
                 secondary=False, supplementary=False):
        self.query_name = name
        self.flag = 99 if first else 147
        if secondary:
            self.flag |= 0x100
        if supplementary:
            self.flag |= 0x800
        self.reference_id = 0
        self.reference_start = 100 if first else 200
        self.reference_end = self.reference_start + 50
        self.mapping_quality = 60
        self.next_reference_id = 0
        self.next_reference_start = 200 if first else 100
        self.template_length = 150 if first else -150
        self.cigartuples = [(0, 50)]
        self.query_sequence = "A" * 50
        self.query_qualities = [30] * 50
        self.tags = {}

    def has_tag(self, tag):
        return tag in self.tags


def audit(**kwargs):
    return TemplateIdentityAudit(source_identity="pinned-source-sha256", synthetic_group="technical-source-only", **kwargs)


class TemplateIdentityTests(unittest.TestCase):
    def test_reciprocal_primary_mates_are_preserved(self):
        value = audit()
        first, second = Record(), Record(first=False)
        original = copy.deepcopy((first.__dict__, second.__dict__))
        value.observe(first)
        value.observe(second)
        receipt = value.finish()
        self.assertEqual(receipt["physical_records_observed"], 2)
        self.assertEqual(receipt["mapped_primary_pairs_reciprocity_checked"], 1)
        self.assertEqual(receipt["distinct_exact_qnames"], 1)
        self.assertEqual(receipt["physical_records_deduplicated"], 0)
        self.assertEqual((first.__dict__, second.__dict__), original)

    def test_identical_primary_end_signatures_are_counted_without_deduplication(self):
        value = audit()
        for record in (Record(), Record(), Record(first=False), Record(first=False)):
            value.observe(record)
        receipt = value.finish()
        self.assertEqual(receipt["physical_records_observed"], 4)
        self.assertEqual(receipt["primary_records"], 4)
        self.assertEqual(receipt["repeated_primary_end_signatures"], 2)
        self.assertEqual(receipt["unique_primary_ends"], 2)
        self.assertEqual(receipt["templates_with_both_primary_mates"], 1)

    def test_optional_auxiliary_tag_variation_is_not_identity_conflict(self):
        first, repeat = Record(), Record()
        first.tags = {"NM": 1, "AS": 45}
        repeat.tags = {"AS": 45, "NM": 2}
        value = audit()
        value.observe(first)
        value.observe(repeat)
        self.assertEqual(value.finish()["repeated_primary_end_signatures"], 1)

    def test_primary_core_field_conflicts_fail_closed(self):
        changes = {"reference_id": 1, "reference_start": 101, "reference_end": 151,
                   "mapping_quality": 59, "next_reference_id": 1, "next_reference_start": 201,
                   "template_length": 151, "cigartuples": [(0, 49), (1, 1)],
                   "query_sequence": "C" * 50, "query_qualities": [31] * 50, "flag": 97}
        for field, changed in changes.items():
            with self.subTest(field=field):
                value = audit()
                value.observe(Record())
                second = Record()
                setattr(second, field, changed)
                with self.assertRaisesRegex(ValueError, "conflicting primary-end identity"):
                    value.observe(second)
                with self.assertRaises(RuntimeError):
                    value.finish()

    def test_secondary_and_supplementary_same_qname_allow_distinct_payloads(self):
        value = audit()
        value.observe(Record())
        secondary = Record(secondary=True)
        secondary.reference_start = 1000
        secondary.query_sequence = "C" * 20
        supplementary = Record(supplementary=True)
        supplementary.cigartuples = [(5, 40), (0, 10)]
        supplementary.query_sequence = "A" * 10
        supplementary.query_qualities = [30] * 10
        value.observe(secondary)
        value.observe(supplementary)
        value.observe(Record(first=False))
        receipt = value.finish()
        self.assertEqual(receipt["physical_records_observed"], 4)
        self.assertEqual(receipt["nonprimary_records"], 2)
        self.assertEqual(receipt["secondary_records"], 1)
        self.assertEqual(receipt["supplementary_records"], 1)
        self.assertEqual(receipt["unique_primary_ends"], 2)

    def test_single_primary_mate_is_explicitly_not_observed_in_scope(self):
        value = audit()
        value.observe(Record())
        receipt = value.finish()
        self.assertEqual(receipt["templates_with_one_primary_mate_not_observed_in_scope"], 1)
        self.assertEqual(receipt["templates_with_both_primary_mates"], 0)

    def test_secondary_only_template_is_counted_without_claiming_mate_absence(self):
        value = audit()
        value.observe(Record(secondary=True))
        receipt = value.finish()
        self.assertEqual(receipt["templates_with_no_primary_alignment_observed_in_scope"], 1)
        self.assertEqual(receipt["templates_with_one_primary_mate_not_observed_in_scope"], 0)

    def test_mapped_pair_reference_start_and_strand_must_be_reciprocal(self):
        for field, changed in (("next_reference_id", 1), ("next_reference_start", 99),
                               ("reference_start", 201), ("flag", 131)):
            with self.subTest(field=field):
                value = audit()
                value.observe(Record())
                second = Record(first=False)
                setattr(second, field, changed)
                with self.assertRaisesRegex(ValueError, "nonreciprocal|strand"):
                    value.observe(second)

    def test_unmapped_mate_flags_must_match_observed_primary_mate(self):
        value = audit()
        value.observe(Record())
        second = Record(first=False)
        second.flag |= 0x4
        with self.assertRaisesRegex(ValueError, "unmapped flags disagree"):
            value.observe(second)

    def test_proper_pair_and_template_length_symmetry_are_not_required(self):
        first, second = Record(), Record(first=False)
        first.flag &= ~0x2
        second.flag &= ~0x2
        first.template_length, second.template_length = 150, -151
        value = audit()
        value.observe(second)
        value.observe(first)
        self.assertEqual(value.finish()["mapped_primary_pairs_reciprocity_checked"], 1)

    def test_primary_read_end_flags_and_pair_status_are_required(self):
        for flag in (0, 1, 1 | 0x40 | 0x80, 0x40):
            with self.subTest(flag=flag):
                record = Record()
                record.flag = flag
                with self.assertRaisesRegex(ValueError, "paired status"):
                    audit().observe(record)

    def test_source_rg_tag_cannot_be_overwritten_or_ignored(self):
        record = Record()
        record.tags = {"RG": "undeclared"}
        with self.assertRaisesRegex(ValueError, "original RG tag"):
            audit().observe(record)

    def test_full_qname_is_not_truncated_or_normalized(self):
        value = audit()
        value.observe(Record(name="HISEQ:22:FC1:1:1101:100:200"))
        value.observe(Record(name="HISEQ:23:FC1:1:1101:100:200"))
        self.assertEqual(value.finish()["distinct_exact_qnames"], 2)
        with self.assertRaises(ValueError):
            audit().observe(Record(name="HISEQ:22:FC1:1:1101:100:200/1"))

    def test_qname_validator_cannot_silently_rewrite_names(self):
        value = audit(validate_qname=lambda name: name.lower())
        with self.assertRaisesRegex(ValueError, "preserve"):
            value.observe(Record())

    def test_memory_entry_limit_fails_without_pruning(self):
        value = audit(max_templates=1)
        value.observe(Record())
        with self.assertRaisesRegex(ValueError, "max_templates exceeded"):
            value.observe(Record(name="HISEQ:22:FC1:1:1101:101:200"))
        with self.assertRaises(RuntimeError):
            value.finish()

    def test_finish_is_idempotent_and_closes_observation(self):
        value = audit()
        value.observe(Record())
        first = value.finish()
        self.assertEqual(first, value.finish())
        first["passed"] = False
        self.assertTrue(value.finish()["passed"])
        with self.assertRaises(RuntimeError):
            value.observe(Record())

    def test_empty_audit_cannot_pass(self):
        with self.assertRaisesRegex(ValueError, "empty"):
            audit().finish()

    def test_signature_does_not_retain_sequence_or_quality_payloads(self):
        value = audit()
        record = Record()
        value.observe(record)
        state = value._templates[record.query_name]
        self.assertEqual(len(state.first.signature), 32)
        self.assertFalse(hasattr(state.first, "query_sequence"))
        self.assertFalse(hasattr(state.first, "query_qualities"))


if __name__ == "__main__":
    unittest.main()
