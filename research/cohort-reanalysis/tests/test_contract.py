"""Pure Python checks; no native build, input decoding, or source access."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from contract import (ContractError, Member, comparison_profile, load_members,
                      normalize_thresholds, normalized_candidates,
                      source_scope_identity, validate_comparison)


def descriptor(source="a", sample="sample-a", fields=3):
    return {
        "version": 1, "semantics": "exact-read-summary-v1", "counting_unit": "read",
        "sampling": "none", "has_reference": True,
        "fields": fields, "schema_version": 1 if fields == 63 else 2,
        "fields_version": 2 if fields & 64 else 1,
        "profile": {"min_mapq": 20, "min_base_quality": 20,
                    "exclude_secondary": True, "exclude_supplementary": True,
                    "exclude_qc_fail": True, "exclude_duplicates": True},
        "sources": [{"role": "alignments", "path": f"/unavailable/{source}.bam",
                     "blake3": source * 64, "bytes": 123},
                    {"role": "reference", "path": "/unavailable/ref.fa",
                     "blake3": "f" * 64, "bytes": 101}],
        "contigs": [{"id": 0, "name": "chr2", "length": 100},
                    {"id": 1, "name": "chr1", "length": 200}],
        "sample_scope": {"version": 1, "mode": "named", "selected_sample": sample,
                         "declared_samples": [sample], "read_groups": [{"id": "rg", "sample": sample}],
                         "allows_unassigned": False},
    }


def member(name):
    return Member(name, Path(f"/unavailable/{name}/evidence-dataset.manifest.json"))


class ComparisonTests(unittest.TestCase):
    def test_sources_and_sample_names_may_differ_and_richer_masks_are_compatible(self):
        first, second = descriptor(), descriptor("b", "different", 127)
        contract = validate_comparison([(member("a"), first), (member("b"), second)])
        self.assertEqual(contract["required_field_mask"], 3)
        self.assertEqual(contract["effective_reference"]["representation_policy"],
                         "exact_source_bytes_no_equivalence")
        self.assertNotIn("sample_scope", contract)
        self.assertNotIn("sources", contract)
        self.assertNotIn("fields", contract)
        third = descriptor("c", "another", 63)
        self.assertEqual(contract, validate_comparison([(member("a"), first), (member("c"), third)]))

    def test_relocation_and_filename_suffix_do_not_establish_reference_format(self):
        first, second = descriptor(), descriptor("b", "b")
        second["sources"][1]["path"] = "/elsewhere/renamed.rref"
        self.assertEqual(comparison_profile(first), comparison_profile(second))
        # The identical content hash establishes byte/representation identity;
        # suffix changes alone are neither format detection nor sequence equivalence.
        second["sources"][1]["blake3"] = "e" * 64
        with self.assertRaisesRegex(ContractError, "incompatible"):
            validate_comparison([(member("a"), first), (member("b"), second)])

    def test_effective_reference_precedence_and_role_are_explicit(self):
        first, second = descriptor(), descriptor("b", "b")
        decoder_reference = copy.deepcopy(second["sources"][1])
        decoder_reference.update(role="cram-reference", blake3="d" * 64)
        second["sources"].append(decoder_reference)
        self.assertEqual(comparison_profile(first), comparison_profile(second))
        second["sources"] = [source for source in second["sources"] if source["role"] != "reference"]
        second["sources"][1]["blake3"] = "f" * 64
        self.assertEqual(comparison_profile(second)["effective_reference"]["role"], "cram-reference")
        with self.assertRaisesRegex(ContractError, "incompatible"):
            validate_comparison([(member("a"), first), (member("b"), second)])

    def test_duplicate_source_scope_alias_cannot_hide_behind_different_index(self):
        first, second = descriptor(), descriptor()
        second["sources"][0]["path"] = "/relocated/same-bytes.bam"
        second["sources"].append({"role": "alignment-index", "path": "/different.bai",
                                  "blake3": "b" * 64, "bytes": 9})
        self.assertEqual(source_scope_identity(first), source_scope_identity(second))
        with self.assertRaisesRegex(ContractError, "multiple specimen aliases"):
            validate_comparison([(member("a"), first), (member("b"), second)])
        different_scope = descriptor("a", "other-sample")
        validate_comparison([(member("a"), first), (member("b"), different_scope)])

    def test_rejects_incompatible_profiles_dictionaries_and_missing_fields(self):
        for mutate in (
            lambda d: d["profile"].update(min_mapq=19),
            lambda d: d["profile"].update(exclude_duplicates=False),
            lambda d: d["contigs"][0].update(name="different"),
            lambda d: d.update(sampling="sampled"),
            lambda d: d.update(counting_unit="fragment"),
            lambda d: d.update(has_reference=False),
            lambda d: d.update(fields=1),
            lambda d: d.update(schema_version=1),
            lambda d: d["sample_scope"].update(mode="pooled"),
            lambda d: d["sample_scope"].update(mode="unknown"),
        ):
            with self.subTest(mutate=mutate):
                other = descriptor("b", "b")
                mutate(other)
                with self.assertRaises(ContractError):
                    validate_comparison([(member("a"), descriptor()), (member("b"), other)])

    def test_no_guessed_filter_default(self):
        value = descriptor()
        del value["profile"]["min_mapq"]
        with self.assertRaisesRegex(ContractError, "explicit complete"):
            comparison_profile(value)


class MemberManifestTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root / "leaf.json").write_text("{}")
        self.path = self.root / "members.json"

    def write(self, members):
        self.path.write_text(json.dumps({"version": 1, "members": members}))

    def test_relative_paths_and_lexical_order(self):
        self.write([{"specimen_id": "z", "dataset": "leaf.json", "subject": "asserted"},
                    {"specimen_id": "a", "dataset": "leaf.json"}])
        members = load_members(self.path)
        self.assertEqual([entry.specimen_id for entry in members], ["a", "z"])
        self.assertEqual(members[0].dataset, (self.root / "leaf.json").resolve())
        self.assertEqual(members[1].metadata, {"subject": "asserted"})

    def test_input_limits_and_invalid_ids(self):
        self.write([{"specimen_id": "a", "dataset": "leaf.json"}])
        with self.assertRaisesRegex(ContractError, "max_manifest_bytes"):
            load_members(self.path, max_manifest_bytes=5)
        for invalid in ("", " ", "bad\nname", "a" * 257):
            self.write([{"specimen_id": invalid, "dataset": "leaf.json"}])
            with self.assertRaises(ContractError):
                load_members(self.path)
        self.write([{"specimen_id": "a", "dataset": "leaf.json"}] * 2)
        with self.assertRaisesRegex(ContractError, "duplicate specimen"):
            load_members(self.path)
        with self.assertRaisesRegex(ContractError, "1..1 members"):
            load_members(self.path, max_members=1)

    def test_duplicate_json_keys_and_implicit_discovery_are_rejected(self):
        self.path.write_text('{"version":1,"version":1,"members":[]}')
        with self.assertRaisesRegex(ContractError, "duplicate JSON key"):
            load_members(self.path)
        self.write([{"specimen_id": "a", "dataset": "."}])
        with self.assertRaisesRegex(ContractError, "leaf manifest file"):
            load_members(self.path)
        self.write([{"specimen_id": "a", "dataset": "leaf.json", "leaves": []}])
        with self.assertRaisesRegex(ContractError, "exactly one"):
            load_members(self.path)


class NativeCandidateTests(unittest.TestCase):
    def test_native_coordinates_alt_union_and_reference_dictionary_order(self):
        sites = [{"contig": 1, "position": 0, "reference": 65, "alternates": [84, 67]},
                 {"contig": 0, "position": 9, "reference": 65, "alternates": [71]},
                 {"contig": 1, "position": 0, "reference": 65, "alternates": [67]}]
        values = normalized_candidates({"sites": sites}, descriptor()["contigs"])
        self.assertEqual([entry.as_dict() for entry in values], [
            {"contig": "chr2", "pos": 10, "ref": "A", "alt": "G"},
            {"contig": "chr1", "pos": 1, "ref": "A", "alt": "C"},
            {"contig": "chr1", "pos": 1, "ref": "A", "alt": "T"},
        ])
        self.assertEqual(normalized_candidates({"sites": []}, descriptor()["contigs"]), [])

    def test_native_candidate_validation(self):
        base = {"contig": 0, "position": 0, "reference": 65, "alternates": [67]}
        for change in ({"alternates": [65]}, {"alternates": [78]}, {"reference": 78},
                       {"position": -1}, {"position": 100}, {"contig": 3}):
            with self.subTest(change=change), self.assertRaises(ContractError):
                normalized_candidates({"sites": [{**base, **change}]}, descriptor()["contigs"])
        with self.assertRaisesRegex(ContractError, "conflicting REF"):
            normalized_candidates({"sites": [base, {**base, "reference": 71}]}, descriptor()["contigs"])
        with self.assertRaisesRegex(ContractError, "ALT rows"):
            normalized_candidates({"sites": [{**base, "alternates": [67, 71]}]},
                                  descriptor()["contigs"], max_candidates=1)

    def test_thresholds_are_exact_positive_and_bounded(self):
        self.assertEqual(normalize_thresholds([20, 10, 5, 10]), (5, 10, 20))
        for values in ([], [0], [-1], [True], [10.0], [(1 << 64)], list(range(1, 10))):
            with self.subTest(values=values), self.assertRaises(ContractError):
                normalize_thresholds(values)


if __name__ == "__main__":
    unittest.main()
