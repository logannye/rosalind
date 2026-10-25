"""Representation, portability and refusal tests; never builds or runs native code."""
from __future__ import annotations

import importlib.util
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import checked_parquet as checked

HAS_ARROW = importlib.util.find_spec("pyarrow") is not None


class InventoryBoundaryTests(unittest.TestCase):
    def test_paths_reject_escapes_and_symlinks(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "payload").write_bytes(b"trusted")
            self.assertEqual(checked._relative_file(root, "payload"), root / "payload")
            for name in ("/payload", "../payload", "x/../payload", "./payload", "x//payload",
                         "x\\payload", "C:payload", "", "x\0y"):
                with self.subTest(name=name), self.assertRaises(ValueError):
                    checked._relative_file(root, name)
            (root / "linked").symlink_to(root / "payload")
            with self.assertRaisesRegex(ValueError, "symlink"):
                checked._relative_file(root, "linked")

    def test_integer_and_locus_do_not_coerce_floats_bools_or_zero(self):
        for value in (True, 1.0, -1, 1 << 64, None):
            with self.subTest(value=value), self.assertRaises(ValueError):
                checked._integer(value)
        self.assertEqual(checked._integer((1 << 64) - 1), (1 << 64) - 1)
        for reference in (None, "", "AA", "N"):
            with self.subTest(reference=reference), self.assertRaises(ValueError):
                checked._locus({"contig": "chr20", "pos": 1, "ref": reference}, {"chr20": 0})
        with self.assertRaisesRegex(ValueError, "one-based"):
            checked._locus({"contig": "chr20", "pos": 0, "ref": "A"}, {"chr20": 0})

    def test_inventory_anchor_and_duplicate_json_keys_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "inventory.json").write_text('{"format": "a", "format": "b"}')
            with self.assertRaisesRegex(ValueError, "anchored digest"):
                checked._load_inventory(root, "0" * 64)
            with self.assertRaisesRegex(ValueError, "duplicate JSON"):
                checked._load_inventory(root)


@unittest.skipUnless(HAS_ARROW, "PyArrow required; hosted pinned-wheel gate must run these cases")
class CheckedParquetTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)

    @staticmethod
    def row(pos=1, ref="A", alt="C", depth=0):
        row = {"contig": "chr20", "pos": pos, "ref": ref, "requested_alts": alt}
        row.update({name: 0 for name in checked.SCALARS})
        row["prefilter_depth"] = row["aligned_depth"] = row["callable_depth"] = depth
        row[ref.lower()] = row[f"{ref.lower()}_fwd"] = depth
        row.update(base_quality_sum=30 * depth, mapping_quality_sum=60 * depth,
                   read_position_sum=10 * depth, read_length_sum=151 * depth)
        row.update({name: [0] * size for name, size in checked.HISTOGRAMS.items()})
        row["base_quality_histogram"][30] = depth
        row["mapping_quality_histogram"][60] = depth
        return row

    def arrow(self, rows, *, batch_rows=1024, schema=None):
        import pyarrow as pa
        source = self.root / "native.arrow"
        schema = schema or checked.legacy_schema()
        with source.open("wb") as handle, pa.ipc.new_stream(handle, schema) as writer:
            for start in range(0, len(rows), batch_rows):
                writer.write_batch(pa.RecordBatch.from_pylist(rows[start:start + batch_rows], schema=schema))
        return source

    def convert(self, rows, **kwargs):
        source = self.arrow(rows, **kwargs)
        destination = self.root / "cache"
        inventory = checked.convert_arrow(source, destination,
            semantics={"sample": "HG002", "filters": {"min_mapq": 20, "min_bq": 20},
                       "selection": "authored synthetic fixture"},
            contig_order=["chr20"], expected_rows=len(rows))
        return destination, inventory

    def test_all_uint64_fields_and_fixed_arrays_survive_complete_readback(self):
        import pyarrow.parquet as pq
        row = self.row(depth=(1 << 64) - 1)
        for offset, name in enumerate(checked.SCALARS):
            row[name] = (1 << 64) - 1 - offset
        row["base_quality_histogram"][17] = (1 << 63) + 1
        row["mapping_quality_histogram"][254] = (1 << 64) - 1
        cache, inventory = self.convert([row])
        actual = pq.ParquetFile(cache / "evidence.parquet").read()
        self.assertTrue(actual.schema.equals(checked.legacy_schema(), check_metadata=True))
        self.assertEqual(actual.to_pylist(), [row])
        self.assertTrue(checked.verify(cache)["full_schema_verified"])
        self.assertEqual(inventory["files"][0]["path"], "evidence.parquet")

    def test_mixed_alt_projection_zero_and_threshold_boundaries(self):
        rows = [self.row(1, depth=0), self.row(2, ref="C", alt="A", depth=5),
                self.row(3, ref="G", alt="A", depth=10), self.row(4, ref="T", alt="A", depth=20)]
        rows[1].update(c=2, c_fwd=2, g=3, g_fwd=1, g_rev=2)
        rows[2].update(g=3, g_fwd=3, t=7, t_fwd=5, t_rev=2)
        rows[3].update(t=9, t_fwd=9, c=11, c_fwd=8, c_rev=3)
        cache, _ = self.convert(rows)
        requests = [{"contig": "chr20", "pos": i + 1, "ref": row["ref"], "alt": alt}
                    for i, (row, alt) in enumerate(zip(rows, "CGTC"))]
        answers, metrics = checked.query(cache, requests, thresholds=(5, 10, 20),
            expected_inventory_sha256=checked.sha256(cache / "inventory.json"))
        self.assertEqual([row["alt_count"] for row in answers], [0, 3, 7, 11])
        self.assertEqual([row["alt_fwd"] for row in answers], [0, 1, 5, 8])
        self.assertEqual([row["alt_rev"] for row in answers], [0, 2, 2, 3])
        self.assertEqual([[row[f"depth_ge_{value}"] for value in (5, 10, 20)] for row in answers],
                         [[False, False, False], [True, False, False], [True, True, False], [True, True, True]])
        self.assertNotIn("a", metrics["projected_columns"])
        self.assertEqual(metrics["batch_readahead"], 1)
        self.assertGreater(metrics["verification_seconds"], 0)
        self.assertIn("schema check", metrics["verification_scope"])

    def test_missing_ref_mismatch_and_duplicate_requests_are_not_answers(self):
        cache, _ = self.convert([self.row()])
        request = {"contig": "chr20", "pos": 1, "ref": "A", "alt": "C"}
        for requests, message in [([{**request, "pos": 2}], "missing requested locus"),
                                  ([{**request, "ref": "G"}], "REF mismatch"),
                                  ([request, request], "duplicate requested locus")]:
            with self.subTest(message=message), self.assertRaisesRegex(ValueError, message):
                checked.query(cache, requests)

    def test_bounded_conversion_and_relocated_source_absent_replay(self):
        import pyarrow.parquet as pq
        cache, inventory = self.convert([self.row(pos) for pos in range(1, 1001)], batch_rows=231)
        metadata = pq.ParquetFile(cache / "evidence.parquet").metadata
        self.assertTrue(all(metadata.row_group(i).num_rows <= 1024 for i in range(metadata.num_row_groups)))
        self.assertEqual(len(list(cache.glob("*.parquet"))), 1)
        anchor = checked.sha256(cache / "inventory.json")
        relocated = self.root / "unrelated" / "saved"
        shutil.copytree(cache, relocated)
        shutil.rmtree(cache)
        (self.root / "native.arrow").unlink()
        result = checked.verify(relocated, expected_inventory_sha256=anchor)
        self.assertEqual(result["rows"], 1000)
        answers, _ = checked.query(relocated, [{"contig": "chr20", "pos": 1000, "ref": "A", "alt": "G"}],
                                   expected_inventory_sha256=anchor)
        self.assertEqual(answers[0]["callable_depth"], 0)
        self.assertEqual(inventory["rows"], 1000)

    def test_wrong_schema_duplicate_unsorted_and_excess_rows_refuse_inventory(self):
        import pyarrow as pa
        for rows in ([self.row(), self.row()], [self.row(2), self.row(1)]):
            source = self.arrow(rows)
            destination = self.root / f"bad-{rows[0]['pos']}"
            with self.assertRaisesRegex(ValueError, "sorted and unique"):
                checked.convert_arrow(source, destination, semantics={"sample": "fixture"},
                                      contig_order=["chr20"], expected_rows=2)
            self.assertFalse((destination / "inventory.json").exists())
        schema = checked.legacy_schema()
        wrong = pa.schema([pa.field(f.name, pa.int64(), nullable=False) if f.name == "a" else f
                           for f in schema], metadata=schema.metadata)
        source = self.arrow([self.row()], schema=wrong)
        with self.assertRaisesRegex(ValueError, "mask-63"):
            checked.convert_arrow(source, self.root / "signed", semantics={"sample": "fixture"},
                                  contig_order=["chr20"], expected_rows=1)
        source = self.arrow([self.row(1), self.row(2)])
        with self.assertRaisesRegex(ValueError, "row bound"):
            checked.convert_arrow(source, self.root / "too-many", semantics={"sample": "fixture"},
                                  contig_order=["chr20"], expected_rows=1)
        source = self.arrow([self.row(pos) for pos in range(1, 1026)], batch_rows=1025)
        with self.assertRaisesRegex(ValueError, "1024-row bound"):
            checked.convert_arrow(source, self.root / "huge-batch", semantics={"sample": "fixture"},
                                  contig_order=["chr20"], expected_rows=1000)

    def test_existing_output_and_consumed_corruption_are_refused(self):
        cache, _ = self.convert([self.row()])
        with self.assertRaises(FileExistsError):
            checked.convert_arrow(self.root / "native.arrow", cache, semantics={"sample": "fixture"},
                                  contig_order=["chr20"], expected_rows=1)
        payload = cache / "evidence.parquet"
        contents = bytearray(payload.read_bytes())
        contents[len(contents) // 2] ^= 1
        payload.write_bytes(contents)
        with self.assertRaisesRegex(ValueError, "payload changed"):
            checked.verify(cache)
        with self.assertRaisesRegex(ValueError, "payload changed"):
            checked.query(cache, [{"contig": "chr20", "pos": 1, "ref": "A", "alt": "C"}])

    def test_inventory_change_during_scan_is_caught_before_return(self):
        cache, _ = self.convert([self.row()])
        original = checked._check_after
        def changed(*args):
            path = cache / "inventory.json"
            path.write_bytes(path.read_bytes() + b" ")
            original(*args)
        with patch.object(checked, "_check_after", side_effect=changed):
            with self.assertRaisesRegex(ValueError, "inventory changed before publication"):
                checked.query(cache, [{"contig": "chr20", "pos": 1, "ref": "A", "alt": "C"}])

    def test_consumed_payload_change_after_scan_is_caught_before_return(self):
        cache, _ = self.convert([self.row()])
        original = checked._check_after
        def changed(*args):
            path = cache / "evidence.parquet"
            content = bytearray(path.read_bytes())
            content[len(content) // 2] ^= 1
            path.write_bytes(content)
            original(*args)
        with patch.object(checked, "_check_after", side_effect=changed):
            with self.assertRaisesRegex(ValueError, "payload changed"):
                checked.query(cache, [{"contig": "chr20", "pos": 1, "ref": "A", "alt": "C"}])


if __name__ == "__main__":
    unittest.main()
