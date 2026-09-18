"""Native/Python cohort agreement, immutable extension and bounded lifecycle."""
import json
import os
import shutil
import subprocess
import tempfile
import unittest
from unittest import mock
from pathlib import Path

import pyarrow as pa
import pysam

from rosalind import EvidenceProcessError, materialize_evidence, open_cohort
from rosalind.cohort import _CommandTransport


class CohortTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.binary = Path(os.environ["ROSALIND_BINARY"]).resolve()
        # Source development tests can explicitly opt out while the surrounding
        # environment contains a prior wheel. Wheel validation leaves this false.
        self.mismatch = os.environ.get("ROSALIND_TEST_ALLOW_VERSION_MISMATCH") == "1"
        self.raw = self.root / "raw"
        self.raw.mkdir()
        self.reference = self.raw / "reference.fa"
        self.reference.write_text(">chr1\n" + "A" * 2200 + "\n")
        pysam.faidx(str(self.reference))
        self.regions = self.root / "targets.bed"
        self.regions.write_text("chr1\t0\t1100\n")
        self.manifests = {}
        for member, depth, alts in [("A", 12, 4), ("B", 8, 1)]:
            bam = self.raw / f"{member}.bam"
            header = {"HD": {"VN": "1.6", "SO": "coordinate"},
                      "SQ": [{"SN": "chr1", "LN": 2200}],
                      "RG": [{"ID": "rg", "SM": f"sample-{member}"}]}
            with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
                for position in [1, 1500]:
                    for n in range(depth):
                        read = pysam.AlignedSegment(out.header)
                        read.query_name = f"{member}-{position}-{n}"
                        read.query_sequence = "C" if n < alts else "A"
                        read.query_qualities = [35]
                        read.reference_id = 0
                        read.reference_start = position
                        read.mapping_quality = 60
                        read.flag = 0
                        read.cigarstring = "1M"
                        read.set_tag("RG", "rg")
                        out.write(read)
            pysam.index(str(bam))
            result = materialize_evidence(
                self.reference, bam, self.root / f"{member}.arrow",
                regions=self.regions, fields=["depths", "alleles"],
                cache_dir=self.root / f"cache-{member}", binary=self.binary,
                allow_version_mismatch=self.mismatch,
            )
            receipt = json.loads(result.manifest_path.read_text())
            self.manifests[member] = Path(receipt["measurements"]["execution.evidence_dataset_manifest"])
        table = self.root / "members.tsv"
        table.write_text("id\tmanifest\n" + "".join(f"{m}\t{p}\n" for m, p in self.manifests.items()))
        self.directory = self.root / "cohort"
        created = self.native("create", "--cohort", self.directory, "--members", table)
        self.snapshot = created["snapshot_id"]
        self.cohort = self.open(self.directory, self.snapshot)
        self.sites = self.write_sites("candidates.vcf", [1, 2])
        self.second = self.write_sites("second.vcf", [1, 1500, 1501])

    def native(self, *args):
        result = subprocess.run([str(self.binary), "cohort", *map(str, args)], capture_output=True, check=False)
        self.assertEqual(result.returncode, 0, result.stderr.decode())
        return json.loads(result.stdout)

    def open(self, directory, snapshot):
        return open_cohort(directory, snapshot, binary=self.binary,
                           allow_version_mismatch=self.mismatch)

    def write_sites(self, name, positions):
        path = self.root / name
        path.write_text("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=2200>\n"
                        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                        + "".join(f"chr1\t{p + 1}\t.\tA\tC\t.\tPASS\t.\n" for p in positions))
        return path

    def table(self, artifact):
        with pa.ipc.open_stream(artifact) as reader:
            return reader.read_all()

    def test_extract_matches_cli_and_preserves_null_zero_and_uint64(self):
        self.assertEqual(self.cohort.inspect()["snapshot_id"], self.snapshot)
        self.assertEqual(self.cohort.verify()["status"], "verified")
        plan = self.cohort.plan(sites=self.second)
        self.assertEqual(plan["status"], "blocked")
        self.assertEqual(plan["original_alignment_records_decoded"], 0)
        output = self.cohort.materialize(self.root / "partial.arrow", sites=self.second, missing="partial")
        native_path = self.root / "native.arrow"
        self.native("extract", "--cohort", self.directory, "--snapshot", self.snapshot,
                    "--sites", self.second, "--missing", "partial", "--format", "arrow-ipc", "-o", native_path)
        self.assertEqual(output.artifact_path.read_bytes(), native_path.read_bytes())
        table = self.table(output.artifact_path)
        self.assertEqual(table.schema.field("callable_depth").type, pa.uint64())
        self.assertEqual(table["callable_depth"].to_pylist(), [12, None, None, 8, None, None])
        observed = self.cohort.materialize(self.root / "observed.arrow", sites=self.sites)
        self.assertEqual(self.table(observed.artifact_path)["callable_depth"].to_pylist(), [12, 0, 8, 0])
        self.assertTrue(output.manifest_path.is_file())

    def test_summary_native_oracle_and_cli_byte_agreement(self):
        output = self.cohort.summarize(self.root / "summary.arrow", sites=self.sites)
        native_path = self.root / "summary-native.arrow"
        self.native("summarize", "--cohort", self.directory, "--snapshot", self.snapshot,
                    "--sites", self.sites, "--format", "arrow-ipc", "-o", native_path)
        self.assertEqual(output.artifact_path.read_bytes(), native_path.read_bytes())
        rows = self.table(output.artifact_path).to_pylist()
        self.assertEqual(rows[0]["n_requested"], 2)
        self.assertEqual(rows[0]["n_observed"], 2)
        self.assertEqual(rows[0]["n_depth_eligible"], 1)
        self.assertEqual(rows[0]["n_alt_supported"], 1)
        self.assertEqual(rows[0]["depth_eligible_support_fraction_numerator"], 1)
        self.assertEqual(rows[0]["depth_eligible_support_fraction_denominator"], 1)
        self.assertEqual(rows[1]["n_depth_eligible"], 0)
        self.assertIsNone(rows[1]["depth_eligible_support_fraction_denominator"])

    def test_extension_and_relocation_preserve_parent_and_original_sources_are_optional(self):
        sources = self.root / "sources.tsv"
        lines = ["id\trole\tpath\n"]
        for member, manifest in self.manifests.items():
            descriptor = json.loads((manifest.parent / "dataset.descriptor.json").read_text())
            lines += [f"{member}\t{s['role']}\t{s['path']}\n" for s in descriptor["sources"]
                      if s["role"] not in ("sites", "regions")]
        sources.write_text("".join(lines))
        plan = self.cohort.plan(operation="extend", sites=self.second, sources=sources, workdir=self.root)
        self.assertFalse(plan["raw_sources_opened"])
        self.assertTrue(plan["source_mapping_members_match"])
        extension = self.cohort.extend(sources, sites=self.second, workdir=self.root)
        self.assertTrue(extension.report["changed"])
        self.assertEqual([member["computed_loci"] for member in extension.report["members"]], [2, 2])
        self.assertEqual(self.cohort.snapshot_id, self.snapshot)
        shutil.rmtree(self.raw)
        moved = self.root / "relocated"
        shutil.move(self.directory, moved)
        child = self.open(moved, extension.cohort.snapshot_id)
        parent = self.open(moved, self.snapshot)
        output = child.materialize(self.root / "extended.arrow", sites=self.second)
        self.assertEqual(self.table(output.artifact_path)["callable_depth"].to_pylist(), [12, 12, 0, 8, 8, 0])
        self.assertEqual(parent.plan(sites=self.second)["status"], "blocked")
        child.verify()
        parent.verify()
        sources.write_text("id\trole\tpath\n")
        noop = child.extend(sources, sites=self.second, workdir=self.root / "unopened")
        self.assertFalse(noop.report["changed"])
        self.assertEqual(noop.cohort.snapshot_id, child.snapshot_id)

    def test_batches_are_lazy_bounded_single_use_and_empty_selection_is_distinct(self):
        sites = self.write_sites("many.vcf", range(1050))
        with self.cohort.batches(sites=sites, workdir=self.root / "batches") as run:
            self.assertIsNone(run._process)
            self.assertFalse(run.artifact_path.exists())
            sizes = [batch.num_rows for batch in run]
            self.assertEqual(sum(sizes), 2100)
            self.assertLessEqual(max(sizes), 1024)
            self.assertTrue(run.result.artifact_path.is_file())
            self.assertTrue(run.result.manifest_path.is_file())
        with self.assertRaises(RuntimeError):
            list(run)
        with self.cohort.batches(sites=sites, workdir=self.root / "early") as early:
            iterator = iter(early)
            next(iterator)
        iterator.close()
        self.assertIsNone(early.result)
        self.assertIsNotNone(early._process.poll())
        empty = self.cohort.materialize(self.root / "empty.arrow", sites=sites, members=[])
        self.assertEqual(self.table(empty.artifact_path).num_rows, 0)

    def test_native_refusal_preserves_existing_outputs_and_exposes_exit_category(self):
        output = self.cohort.materialize(self.root / "preserved.tsv", sites=self.sites, format="tsv")
        before = output.artifact_path.read_bytes()
        with self.assertRaises(EvidenceProcessError):
            self.cohort.materialize(output.artifact_path, sites=self.second, format="tsv", force=True)
        self.assertEqual(output.artifact_path.read_bytes(), before)
        with self.assertRaises(EvidenceProcessError):
            self.cohort.materialize(self.root / "fields.arrow", sites=self.sites,
                                    fields=["depths", "alleles", "strands"])
        self.assertFalse((self.root / "fields.arrow").exists())
        with self.assertRaises(EvidenceProcessError) as refused:
            self.cohort.materialize(self.root / "budget.arrow", sites=self.sites, memory_budget_mb=1)
        self.assertIn(refused.exception.returncode, (3, 4))
        self.assertFalse((self.root / "budget.arrow").exists())
        with self.cohort.batches(sites=self.second, workdir=self.root / "failed") as run:
            with self.assertRaises(EvidenceProcessError):
                list(run)
        self.assertIsNone(run.result)

    def test_large_query_argument_transport_preserves_native_results_and_cleans_up(self):
        scratch = self.root / "requests"
        scratch.mkdir()
        reference = self.cohort.materialize(self.root / "direct.arrow", sites=self.sites)
        with mock.patch.object(_CommandTransport, "argv_bytes", 1), \
                mock.patch("rosalind.cohort.tempfile.tempdir", str(scratch)):
            planned = self.cohort.plan(sites=self.sites, memory_budget_mb=512)
            self.assertEqual(planned["status"], "ready")
            transported = self.cohort.materialize(self.root / "transport.arrow", sites=self.sites)
            self.assertEqual(reference.artifact_path.read_bytes(), transported.artifact_path.read_bytes())
            self.cohort.summarize(self.root / "transport-summary.arrow", sites=self.sites)
            self.assertEqual(list(scratch.iterdir()), [])
            with self.cohort.batches(sites=self.sites, workdir=self.root / "transport-batches") as run:
                iterator = iter(run)
                next(iterator)
                self.assertEqual(len(list(scratch.glob("rosalind-cohort-request-*.json"))), 1)
            iterator.close()
            self.assertEqual(list(scratch.iterdir()), [])
            with self.assertRaises(EvidenceProcessError):
                self.cohort.materialize(self.root / "refused.arrow", sites=self.second)
            self.assertEqual(list(scratch.iterdir()), [])
        # A genuinely large member selection reaches native validation, not E2BIG.
        ids = [f"absent-{'x' * 100}-{index}" for index in range(1000)]
        with mock.patch("rosalind.cohort.tempfile.tempdir", str(scratch)):
            with self.assertRaises(EvidenceProcessError):
                self.cohort.plan(sites=self.sites, members=ids)
        self.assertEqual(list(scratch.iterdir()), [])


class CohortArgumentsTest(unittest.TestCase):
    def test_open_is_lazy_and_validates_only_local_argument_shapes(self):
        cohort = open_cohort("/nonexistent", "a" * 64, binary="/not-a-binary")
        self.assertEqual(cohort.snapshot_id, "a" * 64)
        for snapshot in ["", "A" * 64, "z" * 64, None]:
            with self.assertRaises(ValueError):
                open_cohort("unused", snapshot)
        for value in [0, -1, True, 1.5]:
            with self.assertRaises(ValueError):
                open_cohort("unused", "a" * 64, max_snapshot_bytes=value)
        with self.assertRaises(ValueError):
            cohort.plan(sites="unused", operation="unsupported")
        with self.assertRaises(ValueError):
            cohort.materialize("unused", sites="unused", format="parquet")
        with self.assertRaises(ValueError):
            cohort.materialize("unused", sites="unused", min_callable_depth=0)


if __name__ == "__main__":
    unittest.main()
