#!/usr/bin/env python3
"""Small oracle and harness regressions; no network or representative-data downloads."""
import copy
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import pysam

import representative as harness
import streaming_pysam as oracle


class Fixture:
    def __init__(self, root):
        self.root = Path(root)
        self.reference = self.root / "ref.fa"
        self.reference.write_text(">chr1\n" + "A" * 40000 + "\n")
        pysam.faidx(str(self.reference))
        self.bam = self.root / "reads.bam"
        header = {"HD": {"VN": "1.6", "SO": "coordinate"},
                  "SQ": [{"SN": "chr1", "LN": 40000}], "RG": [{"ID": "rg1", "SM": "S1"}]}
        with pysam.AlignmentFile(str(self.bam), "wb", header=header) as output:
            rows = [(0, "2S4M1I2M1D2M2N2M", "TTACGTCAACCGA", 16, 60, 30),
                    (0, "10M", "A" * 10, 0, 60, 30),
                    (0, "10M", "C" * 10, 1024, 60, 30),
                    (0, "10M", "G" * 10, 256, 60, 30),
                    (0, "10M", "T" * 10, 0, 10, 30),
                    (0, "10M", "A" * 10, 0, 255, 30),
                    (0, "10M", "A" * 10, 0, 60, None),
                    (0, "10M", "A" * 10, 0, 60, 10),
                    (0, "10M", "N" * 10, 0, 60, 30),
                    (16380, "20M", "C" * 20, 0, 40, 25),
                    (32760, "20M", "A" * 20, 16, 50, 35)]
            for number, (start, cigar, sequence, flag, mapq, quality) in enumerate(rows):
                read = pysam.AlignedSegment(output.header)
                read.query_name = f"r{number}"
                read.reference_id = 0
                read.reference_start = start
                read.flag, read.mapping_quality = flag, mapq
                read.cigarstring, read.query_sequence = cigar, sequence
                if quality is not None:
                    read.query_qualities = [quality] * len(sequence)
                read.set_tag("RG", "rg1")
                output.write(read)
        pysam.index(str(self.bam))
        self.cram = self.root / "reads.cram"
        with pysam.AlignmentFile(str(self.bam), "rb") as source, pysam.AlignmentFile(
                str(self.cram), "wc", header=source.header, reference_filename=str(self.reference)) as output:
            for read in source:
                output.write(read)
        pysam.index(str(self.cram))
        self.bed = self.root / "panel.bed"
        self.bed.write_text("chr1\t0\t20\tfirst\nchr1\t10\t30\toverlap\nchr1\t16380\t16400\nchr1\t32760\t32785\n")
        self.vcf = self.root / "sites.vcf"
        self.vcf.write_text("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=40000>\n"
                            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                            "chr1\t1\t.\tA\tG\t.\tPASS\t.\nchr1\t1\t.\tA\tC\t.\tPASS\t.\n"
                            "chr1\t16385\t.\tA\tC\t.\tPASS\t.\nchr1\t32780\t.\tA\tT\t.\tPASS\t.\n")

    def file(self, path):
        return dict(path=str(path), sha256=harness.sha256(path))

    def manifest(self):
        workloads = []
        for alignment, index in ((self.bam, Path(str(self.bam) + ".bai")),
                                 (self.cram, Path(str(self.cram) + ".crai"))):
            workloads.append(dict(id=alignment.suffix[1:] + "-panel", equivalence_group="panel",
                alignments=self.file(alignment), alignment_index=self.file(index),
                reference=self.file(self.reference), reference_fai=self.file(Path(str(self.reference) + ".fai")),
                selection=dict(kind="regions", **self.file(self.bed)), sample="S1", expected_rows=75,
                gates=dict(min_admitted_budgets=1, min_effective_tiles=1, workers=[1])))
        return dict(schema=1, label="synthetic-regression-only", repeats=3, seed=7,
                    cases=[dict(id="base", budget_mib=128, tile_bases=1024, workers=1, cache="none"),
                           dict(id="cache", budget_mib=256, tile_bases=16384, workers=2, cache="cold-resume-dataset")],
                    workloads=workloads)

    def rows(self, alignment=None, kind="regions", width=1024):
        stats = {}
        with pysam.FastaFile(str(self.reference)) as fasta, pysam.AlignmentFile(
                str(alignment or self.bam), "r", reference_filename=str(self.reference)) as bam:
            loci = oracle.selection_loci(kind, self.bed if kind == "regions" else self.vcf,
                                         {"chr1": 0}, {"chr1": 40000})
            rows = list(oracle.evidence_rows(bam, fasta, loci, tile_bases=width, stats=stats))
        return rows, stats


class OracleTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.fixture = Fixture(self.directory.name)

    def tearDown(self):
        self.directory.cleanup()

    def test_cigar_filters_zero_denominators_and_tile_bounds(self):
        baseline, stats = self.fixture.rows(width=1)
        self.assertEqual(len(baseline), 75)
        self.assertEqual(stats["max_retained_loci"], 1)
        row = dict(zip(oracle.HEADER, baseline[0]))
        self.assertEqual((row["prefilter_depth"], row["aligned_depth"], row["callable_depth"]), (9, 5, 2))
        self.assertEqual((row["a"], row["a_fwd"], row["a_rev"]), (2, 1, 1))
        self.assertEqual(row["read_position_sum"], 10)
        self.assertEqual(row["base_quality_histogram"], "30:2")
        for name in ("secondary", "duplicate", "unavailable_mapq", "low_mapq", "unavailable_base_quality", "low_base_quality", "ambiguous_base"):
            self.assertEqual(row["filtered_" + name], 1)
        zero = dict(zip(oracle.HEADER, baseline[20]))
        self.assertEqual(zero["callable_depth"], 0)
        self.assertEqual(zero["requested_alts"], ".")
        for width in (7, 1024, 16384):
            observed, stats = self.fixture.rows(width=width)
            self.assertEqual(observed, baseline)
            self.assertLessEqual(stats["max_retained_loci"], width)
            self.assertLessEqual(stats["max_numeric_array_bytes"], width * 377 * 8)

    def test_bam_cram_and_multiallelic_snv_equality(self):
        for kind in ("sites", "regions"):
            bam, _ = self.fixture.rows(kind=kind)
            cram, _ = self.fixture.rows(self.fixture.cram, kind=kind, width=3)
            self.assertEqual(bam, cram)
        rows, _ = self.fixture.rows(kind="sites")
        self.assertEqual(len(rows), 3)
        self.assertEqual(rows[0][3], "C,G")

    def test_unsorted_selection_and_ambiguous_samples_fail_explicitly(self):
        self.fixture.bed.write_text("chr1\t20\t30\nchr1\t0\t10\n")
        with self.assertRaisesRegex(ValueError, "sorted"):
            self.fixture.rows()
        with self.assertRaisesRegex(ValueError, "ambiguous"):
            oracle.sample_filter({"RG": [{"ID": "a", "SM": "A"}, {"ID": "b", "SM": "B"}]})
        with self.assertRaisesRegex(ValueError, "1..16384"):
            list(oracle.selection_windows([], 0))

    @unittest.skipUnless(os.environ.get("ROSALIND_BENCH_BINARY"), "set ROSALIND_BENCH_BINARY for actual native comparison")
    def test_native_exact_oracle_both_formats_and_selections(self):
        binary = Path(os.environ["ROSALIND_BENCH_BINARY"]).resolve()
        for alignment in (self.fixture.bam, self.fixture.cram):
            for kind in ("sites", "regions"):
                output = self.fixture.root / f"{alignment.suffix[1:]}-{kind}.tsv"
                result = subprocess.run([str(binary), "analyze", "evidence", "--alignments", str(alignment),
                    "--reference", str(self.fixture.reference), "--" + kind,
                    str(self.fixture.vcf if kind == "sites" else self.fixture.bed),
                    "--fields", "all", "--memory-budget-mb", "128", "--output", str(output)],
                    capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                expected = io.StringIO()
                import csv
                writer = csv.writer(expected, delimiter="\t", lineterminator="\n")
                writer.writerow(oracle.HEADER)
                writer.writerows(self.fixture.rows(alignment, kind=kind)[0])
                self.assertEqual(output.read_text(), expected.getvalue())


class HarnessTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.fixture = Fixture(self.directory.name)
        self.path = self.fixture.root / "workload.json"
        self.path.write_text(json.dumps(self.fixture.manifest()))
        self.manifest = harness.load_manifest(self.path)

    def tearDown(self):
        self.directory.cleanup()

    def test_case_inventory_is_repeatable_and_keeps_resume_dependencies(self):
        groups = harness.job_groups(self.manifest)
        self.assertEqual(groups, harness.job_groups(self.manifest))
        self.assertEqual(sum(map(len, groups)), 30)
        for group in groups:
            if len(group) == 3:
                self.assertEqual([job["phase"] for job in group], ["cold", "resumed", "persisted"])
                self.assertEqual(group[0]["cache_key"], group[1]["cache_key"])
                self.assertEqual(group[2]["kind"], "dataset")
        labels = [job["label"] for group in groups for job in group]
        self.assertEqual(len(labels), len(set(labels)))

    def test_fingerprints_and_invalid_matrix_retain_failure_report(self):
        identities, seconds = harness.verify_inputs(self.manifest)
        self.assertEqual(len(identities), 10)
        self.assertGreaterEqual(seconds, 0)
        self.fixture.bed.write_text("chr1\t0\t1\n")
        root = self.fixture.root / "failed"
        report = harness.execute(self.manifest, Path(sys.executable), root, self.path)
        self.assertEqual(report["status"], "failed")
        self.assertIn("fingerprint mismatch", report["error"])
        self.assertTrue((root / "report.json").is_file())
        self.assertEqual(len(report["job_inventory"]), 30)
        changed = self.fixture.manifest()
        changed["repeats"] = 1
        self.path.write_text(json.dumps(changed))
        with self.assertRaisesRegex(ValueError, "three repeats"):
            harness.load_manifest(self.path)
        self.assertEqual(harness.load_manifest(self.path, smoke=True)["repeats"], 1)

    def test_effective_scheduling_and_refusals_are_not_hidden(self):
        workload = self.manifest["workloads"][0]
        workload["gates"] = dict(min_admitted_budgets=3, min_effective_tiles=2, workers=[1, 2, 8])
        measurements = [dict(workload=workload["id"], kind="rosalind", valid=True,
            declared_budget_mib=budget, effective=dict(workers=worker, tile_bases=tile, microtiles=worker, record_visits=worker * 10))
            for budget, worker, tile in ((64, 1, 128), (128, 2, 1024), (256, 8, 1024))]
        self.assertTrue(harness.assess_gates(self.manifest, measurements)[workload["id"]]["passed"])
        measurements[-1]["effective"]["workers"] = 2
        self.assertFalse(harness.assess_gates(self.manifest, measurements)[workload["id"]]["passed"])
        measurements[-1]["valid"] = False
        gate = harness.assess_gates(self.manifest, measurements)[workload["id"]]
        self.assertEqual(gate["observed_budgets_mib"], [64, 128])
        self.assertFalse(gate["checks"]["all_requested_runs_completed_verified_and_equal"])

    def fake_measure(self, mode=None):
        def measure(root, label, argv, stdout):
            content = "#contig\tpos\n" + "".join(f"chr1\t{position + 1}\n" for position in range(75))
            result = dict(argv=list(map(str, argv)), exit_code=0, wall_seconds=0.01, peak_rss_bytes=100,
                          stdout_bytes=0, stdout_sha256="retained-test-value")
            if "--stats" in argv:
                stdout.write_text(content)
                stats = Path(argv[argv.index("--stats") + 1])
                stats.write_text("broken JSON" if mode == "oracle-metadata" else json.dumps(
                    dict(selected_loci=75, record_visits=10, windows=1)))
            elif "verify" in argv:
                stdout.write_text("{}\n")
            else:
                stdout.write_text("")
                Path(argv[argv.index("--output") + 1]).write_text(content)
                receipt = Path(argv[argv.index("--manifest") + 1])
                resumed = "--resume" in argv
                budget = int(argv[argv.index("--memory-budget-mb") + 1])
                values = {"peak_rss_bytes": "100", "predicted_peak_rss_bytes": "90", "contract_verdict": "within",
                    "execution.emitted_loci": "75", "execution.microtile_bases": "1024",
                    "execution.microtiles": "0" if resumed else "1", "execution.record_visits": "0" if resumed else "10",
                    "execution.computed_partitions": "0" if resumed else "1", "execution.reused_partitions": "1" if resumed else "0",
                    "execution.alignment_record_visits": "0", "original_sources_rehashed": "false",
                    "execution.evidence_dataset_manifest": str(root / "portable.manifest.json")}
                if mode == "wrong-rows":
                    values["execution.emitted_loci"] = "74"
                params = dict(memory_budget_mb=str(budget + (1 if mode == "wrong-budget" else 0)), run_status="completed")
                receipt.write_text("broken JSON" if mode == "native-metadata" else json.dumps(dict(params=params, measurements=values)))
                if mode == "rss-breach":
                    result["peak_rss_bytes"] = budget * 1048576 + 1
            if mode == "input-mutation":
                self.fixture.bed.write_text("chr1\t0\t1\n")
            return result
        return measure

    def test_completed_process_metadata_failure_keeps_timings_and_all_jobs(self):
        for mode in ("oracle-metadata", "native-metadata"):
            with self.subTest(mode=mode):
                root = self.fixture.root / mode
                report = harness.execute(self.manifest, Path(sys.executable), root, self.path,
                                         measure=self.fake_measure(mode))
                self.assertEqual(len(report["measurements"]), len(report["job_inventory"]))
                failed = [row for row in report["measurements"] if row.get("status") == "metadata-or-harness-failed"]
                self.assertTrue(failed)
                self.assertTrue(all(row["exit_code"] == 0 and row["wall_seconds"] == 0.01 for row in failed))
                self.assertTrue(all("output_sha256" in row for row in failed))
                self.assertNotEqual(report["status"], "passed")

    def test_budget_rss_and_denominator_gates_reject_false_success(self):
        for mode in (None, "wrong-budget", "rss-breach", "wrong-rows"):
            with self.subTest(mode=mode):
                report = harness.execute(self.manifest, Path(sys.executable), self.fixture.root / str(mode), self.path,
                                         measure=self.fake_measure(mode))
                self.assertEqual(report["status"], "passed" if mode is None else "failed-or-refused", report.get("error"))
                native = [row for row in report["measurements"] if row["kind"] == "rosalind"]
                if mode is None:
                    self.assertTrue(all(row["prediction_validation"]["underestimated"] for row in native))
                    self.assertTrue(all(row["resource_validation"]["passed"] for row in native))
                if mode == "wrong-budget":
                    self.assertTrue(all(not row["resource_validation"]["checks"]["receipt_budget_matches"] for row in native))
                if mode == "rss-breach":
                    self.assertTrue(all(not row["resource_validation"]["checks"]["observed_rss_within_budget"] for row in native))

    def test_midrun_mutation_aborts_and_records_remaining_inventory(self):
        report = harness.execute(self.manifest, Path(sys.executable), self.fixture.root / "changed", self.path,
                                 measure=self.fake_measure("input-mutation"))
        self.assertEqual(report["status"], "failed")
        self.assertEqual(len(report["measurements"]), len(report["job_inventory"]))
        self.assertEqual(report["measurements"][0]["wall_seconds"], 0.01)
        self.assertTrue(all(row["status"] == "not-run" for row in report["measurements"][1:]))
        self.assertTrue(any(not row["valid"] for row in report["final_input_checks"]))

    def test_final_checks_rehash_even_when_snapshot_check_is_not_used(self):
        identities, _ = harness.verify_inputs(self.manifest)
        report = dict(inputs=identities)
        self.assertTrue(harness.final_identity_check(report))
        self.fixture.bed.write_text("chr1\t0\t1\n")
        self.assertFalse(harness.final_identity_check(report))
        self.assertTrue(any(row["sha256"] != row["expected_sha256"] for row in report["final_input_checks"]))

    def test_summaries_exclude_failed_measurements_but_retain_requested_count(self):
        rows = [dict(workload="a", kind="rosalind", case="base", phase="fresh", valid=True,
                     wall_seconds=number, peak_rss_bytes=100) for number in (1, 2, 3)]
        rows.append(dict(rows[0], valid=False, wall_seconds=1000))
        summary = harness.summaries(rows)[0]
        self.assertEqual(summary["requested_repeats"], 4)
        self.assertEqual(summary["completed_repeats"], 3)
        self.assertEqual(summary["wall_seconds"]["median"], 2)


if __name__ == "__main__":
    unittest.main()
