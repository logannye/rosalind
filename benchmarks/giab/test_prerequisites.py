"""Offline regressions for evaluator versioning and prepared read-only inputs."""
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from preflight import REQUIRED_ARTIFACTS, sha256, verify_prepared

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("happy_pin_version", HERE / "happy/pin-version.py")
PIN = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PIN)


class VersionTests(unittest.TestCase):
    def test_source_drift_is_rejected_before_replacement(self):
        with self.assertRaisesRegex(ValueError, "differs from the locked"):
            PIN.pinned_cmake(PIN.UPSTREAM_BLOCK)

    def test_one_guarded_assignment_populates_both_generated_versions(self):
        source = b"# fixture\n" + PIN.UPSTREAM_BLOCK + b"\nconfigure_file(native python)\n"
        with patch.object(PIN, "EXPECTED_CMAKE_SHA256", hashlib.sha256(source).hexdigest()):
            result = PIN.pinned_cmake(source)
            self.assertEqual(result.count(b'set(HAPLOTYPES_VERSION "0.3.15")'), 1)
            self.assertNotIn(b"git describe", result)
            self.assertIn(b"configure_file(native python)", result)
            with self.assertRaises(ValueError):
                PIN.pinned_cmake(result)
        duplicate = PIN.UPSTREAM_BLOCK + b"\n" + PIN.UPSTREAM_BLOCK
        with patch.object(PIN, "EXPECTED_CMAKE_SHA256", hashlib.sha256(duplicate).hexdigest()):
            with self.assertRaisesRegex(ValueError, "exactly one"):
                PIN.pinned_cmake(duplicate)

    def test_candidate_smoke_rejects_missing_python_or_wrong_native_version(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            docker = root / "docker"
            docker.write_text("#!/usr/bin/env python3\n" + """
import os, sys
args = sys.argv[1:]
if '/opt/hap.py/bin/xcmp' in args:
    print('xcmp version ' + ('0.3.14' if os.environ['BAD_VERSION'] == 'native' else '0.3.15'))
elif '/opt/rtg-tools/rtg' in args:
    print('RTG Tools 3.12.1')
else:
    print('Hap.py ' + ('' if os.environ['BAD_VERSION'] == 'python' else '0.3.15'))
""")
            docker.chmod(0o755)
            for mode, message in [("python", "hap.py Python version"), ("native", "hap.py C++ version")]:
                result = subprocess.run(
                    ["bash", str(HERE / "happy/smoke.sh"), "local:test", str(root / mode)],
                    env={**os.environ, "PATH": str(root) + os.pathsep + os.environ["PATH"], "BAD_VERSION": mode},
                    text=True, capture_output=True, check=False,
                )
                self.assertEqual(result.returncode, 2, result.stderr)
                self.assertIn(message, result.stderr)


class PreparedInputTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.prepared = self.root / "prepared"
        self.prepared.mkdir()
        for name in REQUIRED_ARTIFACTS:
            (self.prepared / name).write_bytes(b"fixture\n")
        (self.prepared / "GRCh38.chr20.fa").write_bytes(b">chr20\nAAAA\nAA\n")
        (self.prepared / "GRCh38.chr20.fa.fai").write_text("chr20\t6\t7\t4\t5\n")
        self.write_manifest()

    def write_manifest(self):
        artifacts = {
            path.name: {"bytes": path.stat().st_size, "sha256": sha256(path)}
            for path in self.prepared.iterdir()
        }
        (self.root / "data-manifest.json").write_text(json.dumps({"schema": 1, "prepared_artifacts": artifacts}))

    def test_complete_prepared_inputs_verify_without_writes(self):
        before = {path.name: path.read_bytes() for path in self.prepared.iterdir()}
        self.assertEqual(verify_prepared(self.root), len(REQUIRED_ARTIFACTS))
        self.assertEqual(before, {path.name: path.read_bytes() for path in self.prepared.iterdir()})

    def test_missing_fai_is_rejected_even_when_omitted_from_manifest(self):
        (self.prepared / "GRCh38.chr20.fa.fai").unlink()
        with self.assertRaisesRegex(ValueError, "missing prepared artifact"):
            verify_prepared(self.root)
        self.write_manifest()
        with self.assertRaisesRegex(ValueError, "missing required artifacts.*fa.fai"):
            verify_prepared(self.root)

    def test_changed_or_stale_fai_cannot_pass_preflight(self):
        (self.prepared / "GRCh38.chr20.fa.fai").write_text("chr20\t6\t99\t4\t5\n")
        with self.assertRaisesRegex(ValueError, "differs from data manifest"):
            verify_prepared(self.root)
        self.write_manifest()
        with self.assertRaisesRegex(ValueError, "FAI does not match"):
            verify_prepared(self.root)

    def test_other_prepared_artifact_changes_are_also_rejected(self):
        (self.prepared / "HG002.chr20.bam").write_bytes(b"changed alignment")
        with self.assertRaisesRegex(ValueError, "differs from data manifest.*bam"):
            verify_prepared(self.root)


class PreparationScriptTests(unittest.TestCase):
    def test_preparation_builds_and_manifests_the_extracted_reference_fai(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            harness, data, tools = root / "harness", root / "data", root / "tools"
            harness.mkdir()
            tools.mkdir()
            downloads = data / "downloads"
            downloads.mkdir(parents=True)
            for name in ["prepare.sh", "preflight.py"]:
                shutil.copyfile(HERE / name, harness / name)
            resources = []
            for line in (HERE / "resources.tsv").read_text().splitlines():
                if not line or line.startswith("#"):
                    continue
                identity, filename, _, url = line.split("\t")
                content = b"chr20\t0\t8\n"
                if filename.endswith(".vcf.gz"):
                    content = b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\nchr20\t1\t.\tA\tC\n"
                if filename.endswith(".gz"):
                    content = gzip.compress(content, mtime=0)
                (downloads / filename).write_bytes(content)
                resources.append("\t".join([identity, filename, hashlib.sha256(content).hexdigest(), url]))
            (harness / "resources.tsv").write_text("\n".join(resources) + "\n")
            samtools = tools / "samtools"
            samtools.write_text("#!" + sys.executable + "\n" + """
import os, sys
from pathlib import Path
args = sys.argv[1:]
with open(os.environ['SAMTOOLS_CALLS'], 'a') as log:
    log.write(repr(args) + '\\n')
if args[0] == 'quickcheck':
    pass
elif args[:2] == ['view', '-H']:
    print('@SQ\\tSN:chr20\\tLN:8')
elif args[0] == 'view':
    Path(args[args.index('-o') + 1]).write_bytes(b'prepared alignment')
elif args[0] == 'index':
    Path(args[1] + '.bai').write_bytes(b'prepared index')
elif args[0] == 'faidx' and len(args) == 3:
    print('>chr20\\nACGTACGT')
elif args[0] == 'faidx' and len(args) == 2:
    assert Path(args[1]).read_bytes() == b'>chr20\\nACGTACGT\\n'
    Path(args[1] + '.fai').write_text('chr20\\t8\\t7\\t8\\t9\\n')
elif args[0] == '--version':
    print('samtools fixture')
else:
    raise SystemExit('unexpected samtools invocation: ' + repr(args))
""")
            samtools.chmod(0o755)
            calls = root / "samtools.calls"
            result = subprocess.run(
                ["bash", str(harness / "prepare.sh"), str(data)],
                env={**os.environ, "PATH": str(tools) + os.pathsep + os.environ["PATH"], "SAMTOOLS_CALLS": str(calls)},
                text=True, capture_output=True, check=False,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertIn("['faidx', '{}']".format(data / "prepared/GRCh38.chr20.fa"), calls.read_text())
            manifest = json.loads((data / "data-manifest.json").read_text())
            self.assertIn("GRCh38.chr20.fa.fai", manifest["prepared_artifacts"])
            self.assertEqual(verify_prepared(data), len(REQUIRED_ARTIFACTS))


if __name__ == "__main__":
    unittest.main()
