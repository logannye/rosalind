import importlib.util
import json
import hashlib
import subprocess
import tempfile
import unittest
import zipfile
from pathlib import Path


def module(name):
    spec = importlib.util.spec_from_file_location(name.replace('-', '_'), Path(__file__).with_name(name + '.py'))
    loaded = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(loaded)
    return loaded


prepare = module('prepare-wheel-version').prepare
readiness = module('giab-readiness').readiness
verify = module('verify-wheel-upload').verify


class ReleaseHelpers(unittest.TestCase):
    def fixture(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        root = Path(temporary.name)
        (root / 'Cargo.toml').write_text('[package]\nname = "rosalind-bio"\nversion = "0.4.0"\n')
        (root / 'Cargo.lock').write_text('version = 3\n[[package]]\nname = "rosalind-bio"\nversion = "0.4.0"\n')
        for args in [('init', '-q'), ('config', 'user.email', 'test@example.invalid'), ('config', 'user.name', 'Test'), ('add', '.'), ('commit', '-qm', 'fixture')]:
            subprocess.run(['git', *args], cwd=root, check=True)
        sha = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=root, text=True).strip()
        return root, sha

    def test_stable_wheel_leaves_manifests_unchanged(self):
        root, sha = self.fixture()
        before = [(root / name).read_bytes() for name in ('Cargo.toml', 'Cargo.lock')]
        report = prepare(root, '0.4.0', sha)
        self.assertFalse(report['derived_manifest'])
        self.assertEqual(before, [(root / name).read_bytes() for name in ('Cargo.toml', 'Cargo.lock')])

    def test_rcs_have_distinct_native_and_pep440_versions(self):
        for number in (1, 2):
            root, sha = self.fixture()
            report = prepare(root, f'0.4.0-rc.{number}', sha)
            self.assertEqual(report['python_version'], f'0.4.0rc{number}')
            self.assertIn(f'0.4.0-rc.{number}', (root / 'Cargo.toml').read_text())
            self.assertIn(f'0.4.0-rc.{number}', (root / 'Cargo.lock').read_text())

    def test_mismatched_commit_or_version_refuses_without_edits(self):
        for version, candidate in [('0.5.0', None), ('0.4.0-rc.0', None), ('0.4.0', 'a' * 40)]:
            root, sha = self.fixture()
            before = (root / 'Cargo.toml').read_bytes()
            with self.assertRaises(ValueError):
                prepare(root, version, candidate or sha)
            self.assertEqual(before, (root / 'Cargo.toml').read_bytes())

    def test_missing_evaluator_is_blocked_without_scientific_result(self):
        report = readiness({'generated_image': {'digest': None}})
        self.assertEqual(report['status'], 'blocked')
        self.assertFalse(report['evaluation_performed'])
        self.assertIsNone(report['image'])

    def test_digest_and_source_commit_are_required_for_locked_image(self):
        image = {'repository': 'ghcr.io/example/evaluator', 'platform': 'linux/amd64', 'digest': 'sha256:' + 'a' * 64, 'built_from_commit': 'b' * 40}
        self.assertEqual(readiness({'generated_image': image})['status'], 'ready')
        for field in ('digest', 'built_from_commit'):
            self.assertEqual(readiness({'generated_image': {**image, field: None}})['status'], 'blocked')
        self.assertEqual(readiness({}, 'image:latest')['status'], 'blocked')

    def test_existing_wheel_is_skipped_only_when_bytes_match(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            wheel = root / 'rosalind_bio-0.4.0rc1-py3-none-test.whl'
            with zipfile.ZipFile(wheel, 'w') as archive:
                archive.writestr('rosalind_bio-0.4.0rc1.dist-info/METADATA', 'Name: rosalind-bio\nVersion: 0.4.0rc1\n')
            exact = {'urls': [{'filename': wheel.name, 'digests': {'sha256': hashlib.sha256(wheel.read_bytes()).hexdigest()}}]}
            self.assertTrue(verify(root, 'testpypi', lambda _: exact)[0]['existing_exact'])
            self.assertFalse(verify(root, 'testpypi', lambda _: None)[0]['existing_exact'])
            exact['urls'][0]['digests']['sha256'] = 'a' * 64
            with self.assertRaises(ValueError):
                verify(root, 'testpypi', lambda _: exact)


if __name__ == '__main__':
    unittest.main()
