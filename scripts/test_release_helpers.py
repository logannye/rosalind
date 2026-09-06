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
wheel_artifacts = module('wheel-artifacts')
release_policy = module('verify-release-automation')


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


class WheelArtifacts(unittest.TestCase):
    def fixture(self, version='0.5.0-rc.1'):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        root = Path(temporary.name)
        artifacts = root / 'artifacts'
        artifacts.mkdir()
        candidate = 'b' * 40
        for target, platform in [
            ('aarch64-apple-darwin', 'macosx_11_0_arm64'),
            ('x86_64-apple-darwin', 'macosx_10_12_x86_64'),
            ('x86_64-unknown-linux-gnu', 'manylinux_2_28_x86_64'),
        ]:
            directory = root / target
            directory.mkdir()
            pyversion = version.replace('-rc.', 'rc')
            wheel = directory / f'rosalind_bio-{pyversion}-py3-none-{platform}.whl'
            with zipfile.ZipFile(wheel, 'w') as archive:
                archive.writestr(f'rosalind_bio-{pyversion}.dist-info/METADATA', f'Name: rosalind-bio\nVersion: {pyversion}\n')
            build = directory / 'build.json'
            build.write_text(json.dumps({'source_commit': candidate, 'rust_version': version, 'python_version': pyversion}))
            report = wheel_artifacts.record(directory, build, target)
            (artifacts / wheel.name).write_bytes(wheel.read_bytes())
            (artifacts / f'wheel-build-{target}.json').write_text(json.dumps(report))
        return root, artifacts, candidate

    def test_complete_stable_and_rc_sets_stage_only_the_recorded_bytes(self):
        for version, index in [('0.5.0', 'pypi'), ('0.5.0-rc.1', 'testpypi')]:
            root, artifacts, candidate = self.fixture(version)
            output = root / 'dist'
            identities = wheel_artifacts.stage(artifacts, output, candidate, version, index)
            self.assertEqual(len(identities), 3)
            self.assertEqual(len(list(output.iterdir())), 3)
            for identity in identities:
                self.assertEqual((output / identity['filename']).read_bytes(), (artifacts / identity['filename']).read_bytes())

    def test_wrong_candidate_version_or_index_refuses_before_staging(self):
        for candidate, version, index in [('a' * 40, '0.5.0-rc.1', 'testpypi'), ('b' * 40, '0.5.0-rc.2', 'testpypi'), ('b' * 40, '0.5.0-rc.1', 'pypi')]:
            root, artifacts, _ = self.fixture()
            with self.assertRaises(ValueError):
                wheel_artifacts.stage(artifacts, root / 'dist', candidate, version, index)
            self.assertFalse((root / 'dist').exists())

    def test_missing_extra_or_tampered_artifacts_refuse_before_staging(self):
        for mutation in ('missing-report', 'extra-wheel', 'changed-wheel', 'wrong-target', 'path-escape'):
            root, artifacts, candidate = self.fixture()
            report_path = next(artifacts.glob('wheel-build-*.json'))
            report = json.loads(report_path.read_text())
            if mutation == 'missing-report':
                report_path.unlink()
            elif mutation == 'extra-wheel':
                (artifacts / 'unrecorded.whl').write_bytes(b'extra')
            elif mutation == 'changed-wheel':
                with (artifacts / report['wheel']['filename']).open('ab') as stream:
                    stream.write(b'tampered')
            else:
                if mutation == 'wrong-target':
                    report['target'] = 'x86_64-unknown-linux-gnu'
                    if report_path.name == 'wheel-build-x86_64-unknown-linux-gnu.json':
                        report['target'] = 'aarch64-apple-darwin'
                else:
                    report['wheel']['filename'] = '../escaped.whl'
                report_path.write_text(json.dumps(report))
            with self.assertRaises(ValueError):
                wheel_artifacts.stage(artifacts, root / 'dist', candidate, '0.5.0-rc.1', 'testpypi')
            self.assertFalse((root / 'dist').exists(), mutation)

    def test_wheel_metadata_and_platform_must_match_report(self):
        root, artifacts, _ = self.fixture()
        report = json.loads((artifacts / 'wheel-build-aarch64-apple-darwin.json').read_text())
        wheel = artifacts / report['wheel']['filename']
        with self.assertRaises(ValueError):
            wheel_artifacts.wheel_identity(wheel, 'x86_64-apple-darwin', '0.5.0-rc.1')
        with zipfile.ZipFile(wheel, 'w') as archive:
            archive.writestr('rosalind_bio-0.5.0rc1.dist-info/METADATA', 'Name: rosalind-bio\nVersion: 0.5.0rc2\n')
        with self.assertRaises(ValueError):
            wheel_artifacts.wheel_identity(wheel, 'aarch64-apple-darwin', '0.5.0-rc.1')


class PublisherBoundaries(unittest.TestCase):
    def workflows(self):
        root = Path(__file__).resolve().parents[1] / '.github/workflows'
        return {name: (root / name).read_text() for name in ('rc.yml', 'release.yml', 'wheels.yml')}

    def test_actual_workflows_preserve_top_level_protected_publishers(self):
        release_policy.verify_pypi_boundary(self.workflows())

    def test_reusable_publishing_or_missing_protection_is_rejected(self):
        for filename, old, new in [
            ('wheels.yml', 'contents: read', 'contents: read\n  id-token: write'),
            ('rc.yml', '    environment: release', '    environment: unprotected'),
            ('release.yml', 'https://upload.pypi.org/legacy/', 'https://test.pypi.org/legacy/'),
            ('rc.yml', 'python3 scripts/wheel-artifacts.py stage', 'python3 scripts/unverified-stage.py'),
            ('release.yml', 'needs: [authorize, build-wheels]', 'needs: authorize'),
        ]:
            workflows = self.workflows()
            self.assertIn(old, workflows[filename])
            workflows[filename] = workflows[filename].replace(old, new)
            with self.assertRaises(SystemExit, msg=(filename, old)):
                release_policy.verify_pypi_boundary(workflows)


if __name__ == '__main__':
    unittest.main()
