"""The completion probe must work with historical and current CLI interfaces."""
import hashlib
import os
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class InstallerCompatibilityTest(unittest.TestCase):
    def test_successful_install_does_not_require_historical_version_flag(self):
        for current in (False, True):
            with self.subTest(current=current), tempfile.TemporaryDirectory() as temp:
                root = Path(temp)
                tools, assets, work = (root / name for name in ('tools', 'assets', 'work'))
                for path in (tools, assets, work):
                    path.mkdir()
                bundle = root / 'rosalind-aarch64-apple-darwin'
                bundle.mkdir()
                binary = bundle / 'rosalind'
                binary.write_text('#!/bin/sh\n' + (
                    'if [ "$1" = --version ]; then echo "rosalind 0.5.0"; exit 0; fi\n'
                    'if [ "$1 $2" = "analyze evidence" ]; then exit 0; fi\n'
                    if current else
                    'if [ "$1" = --version ]; then echo "unsupported flag" >&2; fi\n'
                ) + 'exit 2\n')
                binary.chmod(0o755)
                archive = assets / (bundle.name + '.tar.gz')
                with tarfile.open(archive, 'w:gz') as stream:
                    stream.add(bundle, arcname=bundle.name)
                digest = hashlib.sha256(archive.read_bytes()).hexdigest()
                archive.with_name(archive.name + '.sha256').write_text(f'{digest}  {archive.name}\n')
                uname = tools / 'uname'
                uname.write_text('#!/bin/sh\ncase "$1" in -s) echo Darwin;; -m) echo arm64;; esac\n')
                gh = tools / 'gh'
                gh.write_text('#!' + sys.executable + '\n' +
                              'import os, pathlib, shutil, sys\n'
                              'assert sys.argv[1:3] == ["release", "download"]\n'
                              'name = sys.argv[sys.argv.index("--pattern") + 1]\n'
                              'shutil.copyfile(pathlib.Path(os.environ["INSTALLER_FIXTURE_ASSETS"]) / name, name)\n')
                for command in (uname, gh):
                    command.chmod(0o755)
                env = dict(os.environ, PATH=str(tools) + os.pathsep + os.environ['PATH'],
                           INSTALLER_FIXTURE_ASSETS=str(assets), ROSALIND_VERSION='v-test')
                result = subprocess.run(['sh', str(ROOT / 'install.sh')], cwd=work,
                                        env=env, text=True, capture_output=True)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertTrue((work / bundle.name / 'rosalind').is_file())
                if current:
                    self.assertIn('rosalind 0.5.0', result.stdout)
                    self.assertIn('supports exact evidence', result.stdout)
                else:
                    self.assertIn('does not expose --version', result.stdout)
                    self.assertIn('This is a legacy release', result.stdout)
                    self.assertIn('source installation guide', result.stdout)


if __name__ == '__main__':
    unittest.main()
