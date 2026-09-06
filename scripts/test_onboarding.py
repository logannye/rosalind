import argparse
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch
import venv


SPEC = importlib.util.spec_from_file_location("onboarding", Path(__file__).with_name("onboarding.py"))
ONBOARDING = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ONBOARDING)
REPOSITORY = Path(__file__).resolve().parents[1]


class OnboardingTests(unittest.TestCase):
    def fixture(self, root):
        for guide in ONBOARDING.GUIDES:
            path = root / guide
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("# Guide\n")
        (root / "README.md").write_text("[SDK](docs/analyzer-sdk.md)\n[Python](python/README.md)\n")
        (root / "docs/analyzer-sdk.md").write_text("[Schema](nested/schema.json)\n[Example](../example/)\n")
        (root / "docs/nested").mkdir()
        (root / "docs/nested/schema.json").write_text("{}")
        (root / "example").mkdir()
        (root / "example/README.md").write_text("[Semantics](../docs/semantics.md)\n[Script](run.py)\n")
        (root / "example/run.py").write_text("print('example')\n")
        (root / "docs/semantics.md").write_text("[SDK](analyzer-sdk.md)\n")

    def test_staged_guides_keep_transitive_local_links_after_source_is_removed(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source"
            source.mkdir()
            self.fixture(source)
            bundle = root / "bundle"
            self.assertEqual(ONBOARDING.stage_guides(source, bundle), 6)
            source.rename(root / "unavailable-source")
            self.assertEqual(len(ONBOARDING.guide_closure(bundle)[0]), 6)
            self.assertTrue((bundle / "example/run.py").is_file())
            (bundle / "docs/nested/schema.json").unlink()
            with self.assertRaisesRegex(ValueError, "missing local link"):
                ONBOARDING.guide_closure(bundle)

    def test_checker_ignores_commands_remote_links_and_anchors_but_rejects_escape(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.fixture(root)
            document = root / "README.md"
            document.write_text("[web](https://example.invalid/missing) [section](#section)\n"
                                "```sh\n[example](missing-command-path)\n```\n")
            self.assertEqual(list(ONBOARDING.local_links(document, root)), [])
            document.write_text("[outside](../outside.md)\n")
            with self.assertRaisesRegex(ValueError, "escapes bundle"):
                list(ONBOARDING.local_links(document, root))

    def test_staging_does_not_copy_unlinked_files_or_recurse_into_source_root(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source"
            source.mkdir()
            self.fixture(source)
            (source / "example/private-local-data.bam").write_bytes(b"private")
            bundle = root / "bundle"
            ONBOARDING.stage_guides(source, bundle)
            self.assertFalse((bundle / "example/private-local-data.bam").exists())
            (source / "docs/analyzer-sdk.md").write_text("[Repository](../)\n")
            with self.assertRaisesRegex(ValueError, "source/private/build directory"):
                ONBOARDING.stage_guides(source, source / "dist/stage")

    def test_staging_rejects_private_paths_and_ancestor_of_destination(self):
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary)
            self.fixture(source)
            (source / "private").mkdir()
            (source / "private/password.txt").write_text("not for distribution")
            (source / "README.md").write_text("[Private](private/password.txt)\n")
            with self.assertRaisesRegex(ValueError, "source/private/build directory"):
                ONBOARDING.stage_guides(source, source / "dist/stage")
            (source / "README.md").write_text("[Example](example/)\n")
            with self.assertRaisesRegex(ValueError, "contains staging destination"):
                ONBOARDING.stage_guides(source, source / "example/stage")

    def test_documented_sdk_commands_create_the_binary_consumed_by_conformance(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bundle = root / "bundle"
            bundle.mkdir()
            self.fixture(bundle)
            (bundle / "docs/analyzer-sdk.md").write_text(
                (REPOSITORY / "docs/analyzer-sdk.md").read_text()
                .replace("../examples/evidence-analyzer/", "../example/")
                .replace("SEMANTICS.md", "semantics.md")
            )
            commands = root / "bin"
            commands.mkdir()
            binary = commands / "rosalind"
            binary.write_text("#!/usr/bin/env python3\n" + r'''
import json
import pathlib
import sys
if sys.argv[1:] == ["--version"]:
    print("rosalind 0.5.0")
elif sys.argv[1:3] == ["new", "analyzer"]:
    path = pathlib.Path(sys.argv[sys.argv.index("--output") + 1])
    path.mkdir()
    (path / "Cargo.toml").write_text('rosalind-bio = { version = "=0.5.0" }\n')
elif sys.argv[1:3] == ["conformance", "analyzer"]:
    path = pathlib.Path(sys.argv[sys.argv.index("--binary") + 1])
    assert path.is_file(), f"documented binary was never built: {path}"
    print(json.dumps({"passed": True}))
else:
    raise AssertionError(sys.argv)
''')
            cargo = commands / "cargo"
            cargo.write_text("#!/usr/bin/env python3\n" + r'''
import os
import pathlib
import sys
assert pathlib.Path.cwd().name == "locus-qc"
assert pathlib.Path(os.environ["CARGO_HOME"]).name == "cargo-home"
assert "CARGO_TARGET_DIR" not in os.environ
assert not (pathlib.Path.cwd().parent / ".cargo/config.toml").exists()
pathlib.Path("Cargo.lock").write_text("lockfile")
if sys.argv[1] == "build":
    assert "--release" in sys.argv and "--locked" in sys.argv and "--offline" in sys.argv
    path = pathlib.Path("target/release/locus-qc")
    path.parent.mkdir(parents=True)
    path.write_text("built binary")
''')
            binary.chmod(0o755)
            cargo.chmod(0o755)
            args = argparse.Namespace(bundle=bundle, binary=binary, python=None,
                                      candidate_source=None, registry_sdk=True)
            with patch.dict(os.environ, {"CARGO_TARGET_DIR": "must-not-leak"}):
                ONBOARDING.smoke(args)

    def test_candidate_source_version_must_match_packaged_cli(self):
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary)
            (source / "Cargo.toml").write_text('[package]\nversion = "0.5.0-rc.1"\n')
            with self.assertRaisesRegex(ValueError, "does not match"):
                ONBOARDING.candidate_config(source, "0.5.0")
            config = ONBOARDING.candidate_config(source, "0.5.0-rc.1")
            self.assertIn("[patch.crates-io]", config)
            self.assertIn(str(source / "crates/build-info"), config)

    def test_python_launcher_preserves_the_installed_virtual_environment(self):
        with tempfile.TemporaryDirectory() as temporary:
            environment = Path(temporary) / "venv"
            venv.EnvBuilder(with_pip=False, symlinks=True).create(environment)
            executable = ONBOARDING.python_executable(environment / "bin/python")
            prefix = subprocess.check_output(
                [str(executable), "-c", "import sys; print(sys.prefix)"], text=True
            ).strip()
            self.assertEqual(Path(prefix).resolve(), environment.resolve())

    def test_real_bundle_contains_documented_code_assets_and_registry_example_manifest(self):
        with tempfile.TemporaryDirectory() as temporary:
            bundle = Path(temporary)
            ONBOARDING.stage_guides(REPOSITORY, bundle)
            for path in ONBOARDING.CODE_FILES + (
                "integrations/nextflow/examples/evidence/main.nf",
                "integrations/nextflow/modules/rosalind/evidence.nf",
                "integrations/snakemake/Evidence.smk",
            ):
                self.assertTrue((bundle / path).is_file(), path)
            metadata = json.loads((bundle / "ONBOARDING-BUNDLE.json").read_text())
            manifest = (bundle / metadata["rendered_manifest"]).read_text()
            self.assertIn(f'version = "={metadata["sdk_registry_version"]}"', manifest)
            self.assertNotIn('path = "../.."', manifest)
            self.assertEqual((bundle / "examples/evidence-analyzer/src/main.rs").read_bytes(),
                             (REPOSITORY / "examples/evidence-analyzer/src/main.rs").read_bytes())
            tutorial = (bundle / "examples/research-filter/README.md").read_text()
            self.assertNotIn("target/debug/rosalind", tutorial)
            self.assertIn("./rosalind analyze evidence", tutorial)

    def test_maintained_python_snippet_has_defined_row_consumer(self):
        readme = (REPOSITORY / "python/README.md").read_text()
        program = ONBOARDING.snippet(readme, "python-evidence", "python")
        compile(program, "python-README", "exec")
        with self.assertRaisesRegex(ValueError, "exactly one"):
            ONBOARDING.snippet("missing block", "python-evidence", "python")


if __name__ == "__main__":
    unittest.main()
