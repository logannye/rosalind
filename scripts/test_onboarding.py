import argparse
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
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

    def test_adoption_public_schema_exception_does_not_include_partner_records(self):
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "source"
            source.mkdir()
            self.fixture(source)
            public = source / "release/schemas/design-partner-v1.schema.json"
            public.parent.mkdir(parents=True)
            public.write_text('{}')
            private = source / "release/private-design-partners/session.json"
            private.parent.mkdir()
            private.write_text('not for distribution')
            (source / "README.md").write_text(f"[Schema]({public.relative_to(source)})\n")
            bundle = Path(temporary) / "bundle"
            ONBOARDING.stage_guides(source, bundle)
            self.assertEqual((bundle / public.relative_to(source)).read_text(), '{}')
            self.assertFalse((bundle / private.relative_to(source)).exists())
            (source / "README.md").write_text(f"[Session]({private.relative_to(source)})\n")
            with self.assertRaisesRegex(ValueError, "source/private/build directory"):
                ONBOARDING.stage_guides(source, bundle)

    def test_documented_sdk_commands_create_the_binary_consumed_by_conformance(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bundle = root / "bundle"
            # Exercise the full current guide closure; new navigation links must
            # not require a hand-maintained subset or rewritten guide prose.
            ONBOARDING.stage_guides(REPOSITORY, bundle)
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
project_path = (pathlib.Path(sys.argv[sys.argv.index("--manifest-path") + 1]).parent
                if "--manifest-path" in sys.argv else pathlib.Path.cwd())
project = project_path.name
if project == "evidence-analyzer":
    assert sys.argv[1] in ("fetch", "test")
    sys.exit(0)
assert project in ("locus-qc", "candidate-qc")
assert pathlib.Path(os.environ["CARGO_HOME"]).name == "cargo-home"
assert "CARGO_TARGET_DIR" not in os.environ
assert not (pathlib.Path.cwd().parent / ".cargo/config.toml").exists()
pathlib.Path("Cargo.lock").write_text("lockfile")
if sys.argv[1] == "build":
    assert "--release" in sys.argv and "--locked" in sys.argv and "--offline" in sys.argv
    path = pathlib.Path("target/release") / project
    path.parent.mkdir(parents=True)
    path.write_text("built binary")
''')
            binary.chmod(0o755)
            cargo.chmod(0o755)
            args = argparse.Namespace(bundle=bundle, binary=binary, python=None,
                                      candidate_source=None, registry_sdk=True)
            with patch.dict(os.environ, {"CARGO_TARGET_DIR": "must-not-leak"}), \
                    patch.object(ONBOARDING, "smoke_research_workflows") as research:
                ONBOARDING.smoke(args)
            research.assert_called_once()
            self.assertEqual(research.call_args.args[:3], (bundle.resolve(), binary.resolve(), None))

    def test_researcher_example_uses_exact_launchers_and_quoted_private_paths(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            work = root / "work space ' $(touch quoting-escaped)"
            work.mkdir()
            log = root / "launches.jsonl"
            launchers = root / "chosen launchers"
            launchers.mkdir()
            logger = "#!" + sys.executable + "\n" + r'''
import json
import os
from pathlib import Path
import sys
with Path(os.environ["ONBOARDING_TEST_LOG"]).open("a") as output:
    output.write(json.dumps([sys.argv[0], *sys.argv[1:]]) + "\n")
if len(sys.argv) > 2 and sys.argv[1].endswith("prepare.py"):
    Path(sys.argv[2]).mkdir()
'''
            python = launchers / "selected-python"
            binary = launchers / "selected-native"
            for executable in (python, binary):
                executable.write_text(logger)
                executable.chmod(0o755)
            markdown = (REPOSITORY / "examples/research-filter/README.md").read_text()
            program = ONBOARDING.research_snippet(markdown, "researcher-quickstart",
                                                  python, binary, work)
            environment = dict(os.environ, ONBOARDING_TEST_LOG=str(log))
            subprocess.run(["bash", "-euo", "pipefail", "-c", program],
                           cwd=work, env=environment, check=True)
            launches = [json.loads(line) for line in log.read_text().splitlines()]
            self.assertEqual([Path(args[0]) for args in launches],
                             [python, binary, python, binary, binary])
            self.assertEqual(launches[0][-1], str(work / "inputs"))
            self.assertIn(str(work / "inputs/candidates.evidence.vcf"), launches[1])
            self.assertTrue((work / "inputs/research-review.tsv").is_file())
            self.assertFalse((work / "quoting-escaped").exists())
            self.assertNotIn("/tmp/research-filter", program)
            self.assertNotIn("/tmp/rosalind-tutorial-env", program)
            bundled = markdown.replace("target/debug/rosalind", "./rosalind")
            self.assertEqual(ONBOARDING.research_snippet(bundled, "researcher-quickstart",
                                                       python, binary, work), program)
            # Linux commonly places the smoke's own unique root under /tmp.
            private = Path("/tmp/rosalind-onboarding-unique/research-workflows/researcher")
            self.assertIn(str(private / "inputs"), ONBOARDING.research_snippet(
                markdown, "researcher-quickstart", python, binary, private))
            # A caller may explicitly select an existing tutorial environment;
            # it must remain the selected interpreter rather than get replaced.
            selected = Path("/tmp/rosalind-reuse-env/bin/python")
            reuse = (REPOSITORY / "docs/reuse-quickstart.md").read_text()
            self.assertTrue(ONBOARDING.research_snippet(
                reuse, "reuse-quickstart", selected, binary, private).startswith(str(selected)))

    def test_research_smoke_isolates_dependencies_tempfiles_and_exact_binary(self):
        for supplied_python in (True, False):
            with self.subTest(supplied_python=supplied_python), \
                    tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                bundle = root / "bundle"
                ONBOARDING.stage_guides(REPOSITORY, bundle)
                work = root / "smoke"
                work.mkdir()
                binary = root / "exact-selected-native"
                binary.write_text("not executed by this orchestration test")
                python = ONBOARDING.python_executable(sys.executable) if supplied_python else None
                calls = []

                def record(args, **kwargs):
                    args = list(map(str, args))
                    calls.append((args, kwargs))
                    if args[1:3] == ["-m", "venv"]:
                        interpreter = Path(args[3]) / "bin/python"
                        interpreter.parent.mkdir(parents=True)
                        interpreter.symlink_to(sys.executable)

                with patch.object(ONBOARDING, "command", side_effect=record):
                    ONBOARDING.smoke_research_workflows(bundle, binary, python, work,
                                                       {"PATH": "/intentionally-unrelated"})
                research_root = work / "research-workflows"
                expected_python = python or research_root / "venv/bin/python"
                launched = [entry for entry in calls if entry[0][0] == "bash"]
                self.assertEqual(len(launched), 2)
                for _, kwargs in calls:
                    for name in ("TMPDIR", "TMP", "TEMP"):
                        self.assertEqual(Path(kwargs["env"][name]), research_root / "temporary")
                for args, kwargs in launched:
                    self.assertTrue(Path(kwargs["cwd"]).is_relative_to(work))
                    self.assertTrue((Path(kwargs["cwd"]) / "examples/research-filter/prepare.py").is_file())
                    self.assertIn(str(expected_python), args[-1])
                    self.assertNotIn("/tmp/rosalind-", args[-1])
                    selected = Path(kwargs["env"]["PATH"].split(os.pathsep)[0]) / "rosalind"
                    self.assertEqual(selected.resolve(), binary.resolve())
                self.assertEqual(len([args for args, _ in calls if args[1:3] == ["-m", "venv"]]),
                                 0 if supplied_python else 1)
                installs = [args for args, _ in calls if args[1:4] == ["-m", "pip", "install"]]
                self.assertEqual(installs, [] if supplied_python else [
                    [str(expected_python), "-m", "pip", "install", "pysam==0.23.3"]])

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
            self.assertNotIn('path = "../../crates/build-info"', manifest)
            self.assertIn('rosalind-build-info = "=0.1.0"', manifest)
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
