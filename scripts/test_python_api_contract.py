import ast
import importlib.util
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).with_name("python-api-contract.py")
SPEC = importlib.util.spec_from_file_location("python_api_contract", SCRIPT)
CONTRACT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CONTRACT)


class PythonApiContract(unittest.TestCase):
    def fixture(self, source):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        root = Path(temporary.name)
        package = root / "python" / "rosalind"
        package.mkdir(parents=True)
        (package / "__init__.py").write_text(source)
        return root

    def test_comments_docstrings_formatting_and_locations_do_not_change_contract(self):
        first = self.fixture('''"""Module docs."""
class Example:
    """Class docs."""
    async def run(self, value: int = 3):
        """Method docs."""
        return value + 1
''')
        second = self.fixture('''# New comment and leading blank lines.

"""Rewritten module docs."""
class Example:
    """Rewritten class docs."""

    async def run(self, value: int=3): # More documentation.
        """Rewritten method docs."""
        return (value + 1)
''')
        self.assertEqual(CONTRACT.source_contract(first), CONTRACT.source_contract(second))

    def test_signatures_defaults_annotations_and_implementation_are_frozen(self):
        original = "def exported(value: int = 3):\n    return value + 1\n"
        root = self.fixture(original)
        expected = CONTRACT.source_contract(root)
        for changed in [original.replace("value: int", "other: int"),
                        original.replace("= 3", "= 4"),
                        original.replace(": int", ": float"),
                        original.replace("+ 1", "- 1"),
                        original + "\ndef _private():\n    return 2\n"]:
            with self.subTest(changed=changed):
                (root / "python/rosalind/__init__.py").write_text(changed)
                self.assertNotEqual(expected, CONTRACT.source_contract(root))

    def test_module_additions_and_removals_are_frozen(self):
        root = self.fixture("value = 1\n")
        expected = CONTRACT.source_contract(root)
        added = root / "python/rosalind/nested/extra.py"
        added.parent.mkdir()
        added.write_text("value = 2\n")
        self.assertNotEqual(expected, CONTRACT.source_contract(root))
        added.unlink()
        self.assertEqual(expected, CONTRACT.source_contract(root))

    def test_package_is_never_imported_or_executed(self):
        root = self.fixture("raise RuntimeError('must never execute')\nimport unavailable_dependency\n")
        marker = root / "unexpected-import"
        with (root / "python/rosalind/__init__.py").open("a") as source:
            source.write(f"open({str(marker)!r}, 'w').write('executed')\n")
        result = subprocess.run([sys.executable, "-I", str(SCRIPT), "--root", str(root)],
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertFalse(marker.exists())

    def test_new_empty_ast_fields_are_ignored_but_nonempty_fields_are_frozen(self):
        node = ast.parse("def f():\n    return 1\n", feature_version=(3, 9)).body[0]
        expected = CONTRACT.canonical(node)
        node._fields = tuple(name for name in node._fields if name != "type_params") + ("type_params",)
        node.type_params = []
        self.assertEqual(expected, CONTRACT.canonical(node))
        node.type_params = [ast.Name(id="T", ctx=ast.Load())]
        self.assertNotEqual(expected, CONTRACT.canonical(node))

    def test_future_syntax_invalid_sources_and_missing_packages_fail_closed(self):
        for source in ["match value:\n    case 1: pass\n", "def broken(:\n"]:
            with self.subTest(source=source), self.assertRaises(SyntaxError):
                CONTRACT.source_contract(self.fixture(source))
        root = self.fixture("value = 1\n")
        (root / "python/rosalind/__init__.py").unlink()
        with self.assertRaises(ValueError):
            CONTRACT.source_contract(root)

    def test_constant_kinds_remain_distinct_and_json_serializable(self):
        values = ["None", "True", "1", "1.0", "1j", "b'x'", "'x'", "..."]
        observed = set()
        for value in values:
            root = self.fixture("value = " + value + "\n")
            result = subprocess.run([sys.executable, "-I", str(SCRIPT), "--root", str(root)],
                                    capture_output=True, text=True, check=True)
            observed.add(result.stdout)
        self.assertEqual(len(observed), len(values))


if __name__ == "__main__":
    unittest.main()
