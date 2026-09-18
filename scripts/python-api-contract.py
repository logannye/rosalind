#!/usr/bin/env python3
"""Canonical Python source contract; parses files without importing the package.

This intentionally freezes implementation as well as API shape. The wire format
uses Python 3.9 syntax and ignores documentation, locations and empty optional AST
fields, including fields introduced by newer interpreters (such as type_params).
Packaging metadata is parsed and added by xtask, without a Python TOML dependency.
"""

import argparse
import ast
import json
from pathlib import Path
import tokenize


DOCSTRING_SCOPES = (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)


def canonical(node):
    if isinstance(node, ast.AST):
        fields = {}
        for name, value in ast.iter_fields(node):
            if name == "body" and isinstance(node, DOCSTRING_SCOPES):
                if value and isinstance(value[0], ast.Expr):
                    first = value[0].value
                    if isinstance(first, ast.Constant) and isinstance(first.value, str):
                        value = value[1:]
            # Empty optional fields have no semantics in the admitted grammar.
            # Omitting them stabilizes 3.9 ASTs on 3.11/3.12/3.14 interpreters.
            if value is None or (isinstance(value, list) and not value):
                continue
            fields[name] = canonical(value)
        return {"node": type(node).__name__, "fields": fields}
    if isinstance(node, list):
        return [canonical(value) for value in node]
    if isinstance(node, bytes):
        return {"bytes": node.hex()}
    if isinstance(node, complex):
        return {"complex": [node.real.hex(), node.imag.hex()]}
    if isinstance(node, float):
        return {"float": node.hex()}
    if node is Ellipsis:
        return {"ellipsis": True}
    if node is None or isinstance(node, (str, int, bool)):
        return node
    raise TypeError("unsupported AST field value: " + type(node).__name__)


def source_contract(root):
    package = Path(root) / "python" / "rosalind"
    modules = {}
    for path in sorted(package.rglob("*.py")):
        relative = path.relative_to(package).as_posix()
        with tokenize.open(path) as source:
            tree = ast.parse(source.read(), filename=relative, feature_version=(3, 9))
        modules[relative] = canonical(tree)
    if not modules:
        raise ValueError("Python contract requires python/rosalind/*.py sources")
    return {"schema": 1, "syntax": "python-3.9", "modules": modules}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    args = parser.parse_args()
    print(json.dumps(source_contract(args.root), sort_keys=True, separators=(",", ":")))


if __name__ == "__main__":
    main()
