#!/usr/bin/env python3
"""Stage/check offline guide links and execute maintained onboarding snippets."""

import argparse
import json
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from urllib.parse import unquote, urlsplit


GUIDES = ("README.md", "CONTRACT.md", "docs/analyzer-sdk.md", "python/README.md")
EXCLUDED_PARTS = {".git", "target", "dist", "release", "private", "secrets", ".env",
                  "__pycache__", ".venv", "venv", "node_modules"}
MAX_GUIDE_BYTES = 8 << 20
# Files used by maintained code blocks rather than Markdown hyperlinks. These
# narrow reviewed paths are not permission to copy arbitrary example data.
CODE_FILES = (
    "examples/research-filter/prepare.py", "examples/research-filter/join.py",
    "examples/research-filter/sources.json", "examples/evidence-analyzer/Cargo.toml",
    "examples/evidence-analyzer/Cargo.lock", "examples/evidence-analyzer/src/main.rs",
)
CODE_DIRECTORIES = ("integrations/nextflow", "integrations/snakemake")


def local_links(document, root):
    """Resolve local Markdown destinations; anchors and remote URLs need no files."""
    content = document.read_text(encoding="utf-8")
    content = re.sub(r"^```[^\n]*\n.*?^```\s*$", "", content, flags=re.M | re.S)
    for target in re.findall(r"\]\((<[^>]+>|[^\s)]+)(?:\s+\"[^\"]*\")?\)", content):
        parsed = urlsplit(target.strip("<>"))
        if parsed.scheme or parsed.netloc or not parsed.path:
            continue
        destination = (document.parent / unquote(parsed.path)).resolve()
        try:
            destination.relative_to(root)
        except ValueError as error:
            raise ValueError(f"{document}: link escapes bundle: {target}") from error
        if not destination.exists():
            raise ValueError(f"{document}: missing local link: {target}")
        yield destination


def guide_closure(root):
    """Follow linked Markdown and directory READMEs, including their local assets."""
    root = Path(root).resolve()
    pending = [root / name for name in GUIDES]
    documents = set()
    assets = set(pending)
    while pending:
        document = pending.pop()
        if document in documents:
            continue
        if not document.is_file():
            raise ValueError(f"missing onboarding guide: {document}")
        documents.add(document)
        for destination in local_links(document, root):
            assets.add(destination)
            if destination.is_dir():
                readme = destination / "README.md"
                if readme.is_file():
                    pending.append(readme)
            elif destination.suffix.lower() == ".md":
                pending.append(destination)
    return documents, assets


def stage_guides(source, destination):
    source = Path(source).resolve()
    destination = Path(destination).resolve()
    _, assets = guide_closure(source)
    files = set()
    directories = set()

    def validate(asset):
        relative = asset.resolve().relative_to(source)
        if not relative.parts or any(part in EXCLUDED_PARTS for part in relative.parts):
            raise ValueError(f"onboarding link cannot package source/private/build directory: {asset}")
        if asset.is_dir() and destination.is_relative_to(asset.resolve()):
            raise ValueError(f"onboarding directory contains staging destination: {asset}")

    for asset in assets:
        validate(asset)
        if not asset.is_dir():
            files.add(asset)
            continue
        directories.add(asset)
        if (asset / "README.md").is_file():
            # A directory link navigates to its guide. Only explicitly linked
            # assets are included; generated files and ignored inputs are not.
            files.add(asset / "README.md")
            continue
        # Historical result/fixture directory links intentionally identify raw
        # artifacts. Include only Git-tracked files, never arbitrary local files.
        relative = asset.relative_to(source)
        tracked = subprocess.check_output(
            ["git", "-C", str(source), "ls-files", "-z", "--", str(relative) + "/"]
        ).decode().split("\0")
        if not any(tracked):
            raise ValueError(f"linked artifact directory has no tracked files: {asset}")
        for name in filter(None, tracked):
            path = source / name
            validate(path)
            if path.is_file():
                files.add(path)
    for name in CODE_FILES:
        path = source / name
        if path.is_file():
            validate(path)
            files.add(path)
    for name in CODE_DIRECTORIES:
        if not (source / name).is_dir():
            continue
        tracked = subprocess.check_output(
            ["git", "-C", str(source), "ls-files", "-z", "--", name + "/"]
        ).decode().split("\0")
        for name in filter(None, tracked):
            path = source / name
            validate(path)
            if path.is_file():
                files.add(path)
    total = sum(path.stat().st_size for path in files)
    if total > MAX_GUIDE_BYTES:
        raise ValueError(f"offline guides total {total} bytes; link large results remotely instead of bundling them")
    for directory in directories:
        (destination / directory.relative_to(source)).mkdir(parents=True, exist_ok=True)
    for asset in sorted(files):
        output = destination / asset.relative_to(source)
        output.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(asset, output)
    example_manifest = destination / "examples/evidence-analyzer/Cargo.toml"
    if example_manifest.is_file():
        source_manifest = (source / "Cargo.toml").read_text()
        version = re.search(r'^version\s*=\s*"([^"]+)"', source_manifest, flags=re.M).group(1)
        old = 'rosalind-bio = { path = "../.." }'
        new = f'rosalind-bio = {{ version = "={version}" }}'
        manifest = example_manifest.read_text()
        if manifest.count(old) != 1:
            raise ValueError("standalone example path dependency changed; review bundle rendering")
        example_manifest.write_text(manifest.replace(old, new))
        metadata = {
            "source_commit": subprocess.check_output(
                ["git", "-C", str(source), "rev-parse", "HEAD"], text=True
            ).strip(),
            "source_tree_dirty": bool(subprocess.check_output(
                ["git", "-C", str(source), "status", "--porcelain"], text=True
            ).strip()),
            "sdk_registry_version": version,
            "rendered_manifest": "examples/evidence-analyzer/Cargo.toml",
            "original_dependency": old,
            "bundled_dependency": new,
            "note": "Source Cargo.lock is retained; the first registry resolution may update it. Candidate SDK validation requires explicit source patches.",
        }
        (destination / "ONBOARDING-BUNDLE.json").write_text(json.dumps(metadata, indent=2) + "\n")
        example_readme = example_manifest.with_name("README.md")
        if example_readme.is_file():
            example_readme.write_text(
                "> Bundle note: Cargo.toml now pins the exact registry SDK version recorded in "
                "`../../ONBOARDING-BUNDLE.json`. For an unpublished candidate, configure the "
                "explicit source patches in the SDK guide before invoking Cargo.\n\n"
                + example_readme.read_text()
            )
        tutorial = destination / "examples/research-filter/README.md"
        if tutorial.is_file():
            text = tutorial.read_text().replace("target/debug/rosalind", "./rosalind")
            text = text.replace(
                "Build the current development CLI first (`cargo build --locked --bin rosalind`).\nFrom the repository root,",
                "Use the executable included in this native bundle.\nFrom the extracted bundle root,",
            )
            tutorial.write_text(text)
            metadata["rendered_tutorial"] = "examples/research-filter/README.md"
            metadata["tutorial_binary"] = "./rosalind"
            (destination / "ONBOARDING-BUNDLE.json").write_text(json.dumps(metadata, indent=2) + "\n")
    documents, _ = guide_closure(destination)
    return len(documents)


def snippet(markdown, name, language):
    expression = r"<!-- smoke:" + re.escape(name) + r" -->\s*```" + language + r"\n(.*?)\n```"
    matches = re.findall(expression, markdown, flags=re.S)
    if len(matches) != 1:
        raise ValueError(f"expected exactly one {name} onboarding snippet, found {len(matches)}")
    return matches[0]


def command(args, **kwargs):
    subprocess.run([str(arg) for arg in args], check=True, **kwargs)


def python_executable(path):
    # Resolving the venv/bin/python symlink executes its base interpreter and
    # loses the installed environment. Normalize cwd without following symlinks.
    path = Path(os.path.abspath(path))
    if not path.is_file():
        raise ValueError(f"Python executable is missing: {path}")
    return path


def candidate_config(source, native_version):
    """Return explicit Cargo patches after checking this is the matching source."""
    manifest = (source / "Cargo.toml").read_text()
    version = re.search(r'^version\s*=\s*"([^"]+)"', manifest, flags=re.M)
    if version is None or version.group(1) != native_version:
        raise ValueError("candidate source Cargo version does not match packaged executable")
    return "[patch.crates-io]\n" + "".join(
        f"{name} = {{ path = {json.dumps(str(path))} }}\n"
        for name, path in (
            ("rosalind-bio", source),
            ("rosalind-build-info", source / "crates/build-info"),
        )
    )


def smoke(args):
    bundle = args.bundle.resolve()
    binary = args.binary.resolve(strict=True)
    documents, _ = guide_closure(bundle)
    print(f"offline onboarding links: {len(documents)} reachable guides", flush=True)
    native_version = subprocess.check_output([str(binary), "--version"], text=True).strip().split()[-1]
    source = args.candidate_source.resolve(strict=True) if args.candidate_source else None
    with tempfile.TemporaryDirectory(prefix="rosalind-onboarding-") as temporary:
        work = Path(temporary).resolve()
        environment = os.environ.copy()
        environment.pop("CARGO_TARGET_DIR", None)
        environment.pop("PYTHONPATH", None)
        environment["PYTHONNOUSERSITE"] = "1"
        environment["PATH"] = str(binary.parent) + os.pathsep + environment.get("PATH", "")
        # Avoid inheriting user Cargo config/patches, including in registry mode.
        environment["CARGO_HOME"] = str(work / "cargo-home")
        if source is not None:
            config = work / ".cargo/config.toml"
            config.parent.mkdir()
            config.write_text(candidate_config(source, native_version))
            print(f"SDK origin: explicit candidate source {source}; not registry validation", flush=True)
        else:
            print(f"SDK origin: registry-only exact version {native_version}", flush=True)
        sdk = snippet((bundle / "docs/analyzer-sdk.md").read_text(), "legacy-scaffold", "sh")
        command(["bash", "-euo", "pipefail", "-c", sdk], cwd=work, env=environment)
        report = json.loads((work / "locus-qc/conformance.json").read_text())
        if report.get("passed") is not True:
            raise ValueError(f"packaged scaffold conformance failed: {report}")
        generated = (work / "locus-qc/Cargo.toml").read_text()
        if f'version = "={native_version}"' not in generated:
            raise ValueError("scaffold SDK version differs from packaged executable")
        if not (work / "locus-qc/Cargo.lock").is_file():
            raise ValueError("documented analyzer build did not retain a lockfile")
        example = bundle / "examples/evidence-analyzer"
        if example.is_dir():
            for name in ("Cargo.toml", "Cargo.lock", "src/main.rs"):
                if not (example / name).is_file():
                    raise ValueError(f"bundled standalone analyzer is missing {name}")
            standalone = work / "evidence-analyzer"
            shutil.copytree(example, standalone)
            # Resolve the retained source lockfile against the explicitly chosen
            # SDK origin, then compile/test the actual packaged Rust example.
            command(["cargo", "fetch", "--manifest-path", standalone / "Cargo.toml"],
                    cwd=work, env=environment)
            command(["cargo", "test", "--manifest-path", standalone / "Cargo.toml",
                     "--locked", "--offline", "--target-dir", work / "locus-qc/target"],
                    cwd=work, env=environment)
        if args.python:
            python = python_executable(args.python)
            fixture = work / "python-example"
            fixture.mkdir()
            shutil.copy2(bundle / "examples/data/illumina_toy/reference.fa", fixture / "genome.fa")
            command([binary, "sort", "--input", bundle / "examples/data/illumina_toy/alignments.bam",
                     "--output", fixture / "sample.bam"], env=environment)
            # Execute the README embedded in the installed wheel's metadata,
            # not an import or example accidentally taken from the checkout.
            program = r'''
import importlib.metadata
import json
import pathlib
import re
import subprocess
import sys
import sysconfig
import pysam
import rosalind
from rosalind.features import _normalized_version

binary = pathlib.Path(sysconfig.get_path("scripts")) / "rosalind"
assert binary.resolve() == pathlib.Path(sys.argv[1]).resolve()
native = subprocess.check_output([str(binary), "--version"], text=True).strip().split()[-1]
assert _normalized_version(native) == _normalized_version(rosalind.__version__)
assert pathlib.Path(rosalind.__file__).resolve().is_relative_to(pathlib.Path(sys.prefix).resolve())
pysam.faidx("genome.fa")
pysam.index("sample.bam")
with pysam.FastaFile("genome.fa") as reference:
    contig = reference.references[0]
    length = min(100, reference.lengths[0])
    base = reference.fetch(contig, 0, 1).upper()
alternate = next(nucleotide for nucleotide in "ACGT" if nucleotide != base)
pathlib.Path("targets.bed").write_text(f"{contig}\t0\t{length}\tfirst\n")
pathlib.Path("candidates.vcf").write_text(
    "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    f"{contig}\t1\t.\t{base}\t{alternate}\t.\tPASS\t.\n")
readme = importlib.metadata.metadata("rosalind-bio").get_payload()
blocks = re.findall(r"<!-- smoke:python-evidence -->\s*```python\n(.*?)\n```", readme, re.S)
assert len(blocks) == 1, "installed wheel README lacks its executable onboarding example"
namespace = {}
exec(compile(blocks[0], "installed-wheel-README", "exec"), namespace)
assert namespace["rows"] == 1
assert namespace["result"].returncode == 0
assert namespace["summary"].returncode == 0
for result in (namespace["result"], namespace["summary"]):
    subprocess.run([str(binary), "verify", "--manifest", str(result.manifest_path)], check=True)
with pathlib.Path("replay.json").open("w") as output:
    subprocess.run([str(binary), "reproduce", "--manifest", str(namespace["result"].manifest_path),
                    "--inputs", str(pathlib.Path.cwd()), "--binary", str(binary), "--json"],
                   stdout=output, check=True)
print(f"installed wheel README passed with native/Python identity {native}")
'''
            command([python, "-c", program, binary], cwd=fixture, env=environment)
    print("packaged onboarding snippets and conformance passed outside the checkout", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="action", required=True)
    stage = commands.add_parser("stage")
    stage.add_argument("source", type=Path)
    stage.add_argument("destination", type=Path)
    check = commands.add_parser("check")
    check.add_argument("bundle", type=Path)
    run = commands.add_parser("smoke")
    run.add_argument("--bundle", type=Path, required=True)
    run.add_argument("--binary", type=Path, required=True)
    run.add_argument("--python", type=Path)
    origin = run.add_mutually_exclusive_group(required=True)
    origin.add_argument("--candidate-source", type=Path)
    origin.add_argument("--registry-sdk", action="store_true")
    args = parser.parse_args()
    if args.action == "stage":
        print(f"staged {stage_guides(args.source, args.destination)} linked onboarding guides")
    elif args.action == "check":
        documents, _ = guide_closure(args.bundle)
        print(f"offline onboarding links: {len(documents)} reachable guides")
    else:
        smoke(args)


if __name__ == "__main__":
    main()
