"""Research-run integrity, bounded file identities, and create-new publication.

This format is deliberately separate from Rosalind's native receipt format.
Hashes bind bytes and recorded provenance; they do not authenticate an author.
"""
from __future__ import annotations

import ctypes
import errno
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import stat
import subprocess
import sys
import threading
import time

SCHEMA = "rosalind-cohort-research-run-v1"
MAX_MANIFEST_BYTES = 8 << 20


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False)


def digest_json(value):
    return hashlib.sha256(canonical(value).encode()).hexdigest()


def load_json(path, limit=MAX_MANIFEST_BYTES):
    with Path(path).open("rb") as stream:
        raw = stream.read(limit + 1)
    if len(raw) > limit:
        raise ValueError(f"metadata exceeds {limit} bytes: {path}")
    return json.loads(raw)


def write_json(path, value):
    with Path(path).open("x", encoding="utf-8") as stream:
        stream.write(canonical(value) + "\n")


def file_digest(path, algorithm="sha256"):
    if algorithm == "blake3":
        try:
            import blake3
        except ImportError as error:
            raise RuntimeError("native identity checks require the blake3 Python package") from error
        digest = blake3.blake3()
    elif algorithm == "sha256":
        digest = hashlib.sha256()
    else:
        raise ValueError("unsupported digest")
    path = Path(path)
    before = path.stat()
    if not stat.S_ISREG(before.st_mode):
        raise ValueError(f"identity input must be a regular file: {path}")
    count = 0
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
            count += len(block)
    after = path.stat()
    key = lambda s: (s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns)
    if key(before) != key(after) or count != before.st_size:
        raise ValueError(f"input changed while hashing: {path}")
    return digest.hexdigest(), count


class InputLedger:
    """Union of successful children's actual verified dependencies, with a cap."""

    def __init__(self, max_files=16384):
        self.max_files = max_files
        self.entries = {}
        self.bytes_hashed = 0

    def pin(self, path, expected=None, algorithm="sha256"):
        access = Path(path).absolute()
        target = access.resolve(strict=True)
        key = (str(access), algorithm)
        if key not in self.entries and len(self.entries) >= self.max_files:
            raise ValueError("consumed file inventory limit exceeded")
        # Hash now as well as at finalization: successful child identities must
        # agree with earlier uses, not just with whatever happens to exist later.
        actual, size = file_digest(access, algorithm)
        self.bytes_hashed += size
        if expected is not None and actual != expected:
            raise ValueError(f"verified input identity changed: {access}")
        value = {"path": str(access), "resolved_path": str(target),
                 "algorithm": algorithm, "digest": actual, "bytes": size}
        if key in self.entries and self.entries[key] != value:
            raise ValueError(f"repeated input identity changed: {access}")
        self.entries[key] = value
        return value

    def recheck(self):
        for value in self.entries.values():
            path = Path(value["path"])
            if str(path.resolve(strict=True)) != value["resolved_path"]:
                raise ValueError(f"input path retargeted: {path}")
            actual, size = file_digest(path, value["algorithm"])
            self.bytes_hashed += size
            if actual != value["digest"] or size != value["bytes"]:
                raise ValueError(f"input changed before publication: {path}")

    def as_list(self):
        return [self.entries[key] for key in sorted(self.entries)]


def publish_new(staging, destination):
    """Atomic no-replace directory rename on supported benchmark platforms."""
    staging, destination = Path(staging), Path(destination)
    libc = ctypes.CDLL(None, use_errno=True)
    if sys.platform == "linux" and hasattr(libc, "renameat2"):
        rename = libc.renameat2
        rename.argtypes = [ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p, ctypes.c_uint]
        code = rename(-100, os.fsencode(staging), -100, os.fsencode(destination), 1)
    elif sys.platform == "darwin" and hasattr(libc, "renamex_np"):
        rename = libc.renamex_np
        rename.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_uint]
        code = rename(os.fsencode(staging), os.fsencode(destination), 4)  # RENAME_EXCL
    else:
        raise RuntimeError("atomic create-new directory publication requires Linux or macOS")
    if code:
        error = ctypes.get_errno()
        if error in (errno.EEXIST, errno.ENOTEMPTY):
            raise FileExistsError(f"output already exists: {destination}")
        raise OSError(error, os.strerror(error), str(destination))


def source_identity(root):
    root = Path(root)
    def git(*args):
        result = subprocess.run(["git", "-C", str(root), *args], capture_output=True, text=True)
        return result.stdout.strip() if result.returncode == 0 else None
    files = sorted((root / "research/cohort-reanalysis").glob("*.py"))
    files += sorted((root / "python/rosalind").glob("*.py"))
    return {"git_revision": git("rev-parse", "HEAD"),
            "git_dirty": bool(git("status", "--porcelain")),
            "python_sources_sha256": {str(p.relative_to(root)): file_digest(p)[0] for p in files}}


def dependency_versions():
    versions = {"python": platform.python_version()}
    for name in ("pyarrow", "blake3", "psutil"):
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            versions[name] = None
    return versions


class ProcessTreeSampler:
    """Simultaneous RSS samples, not a sum of lifetime process high-water marks."""

    def __init__(self, interval=0.01):
        self.interval = interval
        self.peak = None
        self.samples = 0
        self.stop = threading.Event()
        self.thread = None

    def __enter__(self):
        try:
            import psutil
        except ImportError:
            return self
        root = psutil.Process()
        def sample():
            while not self.stop.is_set():
                total = 0
                try:
                    processes = [root, *root.children(recursive=True)]
                except psutil.Error:
                    processes = [root]
                for process in processes:
                    try:
                        total += process.memory_info().rss
                    except psutil.Error:
                        pass
                self.peak = max(self.peak or 0, total)
                self.samples += 1
                self.stop.wait(self.interval)
        self.thread = threading.Thread(target=sample, daemon=True)
        self.thread.start()
        return self

    def __exit__(self, *_):
        self.stop.set()
        if self.thread:
            self.thread.join()

    def report(self):
        return {"sampled_process_tree_peak_rss_bytes": self.peak,
                "rss_sample_interval_seconds": self.interval, "rss_samples": self.samples,
                "rss_scope": "simultaneous parent and descendants; sampled, not an OS cap"}


def output_inventory(staging, max_files, max_bytes):
    inventory = []
    total = 0
    # All run-created directories have admitted fixed/count-bounded shapes.
    for base, dirs, files in os.walk(staging):
        dirs.sort()
        if any((Path(base) / name).is_symlink() for name in dirs):
            raise ValueError("run output cannot contain symlinks to directories")
        for name in sorted(files):
            path = Path(base) / name
            if path.name == "run-manifest.json":
                continue
            if path.is_symlink():
                raise ValueError("run output cannot contain symlinks")
            sha256, size = file_digest(path)
            total += size
            if len(inventory) >= max_files or total > max_bytes:
                raise ValueError("retained output envelope exceeded")
            inventory.append({"path": str(path.relative_to(staging)), "sha256": sha256, "bytes": size})
    return inventory, total


def verify_run(directory):
    """Verify the research manifest and retained outputs without reopening raw data."""
    directory = Path(directory).resolve(strict=True)
    manifest = load_json(directory / "run-manifest.json")
    if manifest.get("schema") != SCHEMA or manifest.get("status") != "completed":
        raise ValueError("not a completed supported research run")
    claimed = manifest.pop("manifest_sha256", None)
    if claimed != digest_json(manifest):
        raise ValueError("research manifest identity mismatch")
    if not isinstance(manifest.get("outputs"), list) or len(manifest["outputs"]) > 20000:
        raise ValueError("invalid output inventory")
    names = set()
    for output in manifest["outputs"]:
        relative = Path(output["path"])
        if relative.is_absolute() or ".." in relative.parts or str(relative) in names:
            raise ValueError("invalid output path")
        names.add(str(relative))
        path = directory / relative
        if path.is_symlink() or not path.resolve(strict=True).is_relative_to(directory):
            raise ValueError("output escapes run directory")
        actual, size = file_digest(path)
        if actual != output["sha256"] or size != output["bytes"]:
            raise ValueError(f"output identity mismatch: {relative}")
    if not {"summary.jsonl", "scientific.json", "plan.json", "coverage.jsonl"} <= names:
        raise ValueError("research run lacks required outputs")
    if digest_json(load_json(directory / "scientific.json")) != manifest.get("scientific_sha256"):
        raise ValueError("scientific metadata identity mismatch")
    manifest["manifest_sha256"] = claimed
    return manifest
