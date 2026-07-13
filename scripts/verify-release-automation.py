#!/usr/bin/env python3
"""Static, non-mutating security checks for Rosalind release workflows."""

from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
WORKFLOWS = ROOT / ".github" / "workflows"
ACTION_FILES = sorted(WORKFLOWS.glob("*.yml")) + [
    ROOT / "assets" / "scaffold" / ".github" / "workflows" / "ci.yml"
]


def fail(message: str) -> None:
    raise SystemExit(message)


for path in ACTION_FILES:
    text = path.read_text()
    for line_number, line in enumerate(text.splitlines(), 1):
        match = re.search(r"\buses:\s*([^\s#]+)", line)
        if not match or match.group(1).startswith("./"):
            continue
        reference = match.group(1).rsplit("@", 1)[-1]
        if not re.fullmatch(r"[0-9a-f]{40}", reference):
            fail(f"{path}:{line_number}: Action is not pinned by full commit SHA")

release = (WORKFLOWS / "release.yml").read_text()
if re.search(r"^\s+push:\s*$", release, re.MULTILINE) or "tags:" in release:
    fail("stable release workflow must never publish from a pushed tag")
for marker in (
    "workflow_dispatch:",
    "plan_id:",
    "environment: release",
    "release-publish.sh",
    "cancel-in-progress: false",
    "needs: [authorize, build-assets, publish-crates, published-install-smoke]",
    "Create stable tag and GitHub Release last",
):
    if marker not in release:
        fail(f"stable release workflow lacks {marker}")

rc = (WORKFLOWS / "rc.yml").read_text()
for marker in (
    "workflow_dispatch:",
    "plan_id:",
    "environment: rc",
    "contract-snapshot.json",
    "cancel-in-progress: false",
    "needs: [gates, build]",
):
    if marker not in rc:
        fail(f"RC workflow lacks {marker}")

for path in sorted(WORKFLOWS.glob("*.yml")):
    text = path.read_text()
    image_publishers = {"container.yml", "happy-image.yml"}
    if "packages: write" in text and path.name not in image_publishers:
        fail(f"{path}: packages:write is reserved for controlled image workflows")
    if "CARGO_REGISTRY_TOKEN" in text and path.name != "release.yml":
        fail(f"{path}: crates.io token is reserved for release.yml")
    if "workflow_dispatch:" in text and any(
        mutation in text
        for mutation in ("packages: write", "git push origin", "create-pull-request")
    ):
        if "plan_id:" not in text:
            fail(f"{path}: mutating manual workflow lacks authenticated plan_id")

image = (WORKFLOWS / "happy-image.yml").read_text()
for marker in ("platforms: linux/amd64", "push: true", "environment: release", "create-pull-request"):
    if marker not in image:
        fail(f"hap.py image workflow lacks {marker}")

giab = (WORKFLOWS / "giab.yml").read_text()
if "baseline.json" in giab:
    fail("routine GIAB workflow must not mutate baseline.json")

print("release workflow dry-run policy: OK (no remote operation performed)")
