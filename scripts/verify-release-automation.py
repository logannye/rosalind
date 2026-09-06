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


def job(text: str, name: str) -> str:
    match = re.search(r"^  " + re.escape(name) + r":\s*\n(.*?)(?=^  [a-zA-Z0-9_-]+:\s*$|\Z)", text, re.MULTILINE | re.DOTALL)
    if not match:
        fail(f"workflow lacks job {name}")
    return match.group(1)


def verify_pypi_boundary(workflows: dict[str, str]) -> None:
    wheels = workflows["wheels.yml"]
    for forbidden in ("id-token: write", "pypa/gh-action-pypi-publish@", "environment:", "destination:"):
        if forbidden in wheels:
            fail(f"reusable wheel builds must not publish or request publishing privileges: {forbidden}")
    if "wheel-artifacts.py record" not in wheels:
        fail("wheel builds must record the tested candidate wheel bytes")
    for filename, gate, index, url in (
        ("rc.yml", "gates", "testpypi", "https://test.pypi.org/legacy/"),
        ("release.yml", "authorize", "pypi", "https://upload.pypi.org/legacy/"),
    ):
        workflow = workflows[filename]
        if "workflow_call:" in workflow:
            fail(f"{filename}: trusted publisher must be a top-level workflow")
        build, publish = job(workflow, "build-wheels"), job(workflow, "publish-wheels")
        if "uses: ./.github/workflows/wheels.yml" not in build or "id-token: write" in build:
            fail(f"{filename}: reuse candidate-bound wheel builds without OIDC permissions")
        candidate = "${{ needs." + gate + ".outputs.commit }}"
        version = "${{ inputs.version }}" + ("-rc.${{ inputs.number }}" if index == "testpypi" else "")
        for marker in ("candidate_sha: " + candidate, "package_version: " + version):
            if marker not in build:
                fail(f"{filename}: wheel build lacks {marker}")
        for marker in (
            f"needs: [{gate}, build-wheels]", "runs-on: ubuntu-latest",
            "environment: release", "id-token: write", "contents: read",
            "ref: " + candidate, "CANDIDATE_SHA: " + candidate,
            "PACKAGE_VERSION: " + version, "pattern: wheel-*",
            "wheel-artifacts.py stage", "--candidate-sha", "--version",
            "verify-wheel-upload.py --index " + index,
            "repository-url: " + url, "packages-dir: dist", "skip-existing: true",
            "pypa/gh-action-pypi-publish@",
        ):
            if marker not in publish:
                fail(f"{filename}: protected publisher lacks {marker}")
        if publish.index("wheel-artifacts.py stage") > publish.index("pypa/gh-action-pypi-publish@"):
            fail(f"{filename}: validate artifacts before publishing")
        if publish.index("verify-wheel-upload.py") > publish.index("pypa/gh-action-pypi-publish@"):
            fail(f"{filename}: validate registry retries before publishing")
        if "uses: ./.github/workflows/" in publish:
            fail(f"{filename}: OIDC upload cannot be delegated to a reusable workflow")


def main() -> None:
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
        "needs: [authorize, build-assets, publish-crates, published-install-smoke, publish-wheels, publish-container]",
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
        "needs: [gates, build, publish-wheels, publish-container]",
    ):
        if marker not in rc:
            fail(f"RC workflow lacks {marker}")

    for path in sorted(WORKFLOWS.glob("*.yml")):
        text = path.read_text()
        image_publishers = {"container.yml", "happy-image.yml", "rc.yml", "release.yml"}
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

    for parent in (release, rc):
        for child in ("wheels.yml", "container.yml"):
            if f"uses: ./.github/workflows/{child}" not in parent:
                fail(f"release parent must explicitly invoke {child}")
    for name in ("wheels.yml", "container.yml"):
        text = (WORKFLOWS / name).read_text()
        if "workflow_call:" not in text or "candidate_sha:" not in text:
            fail(f"{name}: reusable publication must bind an exact candidate")
        if re.search(r"^  release:", text, re.MULTILINE):
            fail(f"{name}: cannot rely on a GITHUB_TOKEN-created release event")
        if "ref: ${{ inputs.candidate_sha" not in text:
            fail(f"{name}: candidate must control checkout")

    image = (WORKFLOWS / "happy-image.yml").read_text()
    for marker in ("platforms: linux/amd64", "push: true", "environment: release", "create-pull-request", "uses: ./.github/workflows/happy-candidate.yml", "needs: [resolve, candidate]"):
        if marker not in image:
            fail(f"hap.py image workflow lacks {marker}")

    candidate = (WORKFLOWS / "happy-candidate.yml").read_text()
    for marker in ("workflow_call:", "candidate_sha:", "--platform linux/amd64", "happy/smoke.sh", "if: always()"):
        if marker not in candidate:
            fail(f"evaluator candidate validation lacks {marker}")
    if "packages: write" in candidate or "push: true" in candidate or "environment: release" in candidate:
        fail("evaluator candidate builds must run without publication permissions")

    giab = (WORKFLOWS / "giab.yml").read_text()
    if "baseline.json" in giab:
        fail("routine GIAB workflow must not mutate baseline.json")

    verify_pypi_boundary({name: (WORKFLOWS / name).read_text() for name in ("rc.yml", "release.yml", "wheels.yml")})
    print("release workflow dry-run policy: OK (no remote operation performed)")


if __name__ == "__main__":
    main()
