#!/usr/bin/env python3
"""Report whether the opt-in scientific evaluator exists; never fabricate a score."""
import argparse
import json
import re
from pathlib import Path


def readiness(lock: dict, override: str = "") -> dict:
    image = lock.get("generated_image", {})
    selected = override or (f"{image.get('repository')}@{image.get('digest')}" if image.get("digest") else "")
    ready = bool(re.fullmatch(r"[^\s@]+@sha256:[0-9a-f]{64}", selected))
    if not override:
        ready = ready and image.get("platform") == "linux/amd64" and bool(
            re.fullmatch(r"[0-9a-f]{40}", image.get("built_from_commit") or "")
        )
    return {
        "schema": 1,
        "status": "ready" if ready else "blocked",
        "evaluation_performed": False,
        "image": selected if ready else None,
        "reason": "evaluator available; scientific run still required" if ready else
            "publish and review the pinned hap.py image lock before running GIAB",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lock", default="benchmarks/giab/happy/lock.json")
    parser.add_argument("--image", default="")
    parser.add_argument("--output", default="giab-readiness.json")
    args = parser.parse_args()
    report = readiness(json.loads(Path(args.lock).read_text()), args.image)
    Path(args.output).write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(report["status"])


if __name__ == "__main__":
    main()
