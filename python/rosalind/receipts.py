"""Typed offline receipt inspection."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Union


@dataclass(frozen=True)
class Receipt:
    """Portable claims and machine-local measurements from one run."""

    path: Path
    tool_version: str
    subcommand: str
    inputs: tuple[dict[str, str], ...]
    outputs: tuple[dict[str, str], ...]
    params: dict[str, str]
    measurements: dict[str, str]
    raw: dict[str, Any]


def inspect_receipt(path: Union[str, Path]) -> Receipt:
    """Parse a canonical Rosalind receipt without uploading any artifact."""
    path = Path(path)
    raw = json.loads(path.read_text(encoding="utf-8"))
    return Receipt(
        path=path,
        tool_version=raw["tool_version"],
        subcommand=raw["subcommand"],
        inputs=tuple(raw.get("inputs", [])),
        outputs=tuple(raw.get("outputs", [])),
        params=dict(raw.get("params", {})),
        measurements=dict(raw.get("measurements", {})),
        raw=raw,
    )
