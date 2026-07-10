#!/usr/bin/env bash
set -euo pipefail

VERSION="${1:?usage: scripts/post-publish-smoke.sh VERSION}"
case "$VERSION" in *[!0-9A-Za-z.+-]*) echo "invalid version: $VERSION" >&2; exit 2 ;; esac

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
export CARGO_HOME="$WORK/cargo-home"
export CARGO_TARGET_DIR="$WORK/cargo-target"
mkdir -p "$CARGO_HOME" "$CARGO_TARGET_DIR"

cargo install rosalind-bio --version "=$VERSION" --locked --root "$WORK/install"
ROSALIND="$WORK/install/bin/rosalind"
"$ROSALIND" --version
"$ROSALIND" demo --output-dir "$WORK/demo" --json > "$WORK/demo.json"
python3 - "$WORK/demo.json" <<'PY'
import json, sys
report = json.load(open(sys.argv[1]))
if report.get("ok") is not True or report.get("trust") != "reproduced":
    raise SystemExit(f"unexpected demo report: {report}")
PY

"$ROSALIND" new analyzer release-smoke --output "$WORK/release-smoke"
(
  cd "$WORK/release-smoke"
  cargo generate-lockfile
  cargo build --release --locked
)
ANALYZER="$CARGO_TARGET_DIR/release/release-smoke"
test -x "$ANALYZER"
"$ROSALIND" conformance analyzer --binary "$ANALYZER" --json > "$WORK/conformance.json"
python3 - "$WORK/conformance.json" <<'PY'
import json, sys
report = json.load(open(sys.argv[1]))
if report.get("passed") is not True or report.get("kind") != "rosalind-analyzer-conformance":
    raise SystemExit(f"unexpected conformance report: {report}")
PY

echo "published rosalind-bio $VERSION passed fresh-cache demo, scaffold, and conformance smoke tests"
