#!/usr/bin/env bash
# Check a local image; no registry credentials, publication or network access.
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "usage: $0 IMAGE EXPECTED_VERSION NEW_OUTPUT_DIRECTORY" >&2
  exit 2
fi
image=$1
expected_version=$2
output=$3
mkdir "$output"

docker image inspect "$image" > "$output/image.json"
runtime=(--rm --network none --cap-drop ALL --security-opt no-new-privileges)
docker run "${runtime[@]}" --entrypoint id "$image" -u > "$output/uid.txt"
test "$(cat "$output/uid.txt")" = 10001
docker run "${runtime[@]}" "$image" --version > "$output/version.txt"
test "$(cat "$output/version.txt")" = "rosalind $expected_version"
docker run "${runtime[@]}" "$image" analyze evidence --help > "$output/evidence-help.txt"
docker run "${runtime[@]}" "$image" dataset --help > "$output/dataset-help.txt"
# The embedded demo exercises native decoding, managed artifacts, verification,
# and byte replay. It is a packaging smoke, not new biological validation.
docker run "${runtime[@]}" "$image" demo --output-dir /tmp/rosalind-demo --json \
  > "$output/demo.json" 2> "$output/demo.stderr.txt"
python3 - "$output" <<'PY'
import json
import pathlib
import sys

directory = pathlib.Path(sys.argv[1])
image = json.loads((directory / "image.json").read_text())[0]
assert image["Config"]["User"] == "10001:10001", image["Config"]["User"]
assert image["Os"] == "linux" and image["Architecture"] == "amd64"
assert json.loads((directory / "demo.json").read_text())["ok"] is True
print("container smoke: linux/amd64; default UID10001; offline demo verified and replayed")
PY
