#!/usr/bin/env bash
set -euo pipefail

# Resumable crates.io publication. This script deliberately owns no tag or
# GitHub Release mutation; the protected workflow performs those only after the
# published packages pass fresh-cache installation tests.

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
STATE="${RELEASE_PUBLISH_STATE:-$ROOT/publish-state.json}"
REGISTRY="${CARGO_REGISTRY_API:-https://crates.io/api/v1}"
USER_AGENT="rosalind-release-automation/1 (https://github.com/logannye/rosalind)"
PACKAGES=(rosalind-build-info rosalind-receipt rosalind-bio)

test -n "${CARGO_REGISTRY_TOKEN:-}" || {
  echo "CARGO_REGISTRY_TOKEN is required" >&2
  exit 3
}
for tool in cargo curl python3; do
  command -v "$tool" >/dev/null || { echo "missing required tool: $tool" >&2; exit 2; }
done

cd "$ROOT"
mkdir -p "$(dirname "$STATE")"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT
printf '{"schema":1,"packages":{}}\n' > "$work/state.json"

package_version() {
  cargo metadata --format-version 1 --no-deps | python3 -c '
import json, sys
name = sys.argv[1]
for package in json.load(sys.stdin)["packages"]:
    if package["name"] == name:
        print(package["version"])
        raise SystemExit(0)
raise SystemExit(f"workspace package not found: {name}")
' "$1"
}

remote_checksum() {
  local package="$1" version="$2"
  curl --fail --silent --show-error \
    -H "User-Agent: $USER_AGENT" \
    "$REGISTRY/crates/$package/$version" |
    python3 -c 'import json,sys; print(json.load(sys.stdin)["version"]["checksum"])'
}

downloaded_checksum() {
  local package="$1" version="$2" destination="$3"
  curl --fail --silent --show-error --location \
    -H "User-Agent: $USER_AGENT" \
    "$REGISTRY/crates/$package/$version/download" -o "$destination"
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$destination" | awk '{print $1}'
  else
    shasum -a 256 "$destination" | awk '{print $1}'
  fi
}

local_checksum() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$1" | awk '{print $1}'
  else
    shasum -a 256 "$1" | awk '{print $1}'
  fi
}

record_state() {
  local package="$1" version="$2" checksum="$3" disposition="$4"
  python3 - "$work/state.json" "$package" "$version" "$checksum" "$disposition" <<'PY'
import json, os, sys
path, package, version, checksum, disposition = sys.argv[1:]
state = json.load(open(path))
state["packages"][package] = {
    "version": version,
    "sha256": checksum,
    "disposition": disposition,
}
tmp = path + ".new"
with open(tmp, "w") as handle:
    json.dump(state, handle, indent=2, sort_keys=True)
    handle.write("\n")
os.replace(tmp, path)
PY
  cp "$work/state.json" "$STATE.tmp"
  mv "$STATE.tmp" "$STATE"
}

wait_for_registry() {
  local package="$1" version="$2" expected="$3"
  local delays=(5 10 20 40 60 90 120 180)
  local delay observed
  for delay in "${delays[@]}"; do
    if observed="$(remote_checksum "$package" "$version" 2>/dev/null)"; then
      test "$observed" = "$expected" || {
        echo "published checksum conflict for $package $version: expected $expected, registry has $observed" >&2
        exit 5
      }
      return 0
    fi
    sleep "$delay"
  done
  echo "$package $version did not become visible before the bounded timeout" >&2
  exit 4
}

for package in "${PACKAGES[@]}"; do
  version="$(package_version "$package")"
  cargo package -p "$package" --locked
  archive="target/package/$package-$version.crate"
  test -f "$archive" || { echo "cargo package did not create $archive" >&2; exit 5; }
  expected="$(local_checksum "$archive")"

  if registry_sum="$(remote_checksum "$package" "$version" 2>/dev/null)"; then
    remote_archive="$work/$package-$version.crate"
    downloaded="$(downloaded_checksum "$package" "$version" "$remote_archive")"
    if [ "$registry_sum" != "$downloaded" ] || [ "$expected" != "$downloaded" ]; then
      echo "published crate conflict for $package $version: local=$expected api=$registry_sum downloaded=$downloaded" >&2
      exit 5
    fi
    echo "$package $version already exists with the exact package bytes; resuming"
    record_state "$package" "$version" "$expected" "verified-existing"
    continue
  fi

  echo "publishing $package $version"
  if ! cargo publish -p "$package" --locked; then
    # A concurrent or interrupted run may have completed the upload. Re-query
    # and accept only byte-identical package contents.
    registry_sum="$(remote_checksum "$package" "$version" 2>/dev/null || true)"
    test "$registry_sum" = "$expected" || exit 4
  fi
  wait_for_registry "$package" "$version" "$expected"
  downloaded="$(downloaded_checksum "$package" "$version" "$work/published-$package.crate")"
  test "$downloaded" = "$expected" || {
    echo "downloaded bytes differ after publishing $package $version" >&2
    exit 5
  }
  record_state "$package" "$version" "$expected" "published"
done

echo "all workspace crates are published and byte-verified"
