#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
TAG="${1:-rosalind-happy:0.3.15}"
command -v docker >/dev/null || { echo "docker is required" >&2; exit 2; }
docker build --platform linux/amd64 --pull=false -t "$TAG" -f "$HERE/Dockerfile" "$HERE"
image_id="$(docker image inspect "$TAG" --format '{{.Id}}')"
case "$image_id" in sha256:????????????????????????????????????????????????????????????????) ;; *) echo "unexpected image digest: $image_id" >&2; exit 2 ;; esac
echo "local_image_config_id=$image_id"
echo "This is a local image config ID, not a published OCI manifest digest."
"$HERE/smoke.sh" "$TAG" "${2:-$(mktemp -d)}"
echo "Publication must record the registry digest and attestation before GIAB evaluation."
