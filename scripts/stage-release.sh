#!/usr/bin/env bash
set -euo pipefail

TARGET="${1:?usage: scripts/stage-release.sh TARGET}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

NAME="rosalind-$TARGET"
STAGE="dist/$NAME"
rm -rf "$STAGE"
mkdir -p "$STAGE/docs/schema" "$STAGE/examples/data"
cp "target/$TARGET/release/rosalind" "$STAGE/"
cp README.md CONTRACT.md CHANGELOG.md LICENSE-APACHE LICENSE-MIT "$STAGE/"
cp docs/receipt-schema.md docs/receipt-trust.md docs/v0.4-migration.md "$STAGE/docs/"
cp docs/schema/receipt-v5.schema.json docs/schema/trust-report-v2.schema.json "$STAGE/docs/schema/"
cp web/verify/schema/run-attestation-v1.schema.json "$STAGE/docs/schema/"
cp -R examples/data/illumina_toy "$STAGE/examples/data/"
# Keep the top-level README and its linked SDK, Python, semantics, workflow and
# validation guides usable after extracting the archive outside the checkout.
python3 scripts/onboarding.py stage "$ROOT" "$STAGE"
tar -C dist -czf "dist/$NAME.tar.gz" "$NAME"
if command -v sha256sum >/dev/null 2>&1; then
  (cd dist && sha256sum "$NAME.tar.gz" > "$NAME.tar.gz.sha256")
else
  (cd dist && shasum -a 256 "$NAME.tar.gz" > "$NAME.tar.gz.sha256")
fi
