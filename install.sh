#!/bin/sh
# Rosalind installer — download a prebuilt binary for this OS/arch and unpack it.
#
#   curl -fsSL https://raw.githubusercontent.com/logannye/rosalind/main/install.sh | sh
#
# Installs into ./rosalind-<target>/ in the current directory (no sudo, no global
# state). Prefers `gh release download` when the GitHub CLI is available (works for
# private repos); falls back to `curl` against the public release assets.
set -eu

REPO="logannye/rosalind"
VERSION="${ROSALIND_VERSION:-latest}"

say() { printf '%s\n' "$*"; }
die() { printf 'install.sh: %s\n' "$*" >&2; exit 1; }

# --- detect OS/arch -> release target triple ---------------------------------
os="$(uname -s)"
arch="$(uname -m)"
case "$os" in
  Linux)
    case "$arch" in
      x86_64 | amd64) target="x86_64-unknown-linux-musl" ;;
      aarch64 | arm64) target="aarch64-unknown-linux-musl" ;;
      *) die "unsupported Linux arch '$arch' (prebuilt: x86_64, aarch64; build from source for others)" ;;
    esac
    ;;
  Darwin)
    case "$arch" in
      arm64 | aarch64) target="aarch64-apple-darwin" ;;
      x86_64) target="x86_64-apple-darwin" ;;
      *) die "unsupported macOS arch '$arch'" ;;
    esac
    ;;
  *)
    die "unsupported OS '$os' (build from source: see the README)"
    ;;
esac

tarball="rosalind-${target}.tar.gz"
say "Rosalind installer: detected ${os}/${arch} -> ${target} (version: ${VERSION})"

# --- download ----------------------------------------------------------------
if command -v gh >/dev/null 2>&1; then
  say "Downloading via gh (works for private repos)…"
  if [ "$VERSION" = "latest" ]; then
    gh release download --repo "$REPO" --pattern "$tarball" --clobber \
      || die "gh release download failed (is a release published?)"
  else
    gh release download "$VERSION" --repo "$REPO" --pattern "$tarball" --clobber \
      || die "gh release download failed for $VERSION"
  fi
elif command -v curl >/dev/null 2>&1; then
  if [ "$VERSION" = "latest" ]; then
    url="https://github.com/${REPO}/releases/latest/download/${tarball}"
  else
    url="https://github.com/${REPO}/releases/download/${VERSION}/${tarball}"
  fi
  say "Downloading ${url} …"
  curl -fsSL -o "$tarball" "$url" \
    || die "download failed — no release yet, or the repo is private (install gh and retry)"
else
  die "need either 'gh' or 'curl' to download"
fi

# --- verify checksum ---------------------------------------------------------
# A "verifiable" tool must verify its own download. The release ships a
# `<tarball>.sha256` sidecar; fetch it and check before unpacking. Abort on a
# mismatch (corrupt or tampered download); only skip when neither a hash tool
# nor the sidecar is available, and say so loudly.
if command -v sha256sum >/dev/null 2>&1; then sha="sha256sum"
elif command -v shasum >/dev/null 2>&1; then sha="shasum -a 256"
else sha=""; fi

if [ -n "$sha" ]; then
  sumfile="${tarball}.sha256"
  if command -v gh >/dev/null 2>&1; then
    if [ "$VERSION" = "latest" ]; then
      gh release download --repo "$REPO" --pattern "$sumfile" --clobber >/dev/null 2>&1 || sumfile=""
    else
      gh release download "$VERSION" --repo "$REPO" --pattern "$sumfile" --clobber >/dev/null 2>&1 || sumfile=""
    fi
  elif command -v curl >/dev/null 2>&1; then
    if [ "$VERSION" = "latest" ]; then
      sumurl="https://github.com/${REPO}/releases/latest/download/${sumfile}"
    else
      sumurl="https://github.com/${REPO}/releases/download/${VERSION}/${sumfile}"
    fi
    curl -fsSL -o "$sumfile" "$sumurl" 2>/dev/null || sumfile=""
  fi
  if [ -n "$sumfile" ] && [ -s "$sumfile" ]; then
    say "Verifying checksum…"
    $sha -c "$sumfile" \
      || die "checksum verification FAILED — refusing to install a corrupt or tampered download"
    rm -f "$sumfile"
  else
    say "warning: could not fetch the .sha256 sidecar — skipping checksum verification"
  fi
else
  say "warning: no sha256 tool found (sha256sum/shasum) — skipping checksum verification"
fi

# --- unpack ------------------------------------------------------------------
tar -xzf "$tarball"
rm -f "$tarball"
dir="rosalind-${target}"
[ -x "$dir/rosalind" ] || die "unpack failed: $dir/rosalind not found"

say ""
say "Installed: ./$dir/rosalind"
say ""
say "Try the memory contract on the bundled data (≈60 seconds):"
say "  cd $dir"
say "  ./rosalind index --reference examples/data/illumina_toy/reference.fa --output ref.idx"
say "  ./rosalind sort  --input examples/data/illumina_toy/alignments.bam --output sorted.bam"
say "  ./rosalind plan  --index ref.idx --budget-mb 512"
say "  ./rosalind variants --index ref.idx --alignments sorted.bam \\"
say "      --memory-budget-mb 512 --enforce -o calls.vcf"
say "  ./rosalind verify --manifest calls.vcf.manifest.json"
