# Adoption On-Ramp (design)

**Status:** DESIGN SPEC — 2026-06-02. **Branch:** `rosalind/adoption-on-ramp` (off `main` `165a6f1`).
From the reflection audit: the differentiator (the live `plan → enforce → verify` contract) is a
60-second demo, but it is **unrunnable by strangers** today — no releases/tags/containers, only a
from-source build behind a C-toolchain gate. Scope (chosen for max tangible value to users +
forkers): **release core + the `rosalind-budget` GitHub Action**; container deferred.

## 1. Goal

Turn "clone → fight the htslib C-toolchain build" into "**run the memory contract in 60 seconds**,"
and give builders a drop-in way to **enforce the contract in their own CI**.

## 2. Honest constraints (stated up front)

- **CI-only verification.** The release workflow and the Action execute on GitHub runners; the Linux
  **musl-static** build in particular cannot be run-verified locally (cross-compiling htslib's C from
  macOS isn't feasible here). This PR delivers *reviewed* release infrastructure + a **locally-verified
  macOS quickstart**; the Linux build is validated when CI first runs.
- **The `v0.1.0` tag is the user's trigger.** `release.yml` fires on a `v*` tag push (plus
  `workflow_dispatch` for a manual dry run). Merging this PR releases nothing — cutting the tag (the
  outward-facing step) stays the user's call.
- **The Action activates after the first release** (it downloads the release binary). Until `v0.1.0`
  exists it is non-functional infrastructure; documented as such.

## 3. Deliverables

### 3a. `.github/workflows/release.yml`

- **Triggers:** `push: tags: ['v*']` and `workflow_dispatch` (manual, for a dry run).
- **Matrix (native runners, no cross-compile):**
  - `x86_64-unknown-linux-musl` on `ubuntu-latest` — **static** (apt `musl-tools`; `rustup target add`;
    the `-sys` crates build htslib/bzip2/lzma/zlib from bundled C, so static linking should succeed).
  - `aarch64-apple-darwin` on `macos-14` (native arm64).
  - `x86_64-apple-darwin` on `macos-15-intel` (native x86_64).
- **Bundle** per target: the `rosalind` binary + `examples/data/illumina_toy/` (the contract-demo
  fixtures) + `README.md` + `CONTRACT.md` + `LICENSE-APACHE` + `LICENSE-MIT`, into
  `rosalind-<target>.tar.gz` (+ a `.sha256`).
- **Publish:** attach the tarballs to the GitHub Release for the tag (via `softprops/action-gh-release@v2`).
  On `workflow_dispatch` (no tag), build + upload as workflow artifacts only (no release).
- Validate YAML structure locally (`python -c yaml.safe_load` or a manual indentation review); cannot
  run the build here.

### 3b. `install.sh` (curl-pipe one-liner)

Detect `uname -s`/`uname -m` → target triple; download
`https://github.com/logannye/rosalind/releases/latest/download/rosalind-<target>.tar.gz`; verify the
`.sha256`; extract to `./rosalind-<version>/`; print the next-step quickstart. Prefer `gh release
download` when `gh` is available (works even if the repo is private); fall back to `curl`. Fail
loudly with a clear message on an unsupported OS/arch. Locally testable: the OS/arch-detection +
error paths (the actual download needs a published release).

### 3c. README "Quickstart (60 seconds)"

A new section near the top: `curl -fsSL …/install.sh | sh` → `cd` → run the bundled
`plan → variants --enforce → verify` demo on `examples/data/illumina_toy/` (the exact commands already
verified locally end-to-end). Keep the existing from-source "Install & build" section below it for
contributors.

### 3d. `rosalind-budget` GitHub Action (`action.yml` at repo root)

A **composite** action: a consumer adds
```yaml
- uses: logannye/rosalind-budget@v1
  with:
    index: ref.idx
    alignments: sorted.bam
    budget-mb: 4096
    max-depth: 1000
    max-read-len: 250
```
Steps: download the Linux release binary (the action runs on the consumer's Linux runner; no
toolchain needed), run `rosalind variants --index … --alignments … --memory-budget-mb … --enforce
--max-depth … --max-read-len … -o calls.vcf`, so a budget breach **fails the consumer's CI** (exit 3
refuse / exit 4 breach), and upload `calls.vcf.manifest.json` as a build artifact. Inputs:
`index`, `alignments`, `budget-mb` (required), `max-depth` (default 1000), `max-read-len`
(default 250), `version` (default `latest`). Ship:
- `action.yml` (composite).
- `.github/workflows/example-rosalind-budget.yml` — a self-contained example that builds the binary
  from source (so the example runs in *this* repo's CI before any release exists), builds the bundled
  toy index, and runs the same plan→enforce→verify gate — proving the pattern end-to-end in CI now.
- A README/`docs` section documenting the Action and the "memory contract in your CI" value.

## 4. Out of scope (follow-up)

- **Container image** (Dockerfile + ghcr publish) — a second channel for the same binary; valuable for
  HPC/Singularity users, deferred.
- **crates.io publish** — a library-consumer channel; separate.
- Cutting the actual `v0.1.0` tag/release (the user's trigger).

## 5. Self-review

- **Coverage:** release.yml (3a), install.sh (3b), README quickstart (3c), Action + example + docs
  (3d). ✓
- **No placeholders:** the musl-static build and the Action's release-download are concrete; the honest
  caveat is that they run on CI, not locally. ✓
- **Verification:** macOS quickstart + install.sh detection + YAML structure are checked locally; the
  Linux build + the release publish + the Action's download are CI-only (stated). ✓
- **Boundary:** no tag cut; the example workflow proves the Action's *pattern* in this repo's CI today
  (build-from-source), decoupled from the not-yet-existent release. ✓
