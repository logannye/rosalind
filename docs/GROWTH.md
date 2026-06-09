# Rosalind — Production & Growth Roadmap

**Status:** Living growth/GTM plan — 2026-06-09. Companion to [`ROADMAP.md`](ROADMAP.md)
(engineering sequence) and [`OPEN_PROBLEMS.md`](OPEN_PROBLEMS.md) (the research thesis). Where
those order work by *capability*, this orders it by *adoption* for the post-Hacker-News moment.

**Objective:** convert attention into adopters and contributors. **Bandwidth:** solo, near-full-time.
**Horizon:** ~1 quarter.

---

## The reframe

The install path **already works end-to-end**: `v0.1.0` is published with checksum-verified binaries
for macOS arm64/x86_64 and Linux x86_64, `install.sh` resolves, the pinned `@v0.1.0` Action resolves,
CI is green. So the highest-ROI work is **presentation + distribution channels + a contributor
on-ramp — not engineering plumbing.**

Optimize for **installs → returning users → first external PRs**, *not* stars. The telling gap today
is 251 stars / 11 forks / **0 external PRs** — closing that is the whole game.

These map onto [`ROADMAP.md`](ROADMAP.md): Wave 1's demos ≈ **P1.2**, Wave 2's benchmark ≈ **P1.3**,
Wave 3 ≈ **P2.1/P2.2**, Wave 4 ≈ **P4.1**. This document sequences them for conversion.

---

## Wave 0 — Stop the leaks (this week · all hours-scale · do first)

Pure ROI on traffic you already have. README + GitHub settings, no code.

- **0.1** Delete the stale "build from source / available once a release is published" hedges — the
  release works; say so. *(done 2026-06-09)*
- **0.2** Above-the-fold makeover: badge row (CI · release · license · byte-reproducible), the
  one-sentence differentiator, the install one-liner, and a 20–30s asciinema of
  `plan → enforce → verify → reproduce`. *(badges + install one-liner done; asciinema = a tracked issue)*
- **0.3** Link `CONTRIBUTING.md` from the README (a Contributing section). *(done 2026-06-09)*
- **0.4** Enable GitHub Discussions *(done)*; seed 8–15 labelled `good first issue` tickets.
- **0.5** Sharpen the repo About description *(done)*; add a "Jump to" nav atop the README *(done)*.
- **0.6** Be present: a 24h ack-SLA on issues/PRs for the next few weeks.

## Wave 1 — Own the install channels + ship the two viral demos (1–2 weeks)

Meet every audience where they live.

- **1.1 crates.io** — resolve the name squat (the bare `rosalind` crate is a dormant 2016 package):
  publish under a distinct package name, keep `[[bin]] name = "rosalind"`. Unlocks `cargo install`,
  auto-built docs.rs, and `cargo binstall`. Also permanently locks the name.
- **1.2 PyPI wheel** for the Python boundary (maturin/abi3) → `pip install`. Highest *end-user* reach:
  genomics is Python/notebook-first (the polars/sourmash playbook).
- **1.3 Demo A — "Caught you" WASM verifier**: drag a receipt into a browser → VERIFIED; flip one byte
  → TAMPERED, naming the drifted field; 100% client-side. `provenance/mod.rs` is already zero-htslib,
  so it compiles to `wasm32` with no surgery. (ROADMAP P1.2.)
- **1.4 Demo B — the reproduce duel**: split-screen asciinema, `rosalind reproduce` → REPRODUCED
  byte-identical offline vs. a non-deterministic caller that would DIVERGE on a *correct* run.
- **1.5** Add `aarch64-unknown-linux-musl` (and optionally Windows) release targets + `install.sh` arch
  detection — an "edge/field" tool must install on ARM Linux.
- **1.6** List the budget Action on the **GitHub Marketplace** + a one-paste "reproducibility receipt
  gate for your genomics CI" recipe — the stickiest distribution surface.

## Wave 2 — The second-act re-launch + credibility (weeks 3–5)

A single spike decays; a cadence compounds.

- **2.1 Reproducible benchmark/claims harness** (`benchmarks/`, one command, pinned inputs), framed on
  **verifiability** (memory ceiling honored under load; receipt reproduces byte-for-byte), **not**
  throughput, and honestly scoped (simulated, SNV-only). (ROADMAP P1.3.)
- **2.2 Technical deep-dive post** on how the verifiable memory contract + content-addressed receipts
  work, plus a differentiator-first **landing page** with the demo GIF and an explicit honest-scope
  section.
- **2.3 Timed re-launch**: a follow-up *"Show HN"* + cross-posts (r/rust, r/bioinformatics, lobste.rs)
  riding the new demos + crates.io/PyPI.

## Wave 3 — Earn the credibility skeptics named (weeks 4–8, parallel)

- **3.1 A real GIAB HG002 number** (chr20 + CMRG), paired with the receipt — the genotype-aware
  comparator is ready; run it on real externally-mapped reads, framed honestly (SNV-competitive,
  indels trailing). Retires the "simulated only" critique. (ROADMAP P2.1.)
- **3.2 MAPQ into the germline likelihood** (today it is ignored). (ROADMAP P2.2.)

## Wave 4 — Begin the moat, de-risked (weeks 6–12, research track, off the critical path)

- **4.1 √t Phase-D spike (D1a)** — the budget-tunable external-memory blocked SA/BWT build, run **with
  the kill criterion from day one**. Payoff: build a T2T human index on a 32 GB box where `bwa-index`
  OOM-kills, **byte-identical** to the full-RAM build. Completes the contract end to end. (ROADMAP P4.1.)

---

## Explicitly NOT now

- ❌ Lead with **speed** (single-threaded; invites a benchmark fight that buries the real moat).
- ❌ Add `--threads` as a perf play — byte-identity is a *correctness* property; only revisit threading
  as a prerequisite for the D1a build, behind a byte-identity CI gate.
- ❌ Spin up a Discord (Discussions first — an empty chat at peak reads as dead).
- ❌ Position methylation/indel/SV tracks as the moat (they widen reach, not defensibility).
- ❌ Chase stars; don't build the economics on a regulated buyer a solo, pre-1.0 repo can't transact with.

## What to measure (not stars)

crates.io + PyPI + release-download counts · `cargo install` / `pip install` · Action adoption ·
**first external PRs merged** · Discussions activity · returning visitors.

---

> **One-paragraph thesis:** the release works but is invisible and single-channel; the fastest value
> creation is to surface the differentiator visually, open every install channel (cargo/pip/ARM), and
> ship the two demos that make "verifiable memory contract + byte-reproducibility" undeniable — then
> re-launch on that, and only then spend bandwidth earning GIAB credibility and starting the √t moat.
