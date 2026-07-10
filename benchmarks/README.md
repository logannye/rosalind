# Rosalind claims harness

**Re-run it — don't trust it.** `bash benchmarks/run.sh` checks five properties Rosalind
claims about itself. Each is an assertion that **fails the run (non-zero exit) if it's
false**, so the result can't be cherry-picked — and it runs in CI as a standing regression
gate, so a broken claim turns the build red.

```sh
bash benchmarks/run.sh     # builds the release binary, runs on the bundled toy data, ~seconds, zero downloads
```

## Scope — read first

This is a **verifiability** harness: the memory contract + byte-reproducibility. It is
**not** a speed benchmark (the engine is single-threaded) and **not** a real-world
accuracy benchmark (calling here is on **bundled, simulated, SNV-only** toy data). The
index *build* is `O(reference)` and out of scope — only the `variants --index` /
`features --index` streaming paths are bounded. Single-threaded, single-process.
Machine-dependent numbers (peak RSS) are *recorded* but the PASS/FAIL is always on the
qualitative property (predicted ≥ realized, exit codes, byte equality, verdicts), never a
fixed MiB figure.

## The seven claims

1. **Predicted peak is a conservative upper bound.** `plan` predicts a peak from the index
   header alone — *before the BAM is read* — and the realized `peak_rss_bytes` in the
   receipt is ≤ it. *Asserts predicted ≥ realized, same `--max-depth`.*
2. **Honor-or-refuse, up front.** `plan`'s verdict flips `fits` → `refuse` as the budget
   shrinks; `variants --enforce` at a too-small budget exits **3** and writes **no** output
   with cooperative assurance; the receipt says exactly which assurance applied.
3. **Byte-identical text output.** Three `variants` runs produce byte-identical VCFs (one
   recorded BLAKE3); two `features` runs produce byte-identical TSVs. *(Deterministic text
   outputs only; BAM/bgzf is reported INCONCLUSIVE by design.)*
4. **Re-derivable + tamper-evident.** `reproduce` re-derives the output byte-for-byte from
   the receipt (exit **0** REPRODUCED); flipping one byte of `manifest_blake3` makes
   `verify` exit **5** (TAMPERED). *(Tamper-evident, not tamper-proof — signing is planned.)*
5. **`pack` is a placement decision, run-free.** It sums additive per-job predicted peaks;
   every node stays within capacity, and an impossible packing refuses (exit **3**).
6. **External analyzers inherit the platform.** A generated analyzer builds, runs through
   the public contract, verifies, reproduces through its explicitly selected binary,
   and localizes a parameter diff.
7. **The offline front door completes.** `rosalind demo --json` finishes the complete
   index → align → sort → call → reproduce provenance journey from embedded assets.

These are *properties an emergent-peak, non-deterministic caller does not provide* — stated
as properties, not as head-to-head speed or accuracy comparisons (no other tool is run).

## Re-run / verify it yourself

Every row prints the exact `rosalind` subcommand it ran. `git clone … && bash
benchmarks/run.sh` regenerates `benchmarks/results.json` (machine-readable; git-ignored
because it records machine-dependent peak numbers — the *claims* it asserts are portable).

Files: `run.sh` (entrypoint) · `claims.py` (stdlib-only driver, no deps) · `results.json`
(generated). Wired into CI as the **Claims harness** job.
