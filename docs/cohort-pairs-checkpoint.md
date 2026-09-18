# Paired comparison preview checkpoint

This is an engineering checkpoint on `codex/cohort-pairs-preview`, stacked above
`codex/cohort-preview`. It does not change 0.5 or establish independent paired
research use. The branch is deliberately a draft, with `[skip ci]` commits to
avoid starting broad platform/wheel jobs at this stopping point.

Implemented: explicit ordered pair TSV parsing and scope validation; saved-only
serial pair/window execution; exact side counts, missingness and technical depth
screens; right-minus-left sign/magnitude fractions with uint128 decimal-string
components; bounded canonical Arrow/TSV batches; CLI/Python plans and reports;
guarded pair input and complete consumed lineage; tokenized relocated replay;
user contract, examples and a recurring synthetic Fraction-oracle demonstration.
The single-sample interfaces remain unchanged.

Completed local checks for this checkpoint:

- All 71 cohort library tests passed, including the managed artifact failure and
  pair-table mutation/overwrite guards. The core implementation commit is
  `1988363`.
- All-target Clippy passed with warnings denied.
- Both new CLI integration tests passed: explicit direction, missing/zero/low
  observations, three admitted budgets and window widths, byte equality in
  Arrow/TSV, invalid tables/scope, relocated replay without original alignments,
  and refusal after changing the pair table.
- The Python paired test passed against the actual development binary in the
  existing isolated 0.6 environment: native/Python Arrow bytes match, uint64
  counts and string differences retain exact arithmetic and true nulls.
- Python syntax, actionlint and whitespace checks passed.

The integration tests reused compiled dependencies and linked only the development
binary/test target. No new release binary, platform build or wheel was produced.
The current source build is a development test artifact; it is not immutable
commit-bound release evidence.

Resume with these bounded remaining tasks:

1. Refresh the stacked base with the reviewed Git-worktree build-info correction
   from PR #147 (or its merged successor). Rebuild from a clean committed tree,
   verify its producer SHA/dirty fields and retain source identity with results.
2. Execute `run_cohort.py` followed by `run_pairs.py` using a stable copy of that
   binary. Retain the passed or failed reports. The new paired demonstration is
   wired into CI but has not yet been executed locally or on Linux CI.
3. Run the affected replay/unit regressions and full workspace/Python suites on
   that final integrated source. No full workspace or full Python rerun was made
   for this stopping checkpoint; the focused paired paths above did pass.
4. Exercise the supported compiler/platform matrix and matching-wheel install.
   The code avoids trait upcasting to preserve the existing MSRV, but a separate
   MSRV build and clean wheel validation have not been run for these additions.
5. Obtain technical/scientific review before removing draft status. Record real
   paired research use if it occurs; absent external participation is advisory,
   and authored examples must never be counted as independent user value.

The [paired contract and tutorial](cohort-pairs-contract.md) describe the exact
interpretation. Nothing in this checkpoint expands Parquet replay, infers pairing
from metadata, introduces a new biological counting model, or declares an OS
memory guarantee without measured Linux evidence.
