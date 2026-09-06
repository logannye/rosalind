# Adoption and reusable evidence implementation

This sequence follows the evidence-engine foundation merged in PR #97. Work remains
ordered; external publication and adoption prerequisites are tracked separately
from code and local verification. The original scientific invariant and historical
artifact verification remain requirements throughout.

| Order | Delivery | Status |
|---|---|---|
| 1 | Shared panel defaults and explicit named/unknown/pooled sample scope | Implemented; 11 sample tests, Rust/CLI threshold parity, Python/native parity, cache/replay and existing evidence regressions pass |
| 2 | Installable evidence-engine RC, executable onboarding, candidate publication gates | Implemented for 0.5.0-rc.1; native bundle and fresh wheel onboarding pass locally; all 16 PR checks pass and PR #98 merged; immutable RC build dispatched at ce6c489; protected publication pending |
| 3 | Physical field projection and bounded coalesced sparse fetches | Implemented; all 64 masks, frozen full-output hashes, cache/replay/workers, and corrected 18-run pressure matrix verified locally; all platform PR gates passed and PR #99 merged |
| 4 | VCF.gz/BCF selection, record-preserving annotation, optional per-allele quality/position sums | Implemented; all 128 masks, legacy byte goldens, independent allele sums, record/genotype preservation and three-budget/worker replay pass locally; 514-test workspace suite plus focused follow-up checks and Python 3.9/3.11 pass locally; all platform PR gates passed and PR #101 merged |
| 5 | Verified subset/superset dataset reuse, persisted analysis, bounded Parquet and Python/R/SQL adapters | Implemented; full workspace/Clippy/MSRV and Python3.9/3.11 pass; Linux native/R and DuckDB checks pass; fresh wheel installs, dataset queries/export and onboarding pass. All 15 platform PR gates passed; PR #102 merged |
| 6 | Evidence-native artifact runner, scaffold, and conformance | Implemented; full workspace, Clippy/MSRV, 19 external conformance checks, focused lifecycle/cancellation tests and fresh Rust 1.83 projects pass; final packaged/platform validation in progress |
| 7 | Representative BAM/CRAM resource measurements and independent user validation | Pending; external users have not been recruited |

## Acceptance

1. Default panel results agree through Rust, CLI, and Python at threshold
   boundaries. Samples cannot silently mix: scope is explicit in receipts,
   scientific identity, and cache identity. Selection/pooling replay correctly.
2. One immutable candidate builds matching binaries, wheels, crates, and OCI.
   Clean installations run the researcher tutorial, verification, and relocated
   replay; external SDK builds use published dependencies without source patches.
   Preserve protected review, credentials, relevant caller evidence, and the
   seven-day RC soak for releases from 0.5 onward.
3. Omitted capabilities allocate/encode no corresponding buffers. Present values
   equal full evidence, with bytes invariant across admitted execution settings.
   Nearby sparse loci share bounded fetches without widening scientific selection.
4. Equivalent VCF encodings yield equal selected evidence; annotation preserves
   records and allele correspondence. Unsupported variants and REF mismatches
   remain explicit. Per-allele quality swaps are distinguishable with oracle
   agreement and a declared memory model.
5. Compatible persisted evidence answers subset requests without alignment
   decoding. Incompatible filters/capabilities never reuse insufficient summaries.
   Each reused artifact retains its producer, original byte identity, and verified
   lineage. Export adapters have bounded buffers.
6. A standalone analyzer implements its reducer/output while inheriting refusal,
   cancellation, atomic publication, receipts, and relocated replay.
7. High-depth panel and chromosome-scale sparse workloads vary effective tiles,
   budgets, workers, input encoding, and cold/resumed cache state. Retain three
   repetitions, correctness comparisons, full time/RSS/I/O costs, and new-engine
   cgroup checks. Three non-authors complete real tasks; two teams return after
   30 days. Record installation time, integration effort, and work avoided.

## Publication prerequisites inspected 2026-09-06

The `rc` and `release` environments retain required review and currently allow
zero deployment refs. No repository/environment secret names are configured.
PyPI/TestPyPI trusted-publisher mappings cannot be confirmed from their public
project APIs. These prerequisites have not been changed or waived.

The latest public release remains v0.1.0. Candidate builds and authored integration
examples are not registry publication or independent adoption.

## Candidate preparation evidence

Source version 0.5.0 is synchronized across Cargo, Python, release policy, and
standalone examples. The local workspace suite passed 483 tests; the Python source
suite passed 11 tests on fresh Python 3.9 and 3.11 environments. A native archive
and Python 3.11 wheel were installed outside the checkout and completed documented
analysis, byte verification, offline replay, scaffold conformance, and the bundled
evidence reducer example. Candidate SDK checks used explicit source patches;
registry-only SDK installation remains unverified until its crates are published.

Publication workflows bind platform reports to the exact tested wheel bytes and
candidate commit, and use top-level protected OIDC upload jobs. GIAB preparation
validates its indexed reference and input hashes; evaluator source version pinning
is regression-tested. The full evaluator container requires its CI smoke run.
These local results are not public release or independent-user evidence.

## Physical projection evidence

[Retained findings](findings/field-projection-2026-09-06/README.md) include the
initial admitted-budget failure and the corrected matrix. Buffers are reused
across execution windows, omitted groups have no corresponding storage, and cache
readers reject incompatible fields before record allocation. Synthetic results
are tracked separately from representative BAM/CRAM performance and external use.
