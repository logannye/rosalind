# Adoption and reusable evidence implementation

This sequence follows the evidence-engine foundation merged in PR #97. Work remains
ordered; external publication and adoption prerequisites are tracked separately
from code and local verification. The original scientific invariant and historical
artifact verification remain requirements throughout.

| Order | Delivery | Status |
|---|---|---|
| 1 | Shared panel defaults and explicit named/unknown/pooled sample scope | Implemented; 11 sample tests, Rust/CLI threshold parity, Python/native parity, cache/replay and existing evidence regressions pass |
| 2 | Installable evidence-engine RC, executable onboarding, candidate publication gates | Prerequisites inspected; publication not performed |
| 3 | Physical field projection and bounded coalesced sparse fetches | Pending |
| 4 | VCF.gz/BCF selection, record-preserving annotation, optional per-allele quality/position sums | Pending |
| 5 | Verified subset/superset dataset reuse, persisted analysis, bounded Parquet and Python/R/SQL adapters | Pending |
| 6 | Evidence-native artifact runner, scaffold, and conformance | Pending |
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
