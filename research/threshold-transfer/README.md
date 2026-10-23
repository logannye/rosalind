# Frozen thresholds across coverage and specimen changes

Status: the [metadata-preflight continuation](CONTINUATION_RESULTS.md) completed
on fresh HG001 data. Input compatibility and every correctness check passed, but
the unchanged scientific screen failed. Primary total retained candidate count fell only 0.44%;
two depth conditions lacked the required support for their recall bounds.
No advancement or selective escalation follows. Frozen models, thresholds and
scientific gates are unchanged; the [original intake stop](RESULTS.md) is preserved.

This study asks whether fixed pooled counts and rich-evidence models can retain
at least 99% of assessed positive candidates under changed coverage while reducing
downstream candidate volume. It preserves the original P0 stop and P1 lock.

The [overnight plan](OVERNIGHT_PLAN.md) records the original time and resource budget;
the separate [input contract](INPUT_CONTRACT.md) bounds the completed continuation. The
machine-readable [protocol](protocol.json) pins source artifacts, data sources,
conditions, finite policy family and gates before acquiring fresh HG001 evidence.
[Data sources](DATA_SOURCES.md) explains sample, library, alignment and truth-release
differences and the limits of partially acquired BAM identity.

## Comparison

The base models remain M10_pooled (22 counts features) and M11_pooled (90 rich
features), from hosted calibration run 39895644473. Neither model is fitted again.
All newly constructed policies are shared across conditions, using raw model
scores rather than a condition-specific probability calibrator.

| Policy | Construction | Selection eligibility |
|---|---|---|
| Existing thresholds, with and without low-depth pass-through | Historical HG002 thresholds | Diagnostic only |
| Conservative global threshold | HG002 calibration positives, fixed recall quantiles | Yes |
| Small depth table | Same quantiles in fixed depth bins; sparse bins pass through | Yes |
| Retain all | No model needed; all candidates remain unresolved | Yes |

Global/table quantiles are 0.99, 0.995, 0.999 and 1.0. Bins are 0–4, 5–9,
10–19, 20–39 and 40+ callable reads. All new policies pass through depth below 5;
depth-table bins with fewer than 300 calibration-positive rows also pass through.
Both derivative depths can contribute calibration rows; these are not independent
people or independent replicates. No per-bin Cartesian threshold search occurs.

The separate HG002 threshold role selects the least total retained candidate
volume, including unknown labels, among eligible policies satisfying the recall
screen in both development conditions. A paired 100 kb spatial bootstrap uses a
simultaneous maximum-deviation correction over this finite selection family.
Nonstructural policies need at least 20 positive-bearing windows and two
error-bearing windows; degenerate zero-error tails are not declared certain.
Retain-all provides structural candidate-list retention only.

These are empirical spatial uncertainty estimates, not distribution-free
guarantees for new people, libraries, or arbitrary coverage shifts. Sparse depth
slices are descriptive and not automatically certified by an aggregate bound.

## Untouched evaluation and stop rules

The selected policies and model identity are published before full fresh acquisition
and evaluation, after the outcome-blind metadata preflight in the continuation.
HG001 chr20:55,000,000–59,000,000 (zero-based half-open) is evaluated once with
35/300, 10/300 and 5/300 template retention from the same nominal 300x source.
Actual callable depth is measured. Candidates are generated separately per
condition, without reference to truth labels. Missing candidate truth SNVs are
reported separately from filtering losses. No fresh outcomes select a policy.

Every declared condition must satisfy the 99% absolute candidate-retention recall
lower-bound screen for both selected policies and the 0.005 rich-minus-counts
noninferiority margin. The primary 10/300 condition additionally requires:

- At least 25% and 10 fewer assessed false candidates, with positive paired
  confidence-bound improvement.
- At least 1% fewer assessed retained candidates and 1% fewer total retained
  candidates including unknowns.
- Complete correctness, identity, extraction, caller and saved-feature cost records.

Unresolved/pass-through candidates remain in workload and denominators. Retention
is not autonomous variant-call correctness. Retained-list size is a workflow
proxy; the study does not measure human review time. A pass only nominates further
workflow validation. It never authorizes selective escalation, P1 or product
promotion. Failure or unsupported uncertainty is not permission to change gates.

## Execution and evidence

Use `.github/workflows/threshold-transfer.yml` on the dedicated development branch.
It runs all scientific tests and native compilation on hosted Ubuntu 24.04 with
the existing exact Python dependency wheels. No local native/scientific execution.
The workflow freezes policies before a separate fresh-evaluation step, retains
failure records, and publishes scientific disposition separately from job status.

Artifacts include frozen policy bytes, source identities, independent count audits,
counts-only/rich feature parity, raw-unavailable replay, every candidate decision,
depth diagnostics, candidate-generation recovery, caller reruns and measured costs.
Fully acquired reference/truth/index files are checksum verified. Remote read
ranges have per-range SHA256 receipts; the full large BAM checksum is not verified.

Counts and rich extraction use identical candidate sites, with extraction order
alternating across conditions. Standalone Parquet feature tables provide a simple
columnar baseline for seven counterbalanced repeated-query measurements. Projected
1/5/20-query totals are arithmetic projections, not observed end-to-end crossovers.
Structural retain-all policies do not need models or saved features; optional
serialization overhead is labeled separately from necessary policy cost.
