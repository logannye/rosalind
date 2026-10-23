# Recall-constrained threshold transfer: bounded overnight study

Authorized by Logan on 2026-10-22 America/Los_Angeles. Started 2026-10-23
04:22 UTC; stop starting work at 12:22 UTC (05:22 Pacific), approximately eight
hours later. This is a new study. The original P0 stop, P1 lock, and failed
calibration advancement screen remain historical facts.

## Objective

Determine whether a conservative, frozen threshold policy can preserve candidate
recall under depth changes on fresh GIAB data, and whether rich evidence reduces
useful workload compared with counts under exactly the same requirements.
No new estimator architectures, base-model refits, or selective-escalation router.

## Execution sequence

1. Inspect and pin prior HG002 model/artifact identities. Freeze a finite threshold
   family, data sources, intervals, depth conditions, selection rules and gates.
2. Implement policies and independent edge-case/provenance tests. Reuse the existing
   evidence extraction and labeling contracts. All executable tests and scientific
   processing run on GitHub-hosted Linux; local work is editing and static checks.
3. Construct thresholds using the HG002 calibration role; select the counts and rich
   policies on the separate threshold role. Use identical finite policy families
   and uncertainty procedures. Freeze selections before fresh evidence is scored.
4. Acquire bounded indexed fresh GIAB Illumina data from an unrelated specimen,
   preferably HG001, with matched GRCh38 reference and benchmark regions. Verify
   exact source compatibility first. Generate candidates independently at each
   declared depth. Evaluate once; no changes informed by fresh outcomes.
5. Audit identities, candidate-generation losses, filtering recall, unresolved
   candidates, false positives, storage and repeated-processing cost. Record a
   completed pass/fail/inconclusive scientific decision separately from CI success.
6. Reserve the last hour for verification, coherent commits and a morning report.
   Stop early if a decisive scientific gate fails; finish independent documentation
   and engineering verification without additional hypothesis searches.

## Scientific boundaries

- Freeze the existing pooled F0 counts and F1 rich models; no fresh-person fitting.
- Compare existing global thresholds, conservative global thresholds, a small
  fixed measured-depth table, and low-depth pass-through/retain-all baselines.
- Target 99% candidate-list retention recall, with one-sided 95% uncertainty.
  Spatial block-bootstrap bounds are empirical uncertainty estimates, not a
  distribution-free guarantee across people, libraries or arbitrary depth shifts.
  Correct development screening for the finite policy search; disclose sparse or
  degenerate support. Retain-all gives deterministic candidate-list retention only.
- Require the absolute recall floor and a rich-minus-counts recall-loss margin of
  0.005. Both arms can fail together; equivalence alone cannot establish reliability.
- For practical screening require at least 25% and 10 fewer false candidates, a
  positive paired confidence bound, at least 1% fewer assessed retained candidates,
  and at least 1% fewer total retained candidates including unknowns. These volume
  measures are workload proxies, not measured manual review time.
- A positive screen nominates further workflow validation only. Full extraction,
  storage and repeated-analysis costs must be recorded before any computational
  advantage claim. Selective escalation, P1, production defaults and clinical
  claims remain unauthorized by this experiment.
- All unknown/deferred candidates remain explicit and incur workload. Report
  candidate-generation truth recovery separately from conditional filtering recall.
- All downsampled derivatives and matching loci remain paired. Fresh data must not
  select thresholds, models, bins, regions, metrics or the winning arm.

## Overnight operational bounds

- Use the current dedicated development branch; preserve existing results and main.
- One active scientific workflow at a time; at most three hosted attempts total,
  at most six hosted runner-hours, and no new paid cloud service or GPU provisioning.
- Per job: hard timeout at most 180 minutes, 14 GiB process-tree RSS, 30 GiB scratch,
  4 GiB acquired source bytes, 2 GiB retained scientific artifact, and sufficient
  disk checks. Use bounded chr20 data, not whole-genome acquisition.
- At most two repairs of execution/provenance defects, preserving failed attempts.
  A statistical failure is not permission to alter gates or consume another test.
- No local Rust/native builds, tests, scientific model execution or data acquisition.
- No merges, releases, outreach or changes to product defaults. Commit/push study
  code, aggregate reports and provenance; keep bulky inputs/models/results ignored.
- Maintain `results/threshold-transfer/overnight-session.json` with source hashes,
  active run IDs, phase, attempt count, fresh-test consumption and next action.
- At the deadline, stop launching work, cancel any remaining study-owned compute,
  preserve partial evidence and mark incomplete validation honestly. Pause the
  heartbeat on completion or deadline. Stay quiet on unchanged status.

## Morning report

Give tested revision/platform and hosted run links; exact data/source identity;
recall and uncertainty by condition; false-candidate and workload differences;
actual costs and scope; stop/advance disposition; failures and unresolved limits.
If no scientifically valid result completes, report that rather than an inferred win.
