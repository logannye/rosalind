# Saved HG002 calibration: ranking benefit, escalation remains stopped

**Completed on hosted compute. Rich evidence improved reduced-depth ranking
after stronger count baselines, but failed the frozen recall-reliability screen.
The result is `inconclusive_validation_screen`; selective escalation and P1
remain disabled.** This is exploratory reuse of previously inspected HG002
evidence, not independent validation or a clinical result.

The successful [hosted run](https://github.com/logannye/rosalind/actions/runs/39895644473)
executed source `5808c19a45a3eda5746df0c65e3a252fb79ab8fd`. Its
[complete report](reports/calibration-v1/report.json),
[model freeze](reports/calibration-v1/frozen-models.json), and
[verification receipt](reports/calibration-v1/hosted-validation.json) are retained.
The [protocol](PROTOCOL.md) was frozen before the experiment. Neither earlier
infrastructure repair changed its dataset, model settings, margins, or selection
rules. The original P0 `stop_no_signal` decision remains unchanged.

## What calibration explained

At reduced depth, the original counts-only logistic regression had 84.575%
candidate recall and 3 false positives using its original-depth threshold.
Selecting a reduced-depth threshold on separate development rows restored
recall to 98.969%, with 417 false positives. The same model's AP stayed
0.996987. Thus much of the dramatic recall loss was a threshold-transfer problem;
recovering recall also exposed a substantial false-positive cost.

Condition-specific positive-slope sigmoid calibration preserved those decisions
exactly while improving Brier score from 0.067241 to 0.036258. All 32
role/condition/base-model monotonic controls passed. Probability correction
improved the scoring rule, while changing the threshold recovered recall.
Neither added ranking information to the frozen model.

Training on both paired depth conditions improved the count baseline further.
The development threshold role nominated pooled-depth histogram gradient
boosting for both counts (`M10_pooled`) and rich features (`M11_pooled`). These
choices were frozen before computing the new test-role metrics. All prespecified
arms remain visible in the complete report; no winner was substituted afterward.

## Matched recall and threshold transfer answer different questions

The test region contained 6,726 assessed candidates per condition: 5,530 positives
and 1,196 negatives, plus 907 unknown labels retained outside accuracy metrics.

At the descriptive empirical 99% recall target, both selected arms attained
exactly 99.0054% recall (5,475 true positives and 55 false negatives):

| Condition | Counts FP | Rich FP | Reduction |
|---|---:|---:|---:|
| Reduced depth | 219 | 141 | 78 fewer; 35.6% |
| Original depth | 6 | 6 | None |

The reduced-depth paired 95% genomic-window bootstrap interval was 23–127 fewer
false positives. All 2,000 resamples were valid across 40 nonempty 100 kb windows.
This supports an exploratory ranking advantage in the reduced condition. These
curve points use evaluation labels to match recall; their thresholds are not a
validated policy and were never copied into the frozen models.

Using the actual thresholds selected on the separate development threshold role:

| Reduced-depth selected arm | TP | FP | FN | Attained recall |
|---|---:|---:|---:|---:|
| Depth-aware counts | 5,492 | 271 | 38 | 99.3128% |
| Rich evidence | 5,472 | 131 | 58 | 98.9512% |

Rich-minus-counts recall was -0.3617 percentage points, with a paired 95% interval
of -0.6577 to -0.0906 points. Its lower bound failed the requirement to exceed
-0.5 points. The point estimate passed that margin, but the uncertainty did not.
Rich evidence also fell slightly below the desired 99% recall on this region.
The 25% relative and 10-candidate absolute FP screens, positive FP-gain interval,
data-support checks, and cost-measurement completeness checks passed. **Failure to establish the
recall margin keeps the overall screen inconclusive.**

The lowest measured callable-depth bin (0–4 reads) illustrates the tradeoff:
counts had 201 TP, 47 FP and 6 FN; rich evidence had 192 TP, 17 FP and 15 FN.
This is a descriptive subgroup of the same candidates, not a new selection rule.
Reliable operating thresholds at low coverage remain unresolved.

## Measured saved-evidence costs

Seven counterbalanced warmed repetitions used the same hosted machine, numerical
thread settings and candidate pool, with each arm applying its frozen decision rule.

| Measurement | Counts | Rich | Rich/counts |
|---|---:|---:|---:|
| Features | 22 | 90 | — |
| Canonical compressed feature table | 1,604,890 B | 7,378,748 B | 4.60x |
| Median matrix construction + calibration + prediction | 0.140 s | 0.271 s | 1.93x |
| Median full-table parsing + selected-row prediction | 0.685 s | 2.195 s | 3.20x |
| Selected serialized model | 361,195 B | 414,615 B | — |

Inference scored 7,633 reduced-depth test candidates, including unknowns. The
parsing measurement scanned each complete 55,660-row compressed feature table,
then selected those same test candidates. It measures that serialization and
access pattern, not native-versus-Parquet performance or indexed random access.
Absolute replay overhead was small in this pilot; the ratios alone do not show
that rich evidence is economically unattractive.

The full run took 93.32 seconds including artifact verification and loading,
with sampled peak RSS of 1,398,820,864 bytes (about 1.30 GiB). It fitted four new
base models and 18 calibrators. Source-artifact downloading, dependency setup,
and hosted tests are outside that timer. No GPU, native build, fresh alignments,
or new truth data were needed. All resource ceilings passed.

Rich extraction had already occurred before this experiment. Initial acquisition,
selective retrieval, raw-read fallback, caller execution and human-review savings
were not measured. **A full workflow cost advantage remains unknown.**

## Decision and retained evidence

Do not implement selective escalation or start P1. Keep counts as the baseline
for further development. The result supports investigating reliable thresholds
and coverage diagnostics; it does not establish equivalence between counts and
rich evidence or erase the observed ranking gain. Any subsequent advancement
needs a separately justified protocol, fresh specimen/library evidence, a
reduced-depth caller comparison and full workflow costs.

All 42 hosted tests passed: 12 existing study tests and 30 calibration/provenance
tests. The source ZIP, every consumed member, original inventory entries, trusted
model/dependency identities and final input recheck passed. The documented P0
telemetry omission from its inventory is handled explicitly; that file still
passes the pinned ZIP and individual SHA256 checks. The local result download
also verified all 11 output-inventory entries, completion/report/manifest hashes,
model/freeze identities and exact executed analysis sources.

The complete artifact is 26,014,461 compressed bytes and remains in Actions for
30 days. A local copy, full predictions, models, feature tables, and logs are under
`results/depth-calibration/39895644473`. Compact receipts and the complete metric
report are committed under `reports/calibration-v1`. The
[validation history](reports/validation-history.json) preserves the two earlier
failures: a degenerate synthetic fixture and an inventory-adapter mismatch.
Neither reached HG002 model fitting.

All findings use the original candidate list, deterministic 10/35 template
retention, and the same person's conservative GIAB SNV assessment domain.
The paired spatial intervals do not establish donor, library, platform or
arbitrary-shift robustness. Candidate-generation recovery is retained separately
in [the source recovery record](reports/calibration-v1/source-truth-recovery.json).
