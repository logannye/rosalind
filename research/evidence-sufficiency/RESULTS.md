# HG002 pilot result: stop before the locked study

**P0 is implemented and completed. The sole permitted expansion passed all
correctness and data-adequacy checks, but failed the predefined predictive-signal
gate. P1 was not started.** This is a feasibility result within HG002, not a
population or clinical result.

The final hosted run is
[39879349288](https://github.com/logannye/rosalind/actions/runs/39879349288), at
source `5a4152a5eb7306662be20bf6867147ce522a5198`. Its
[machine-readable decision](reports/p0-expanded/decision.json),
[full report](reports/p0-expanded/report.json), and
[execution identity](reports/p0-expanded/identity.json) are retained here.
Native builds and scientific execution ran on hosted Ubuntu; local work was
editing, pure tests, report inspection and hash verification.

## Decision against the frozen gate

The original 3.99 Mb run completed, but its calibration and threshold windows
had zero GIAB benchmark-confident bases. It fitted no models. The
[initial report](reports/p0-initial/RESULT.md) records that insufficient-data
result. The single count-driven replacement split was frozen at 15.99 Mb before
running it; its former test was retired into development and a fresh 55–59 Mb
test remained isolated. An early benchmark-mask coverage check now catches
unassessable blocks before expensive extraction.

The expanded development set contained 3,778 assessed negative candidates and
14,021 positives. Both classes occurred in both training folds, calibration and
threshold-selection roles. It generated 27,747 candidate loci across all roles.

| Learner | Counts-only mean blocked-CV AP | Rich-evidence mean blocked-CV AP | Gain | Required gain |
|---|---:|---:|---:|---:|
| Logistic regression | 0.999635916 | 0.999675053 | 0.000039137 | 0.005 |
| Histogram gradient boosting | 0.999446629 | 0.999835287 | 0.000388658 | 0.005 |

Neither learner passed. The recorded outcome is `stop_no_signal`; the test did
not influence this decision. No further P0 expansion or P1 execution is enabled
by this result. Previous cohort performance stops also remain in force.

The counts-only baselines are already close to AP's ceiling of 1. At these
observed levels, even perfect ranking could not add 0.005 AP. The prespecified
margin was therefore unattainable in this setting. We retain the stop rather
than lower the margin after seeing the data. This does not establish that richer
evidence has no practical use; it shows this experiment cannot support the
planned advancement claim.

## Correctness and retained evidence

The hosted run passed all native workspace tests, 29 existing Python integration
tests, 38 independent evidence-oracle tests and all 54 pilot tests. The real-data
run passed independent scalar/strand/filter count checks at 64 loci per condition,
exact native/Parquet observation and feature equality, saved-only replay with the
original input paths unavailable, and frozen-model prediction equivalence.
Per-allele quality sums additionally have native fixture and feature consistency
checks; the 64-locus real-data oracle is not a full independent per-allele audit.

All 3,967 files in the retained inventory verified by SHA256 and byte size:
401,622,465 bytes, including features, models and diagnostics. The frozen model
and model-manifest hashes also match the study report. The model/threshold
manifest was written before test prediction and evaluation. Original inputs and
consumed evidence were rechecked before completion.

The [retained file inventory](reports/p0-expanded/run-manifest.json) and compact
reports are committed. The full 32,008,328-byte compressed artifact is retained
in Actions for 30 days, with a local copy under the ignored
`results/evidence-sufficiency/39879349288` directory. Native and Parquet files
represent the same observations; this experiment does not establish a backend
speed advantage.

## Descriptive test results and practical implications

The frozen original-depth test evaluated 6,726 labeled candidates: 5,530
positives and 1,196 negatives. It retained 907 unknown labels separately. Before
test evaluation, the strongest simple comparator was bcftools confidence (B1),
the rich challenger was M11, and the probability comparator was count LR (M00).
These test results are descriptive, not an alternative advancement gate.

| Frozen arm | True positives | False positives | Candidate recall | Brier score |
|---|---:|---:|---:|---:|
| B0: depth/fraction thresholds | 5,524 | 33 | 99.8915% | Unavailable |
| B1: bcftools confidence | 5,502 | 3 | 99.4937% | 0.003534 |
| M00: counts + logistic regression | 5,498 | 7 | 99.4213% | 0.003209 |
| M01: rich evidence + logistic regression | 5,513 | 6 | 99.6926% | 0.002424 |
| M10: counts + gradient boosting | 5,507 | 10 | 99.5841% | 0.003280 |
| M11: rich evidence + gradient boosting | 5,512 | 9 | 99.6745% | 0.002364 |

The selected rich challenger did not reduce false positives relative to B1
(9 versus 3), although its recall was higher. These are different attained test
operating points, not a matched-recall superiority comparison. Candidate
generation found 5,530 of 5,546 assessed truth SNVs (99.7115%). The
[truth-recovery report](reports/p0-expanded/truth-recovery.json) includes those
16 missed candidates in its denominator. All claims retain the conservative
GIAB mask and complex-event exclusions described in
[label-domain.json](reports/p0-expanded/label-domain.json).

The reduced-depth stress test exposed a more concrete deployment concern.
With original-depth calibration and thresholds held fixed, count LR's candidate
recall fell to 84.5750%; rich LR reached 95.1899%. Their Brier scores were
0.06724 and 0.03133, respectively. The simple depth/fraction rule reached
94.3400% recall with 121 false positives, versus 95.1899% and 37 for rich LR.
However, this uses the original candidate list, a 10/35 template-retention
perturbation of the same specimen, and no rerun of B1. It does not establish
de novo low-depth calling performance or a benefit over a properly calibrated
counts-only alternative.

The useful product lesson is to assess calibration and threshold transfer when
depth changes before selecting a richer model. A future, separately frozen
experiment could test a practical false-positive/recall tradeoff under depth
shift, against depth-aware count calibration, simple thresholds and a rerun
caller. That study remains unstarted; the present plan stops at the P0 gate.

## Observed compute and storage

| Measurement | Expanded P0 |
|---|---:|
| Analysis wall time, excluding native build/test setup | 353.15 s |
| Input preparation | 92.26 s |
| Truth-blind candidate scan | 131.66 s |
| Rich extraction and initial saved reads, both conditions | 10.23 s |
| Caller baseline | 63.38 s |
| Model study phase | 17.04 s |
| Saved-only replay, both conditions | 1.88 s |
| Prediction-equivalence audit | 16.49 s |
| Downloaded input bytes, including retries | 2,023,232,866 |
| Sampled process-tree peak RSS | 2,289,688,576 bytes |
| Sampled analysis scratch peak | 7,662,047,753 bytes |
| Two native plus two Parquet evidence folders | 70,838,263 bytes |

All declared resource ceilings passed. Peaks were sampled, and these phase
measurements are not a controlled native-versus-Parquet speed benchmark. No GPU
was used. The validated harness, exact saved evidence and reproducible negative
advancement decision are the deliverables; no new predictive product model was
deployed.
