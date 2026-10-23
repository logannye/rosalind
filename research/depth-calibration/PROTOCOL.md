# Saved HG002 depth-calibration diagnosis

This is a separately specified, exploratory reanalysis of the saved HG002 P0
evidence. It asks whether depth-aware counts recover the benefit previously
associated with richer evidence. It does not reopen P0, satisfy its failed
advancement gate, authorize P1, or establish independent-specimen performance.
The earlier P0 test results have already been inspected and motivated this
question. Every result in this experiment is therefore exploratory, including
results from rows whose saved role is named `test`.

The numerical criteria below are research engineering screens for nominating a
fresh validation study. They are not clinical requirements or confirmatory
statistical tests. Freeze this protocol and the executable configuration before
the hosted run; retain their hashes and the exact tested source revision.

## Evidence and scope

Use the saved records and trusted, identity-checked model artifacts from expanded
P0 run `39879349288`, executed at source
`5a4152a5eb7306662be20bf6867147ce522a5198`. Verify input sizes and SHA256 hashes
against the retained inventory before analysis. Do not silently replace missing
or incompatible inputs. Loading a serialized model requires the established
trusted artifact provenance as well as its matching digest; a digest alone does
not make an arbitrary pickle safe.

The two conditions are original depth and deterministic template retention at
10/35. The latter is called `downsampled`; it is not a measured guarantee of
exactly 10x coverage. Both conditions use the original-depth candidate list.
This is a candidate-decision experiment, not de novo low-depth variant calling.
Candidates absent from that list cannot become predictions in this experiment.
Carry forward the original candidate-generation recovery report separately.

The assessed domain remains the saved GIAB benchmark-confident SNV domain with
its original complex-event exclusions. Preserve unknown labels and report their
counts by role and condition; never convert them into negatives. Report unknown
prediction counts but exclude unknown outcomes from supervised fitting and
accuracy metrics. Do not add or delete loci based on current model errors.

## Split discipline

Preserve all saved train, calibration, threshold-selection, and test roles. A
locus, every ALT candidate at that locus, and all its depth derivatives remain
in the same role. Verify candidate identity, condition pairing, block membership,
feature schema, and that F0 values are unchanged within F1. Reject duplicate
candidate-condition records or role-crossing loci.

Fit imputers, scalers, and base models on training rows only. Fit recalibration
models on calibration rows only. Select operating thresholds and nominate model
arms on threshold-selection rows only. Publish the frozen model identities,
arm nominations, and operating points before computing the new test-role
metrics. This separation avoids additional fitting leakage; it does not make
the previously inspected test region an untouched validation set.

Pooled training includes both conditions from the training role. Those rows are
correlated observations, not twice as many independent biological examples.
Report both candidate-condition row counts and unique candidate/locus counts.
Condition names, coordinates, roles, labels, and specimen identifiers must not
enter the base feature vectors. Existing F0 already includes callable depth,
counts, log counts, fractions, and base encodings. “Depth-aware” here refers to
training or probability correction across depth conditions, not a newly added
depth field.

## Prespecified model comparisons

Use the existing fixed LR and histogram gradient boosting settings without a
hyperparameter search. Compare the same learners with F0 and F1 so a feature
claim is not confounded with model capacity. Keep these distinct interventions:

1. **Frozen transfer controls:** apply the original fitted models,
   calibrators, and thresholds unchanged to both conditions.
2. **Threshold correction:** keep the original model and score transformation;
   choose a condition-specific operating threshold using only that condition's
   threshold-selection rows. Use the model's raw logit for corrected decisions
   to avoid artificial ties from saturated sigmoid probabilities; retain
   calibrated probabilities separately for proper scoring rules.
3. **Monotonic probability correction:** fit a positive-slope sigmoid to the
   frozen model score using condition-specific calibration rows, then select
   condition-specific thresholds on the separate threshold role. Verify the
   fitted slope and preserve ties. A monotonic transform cannot improve ranking;
   report ranking and acceptance agreement with threshold correction. Numerical
   probability clipping or ties must not be presented as new information.
4. **Pooled-depth learning:** fit F0 and F1 LR/HGB on the paired training
   conditions, with the same calibration and threshold-selection rules.
5. **Depth-conditioned probability prediction:** fit the prespecified LR using
   model score, measured log depth, and their interaction on calibration rows,
   for both F0 and F1. This transformation can reorder candidates and is a new
   depth-conditioned predictor, not pure monotonic recalibration. Keep its
   interpretation separate from intervention 3.

Include the fixed-grid depth/allele-fraction baseline B0. Include saved bcftools
confidence B1 only for the original condition, where caller scores exist.
Never reuse original caller scores as a downsampled caller result. This saved
evidence experiment cannot establish superiority over a caller rerun at reduced
depth.

Nominate one counts arm and one rich arm using downsampled threshold-role
outcomes: among eligible arms, minimize false positives at the 0.99 recall
target, then prefer higher attained recall, then the stable arm identifier.
Freeze the candidate arm lists and tie-breaking order in the executable
configuration. Keep all prespecified arm results visible; do not select a new
winner after reading test-role outcomes.

## Operating points and descriptive matched-target comparison

For deployable operating points, choose the highest score threshold reaching
0.99 recall on the appropriate threshold-selection rows and freeze it. Apply it
unchanged to test-role rows. Report TP, FP, FN, TN, attained recall, precision,
and unknown counts for both conditions. Include Brier score, log loss, average
precision, and reliability-bin counts for probability-producing arms. Brier
score combines discrimination and calibration; improvement alone is not proof
of calibration.

Separately compute a descriptive test-role curve comparison at target recall
0.99; report 0.95 and 0.995 as secondary descriptive targets. For each arm and
target, use the highest observed score threshold that reaches the target, with
all tied scores accepted together. This calculation uses test labels and is
explicitly an empirical comparison, not a newly deployable threshold or an
unbiased validation of threshold transfer. Never copy these thresholds into the
frozen policy.

Report actual attained recall beside every target-recall FP count. Ties can make
the achieved recalls differ; call this a same-target comparison rather than
asserting exact matched recall. Do not break ties using labels, discard tied
candidates, or hide unattainable operating points. An apparent FP reduction
that depends on materially lower recall cannot support the proposed benefit.

## Paired uncertainty and advancement screen

Use 2,000 deterministic paired genomic-window bootstrap resamples with 100 kb
windows and the seed frozen in executable configuration. Resample windows,
not individual rows. Every candidate, ALT, and condition from a sampled window
travels together. Recompute empirical target-recall thresholds within each
resample for empirical curve comparisons; preserve the frozen thresholds for
the deployment comparison. Report the number of nonempty windows, valid and
invalid resamples, and the 2.5th/97.5th percentile intervals. Do not quietly omit
failed resamples or zero-denominator ratios. Absolute FP differences remain
reportable when a ratio is undefined.

Require at least 20 nonempty windows and 1,900 valid resamples for a usable
advancement interval; otherwise the screen remains inconclusive.

These intervals describe variation among windows within one previously
inspected chromosome region of one individual. Their interpretation depends on
the chosen window scale and does not establish independence across people,
libraries, platforms, or arbitrary shifts. They are exploratory uncertainty
summaries, not confirmatory significance claims.

A rich arm can nominate a fresh predictive validation study only if all of the
following hold for the primary downsampled comparison with the frozen counts
nominee:

- At empirical target recall 0.99, FP count is reduced by at least 25% and by at
  least 10 candidates. The counts arm must have nonzero FP count for the ratio.
- The paired 95% interval for absolute FP improvement has a lower bound above
  zero. An unavailable or inadequately supported interval is inconclusive.
- At the frozen operating points, rich-minus-counts recall is no worse than
  -0.005, and its paired 95% interval lower bound is above -0.005. If a usable
  interval cannot be obtained, do not declare this requirement passed.
- The empirical same-target comparison does not obtain its FP advantage by
  lowering attained recall by more than 0.005 relative to counts. Report the
  exact difference even when the criterion passes.
- Input identity, split/feature checks, complete predictions, monotonic-control
  checks, and the execution resource limits all pass.

The original-depth outcomes and proper scoring rules remain visible secondary
results and identify regressions or limits to the nomination. Passing the screen
is `nominate_fresh_validation`, not established benefit. If a practical margin
is missed, report that the prespecified screen was not met. If intervals span
the target, counts are inadequate, or a cost comparison is unavailable, report
the corresponding uncertainty rather than concluding equivalence or no value.
Do not add model arms or tune margins in response to the observed result.

## Costs and selective escalation boundary

Measure saved feature serialization size after gzip compression and selected
model artifact sizes. Measure parsing, matrix preparation, and inference for
the selected counts and rich paths over seven warmed repetitions with the
same candidate set, hosted machine, and thread settings. Preserve each
repetition and summarize medians and ranges. Alternate or counterbalance order
so one arm does not always benefit from earlier cache warming. Report timing
scope, preprocessing included, and exclusions. Record fit/calibration wall and
CPU time separately from repeated inference. Do not infer hosted dollar costs
without an applicable provider rate.

The saved artifact already contains rich evidence. Its acquisition cost is sunk
in this run. Smaller projected features or faster inference do not measure
selective extraction, omitted-field acquisition, read-level fallback, random
I/O, hashing, or end-to-end escalation cost. Historical aggregate extraction
timings are context, not a matched marginal-cost measurement. Report total
acquisition/escalation advantage as **unknown**.

Consequently this run cannot authorize selective escalation, even if the rich
arm passes the exploratory predictive screen and has modest saved inference
overhead. `proceed_to_selective_escalation` and `proceed_to_p1` remain false.
A favorable result may nominate a separately frozen fresh-data and complete-cost
validation. Only a demonstrated practical predictive advantage together with
measured total workflow cost can justify pursuing selective escalation. A
negative or inconclusive result retains the cheaper adequate approach and the
evidence needed to explain the decision.

## Execution and failure handling

Run model fitting, scientific metrics, performance measurements, and any native
build/test on hosted compute. Local work is limited to editing, static
inspection, and non-scientific validation that does not compile native code.
Use saved evidence without fresh genome downloads or candidate regeneration.
Set bounded CPU, memory, elapsed-time, artifact-size, and retry limits in the
executable configuration and workflow. Record the runtime, dependency versions,
source revision, protocol/input identities, command, outputs, and failures.

Missing classes, incompatible models, missing feature fields, non-finite scores,
unavailable caller results, failed identity checks, incomplete predictions,
and invalid uncertainty summaries must remain explicit. Never substitute zero
for absent evidence or silently downgrade a required check. A failed run may
be repaired for an identified implementation error with the change recorded;
the dataset does not regain untouched status. Preserve results and the original
P0 stop regardless of the new outcome.
