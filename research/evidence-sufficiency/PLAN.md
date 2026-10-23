# Evidence sufficiency for reliable SNV-candidate decisions

Analysis and implementation plan, 2026-10-22. Status: **proposed research; not executed**.
Repository inspected at `40773f30a15ed50e38ed8a9cbe2842ed2eb6070c`.

## 1. Decision and scope

Determine whether richer saved read summaries improve calibrated evaluation of
supplied SNV candidates on unseen specimens, and measure the cost of retaining
and using that information. If that succeeds, test whether selective acquisition
of richer evidence can reduce total work under an explicit error budget.

The observation is one **specimen / library / depth condition / locus / ALT**.
The biological unit is the individual, with family relationships retained.
Libraries, genomic blocks, sequencing reads, and downsampled copies are not new
biological replicates. The outcome is benchmark-supported presence of a supplied
SNV allele, not pathogenicity, a clinical diagnosis, or a complete genotype call.

The existing [bounded-cohort result](../../docs/findings/bounded-cohort-reanalysis-2026-10-22/README.md)
establishes exact reduction and replay, but no useful performance advantage for
that wrapper. Preserve its stop decision. The new study operates on individual
evidence rows, not its specimen-collapsed summaries. The historical 1,000-locus
HG002 archive remains an engineering regression fixture, not a representative
training or evaluation set.

**Primary hypothesis:** richer retained summaries reduce false-positive accepted
candidates compared with the strongest prespecified count-only baseline, while
preserving recall and calibration on a locked, unrelated individual's data.

**Separate systems hypothesis:** Rosalind provides useful integrity or resource
behavior relative to an implementation using the same features in Parquet.
Predictive gains from features or models do not establish a storage-engine gain.

This plan does not claim a novel easy/hard inference cascade. DeepVariant already
uses a small summary-feature model to route uncertain cases to its larger model
[S6]. The research opportunity is a measured account of evidence sufficiency,
calibration under shift, replay, and the cost of obtaining additional evidence.

## 2. Smallest useful first experiment: P0

Implement a **four-megabase, one-person feasibility experiment** before building
the larger study. It answers whether labels, feature extraction, model fitting,
and cost accounting work and whether richer summaries contain a promising signal.
It cannot establish transfer between people or laboratories.

| Item | Fixed P0 choice |
|---|---|
| Reads | Pinned HG002 NovaSeq PCR-free 35x chromosome-20 BAM and BAI [S1] |
| Reference | Full matching GRCh38 no-alt FASTA, uncompressed and indexed; compare the entire dictionary |
| Truth | Existing pinned HG002 v5.0q small-variant VCF and benchmark BED [S1, S2] |
| Training | BED intervals `chr20 10000000 11000000` and `chr20 20000000 21000000` |
| Probability calibration | `chr20 30000000 30495000` |
| Decision-threshold selection | `chr20 30505000 31000000` |
| Pilot test | `chr20 40000000 41000000` |
| Conditions | Original reads and a deterministic template-preserving subset with fraction `10/35`; report achieved depth |
| Models | Count-only and rich-feature logistic regression; count-only and rich-feature histogram gradient boosting |
| Limits | 100,000 distinct candidate loci per condition; 2 GiB retained evidence; 30 GiB scratch; 2-hour analysis wall-time ceiling |

Coordinates are zero-based, half-open. The intervals are proposed before outcome
inspection and are not selected from truth-positive loci. Intersect labels with
the applicable benchmark regions; retain excluded base and candidate counts.
The intervals leave a 10 kb boundary guard between calibration-fit and
threshold-selection regions, so the usable total is 3.99 Mb. All depth derivatives of a
locus keep its partition. Freeze the exact BEDs before acquisition.

Generate candidates from original-depth reads with a frozen permissive rule:
at least two callable observations of a non-reference A/C/G/T allele and ALT
fraction at least 0.02, under the explicit Rosalind default filter contract.
Use a bounded, streamed `depths,alleles` scan over the selected intervals and
discard noncandidate rows after tallying denominators. Do not load truth during
candidate generation. Record every ALT passing the rule; never inject truth
variants to improve candidate recall. Charge this first scan to the workflow.

Freeze the original-depth candidate list and extract both original-depth and
downsampled evidence at those same positions, including measured zero rows. This
paired comparison tests information loss at fixed candidates. A separately
labeled secondary analysis may regenerate candidates at reduced depth to study
candidate-generation failures; it must not be substituted for the paired result.

Downsampling uses a seeded hash of source identity, read group and template name
to retain or discard all records belonging to a template together. It creates a
new immutable BAM and evidence identity. Marginal quality histograms and counts
cannot simulate this operation. Hash seed and fraction are part of the protocol.

P0 uses one fixed setting per model, no test-driven tuning, and reports all four
arms. Fit preprocessing on training only, probability calibration on its assigned
region, and thresholds on the separate selection region. Score the pilot test
once. At P0's size, calibration and effect estimates may be unstable; record that
rather than presenting a statistical win.

**P0 deliverables:** compact per-candidate features and predictions, label and
exclusion counts, candidate-generation sensitivity, four-arm metrics, feature
ablations on development data, timing/storage records, exact source identities,
and a replay from saved evidence with original BAM paths unavailable.

**P0 gate:** correctness and native/Parquet equivalence must pass. Require both
classes in every fitting/calibration role and at least 50 negative candidates in
the combined development comparison. Use the two training blocks as two grouped
development-CV folds; the exploratory go criterion is a mean average-precision
gain of at least 0.005 for F1 over F0 within at least one learner. This is a signal
to fund P1, not a statistical discovery. The pilot-test score cannot determine
this gate. If class counts are inadequate, one bounded development expansion to
16 Mb may be planned using unused HG002 blocks. Record the new split and protocol
before running it. The original pilot test becomes development evidence; it never
becomes the final test. Adequate counts with no signal stop this study; a different
mechanism requires a separately stated hypothesis and protocol rather than an
automatic exception. Do not keep expanding to search for a win.

## 3. Open-access data and acquisition readiness

Access checks performed for this plan were public catalog, HEAD, and small header
or range reads. No bulk sequence data was downloaded. Located objects are not yet
validated execution inputs. Each acquired file or derivative needs its own hash.

| Resource | Use | Access/size evidence | Readiness and limits |
|---|---|---|---|
| HG002 35x NovaSeq PCR-free chr20 [S1] | P0 | Anonymous BAM HEAD: 1,093,938,924 bytes; BAI: 196,136 bytes; byte ranges advertised | URLs and SHA256s already in repository; fresh acquisition must verify them |
| Matching no-alt reference [S1] | P0 | Repository historical compressed size 886,344,618 bytes | Keep full dictionary; decompression adds roughly 3 GB, to be measured |
| HG002, HG001, HG005 matched 30x NovaSeq PCR-free release [S3] | P1 independent-person test | Public object listings: approximately 43.7–43.9 GB per BAM, approximately 9.4 MB per BAI; small HG001/HG005 range reads returned HTTP 206 | Slice selected intervals on hosted compute; exact hs38DH reference remains a mandatory unresolved preflight item |
| GIAB truth and difficult-region annotations [S2, S4] | Labels and error strata | Public NCBI release catalogs; repository pins HG002 v5.0q and stratifications v3.1 | NIST now lists stratifications v3.6. Either freeze the dated v3.1 annotation deliberately or pin v3.6 before running; do not silently update |
| Same sequencing release, a separate PCR-positive or instrument condition [S3] | First real library-shift extension | Release catalog describes multiple preparations/instruments | Exact files and provenance still to select; original person/family group remains unchanged |
| SEQC2 HCC1395/HCC1395BL [S7] | Later somatic study of site/preparation/purity effects | Public SRA accession SRP162370 and reference-material data repository | New tumor/normal protocol required; exact selected files not verified; technical repeats are one biological pair |
| Cancer GIAB HG008 [S8] | Later study of pancreatic tumor/normal material and passages/clones | NIST states public release without embargo; draft v0.4 clonal/truncal small-variant benchmark dated 2026-10-06 | Material passage/clone and benchmark applicability must match; draft status retained; not part of germline P0/P1 |

### Exact starting resources

Use the existing [resources.tsv](../../benchmarks/giab/resources.tsv) for P0 BAM,
BAI, reference, HG002 v5.0q truth and annotations. Do not modify those historical
pins. The first BAM is:

`https://storage.googleapis.com/deepvariant/case-study-testdata/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.chr20.bam`

The matched P1 release has this base URL:

`https://storage.googleapis.com/brain-genomics-public/research/sequencing/grch38/bam/novaseq/wgs_pcr_free/30x/`

| Filename under that base | Catalogued BAM bytes |
|---|---:|
| `HG002.novaseq.pcr-free.30x.dedup.grch38.bam` | 43,667,905,004 |
| `HG001.novaseq.pcr-free.30x.dedup.grch38.bam` | 43,820,263,410 |
| `HG005.novaseq.pcr-free.30x.dedup.grch38.bam` | 43,919,189,723 |

Append `.bai` for the indices. Use remote indexed extraction to materialize small
local BAMs before Rosalind runs. Use a normalized BED union and samtools' multi-region
iterator or equivalent physical-record-once traversal. Preserve headers, overlapping
mates and legitimate duplicate records; never deduplicate by QNAME or identical
record contents. Hash the actual derivatives.
Record remote object generation/version where available, length, checksum metadata,
selected intervals, command, and tool version. A derivative hash does not claim
that all approximately 44 GB of its source BAM was downloaded and verified.

The inspected HG005 header contains 3,366 contigs and an aligner command naming
`hs38DH.fa`, including ALT/decoy/HLA sequences. The P0 no-alt FASTA is therefore
not a substitute. Resolve the matching reference through the release/reference
provenance, compare all contig names/order/lengths and available sequence checksums,
and verify selected reference bases. If unresolved, this P1 acquisition route is
blocked. Choose and document another fully matched release before freezing P1;
do not improvise a reference or silently realign the selected regions.

**Truth-version policy:** use the currently recommended, explicitly pinned release
for each specimen. NIST lists HG002 v5.0q and says it supersedes HG002 v4.2.1;
v4.2.1 remains the listed small-variant release for the other GIAB samples [S2].
Do not downgrade HG002 solely to make version strings identical. Record release
and confidence-mask differences as a source of evaluation-domain differences.
A secondary HG002 label-stability analysis can use the old/new mask intersection,
with the deprecated release clearly labeled; it does not replace the primary labels.

HG001 truth directory:
`https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38/`

HG005 truth directory:
`https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/ChineseTrio/HG005_NA24631_son/NISTv4.2.1/GRCh38/`

Each contains `{sample}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz`, its `.tbi`, and
`{sample}_GRCh38_1_22_v4.2.1_benchmark.bed`. These filenames were located in public
directory listings; content hashes remain to be pinned before execution.

## 4. Locked P1 split and evaluation domain

Keep HG002 as development/training because P0 already informed design. Use HG001
for calibration and HG005 as the final test. All three P1 BAMs must come from the
matched 30x release; do not pool P0's 35x/no-alt alignment with these data.

| Partition | Person/family | Initial genomic budget | Permitted use |
|---|---|---:|---|
| Training and internal tuning | HG002; HG003/HG004 belong to the same family | 16 nonoverlapping 1-Mb blocks | Fit models; three-fold blocked development CV; choose model settings |
| Probability calibration | HG001 | 8 distinct 1-Mb blocks | Fit a sigmoid probability calibrator for each frozen model |
| Operating-point selection | HG001 | 8 further 1-Mb blocks | Choose decision thresholds and the baseline challenger |
| Final test | HG005; HG006/HG007 belong to the same family | 32 further 1-Mb blocks | Score frozen predictions and compare once |

Select all 64 blocks with a fixed seeded hash of autosomal reference coordinates,
with a 10 kb guard between roles. Use public reference masks/annotations to allocate
both ordinary and difficult contexts; record the selected composition. Require
at least 50% applicable benchmark coverage per block, based only on public BED
membership, not observed errors or truth-positive density. Coordinates assigned
to a calibration/test role cannot appear in training even in another person.
Exclude the P0 windows from P1 calibration/test roles. Freeze exact BEDs and hashes.

No person/family, library derivative, template, or locus crosses its designated
role. No coordinates, sample IDs, truth-derived annotations, or benchmark membership
enter the learned feature vector. Caller confidence is a separate baseline, not
a hidden rich-model feature. Multiple ALTs at a locus remain in the same block.

Size feasibility from training/calibration event rates and paired-error/block
variance before opening test predictions or labels. A prespecified expansion to
at most 128 Mb total (32 training, 16+16 calibration, 64 test) is allowed at this
stage only. Select the expansion from the frozen coordinate ordering. Once test
results are opened, no expansion, threshold change, or model choice is permitted.
If the design is underpowered within that ceiling, stop as inconclusive.

P1 uses P0's exact candidate rule on nominal 30x reads. Fit models, preprocessing,
probability calibration and thresholds on original-depth evidence only; never
pool depth derivatives into the training population. Freeze the original-depth
candidate universe per person for the paired stress evaluation.

The primary endpoint uses nominal 30x inputs. Apply the same frozen models and
calibrators to a template-hash `1/3` depth subset for a secondary 10x stress test;
measure actual depth. This is a controlled depth perturbation, not a new library
or independent biological replication. Claim real batch transfer only after a
separately acquired library/instrument condition is evaluated.

One held-out person supports a result about that person and domain. Bootstrap
intervals over genomic blocks quantify within-person uncertainty; they do not
create population-level biological replication.

## 5. Labels, denominators and features

Use exact normalized `(contig, position, REF, ALT)` matching for supported SNVs.
A benchmark-supported SNV allele is positive only if the benchmark genotype
actually contains that ALT index; merely appearing in a multiallelic ALT list is
insufficient. A missing or partial benchmark genotype is unknown. An unmatched
supported candidate inside an applicable confident region is a putative negative
only when it does not overlap an excluded benchmark event. For v1, exclude the
reference span of any non-SNV, symbolic, missing-genotype or unresolved benchmark
record plus 10 bp on each side; merge and record this conservative mask. Do not
decompose complex alleles into training labels. Outside-mask candidates,
unsupported alleles, unresolved REF, and excluded-event ambiguity remain unknown.
Apply the truth release's README rules. Independently check a deterministic audit
subset with the existing hap.py/vcfeval route or manual benchmark adjudication.

Always publish counts of evaluated, excluded and unobserved candidates by reason.
Report separately (a) truth SNVs proposed by the candidate generator, (b) recall
conditional on proposed candidates, and (c) end-to-end truth-SNV recovery in the
assessed domain. Class balancing or negative subsampling is prohibited in test
and calibration. Training resampling, if later introduced, requires a protocol
revision and correction before interpreting probabilities.

| Feature set | Contents |
|---|---|
| F0: counts | Callable depth; REF/ALT/other nucleotide counts; corresponding fractions; log1p counts; explicit undefined indicators; common REF/ALT base encodings |
| F1: compact richer evidence | F0 plus per-allele strand fractions, ALT/REF quality means and differences, mapping-quality means, read-position and read-length summaries, and exclusive filter-failure fractions |
| F2: optional information ablation | F1 plus global callable BQ/MAPQ histogram quantiles and tail fractions; secondary, not a second primary hypothesis |

Use the actual native field groups. F1 needs `depths,alleles,strands,quality-sums,
read-position,allele-quality`; F2 adds `quality-histograms`. `all-supported`
includes allele-specific quality; historical `all` does not. Scan counts first,
then extract rich evidence only at candidates with an explicit `cache_dir` and
record the resulting portable `evidence-dataset.manifest.json`. An ordinary
materialized Arrow artifact and receipt alone are not a portable dataset leaf.
Read the leaf from
`open_dataset(...).batches(...)`, require complete stream exhaustion and a successful
receipt, and preserve raw uint64 counts separately from deliberate model floats.

Zero-depth fractions are undefined, not zero; absent saved rows are unavailable
evidence. Training-fitted imputation plus missingness indicators may handle undefined
feature values, but must not fabricate a measured row. Expose unavailable-candidate
rates and refuse a completed primary run if extraction was unexpectedly incomplete.
Use only callable-read semantics; the histograms cannot recover excluded reads,
haplotypes, fragment consensus, or joint read-level relationships. A ratio of
position sums to length sums is not the mean of per-read relative positions.

## 6. Baselines and model budget

Use a feature-by-model factorial so a stronger algorithm cannot masquerade as a
better evidence representation:

| Arm | Representation | Learner/decision rule |
|---|---|---|
| B0 | Depth and ALT fraction | Calibrated grid of simple acceptance thresholds |
| B1 | Same supplied candidate universe | Pinned bcftools call confidence; candidate not emitted gets a defined lowest score, never disappears |
| M00 / M01 | F0 / F1 | Regularized logistic regression |
| M10 / M11 | F0 / F1 | Histogram gradient-boosted trees |
| D | Identical F0/F1 rows and models | Bounded, content-checked Parquet loading |

Pin and record bcftools preprocessing, filters, BAQ/overlap policy and version;
its different statistical assumptions are part of this baseline. Native feature
equivalence is assessed against the same producer/contract, with an independent
small read-count oracle, not by requiring bcftools AD to equal Rosalind read counts.

B0's fixed depth grid is `{0,5,10,20}` and ALT-fraction grid is
`{0.02,0.05,0.10,0.20,0.30}`. Choose the eligible pair with the fewest calibration
false positives subject to 99% calibration recall. Tie-break by higher recall,
then lower depth threshold, then lower ALT-fraction threshold. All score-based
arms choose the highest score threshold meeting that recall target; identical
scores remain tied. An arm unable to meet the target is ineligible, with its
failure reported. If none qualifies, stop as data-limited/inconclusive.

P0 uses LR `C=1` and HGB `max_iter=200`, `learning_rate=0.05`, `max_leaf_nodes=15`,
`l2_regularization=1`, `early_stopping=False`, fixed random seed. P1 allows LR
`C in {0.1,1}` and HGB `max_leaf_nodes in {7,15}`, `l2_regularization in {1,10}`;
other settings stay fixed. This is 12 configurations across F0/F1, at most 36
blocked development-CV fits plus final fits. Choose settings by development
average precision with log loss as a fixed tie-breaker. Preprocessing is fitted
inside each fold. Choose the F1 learner for the primary rich challenger by the same
training-CV objective. No hyperparameter selection from HG001/HG005 errors.

Use HG001's probability-calibration blocks for sigmoid calibration, then select
thresholds on its distinct threshold blocks to target 99% candidate recall.
Choose the strongest eligible simple baseline among B0, B1, M00 and M10 by the
fewest false positives at that operating point on those selection blocks; tie-break
by higher recall, then fixed arm order B0/B1/M00/M10. Independently nominate the
F0 LR/HGB probabilistic comparator using training-CV average precision and the
log-loss tie-break, because B0's two-threshold rule has no probability score.
Freeze these comparators, the rich challenger,
all preprocessing/calibration parameters and thresholds before test access.
Keep each same-learner F0/F1 contrast as an explanatory secondary result.

Parquet must preserve the exact same feature rows, masks and label universe.
Compare canonical raw features exactly; compare deterministic derived features
and fixed-model predictions under a documented numeric tolerance. Use projection,
bounded batches and ordinary content/manifests checks. Charge common extraction
to both paths. Separate native-reader versus Parquet-reader costs from any future
comparison of independently implemented feature producers.

## 7. Success, stop and invalidation criteria

These are research engineering margins, not clinical standards. Freeze them
before P1; revise only before opening its final test and record why.

**G0 — admissibility.** Proceed only with matched reference/sample identity,
resolved truth applicability, deterministic candidate generation, adequate input
access, bounded resources, and a passing exact-source extraction/feature test suite.
Unresolved reference or benchmark identity is blocked, not an ML failure.

**G1 — P0 feasibility.** Native and Parquet paths reproduce identical observations
and equivalent features/predictions; unknowns and zeros remain distinct; labeled
examples satisfy the class-count checks; replay succeeds. Apply the numerical
development-signal rule in Section 2. If it fails after any allowed count-driven
expansion, retain the benchmark and stop feature/model expansion.

**G2 — locked P1 predictive benefit.** The rich challenger must satisfy all of:

- False-positive count ratio versus the frozen strongest simple baseline is at
  most **0.75**, and the paired genomic-block-bootstrap upper 95% bound is below 1.
- The lower 95% bound for the difference in candidate recall is above **-0.005**
  (a 0.5 percentage-point noninferiority margin). Do not retune test thresholds to
  manufacture matched recall; show actual recall and candidate-generation recall.
- The upper 95% bound on calibrated Brier-score degradation against the nominated
  probabilistic F0 comparator is at most **0.005**. This preserves a proper scoring
  rule; it does not prove calibration, because Brier score also reflects
  discrimination. Report reliability plots and log loss, including difficult
  contexts, before making any separate calibration claim.
- At least 2,000 benchmark-positive proposed candidates, 100 baseline false
  positives, and 32 evaluable test blocks support the comparison. These are minimum
  information checks, not a power guarantee. Report difficult-region results and
  explain concentrated failures before extending the scope of any claim.

Use 2,000 paired block-bootstrap resamples with a fixed seed. Report ratio estimates
as unavailable when the baseline has zero false positives. Test class prevalence
stays unchanged and exclusions are visible. Passing nominal 30x does not imply
passing the 10x stress test or a new-library condition.

**G2 outcome taxonomy:**

- `benefit_supported`: all prespecified primary gates pass for the stated domain.
- `practical_effect_ruled_out`: sufficient data and the lower 95% FP-ratio bound
  exceeds 0.75, or the required recall/probability-quality margin is demonstrably
  violated. State which claim is ruled out.
- `inconclusive`: intervals span the target, event counts are inadequate, or the
  resource ceiling prevents a decisive test. Nonsignificance is not no benefit.
- `invalid`: leakage, changed test-driven selection, incompatible sources, failed
  extraction or incorrect labels. Fix engineering on development data; a used test
  cannot silently be reset to untouched status.

All non-success outcomes stop automatic expansion within this plan. Preserve
results and propose a new hypothesis only when the observed failure gives a
specific reason; do not increase model size or seek another favorable test split.

**G3 — systems value.** Reader/backend superiority requires a separate matched
workflow experiment. A useful initial target is at least 2x lower total cost or a
material resource-capability difference under equivalent output and integrity
requirements. Include materialization, hashing, training, calibration and replay;
report one-shot and K=5/20 repeated-analysis totals with three counterbalanced
repetitions of the affordable comparisons. No engine win is required for an ML
research result. If Parquet is equivalent and operationally better, use it for
the model-facing table and retain Rosalind's evidence/contract role.

## 8. Compute, storage and execution policy

All numbers below are initial planning ceilings or estimates, not measured runtime
claims. No GPU is needed for P0/P1. Acquisition and reference handling may dominate.

| Stage | Proposed execution allocation | Initial budget |
|---|---|---|
| Static development / unit fixtures | Local Python and documentation checks using installed compatible packages | No native compilation or package operation that implicitly builds Rust |
| P0 | Hosted Linux, 4 vCPU, 16 GiB RAM, 30 GiB scratch | 2 hours analysis; up to approximately 8 CPU-hours; at most 4 GiB network acquisition; 2 GiB retained evidence |
| P1 | Hosted Linux, 8 vCPU, 32 GiB RAM, 50 GiB scratch | 4 hours per complete analysis attempt; approximately 32 CPU-hours; 20 GiB acquisition ceiling; 5 GiB retained study evidence |
| Native build/validation | Existing hosted exact-revision CI pattern | Separate from scientific timing; one reusable build directory and one native validation process |
| Selective read-level escalation | Allocate only after G2 and an explicit fallback benchmark design | CPU/GPU needs depend on the pinned fallback; no reservation or spend in this plan |

Before execution, verify actual runner capacity and account permissions. These
allocations are not claims that a particular current hosted runner supplies them.
Track transferred bytes including retries, process-tree peak RSS, elapsed time,
CPU use, scratch high-water mark, and retained bytes. No cloud resources are
provisioned by this plan, and no monetary cost estimate is implied without a rate.

Start P0 with the fully pinned approximately 1.1 GB chr20 BAM if practical; an
indexed slice is also acceptable if its derivative provenance is explicit.
Never download all three approximately 44 GB P1 BAMs under the stated budget.
If indexed extraction exceeds the network/storage ceiling, stop and reassess
access rather than treating remote data as free or unbounded.

Full BQ/MAPQ histograms contain 349 uint64 bins, or 2,792 bytes per locus before
other fields. Rich evidence is roughly 3 KiB per locus uncompressed. Extracting
every chr20 position would approach 180 GiB before overhead. P0's candidate-first
100,000-locus cap keeps one rich condition near 300 MiB before overhead. Compact
F1 omits those histograms; F2 is optional. Admission failures do not justify
silently dropping test candidates or selectively choosing easier intervals.

Existing `benchmarks/giab/prepare.sh` makes a chromosome-only FASTA while a sliced
BAM can retain a full header dictionary. Native evidence requires equal dictionary
contig counts and matching names/lengths. Use the full matched reference as the
safe default. Do not reuse that preparation path unchanged for this study.

## 9. Implementation sequence and reviewable outputs

Keep this in `research/evidence-sufficiency/`; no public cohort API is required.
Use the existing Python/native evidence interfaces before proposing Rust changes.
Effort estimates assume dataset/reference access resolves; they are not deadlines.

| Work package | New files / change boundary | Acceptance evidence | Effort estimate |
|---|---|---|---|
| W0: freeze P0 inputs and protocol | `resources.tsv`, `protocol-p0.json`, `requirements.lock` | URLs, hashes, BAM/reference dictionary audit, exact BEDs, candidate rule, seeds and ceilings; no unresolved P0 input identity | 0.5–1 day |
| W1: prepare and label | `prepare.py`, `labels.py`, `splits.py` | Truth-blind candidates; bounded derivatives; explicit unknowns; family/locus/derivative split audit; candidate-recall accounting | 1–2 days |
| W2: feature table | `features.py` | Correct ALT vector indexing and ratios; exact integer source columns; native/Parquet equivalence; complete successful stream receipts | 1 day |
| W3: bounded P0 study | `study.py`, `evaluate.py`, `report.py` | Four frozen model arms, calibration separation, costs, replay, P0 go/stop report | 1–2 days |
| W4: lock and run P1 | `protocol-p1.json`, `locked-splits/`, new immutable dataset registry entries | Reference resolved; development-based sizing; precommitted test identity; one final evaluation with intervals and all outcomes | 2–3 days plus acquisition/compute |
| W5: independent reproduction | `docs/findings/evidence-sufficiency-<date>/`, reproducible workflow | Second clean environment reproduces labels/features and metrics within declared tolerance; one concise report and figure data | 1 day |

Use Python 3.11 and the existing pinned wheel-based NumPy/PyArrow/pysam/blake3/psutil
environment where compatible. Add pinned scikit-learn/SciPy wheels; LR and
`HistGradientBoostingClassifier` avoid a GPU or another training framework.
Pin samtools/bcftools in the hosted environment. Never install Rosalind from source
as a shortcut that triggers a local Rust build.

Add `.github/workflows/evidence-sufficiency.yml`: small fixture and pure-Python
checks on relevant PRs; manually dispatched bounded real-data runs; source, native
binary, dependency, protocol, data and model hashes in every result. Reuse existing
native conformance/replay/CRAM gates when native code changes; an older executable
does not validate changed Rust. Do not add bulk downloads to every PR.

Required focused tests cover multiallelic matching, benchmark boundaries and
complex-event unknowns; missing candidates in recall denominators; no family or
derivative leakage; train-only preprocessing and separate calibration roles;
zero versus absent observations; ALT quality-vector and strand indexing;
feature/backend equivalence; paired block metrics; failed streams/mutation;
and replay from relocated saved evidence. Existing known-answer fixtures remain
unchanged. Use independent expected values, not tests that duplicate the code.

Persist compact protocol, split hashes, all model settings, held-out predictions,
exclusions, metrics, bootstrap summaries, cost measurements and figure data beyond
short CI artifact retention. Large BAMs/evidence stay outside git. Retain useful
failure records. Identify and clean only this study's regenerable inactive outputs
at validation boundaries; preserve reusable inputs and scientific evidence.

## 10. Conditional next study: selective acquisition

After G2, compare three actions: decide from F0, acquire/use F1, or request a
pinned read-level fallback. Begin with simple confidence and depth policies;
train a routing policy only against those anchors. Benchmark against random
escalation at the same workload, simple confidence routing, and always invoking
the richer/fallback method. Include an established small-model cascade where
the supported task and model training provenance make comparison defensible [S6].

Obtaining omitted fields or discarded read information requires explicit fresh
extraction and a new evidence identity. The benchmark harness may collect rich
data in advance to estimate policies, but counterfactual savings are only a
simulation until the selective workflow is actually timed end to end. Charge
random I/O, decode granularity, repeated hashing, context padding and fallback
setup. Masking already-computed features does not save their acquisition cost.

Predeclare an error/abstention budget and a practical target such as 2x lower total
cost than always-rich processing at noninferior recall/calibration. Count every
abstention and fallback in the denominator. A high-quality subset alone is not a
complete workflow. Freeze the policy before another unused test set; P1's used
test can support development but cannot certify the new policy. Audit fallback
pretraining overlap with GIAB specimens before calling them unseen to that model.

The final portfolio artifact should explain one defensible result: which retained
information improved predictions, where it stopped helping, what failed under
shift, and how much work reliable decisions required. A negative result remains
useful when it rules out a concrete approach and preserves a simpler solution.

## Sources and evidence status

- **S1 — observed repository pins:** [GIAB resources](../../benchmarks/giab/resources.tsv)
  and [benchmark README](../../benchmarks/giab/README.md). BAM/index access and size
  were refreshed on 2026-10-22; reference size comes from retained repository
  preparation evidence, not a new full download.
- **S2 — primary benchmark authority:** [NIST Genome in a Bottle](https://www.nist.gov/programs-projects/genome-bottle).
  Current HG002 release, older-release deprecation, other-sample releases and
  public data access checked 2026-10-22. Source applicability still requires the
  per-release README and material/reference audit at execution.
- **S3 — sequencing release:** [Primary sequencing-release preprint](https://www.biorxiv.org/content/10.1101/2020.12.11.422022v1.full),
  [GIAB data indexes](https://github.com/genome-in-a-bottle/giab_data_indexes), and
  the exact public GCS paths above. Anonymous catalog/range/header checks establish
  availability, not full source-byte verification. [BWA reference documentation](https://github.com/lh3/bwa/blob/master/README-alt.md)
  informs the hs38DH investigation; it does not establish the exact reference hash.
- **S4 — interpretation and evaluation:** [GIAB FAQ](https://github.com/genome-in-a-bottle/giab_FAQ),
  [NIST benchmarking guidance](https://www.nist.gov/publications/best-practices-benchmarking-germline-small-variant-calls-human-genomes).
  These support confidence-region restrictions, error stratification and explicit
  uncertainty; numeric project success gates are this plan's proposals.
- **S5 — observed implementation:** [native semantics](../../docs/SEMANTICS.md),
  [Python dataset API](../../python/rosalind/dataset.py),
  [evidence API](../../python/rosalind/evidence.py), and
  [existing hosted workflow](../../.github/workflows/cohort-reanalysis.yml).
- **S6 — prior art:** [DeepVariant 1.8 small-model details](https://github.com/google/deepvariant/blob/r1.8/docs/deepvariant-small-model-details.md).
  Summary-feature routing is established prior work. This is a pinned historical
  architectural reference; execution would require a separately pinned comparator.
- **S7 — optional somatic extension:** [SEQC2 data paper](https://www.nature.com/articles/s41597-021-01077-5),
  [maintainer data guide](https://github.com/bioinform/somaticseq/blob/master/docs/seqc2.md),
  [public reference-material repository](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/).
  Dataset family located; selected-file readiness remains unknown.
- **S8 — optional evolving tumor benchmark:** [NIST Cancer GIAB](https://www.nist.gov/programs-projects/cancer-genome-bottle)
  and [HG008 draft v0.4 small-variant release](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/analysis/NIST_HG008-T_somatic-smvar_DraftBenchmark_V0.4-20261006/).
  Public-release and draft status checked 2026-10-22. Exact study inputs and
  material-to-benchmark applicability remain unselected.
