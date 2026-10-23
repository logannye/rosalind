# HG001 threshold transfer completed: practical screen not met

The metadata preflight and explicit source-level identity adapter resolved the
GIAB input incompatibility. The frozen comparison completed on hosted Linux with
all correctness checks passing, but **richer evidence did not meet the predefined
practical screen**. No workflow-validation nomination, selective escalation, P1
advancement or product promotion follows. Original P0, calibration and intake
failure history remains preserved.

The verified [hosted run](https://github.com/logannye/rosalind/actions/runs/39997188782)
tested source `10096419d00c18a03f64098c8c5f35701560bb59` on Ubuntu 24.04.
All 246 Python tests, Rust workspace tests, formatting and the release build
passed. No native build, executable test or scientific processing ran locally.
The [full hosted report](reports/continuation-attempt-2/report.json) and
[validation receipt](reports/continuation-attempt-2/hosted-validation.json) bind
the findings to source, run, artifacts and test logs.

The pooled models remain M10 (22 count features) and M11 (90 rich features), with
zero refits. Selection on the separate HG002 threshold role still nominates
`M10_pooled::retain_all` and `M11_pooled::depth_table_q0.999`. Exact policy bytes
were published before full fresh acquisition and evaluation, after the outcome-blind
metadata probe. The model, protocol and policy hashes are unchanged:

- Model pickle: `754bf9772105ebedc684076db9e5df51894ce4795e22987737b315d8b55288b9`.
- Protocol: `27cdb262fef4ed2e23f4f78914be766b34cfd89dd5898f4fc270f3bb049e3d1a`.
- Policy: `899ee3f34e9c3fb10d04da23d0637d07780b552f23239724c894b980c4ae0dbd`.

Fresh HG001 chr20:55,000,000–59,000,000 contains 5,319 assessed truth SNVs.
Candidates were generated independently at each frozen template-retention fraction,
using the unchanged two-ALT/2% rule. The depth conditions are paired derivatives
of one specimen and region. Fractions are thinning settings, not measured coverage.

| Template fraction | Candidate-generation recall | Rich positive candidates retained | Rich empirical recall lower bound | Rich truth retained end to end |
|---|---:|---:|---|---:|
| 35/300, original | 5,319/5,319 = 100% | 5,318/5,319 = 99.9812% | 99.9411%, **unsupported** | 99.9812% |
| 10/300, primary | 4,983/5,319 = 93.6830% | 4,981/4,983 = 99.9599% | 99.9197%, supported | 93.6454% |
| 5/300, novel low | 3,732/5,319 = 70.1636% | 3,732/3,732 = 100% | 99.9599%, **unsupported** | 70.1636% |

The original and novel-low bounds fail the frozen requirement for at least two
error-bearing 100 kb windows: they have one and zero, respectively. This is
insufficient bootstrap-tail support, **not observed retention below 99%**.
Numerical bounds cannot override their unsupported status. All conditions contain
40 positive-bearing windows, and all 2,000 bootstrap replicates were valid.
The estimates describe empirical spatial variation; they do not guarantee transfer
to arbitrary people, platforms or depth slices.

Counts retain-all preserves every generated candidate structurally, without
model scoring or stored features. Its end-to-end truth retention therefore equals
candidate-generation recall. Rich-minus-counts recall noninferiority passed in
all three conditions. The primary paired 95% interval was
[-0.09969, 0] percentage points, inside the fixed 0.5-point loss margin.

| Primary practical requirement | Observed counts → rich | Result |
|---|---|---|
| At least 25% fewer assessed false candidates | 69 → 53; 23.1884% fewer | Fail |
| At least 10 fewer false candidates, positive gain interval | 16 fewer; paired 95% interval 8–24 | Pass |
| At least 1% fewer assessed retained candidates | 5,052 → 5,034; 0.3563% fewer | Fail |
| At least 1% fewer total retained candidates | 5,697 → 5,672; 0.4388% fewer | Fail |

The primary rich policy removes 16 false candidates, two true candidates and
seven unknown-label candidates. It retains 638 unknowns and 487 unresolved
deferrals; these remain workload. Original-depth false candidates fall 380→12
and total retention 6,805→6,127, but that secondary result cannot rescue the
failed primary screen or unsupported cross-depth recall bounds. Novel-low
retention is 4,165→4,163, with no reduction in assessed false candidates.

At reduced depth, candidate generation misses 336 truth SNVs before filtering;
at novel low depth it misses 1,587. High conditional filtering recall therefore
does not establish reliable end-to-end recall. The separate caller recovers
5,318, 4,958 and 4,009 truth SNVs across the three conditions. Its caller policy
differs from the candidate rule, so it is a separate candidate/calling pipeline;
both end-to-end evaluations use the same 5,319 assessed-truth denominator.
Retained candidates, including unresolved cases, are not
final accurate variant calls.

Measured candidate callable-depth bins [0,5), [5,10), [10,20), [20,40), and 40+
contain respectively `[38,88,741,5362,576]`, `[487,3187,2007,16,0]`, and
`[1974,2089,102,0,0]` rows. Whole-region mean/median depths were not exported;
these candidate-conditioned distributions must not be presented as whole-region
coverage measurements.

The preflight took **0.961 seconds and 18,559,920 downloaded bytes**, including
the pinned BAI and bounded BAM ranges. It checked the exact retained header,
indexed access and 64 read identities before full reference/truth acquisition.
The full-region audit then reconciled 7,490,699 physical records, 3,805,142 exact
QNAMEs and 3,685,557 observed mate pairs. It rewrote no names and deduplicated no
records. Derived BAMs use an explicitly synthetic source-level read group;
original library provenance remains unavailable. The original files are unchanged.

All 12 external checks passed, including full template identity, coordinate order,
independent count parity, native/Parquet feature parity, standalone-F0 agreement, raw-unavailable
replay and input/preflight rechecks. Saved replay made zero alignment visits.
Full reference, truth and index checksums were verified on hosted compute.
The whole 589.7 GB source BAM checksum was not verified: source metadata and
range receipts bind the partial acquisition. Missing mates are reported as absent
from the selected scope, without inferring their location.

| Primary measured cost | Counts F0 | Rich F1 |
|---|---:|---:|
| Matched-site extraction, verification and first saved read | 2.102 s | 2.821 s |
| Retained native evidence | 4,430,896 bytes | 7,543,178 bytes |
| Standalone feature table, 5,697 rows | 134,772 bytes | 417,540 bytes |
| Median saved-table parse and policy execution, seven rounds | 0.099 s | 0.485 s |
| Median warm policy execution | 0.003 s | 0.202 s |

Counts retain-all does not require the feature tables or inference represented by
these optional F0 costs. Shared candidate scanning took about 40 seconds per
condition; the primary caller rerun took 10.356 seconds. One counterbalanced
extraction per arm/depth does not establish a general speed advantage. Repeated-query
totals are arithmetic projections. Human review time, selective fallback and net
end-user benefit remain unmeasured; stored evidence replay success alone does not
establish a complete workflow advantage.

The successful job used **730 runner-seconds**; its fresh phase took 418.240
seconds, downloaded 1,772,676,465 bytes including preflight, and recorded sampled
peaks of 2,647,212,032 RSS bytes and 4,912,823,912 scratch bytes. Sampled peaks can
miss brief maxima. The first new attempt
[39996264574](https://github.com/logannye/rosalind/actions/runs/39996264574) used
38 seconds and stopped on a synthetic SAM section-order assertion before any
preflight or fresh acquisition. Its [failure receipt](reports/continuation-attempt-1/receipt.json)
and narrow test-only repair are preserved. New usage totals **768 seconds**;
including the original three attempts, **1,386 seconds (23m06s)**. These are
provider job spans, not billed amounts. Both continuation attempts are consumed;
no study compute or overnight automation remains active.

Result artifact `12520577858` is 16,529,974 bytes, SHA256
`0b389958f07e9f12b528e9d290fa996f1fcf5b9fc42a286c1f6bb67b423528af`.
Freeze artifact `12520991213` is 16,378 bytes, SHA256
`9c0411fa173f324a77affb0cc60c8605cbf3171509e6782e98012ea07f0bcbd0`.
Independent review verified every one of 3,063 manifest entries and exact equality
between the separately published freeze and final result. Raw acquisition scratch
and the BAI are intentionally excluded from the artifacts; hosted recheck receipts
remain, and aggregate reports/provenance are committed. Bulky evidence stays ignored.

The appropriate practical conclusion is to retain the simple baseline for this
scope. A future, separately frozen study could target the much larger low-depth
candidate-generation losses and compare total review burden and caller reruns.
This consumed HG001 interval is now diagnostic evidence, not an untouched holdout.
No model, threshold, gate, interval or candidate rule was changed after its results.
