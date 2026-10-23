# Threshold transfer: frozen policies, fresh validation incomplete

**The study stopped after its third permitted hosted attempt. Engineering tests
passed, but fresh HG001 evaluation could not start because the source BAM has no
read-group records. Recall reliability and practical benefit on fresh data remain
unmeasured. No advancement is authorized.** This is an acquisition incompatibility
with Rosalind's audited input contract, not a negative predictive result or evidence
that the GIAB specimen is invalid.

## What completed

The final [hosted run](https://github.com/logannye/rosalind/actions/runs/39991819570)
tested source `47773f64cb2c3b9e87c9ce369b03891b0766a32a` on GitHub-hosted Ubuntu
24.04, Linux x86_64, Python 3.11.17. All 211 Python tests passed: 54 original
evidence-study, 30 calibration, 60 threshold-transfer, 29 native-Python integration
and 38 independent evidence-oracle tests. Rust workspace tests, formatting checks
and the release build passed. No native builds, executable tests or scientific
processing ran locally.

The original pooled models were not refitted. Thresholds were constructed on the
HG002 calibration role and selected under identical recall constraints on the
separate threshold role. The policies were frozen and published before acquisition:

| Arm | Frozen selected policy |
|---|---|
| Counts | `M10_pooled::retain_all` |
| Rich evidence | `M11_pooled::depth_table_q0.999` |

Counts retain-all needs no model scoring or stored features to apply its retention
decision. Its candidate retention is structural; this does not establish candidate
generation sensitivity or final variant-call accuracy.

In the HG002 reduced-depth **development selection** rows, counts retained all
4,302 candidates, while rich retained 3,422. Rich retained 2,303 of 2,306 assessed
positive candidates, with 418 false candidates, 701 unknowns retained and 251
unresolved deferrals. Its empirical simultaneous recall lower bound was 99.0776%,
against 99.8699% attained retention recall. These figures helped select the policy;
they are not independent confirmation. The original-depth rich policy had the same
positive retention and lower bound, 21 false candidates and 2,708 total retained.
The complete [development report](reports/attempt-3/policies/development-report.json)
retains the finite selection family and support checks.

The final workflow required byte equality with the first pre-acquisition freeze:

- Policy SHA256: `899ee3f34e9c3fb10d04da23d0637d07780b552f23239724c894b980c4ae0dbd`.
- Protocol SHA256: `27cdb262fef4ed2e23f4f78914be766b34cfd89dd5898f4fc270f3bb049e3d1a`.
- Model pickle SHA256: `754bf9772105ebedc684076db9e5df51894ce4795e22987737b315d8b55288b9`.

All freeze files in the final result artifact match its separately published freeze
artifact. [Hosted validation](reports/attempt-3/hosted-validation.json) binds test
logs, artifact IDs, hashes, source and job timing. The final result ZIP is artifact
`12517903067`, SHA256
`9f63cde6a482771ddb69b4d99b59549c98d71d4f67e4e7c5bd45a6450b25ec75`;
the freeze ZIP is `12518421880`, SHA256
`0d86dc5cfef0315493aa1afe019d433dadf00e1b869618b55bc56217b4a0f155`.

## Why fresh evaluation stopped

The frozen source is the GIAB HG001/NA12878 nominal 300x novoalign BAM documented
in [Data sources](DATA_SOURCES.md). The hosted adapter verified the full reference,
BAI, benchmark VCF and BED downloads against their published checksums and inspected
the remote BAM header. Its retained [original header](reports/attempt-3/fresh/source-header.sam)
contains 2,580 reference entries and **zero `@RG` entries**. Its SHA256 is
`2a507e6a549c0e6b62bdccb14da7c4b10dec3c7882fb04e357782672ceabc856`.

The [admission receipt](reports/attempt-3/fresh/source-sample-admission.json)
records `accepted: false` and `source header has no read groups`. This occurs before
the reference/index compatibility audit, interval read fetching, downsampling,
candidate generation, candidate truth labeling or model scoring. Reference/truth
files were acquired and benchmark BED coverage checked; it would be inaccurate to
say no fresh data were accessed.

The preceding attempt had rejected the header without retaining its contents. A
source-specific finite library-alias repair was therefore a hypothesis, explicitly
documented as unconfirmed. The final receipt shows that missing read groups, rather
than an unrecognized library alias, prevent admission. Synthesizing groups or
dropping them from the read-template identity calculation would require a separately
reviewed contract; this study did neither. Partial BAM range receipts do not verify
the whole BAM's published checksum or independently authenticate the donor.

## Fresh metrics and gates

| Declared condition | Candidate-generation recall | Counts/rich retention and uncertainty | Workload/cost advantage | Gates |
|---|---|---|---|---|
| 35/300 template retention | Unmeasured | Unmeasured | Unmeasured | Not evaluable |
| 10/300, primary | Unmeasured | Unmeasured | Unmeasured | Not evaluable |
| 5/300, novel low depth | Unmeasured | Unmeasured | Unmeasured | Not evaluable |

These fractions are prescribed thinning conditions, not measured depth claims.
There are no fresh false-positive, end-to-end truth-retention, uncertainty,
extraction, caller or replay-cost results. The 99% absolute recall floor, 0.5
percentage-point noninferiority margin and practical workload screens remain
unevaluated. Missing results are explicit nulls in [study status](reports/study-status.json).
The original P0 stop, P1 lock and failed calibration advancement screen remain in
force. No selective-escalation router, workflow advancement or product claim follows.

## Attempts, interruption and observed cost

| Attempt | Source | Hosted result | Job seconds |
|---|---|---|---:|
| [39901841463](https://github.com/logannye/rosalind/actions/runs/39901841463) | `c42ecf3` | Workflow YAML rejected before jobs | 0 |
| [39902020683](https://github.com/logannye/rosalind/actions/runs/39902020683) | `859a128` | Tests/build and freeze passed; header admission failed | 275 |
| [39991819570](https://github.com/logannye/rosalind/actions/runs/39991819570) | `47773f6` | Tests/build and identical freeze passed; zero read groups confirmed | 343 |

Total observed hosted job time was **618 seconds, or 10 minutes 18 seconds**. This
is the sum of provider job start/completion spans, not a billing receipt. All three
attempts and both allowed repairs were consumed; there is no fourth attempt.

Usage exhaustion interrupted the second repair during the overnight window. After
the user replenished usage and explicitly requested resumption, the last attempt
ran in a separately recorded [operational window](execution-window.json). The
scientific protocol and frozen policy bytes did not change. The expired overnight
automation was deleted at the user's request. No study-owned compute remains active.

The final failed fresh phase took 44.94 seconds. Receipts account for 1,028,908,345
downloaded bytes, including 852,024 BAM bytes across four range receipts. Sampled
peak process-tree RSS was 192,356,352 bytes; sampled scratch was 4,178,937,770 bytes.
No resource-limit failure was recorded; sampled peaks can miss brief maxima.
These are ingestion costs, not a completed
workflow cost comparison; human review time and net end-user benefit are unknown.

## Appropriate next scope

Resolve intake compatibility before another predictive experiment. A cheap hosted
metadata preflight should check real sample/read-group metadata, indexed access and
reference identity before downloading the full reference and truth files. A passing
receipt should be bound to the subsequent acquisition and evaluation. For this BAM,
either establish a defensible source/read-template identity policy with explicit
limits, or select a compatible source under a newly frozen acquisition plan.
Neither path should invent missing library provenance.

Once intake is established, the meaningful scientific question remains whether an
unchanged conservative policy meets the same absolute recall and noninferiority
requirements on fresh data while reducing total downstream workload enough to
justify richer evidence. This run supplies no answer to that question.
