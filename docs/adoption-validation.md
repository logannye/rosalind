# Researcher and builder validation kit

Use this kit to observe three useful tasks with non-author participants: answering
a research question from candidate evidence, adding an external analyzer, and
integrating a workflow. Record the onboarding fixture and the participant's own
permitted task separately. A completed fixture establishes that the setup works;
`real_task_completed` describes the participant's actual work.

The [session template and schema](../examples/adoption/) are blank measurement
forms. No participants, successful sessions, time savings, or return visits have
been recorded by creating this kit. The acceptance target remains
[three non-authors completing tasks and two teams returning after 30 days](ADOPTION_ROADMAP.md#acceptance).
This kit adds no publication gate.

## Before a session

Copy [session.template.json](../examples/adoption/session.template.json) outside
the repository or under ignored `release/private-design-partners/`. Use one copy
per participant/task attempt. Start with `record_kind: "session"`; leave unknown
values `null`. Record a full tested commit, executable/package hash, installation
origin, OS/architecture, and whether the participant is a non-author. Distinguish
published-registry, candidate-wheel, offline-bundle, and source-checkout installs.
A candidate using Rust source patches is a candidate SDK test.

Measure wall time from the first installation action through a working `rosalind
--version`; also record active person-minutes and waiting time. Separately measure
time to the first verified result, integration effort, and maintainer assistance.
Retain failed attempts and their diagnosis. Do not label a blocked attempt as zero
minutes or a successful installation. Follow the selected artifact's
[installation instructions](../python/README.md) and the
[SDK's compiler prerequisites and candidate patches](analyzer-sdk.md).

Put that exact `rosalind` on PATH. The commands below run from this repository's
root and use a new directory. Fixture preparation needs Python with
`pysam==0.23.3`, as in the [pinned researcher tutorial](../examples/research-filter/).
Preparation downloads about 118 KB of content-locked public SAMtools examples;
record that setup time separately from analysis time.

```sh
export ADOPTION_WORK=/tmp/rosalind-adoption-session
python3 examples/research-filter/prepare.py "$ADOPTION_WORK"
rosalind --version
rosalind doctor --json > "$ADOPTION_WORK/doctor.json"
```

## Task 1: researcher candidate evidence and a second question

Ask the researcher to review supplied SNVs using exact read observations, then
answer a changed selection or summary question from saved evidence. The fixture
uses an illustrative review rule, not a calibrated variant caller or clinical
threshold. Read the [counting/filter semantics](SEMANTICS.md) before comparing tools.

<!-- adoption:researcher -->
```sh
rosalind analyze evidence --reference "$ADOPTION_WORK/ex1.fa" \
  --alignments "$ADOPTION_WORK/sample.bam" --sites "$ADOPTION_WORK/candidates.vcf" \
  --memory-budget-mb 128 --output "$ADOPTION_WORK/evidence.tsv"
python3 examples/research-filter/join.py "$ADOPTION_WORK/candidates.vcf" \
  "$ADOPTION_WORK/evidence.tsv" > "$ADOPTION_WORK/research-review.tsv"
rosalind verify --manifest "$ADOPTION_WORK/evidence.tsv.manifest.json"
rosalind reproduce --manifest "$ADOPTION_WORK/evidence.tsv.manifest.json" \
  --inputs "$ADOPTION_WORK"

rosalind analyze evidence --reference "$ADOPTION_WORK/ex1.fa" \
  --alignments "$ADOPTION_WORK/sample.bam" --regions "$ADOPTION_WORK/targets.bed" \
  --fields depths,alleles,strands,quality-sums --memory-budget-mb 256 \
  --cache-dir "$ADOPTION_WORK/cache" --format arrow-ipc -o "$ADOPTION_WORK/panel.arrow"
export ADOPTION_DATASET="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["measurements"]["execution.evidence_dataset_manifest"])' "$ADOPTION_WORK/panel.arrow.manifest.json")"
rosalind dataset verify --dataset "$ADOPTION_DATASET"
rosalind dataset extract --dataset "$ADOPTION_DATASET" \
  --sites "$ADOPTION_WORK/candidates.vcf" --fields depths,alleles,strands,quality-sums \
  --format tsv --memory-budget-mb 128 -o "$ADOPTION_WORK/saved-candidates.tsv"
python3 examples/research-filter/join.py "$ADOPTION_WORK/candidates.vcf" \
  "$ADOPTION_WORK/saved-candidates.tsv" > "$ADOPTION_WORK/saved-review.tsv"
cmp "$ADOPTION_WORK/research-review.tsv" "$ADOPTION_WORK/saved-review.tsv"
rosalind dataset panel-qc --dataset "$ADOPTION_DATASET" \
  --regions "$ADOPTION_WORK/targets.bed" --memory-budget-mb 128 \
  -o "$ADOPTION_WORK/saved-panel.tsv"
rosalind verify --manifest "$ADOPTION_WORK/saved-candidates.tsv.manifest.json"
rosalind verify --manifest "$ADOPTION_WORK/saved-panel.tsv.manifest.json"
```

Completion evidence: a useful reviewed candidate table, verified receipts, an
explicit explanation of zero-depth versus missing loci, and a second query answered
from the saved dataset. Stored fields and filtering must cover the new question;
unsupported filters or missing positions require fresh extraction. The fixture has
four candidate records. Record the participant's own input shape and result checks
without publishing genomic records. [Python/R/SQL examples](../examples/persisted-evidence/)
provide optional downstream interfaces.

## Task 2: builder external analyzer

Ask the builder to add one statistic or output needed by their project, declare its
fields and retained memory, and retain successful output/receipt behavior under a
small budget. Record active implementation/debugging minutes and assistance. First
run the [managed scaffold](analyzer-sdk.md) unchanged to isolate setup problems:

<!-- adoption:builder -->
```sh
rosalind new analyzer adoption-qc --api evidence --output "$ADOPTION_WORK/adoption-qc"
# Before Cargo, apply the SDK guide's explicit patches for an unpublished candidate.
unset CARGO_TARGET_DIR
(cd "$ADOPTION_WORK/adoption-qc" && cargo test && cargo build --release --locked --offline)
export ADOPTION_ANALYZER="$ADOPTION_WORK/adoption-qc/target/release/adoption-qc"
rosalind conformance analyzer --api evidence --binary "$ADOPTION_ANALYZER" --json \
  > "$ADOPTION_WORK/analyzer-conformance.json"
"$ADOPTION_ANALYZER" run --reference "$ADOPTION_WORK/ex1.fa" \
  --alignments "$ADOPTION_WORK/sample.bam" --sites "$ADOPTION_WORK/candidates.vcf" \
  --memory-budget-mb 128 --enforce -o "$ADOPTION_WORK/analyzer-summary.tsv"
rosalind verify --manifest "$ADOPTION_WORK/analyzer-summary.tsv.manifest.json"
rosalind reproduce --manifest "$ADOPTION_WORK/analyzer-summary.tsv.manifest.json" \
  --inputs "$ADOPTION_WORK" --binary "$ADOPTION_ANALYZER"
```

Completion evidence: the builder's own statistic and scientific test, declared
requirements, conformance JSON, and physical replay with the explicitly chosen
binary. Conformance tests lifecycle behavior; the builder still owns the statistic's
scientific oracle. Record lines of integration code only with a defined counting
method; person-minutes and removed orchestration steps are usually more informative.

## Task 3: workflow integration

Ask a workflow maintainer to connect one sample to an existing pipeline, carry its
plan into scheduler memory, produce panel and positional artifacts together, and
verify the receipt. Record executor/tool versions, integration effort, scheduler
request, declared Rosalind budget, and any independently checked OS limit.

The maintained [Nextflow example](../integrations/nextflow/) uses Nextflow 25.04.8
and Java 21. This exact local mode uses the selected executable on PATH and disables
Docker; it does not establish container, Slurm, or cgroup validation.

<!-- adoption:workflow -->
```sh
rosalind reference build --fasta "$ADOPTION_WORK/ex1.fa" --output "$ADOPTION_WORK/reference.rref"
nextflow run integrations/nextflow/examples/evidence/main.nf \
  -c integrations/nextflow/examples/evidence/nextflow.test.config \
  --reference "$ADOPTION_WORK/reference.rref" --bam "$ADOPTION_WORK/sample.bam" \
  --index "$ADOPTION_WORK/sample.bam.bai" --targets "$ADOPTION_WORK/targets.bed" \
  --budget_mb 256 --outdir "$ADOPTION_WORK/workflow-results" \
  -work-dir "$ADOPTION_WORK/nextflow-work" -with-trace "$ADOPTION_WORK/workflow-trace.tsv"
```

[Snakemake](../integrations/snakemake/) is an equivalent supported choice. Use its
`Evidence.smk`, copy `evidence.example.yaml` with absolute input paths, and run the
local profile. For a site executor or container, record a separate real run with
its actual configuration and published immutable image digest. A requested
scheduler amount does not establish an enforced OS cap. Completion evidence is a
verified artifact pair in the maintainer's own workflow and a correctly derived
scheduler request. Record which previous extraction/merge/wrapper steps were removed.

## Measurements and a 30-day return

Keep one `attempts` entry for each measured command or pipeline run, including
failures. Preserve raw logs/receipts privately and use anonymous labels in the form.
Do not copy participant paths, sample names, coordinates, credentials, hostnames,
organizations, contacts, or raw interview notes into a public record.

| Measure | How to record it |
|---|---|
| Input shape | BAM/CRAM, byte/record/contig counts, selected loci/targets, read-length and depth summaries, named/unknown/pooled scope; no sample IDs |
| Costs | Wall seconds, active person-minutes, installation versus preparation versus integration, assistance and verification/export time separately |
| Resource limit | Declared bytes, scheduler request, actual checked OS cap, and native RSS from its receipt; Python/SQL/workflow retention needs its own measurement |
| Repeat work avoided | Compare the same question, input/profile/fields and correctness result; record initial materialization and storage costs plus every reused query and verification |
| Proven reuse | `execution.reused_loci`/`computed_loci`, partition counts, and alignment record visits when emitted; unavailable counters stay null |
| Savings | Three matched repetitions for a timing claim; report raw attempts and medians, cold/warm/unknown filesystem-cache state, and failures. Distinguish measured savings from participant estimates |
| Return | At least 30 days after completion, record an observed new useful task by the same anonymous team, or an explicit negative/unreachable observation. No response remains unknown |

Dataset queries record `execution.alignment_record_visits = 0`; native evidence
uses `execution.record_visits`, and managed analyzers use
`execution.alignment_record_visits`. Visits can include rereading one alignment in more
than one window; they are not unique read counts or bytes read. Count an avoided
full extraction only when a saved query answers an actual additional question.
Include initial extraction, hashing, verification, and export in total-cost
comparisons. A reused task is not automatically faster on a small input.

The existing roadmap's return target is two distinct teams doing meaningful work
after 30 days. CI reruns, maintainer demos, stated intentions, and a scheduled
follow-up are recorded separately from observed return. This kit schedules or
contacts nobody; the follow-up fields start `not-scheduled` and null.

## Existing partner records and factual reporting

The [release partner schema](../release/schemas/design-partner-v1.schema.json)
accepts only its existing fields and persona-specific scenarios. Generate that
packet with `cargo xtask partners init`, as documented in
[maintainer releases](MAINTAINER_RELEASES.md#design-partner-gate), and validate the
completed sanitized feedback with `cargo xtask partners validate --input feedback.json --json`.
These three task categories do not replace the existing analyzer-builder,
workflow-hpc, or constrained-offline scenario lists. Link the same anonymous partner
ID where appropriate; keep supplemental session metrics outside the release gate.
Only consented, reviewed, anonymized records belong in `release/design-partners/`.

Use the [unfilled report template](../examples/adoption/report.template.md):
**who was eligible** (anonymous non-author/team
counts), **what completed or failed**, **measured installation/integration costs**,
**resource fit and work avoided**, **assistance and unresolved blockers**, and
**30-day observation with dates**. List incomplete and unreachable sessions in the
denominator. Never fill absent measurements with estimates without labeling them.

Existing [SDK installation/conformance findings](findings/evidence-sdk-2026-09-06/),
[portable dataset interoperability](findings/portable-evidence-2026-09-06/), and
[synthetic projection resource measurements](findings/field-projection-2026-09-06/)
provide reproducible technical baselines. Their fixtures, machines, package origins,
and limits are recorded there. They do not supply independent participant results,
representative workload performance, or 30-day return evidence for this kit.
