# Inspect candidate evidence and target coverage

Use this workflow when you already have candidate SNVs or a panel BED and want
to inspect the reads supporting them. Rosalind extracts evidence; candidate
generation and biological interpretation remain separate steps.

First [install the source preview](installation.md). If you do not have prepared
inputs, run the [pinned NA18507 tutorial](../examples/research-filter/README.md).
It prepares a tiny public dataset and four candidate SNVs, joins the evidence to
an illustrative screen, and verifies and replays the result. It is not an accuracy
benchmark.

## Bring your own inputs

Use a fresh output directory and provide:

- A coordinate-sorted BAM with BAI/CSI, or a supported CRAM with CRAI and a local decoding FASTA.
- An uncompressed reference FASTA with FAI that matches the alignment's contig names and lengths.
- Candidate SNVs in VCF, VCF.gz, or BCF, or targets in zero-based half-open BED.

Candidate REF and ALT alleles must be single A/C/G/T bases. The default profile
requires MAPQ 20 and base quality 20 and excludes secondary, supplementary,
duplicate-flagged, QC-failed, and unmapped reads. It counts overlapping mates
separately. Read [core concepts](concepts.md) before treating counts as molecules
or comparing them with another tool's output.

CRAM support has a [specific decoder profile](SEMANTICS.md#cram-decoder-admission)
and performs a complete validation pass before extraction. That setup cost also
applies to `--plan` and cache resume. For the first run, BAM is a simple starting
point if it is already available.

## Extract supplied SNV evidence

From the directory containing your inputs, with `rosalind` on PATH:

```sh
rosalind analyze evidence \
  --reference genome.fa --alignments sample.sorted.bam --sites candidates.vcf \
  --memory-budget-mb 256 --output evidence.tsv \
  --annotated-variants candidates.evidence.vcf

rosalind verify --manifest evidence.tsv.manifest.json
```

Open `evidence.tsv` in a text viewer or analysis tool. Each selected locus has
reference and allele counts, depths, strand counts, and quality/filter summaries.
`callable_depth` means the observations passed the declared technical filters.
It does not establish genotype confidence. The
[pinned tutorial's join](../examples/research-filter/README.md) shows how to connect
supplied ALT alleles to these columns without conflating missing loci with zero.

`candidates.evidence.vcf` preserves the supplied records and adds `RSL_DP`
(callable depth), `RSL_AD` (REF/ALT read counts), and strand evidence. The receipt
covers both outputs. These annotations do not replace the original genotype calls.

For Arrow output, add `--format arrow-ipc` and choose an `.arrow` output name.
For per-position evidence throughout targets, replace `--sites candidates.vcf`
with `--regions targets.bed` and remove `--annotated-variants`; annotation requires
supplied variant records. Add `--plan` before running to inspect the declared
resource model. Inputs must remain immutable throughout a run.

## Check the panel in one traversal

```sh
rosalind analyze panel-qc \
  --reference genome.fa --alignments sample.sorted.bam --regions targets.bed \
  --min-callable-depth 10 --memory-budget-mb 256 \
  --position-output positions.arrow --output panel.tsv

rosalind verify --manifest panel.tsv.manifest.json
```

`panel.tsv` summarizes each original BED target. Uncovered positions stay in the
target's full denominator; overlapping targets remain separate. The optional
`positions.arrow` contains evidence from the same traversal. The 10x threshold is
a configurable technical screen, not a clinical criterion.

## Know what to keep

| Artifact | Purpose |
|---|---|
| Evidence or panel output | The values you analyze |
| Annotated candidate VCF, when requested | Supplied records plus read-evidence INFO fields |
| Its `.manifest.json` receipt | Input/content identities, settings, output verification, and replay recipe |
| Immutable prepared inputs and matching producer | Required to reproduce extraction |
| Portable dataset directory, when requested | Required for offline evidence reuse |

To repeat the extraction from its receipt:

```sh
rosalind reproduce --manifest evidence.tsv.manifest.json --inputs .
```

Verification checks the output against its receipt. Replay reruns the recorded
calculation with the matching producer and inputs. Neither proves biological
accuracy. Existing outputs are protected; use new names for a new analysis.

Next, try the [reuse quickstart](reuse-quickstart.md) to query a relocated dataset
without the original alignments, or [annotate your VCF](variant-annotation.md).
Use [troubleshooting](troubleshooting.md) if sample selection, input compatibility,
or resource admission refuses a request.
