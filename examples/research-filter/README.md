# Supplied candidates → exact evidence → research review

Start here to inspect a small set of supplied SNVs without preparing your own
data. You will produce a table of exact read evidence, an illustrative research
screen, and a receipt that can verify and replay the extraction. For the current
source preview, follow [installation](../../docs/installation.md) first; the
public 0.1.0 release does not provide this evidence workflow.

This small tutorial uses real-origin alignments from **NA18507**, distributed in
the official SAMtools example. Upstream documents extraction of two build36
regions (chr2:2043966–2045540 and chr20:67967–69550), renamed as `seq1`/`seq2`, and
MAQ alignments adjusted with fixmate. See the
[pinned upstream provenance](https://github.com/samtools/samtools/blob/da72567097265a61650a081c9f68d4a9f45bd105/examples/00README.txt).
These are tiny historical human data slices, not a representative cohort or a
truth set. We do not infer any person's phenotype or medical status.

`sources.json` locks the upstream commit and SHA256 of the reference, compressed
SAM, and provenance text (about 118KB total). Preparation downloads these public
files, verifies them, constructs local indexes, and uses independent bcftools 1.21
through pysam 0.23.3 to generate candidate SNVs. The candidate calls are supplied
inputs to Rosalind, not evidence of Rosalind caller accuracy.

Build the current development CLI first (`cargo build --locked --bin rosalind`).
From the repository root, using Python 3.9+ and an isolated environment:

<!-- smoke:researcher-quickstart -->
```sh
python3 -m venv /tmp/rosalind-tutorial-env
/tmp/rosalind-tutorial-env/bin/pip install pysam==0.23.3
/tmp/rosalind-tutorial-env/bin/python examples/research-filter/prepare.py /tmp/research-filter

target/debug/rosalind analyze evidence \
  --reference /tmp/research-filter/ex1.fa \
  --alignments /tmp/research-filter/sample.bam \
  --sites /tmp/research-filter/candidates.vcf \
  --memory-budget-mb 256 --output /tmp/research-filter/evidence.tsv \
  --annotated-variants /tmp/research-filter/candidates.evidence.vcf

python3 examples/research-filter/join.py \
  /tmp/research-filter/candidates.vcf /tmp/research-filter/evidence.tsv \
  > /tmp/research-filter/research-review.tsv
target/debug/rosalind verify --manifest /tmp/research-filter/evidence.tsv.manifest.json
target/debug/rosalind reproduce --manifest /tmp/research-filter/evidence.tsv.manifest.json \
  --inputs /tmp/research-filter
```

Use a new output directory on a repeat; preparation and artifacts are create-new.
`preparation.json` records source hashes, tool versions, exact candidate-generation
arguments, and prepared file hashes. The pinned preparation produces four candidate
SNV records. `join.py` retains only the small candidate set while reading evidence
rows sequentially, and rejects missing/duplicate loci or mismatched REF.

The example labels a candidate `review` when callable depth is at least 10, ALT
read count at least 3, and both strands have at least one ALT read. These arbitrary
thresholds demonstrate an auditable research screen; they are not calibrated
confidence, a diagnosis, or a replacement for validated calling/filtering.

## Read the result

`research-review.tsv` reports each supplied REF/ALT pair with its callable depth,
ALT read count, ALT strand counts, and example screen label. `review` means only
that the row met the three stated thresholds. `insufficient_example_support`
means it did not; neither label establishes biological absence or significance.
The counts and labels should be read from your generated table, not inferred
from the four-candidate input count.

| File | What it gives you |
|---|---|
| `preparation.json` | Pinned public source identities and independent candidate-generation provenance |
| `candidates.vcf` | The four supplied SNV records |
| `evidence.tsv` | Exact read evidence at the selected loci |
| `candidates.evidence.vcf` | Original variant records with added read-evidence INFO fields |
| `research-review.tsv` | The supplied candidates joined to the illustrative screen |
| `evidence.tsv.manifest.json` | Extraction settings, content identities, verification, and replay recipe |

The annotated VCF preserves the supplied calls and adds `RSL_DP` for callable
depth, `RSL_AD` for counts in REF/ALT order, and strand counts. Verification and
replay cover both evidence and annotated variants. See the
[annotation reference](../../docs/variant-annotation.md) for the full field contract.

## Add target QC

For target-level denominators and the same per-position evidence in one pass:

```sh
target/debug/rosalind analyze panel-qc \
  --reference /tmp/research-filter/ex1.fa \
  --alignments /tmp/research-filter/sample.bam \
  --regions /tmp/research-filter/targets.bed --memory-budget-mb 128 \
  --position-output /tmp/research-filter/positions.arrow \
  --output /tmp/research-filter/panel.tsv
```

Read [the semantic profile](../../docs/SEMANTICS.md) before comparing counts with
another pileup implementation. In particular, this evidence counts reads,
excludes unavailable qualities, and does not apply BAQ or mate-overlap correction.

Continue with the [reuse quickstart](../../docs/reuse-quickstart.md) to relocate
saved evidence and query it after deleting only generated tutorial inputs, or
the [researcher quickstart](../../docs/researcher-quickstart.md) to use your own
data. See [troubleshooting](../../docs/troubleshooting.md) for input and resource
failures, and the [validation guide](../../docs/adoption-validation.md) to record
what this workflow helped you accomplish.
