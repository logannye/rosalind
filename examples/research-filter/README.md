# Supplied candidates → exact evidence → research review

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

```sh
python3 -m venv /tmp/rosalind-tutorial-env
/tmp/rosalind-tutorial-env/bin/pip install pysam==0.23.3
/tmp/rosalind-tutorial-env/bin/python examples/research-filter/prepare.py /tmp/research-filter

target/debug/rosalind analyze evidence \
  --reference /tmp/research-filter/ex1.fa \
  --alignments /tmp/research-filter/sample.bam \
  --sites /tmp/research-filter/candidates.vcf \
  --memory-budget-mb 128 --output /tmp/research-filter/evidence.tsv

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
