# Annotate candidate SNVs with exact read evidence

The development CLI accepts VCF, VCF.gz, and BCF site selections without a variant
index. Alignment inputs still require BAI/CSI or CRAI. Only A/C/G/T SNVs with
distinct ALT alleles are supported; indels, symbolic alleles, reference mismatches,
malformed records and undeclared header fields fail explicitly.

```sh
rosalind analyze evidence --reference genome.fa --alignments sample.sorted.bam \
  --sites candidates.vcf.gz --sample SAMPLE \
  --fields depths,alleles,strands,allele-quality --memory-budget-mb 256 \
  --format arrow-ipc --output evidence.arrow \
  --annotated-variants candidates.evidence.vcf.gz

rosalind verify --manifest evidence.arrow.manifest.json
rosalind reproduce --manifest evidence.arrow.manifest.json --inputs . --no-attest
```

Use `--plan` on the analysis command to check admission for your selection and
header. Choose `.vcf`, `.vcf.gz`, or `.bcf` for the annotated destination. The
original record order, duplicate records, ALT order, IDs, quality values, FILTER,
existing INFO, sample order, FORMAT and phased genotypes are preserved semantically.
HTSlib may normalize their text representation. Rosalind adds INFO fields; it does
not rewrite calls or infer which VCF sample corresponds to an alignment sample.
The output header records the explicit named, unknown, or pooled alignment scope.

| INFO | Meaning |
|---|---|
| `RSL_DP` | Callable A/C/G/T read depth |
| `RSL_PF` | Prefilter aligned base observations; deletion/reference-skip positions excluded |
| `RSL_ED` | Eligible aligned base observations before base-quality filtering |
| `RSL_AD` | Callable read counts in original REF, then ALT order (`Number=R`) |
| `RSL_ADF`, `RSL_ADR` | Forward/reverse counts, when `strands` is requested |
| `RSL_BQS`, `RSL_MQS` | Exact per-allele base/mapping-quality sums |
| `RSL_RPS`, `RSL_RLS` | Exact per-allele stored-SEQ position/read-length sums |

The last four fields require `allele-quality` and contain decimal unsigned integer
strings, preserving the full `uint64` range. Counts use VCF Integer and refuse
values exceeding its signed 32-bit range; Arrow/TSV retain the exact counts.
Read positions are zero-based offsets, adjusted for reverse orientation, and
cannot reconstruct bases removed before alignment. Divide each sum by the
corresponding count for a mean; a zero count makes the mean undefined. An observed
zero is distinct from an unrequested field, which is absent.

`RSL_AD` can sum to less than `RSL_DP` when reads support a base not listed in the
record's alleles. Overlapping mates count as separate reads. See
[SEMANTICS.md](SEMANTICS.md) for quality, flag, sample and CIGAR rules.

In Arrow, `allele-quality` adds four fixed-size `uint64[4]` columns in A/C/G/T order:
`allele_base_quality_sum`, `allele_mapping_quality_sum`, `allele_read_position_sum`,
and `allele_read_length_sum`. TSV uses four comma-separated values per column.
Historical `--fields all` remains the original six groups; `all-supported` opts
into all seven. Existing schema-1 and field-mask-v1 artifact bytes are unchanged.
Masks containing the new group use schema 2 and field-mask version 2.

```python
from rosalind import materialize_evidence

result = materialize_evidence(
    "genome.fa", "sample.sorted.bam", "evidence.arrow",
    sites="candidates.bcf", annotated_variants="candidates.evidence.bcf",
    fields=["depths", "alleles", "allele-quality"], memory_budget_mb=256,
)
print(result.annotated_variants_path)
```

Annotation uses verified canonical Arrow partitions. It retains at most one
16,384-base partition while restoring original variant order; widely interleaved
records may reread partitions. A temporary cache is cleaned on completion; use
`--cache-dir` and `--resume` for persisted reuse. The planner includes physical
source fields, lookup state, native headers, records, formatting and compression.
Variant header/record envelopes are configurable through
`--max-variant-header-bytes` and `--max-variant-record-bytes`. HTSlib allocations
are checked after decoding, so these remain cooperative limits.

The evidence artifact, annotated variants and receipt publish atomically. An
existing Rosalind INFO tag is a conflict, including with `--force`. Invalid input
or an I/O error preserves previous destinations; resource failures identify any
partial outputs explicitly. The receipt binds the original variant file's bytes
as well as the normalized scientific selection, because record order and duplicate
records matter to annotation. Byte verification and replay cover both outputs.
Compressed annotation replay uses the new versioned annotation contract; a native
codec change can produce a byte divergence. Historical BAM/BGZF replay exclusions
remain unchanged. No semantic comparison substitutes for physical verification.
