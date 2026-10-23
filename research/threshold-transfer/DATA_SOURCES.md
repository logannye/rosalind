# Independent HG001 sources and acquisition boundary

Metadata was inspected on 2026-10-22 before acquiring or evaluating fresh reads.
The adapter freezes the exact URLs, published MD5 values, object sizes and, for
the partially acquired BAM, Last-Modified value in `DEFAULT_SOURCES`. The study
protocol copies this dictionary before hosted execution. No HG001 observations
or truth labels may fit a model, recalibrate probabilities, select thresholds,
choose arms, or change the interval selection.

## Official source files

The reads are GIAB HG001, also named NA12878: NIST Illumina HiSeq data aligned
with novoalign 3.02.07 to GRCh38 with hs38d1 decoys and without alternative loci.
The nominal source depth is 300x. This differs from the HG002 NovaSeq library and
alignment workflow; the experiment tests a combined specimen/library/pipeline
shift and cannot isolate the contribution of each one.

- [GIAB alignment README](https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/NHGRI_Illumina300X_novoalign_bams/README)
- [GIAB alignment checksum manifest](https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/NHGRI_Illumina300X_novoalign_bams/novobam.HG001.md5sum)
- [GIAB alignment index](https://github.com/genome-in-a-bottle/giab_data_indexes/blob/master/NA12878/alignment.index.NA12878_Illumina300X_wgs_novoalign_GRCh37_GRCh38_NHGRI_03082016)

| Resource | Official published MD5 | Bytes |
|---|---|---:|
| HG001.GRCh38_full_plus_hs38d1_analysis_set_minus_alts.300x.bam | `f48bf321c34f32b7c08feb17b9e5b220` | 589,702,020,566 |
| Matching BAM BAI | `00fd8ac75a92daf341e79008e7e9d900` | 12,006,264 |
| GRCh38 no-alt plus hs38d1 reference FASTA gzip | `a056c57649f3c9964c68aead3849bbf8` | 874,637,925 |
| HG001 GRCh38 v4.2.1 benchmark VCF gzip | `ec61d76ad725d808e3a9e1c8acec402a` | 125,932,193 |
| HG001 GRCh38 v4.2.1 benchmark BED | `e1ea90c830f981d29a46043eb618e7bf` | 15,479,939 |

The read BAM Last-Modified value is `Tue, 08 Mar 2016 18:53:34 GMT`.
The benchmark file names are `HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz`
and `HG001_GRCh38_1_22_v4.2.1_benchmark.bed`.

- [HG001 GRCh38 benchmark directory](https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38/)
- [GIAB v4.2.1 README](https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/README_v4.2.1.txt)
- [HG001 benchmark checksum manifest](https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/md5.in)
- [NCBI matching reference directory and checksums](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/md5checksums.txt)

The matching reference is the specific file named in the alignment README:
`GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz`. Its complete
ordered contig names and lengths must match the BAM header. The adapter also
computes the chr20 sequence MD5 and checks BAM chr20 M5 when supplied. If the
protocol supplies `expected_chr20_md5` from the earlier HG002 reference, it must
match. Otherwise cross-specimen chr20 sequence equality is explicitly reported
as not separately verified; a shared assembly label alone is not that check.

The fresh benchmark uses v4.2.1 while the HG002 development evidence used v5.0q.
Each specimen is assessed only within its own benchmark BED and the existing
conservative SNV/complex-event exclusion rules. This is not an identical-domain
cross-release comparison or a clinical benchmark.

## Bounded acquisition

The whole 589.7 GB BAM is never downloaded. The hosted adapter downloads its
12 MB BAI and runs standard pysam/htslib indexed fetches only for the frozen
chr20 intervals. A loopback HTTP relay translates each upstream request into
closed HTTPS byte ranges of at most 8 MiB, requiring HTTP 206 and exact
Content-Range, size and Last-Modified values. It rejects an HTTP 200 whole-object
fallback. All received BAM bytes, including bytes read before a seek closes a
connection, share the same 4 GiB acquisition budget as the reference and truth
files. No range may silently bypass this accounting.

Published MD5 checksums are verified for the fully acquired reference, BAI, VCF
and BED, and SHA256 identities are recorded. The whole BAM's published MD5 is
source metadata only; it cannot be verified from partial acquisition. Every
received range prefix has its own length and SHA256 receipt, and every derived
local BAM is hashed. These identities establish acquired-byte provenance, not
independent authentication of the biological material.

Template retention uses the existing SHA256 function with one frozen seed,
source identity, read group and QNAME. Fractions 35/300, 10/300 and 5/300 are
nested under the same hash and are nominal thinning conditions, not guarantees
of measured 35x, 10x or 5x depth. Every physical alignment overlapping selected
intervals is considered; duplicate QNAMEs are not used to deduplicate records.
Reads crossing two selected intervals are owned by their first overlapping
interval, preserving physical-record multiplicity without double counting.

## Scientific and cost outputs

Candidates are generated independently for every depth using the unchanged
count/fraction rule and no truth mask or truth-label filtering. Every condition
gets its own candidate-generation recall against all eligible benchmark SNVs.
All candidate labels outside the assessed domain remain unknown. The returned
records contain the original feature schema and an explicit specimen field
outside model features.

The adapter extracts counts-only and rich evidence at identical candidate sites
with alternating extraction order across conditions. It derives F0 directly
from the counts-only measurement and requires exact equality with F0 inside the
rich features. It retains independent count audits, native/Parquet feature
equivalence and saved-only replay checks. Candidate scanning, paired extraction,
bytes and a bcftools rerun for every condition are timed separately. The caller's
unthresholded SNV counts and its unknown labels remain visible.

There is one extraction measurement per arm and depth, shared source acquisition
and candidate scanning, and ordinary OS caching. These measurements do not
establish a general speed advantage, selective fallback cost, deployment cost,
or human review savings. No selective escalation or new locked study is
authorized by this acquisition adapter.
