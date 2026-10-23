# P1 reference and bounded-access preflight

Checked 2026-10-22 (Pacific). **Metadata and small range reads only; no P1 study
was run.** Full reference acquisition, complete sequence verification, exact
P1 input pins, and slice validation remain pending the P0 gate. This note does
not authorize or report a P1 result.

## Reference located and dictionary compatibility verified

The [1000 Genomes reference directory](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/)
provides the following exact objects:

| Object | Verified bytes | SHA256 of metadata downloaded in this preflight |
|---|---:|---|
| [GRCh38_full_analysis_set_plus_decoy_hla.fa](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa) | 3,263,683,042, HEAD only | Full FASTA not downloaded or hashed |
| [GRCh38_full_analysis_set_plus_decoy_hla.fa.fai](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai) | 154,196 | `f361f004bdae5ca68d632458b01a3848d02834ac7176f378e177344d148a6a6d` |
| [GRCh38_full_analysis_set_plus_decoy_hla.dict](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.dict) | 756,744 | `91f7a9993b70c9c85b3a51a20b888d866cd55938ce1ccd3a325a38954d2fbdeb` |

The FASTA advertises byte ranges, has ETag `"c287d9e2-510de0417b540"`, and
Last-Modified `Mon, 09 Mar 2015 17:03:41 GMT`. These are HTTP validators, not
cryptographic sequence hashes. Its dictionary contains 3,366 contigs totaling
3,217,346,917 bases, with an `M5` sequence checksum for every contig.

The [directory README](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/README.20150309.GRCh38_full_analysis_set_plus_decoy_hla)
states that this FASTA came from Heng Li's bwakit 0.7.12 and adds HLA sequences
to the GRCh38 full analysis set plus hs38d1. The reference creator's
[version-pinned URL manifest](https://github.com/lh3/ref-gen/blob/08c35a69166f4feb0e33295bc44fa9b0002afa40/URL.txt)
explicitly maps `hs38DH` to the
[NCBI mirror of this FASTA](https://ftp-trace.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa).
The mirror independently returned HEAD 200, the same 3,263,683,042-byte size,
and byte-range support. Equal size is not a full mirror-content comparison.

Two small full-contig range checks on the EBI FASTA verified the FAI offsets and
published sequence M5 values: `chrM` (16,569 bases, MD5
`c68f52674c9fb33aef52dcf399755519`) and `HLA-A*01:01:01:01` (3,503 bases,
MD5 `01cd0df602495b044b2c214d69a60aa2`). Their combined FASTA response bodies
were 20,356 bytes. These spot checks do not validate the remaining sequences.

## All three matched 30x BAM headers checked

Each exact BAM URL below returned HTTP 206 for bytes `0-65535`. Decompressing
the available BGZF members recovered the complete SAM header and binary BAM
dictionary. **All three dictionaries match the reference in all 3,366 names,
lengths, and order**, both in `@SQ` records and the binary dictionary.

| Sample / exact BAM object | Object bytes | GCS generation | Whole-object MD5 reported by source metadata |
|---|---:|---|---|
| [HG001.novaseq.pcr-free.30x.dedup.grch38.bam](https://storage.googleapis.com/brain-genomics-public/research/sequencing/grch38/bam/novaseq/wgs_pcr_free/30x/HG001.novaseq.pcr-free.30x.dedup.grch38.bam) | 43,820,263,410 | `1607665356371538` | `6de1148dc3040bc441cc3361dff64a2b` |
| [HG002.novaseq.pcr-free.30x.dedup.grch38.bam](https://storage.googleapis.com/brain-genomics-public/research/sequencing/grch38/bam/novaseq/wgs_pcr_free/30x/HG002.novaseq.pcr-free.30x.dedup.grch38.bam) | 43,667,905,004 | `1607665356367383` | `7378e27ba9a623ccb5eb1a4353aaebd9` |
| [HG005.novaseq.pcr-free.30x.dedup.grch38.bam](https://storage.googleapis.com/brain-genomics-public/research/sequencing/grch38/bam/novaseq/wgs_pcr_free/30x/HG005.novaseq.pcr-free.30x.dedup.grch38.bam) | 43,919,189,723 | `1607665356396692` | `793562eb86c4fdb2ac17ed21f00200f6` |

Source metadata is available from the public
[GCS release object listing](https://storage.googleapis.com/storage/v1/b/brain-genomics-public/o?prefix=research%2Fsequencing%2Fgrch38%2Fbam%2Fnovaseq%2Fwgs_pcr_free%2F30x%2F&maxResults=100).
Appending `.bai` to each BAM URL identifies its index. Previously checked
catalog sizes are 9,371,552 / 9,390,960 / 9,384,032 bytes for HG001 / HG002 /
HG005 respectively; full index acquisition and hashes remain execution inputs
to pin. Source-reported whole-BAM hashes were **not independently recomputed**.

Each BAM has one read group and the expected sole sample name. Each `@PG`
records BWA `0.7.17-r1194-dirty` using
`references/hs38dh_bwa_index/hs38DH.fa`. **No BAM `@SQ` record contains `M5`.**
The names/lengths fingerprint, computed as SHA256 of concatenated UTF-8
`name\tdecimal_length\n` records in order, is
`dd392f21162b59104a6dd43aad8e33a4010092d0eab30685594f2b8197a57980` for the
reference and all three BAMs. This fingerprint does not include sequence bases.

## Proposed complete reference verification after the gate

The creator's [Zenodo record 8045374](https://zenodo.org/records/8045374) supplies
[hs38DH.genome.tgz](https://zenodo.org/api/records/8045374/files/hs38DH.genome.tgz/content).
The [public record metadata](https://zenodo.org/api/records/8045374) reports
900,079,630 bytes and MD5 `f475410c90d6157aeaef10f45d7e25f5`; an anonymous HEAD
independently confirmed HTTP 200 and that size. No archive body was acquired.
The creator's
[version-pinned packaging recipe](https://github.com/lh3/ref-gen/blob/08c35a69166f4feb0e33295bc44fa9b0002afa40/ref-gen.mak)
includes FASTA, FAI, and dictionary, and rewrites FASTA capitalization/line
wrapping. Therefore an archive-derived FASTA is not assumed byte-identical to
the EBI FASTA, and the EBI FAI offsets must not be reused for it.

Proposed acquisition gate:

1. Freeze the archive URL, size, published MD5, metadata/dictionary hashes above,
   and BAM generations in a new P1 resource lock. Count all received bytes,
   including retries, against the aggregate budget.
2. Download the archive on hosted compute, verify its exact size and published
   MD5, and record its computed SHA256. Extract only expected regular files
   under bounded expanded-byte limits; reject links and path traversal.
3. Recompute every sequence's SAM-compatible M5 from the acquired FASTA and
   compare all 3,366 names, lengths, order, and M5s to the pinned EBI dictionary.
   Regenerate the FAI for the actual FASTA. Record SHA256s of the actual FASTA,
   generated FAI, and sequence dictionary, then recheck identity after use.
4. Recheck each pinned BAM's full binary and SAM dictionary and sole sample,
   and verify reference alleles for all truth and candidate positions used.
   Any mismatch stops acquisition before model fitting.

This establishes a documented `hs38DH` source with independently verifiable
sequence identity and a full dictionary match. The absent BAM `M5` values and
absent aligner-input checksum prevent cryptographic proof that the historical
aligner consumed these exact bases. Preserve that limitation; neither a filename
nor a dictionary match supplies the missing checksum. If the study requires
that stronger historical proof, publisher confirmation or a different fully
documented release is still necessary. Selected reference-allele checks cannot
prove the entire historical alignment reference either.

## Bounded BAM slicing proposal; not executed

After the gate, freeze the P1 BEDs before accessing study reads. Acquire and hash
the three small BAIs, and pin both BAM and BAI generations. Use an indexed
multi-region traversal over each sample's normalized BED union, preserving each
physical source record once, its flags/tags, and the full header. The
[samtools multi-region iterator](https://www.htslib.org/doc/samtools-view.html)
(`-M -L`) avoids repeated emission caused by overlapping requested regions;
it is not QNAME-based removal of legitimate duplicate alignments. Do not add
`--fetch-pairs`: out-of-region mates need not be acquired to retain overlapping
evidence, and that option expands the read domain.

Use an acquisition layer that counts every remote response-body byte across
reads, indices, reference, truth, and retries, enforces deadlines and the
**20 GiB aggregate ceiling**, checks HTTP 206/Content-Range for BAM range reads,
and verifies the expected generation on every response. Bind requests to the
frozen GCS generation using supported
[generation preconditions](https://docs.cloud.google.com/storage/docs/request-preconditions).
An opaque remote `samtools` invocation plus output-file size is not a download
meter. A controlled metering relay or an explicitly metered indexed range reader
is required before claiming the ceiling is enforced. Record derivative BAM
SHA256s, selected BEDs, source object identities, physical record counts, tool
versions, received-byte totals, and actual disk usage. Do not describe a slice
hash as verification of the unread remainder of the approximately 44 GB source.

For planning only, uniformly scaling source size by the initial 64 Mb total
across people gives roughly 0.9 GB of compressed alignment data (approximately
1.8 GB at 128 Mb), before BGZF/index overfetch, coverage variation, and retries.
The reference archive adds 0.90 GB compressed and approximately 3.3 GB expanded;
three BAIs add approximately 28 MB, with truth/annotation downloads additional.
These are estimates, not measured P1 costs or guarantees. The 20 GiB network and
50 GiB scratch limits remain hard stop conditions even if a difficult region or
an access implementation exceeds the estimate.
