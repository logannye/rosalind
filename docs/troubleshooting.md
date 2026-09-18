# Troubleshooting

Start by recording `rosalind --version` and, for a source build,
`git rev-parse HEAD`. Keep the exact command and error text. The guides describe
the 0.5.0 source preview; the public 0.1.0 release has a different command surface.

| Symptom | Check and next step |
|---|---|
| `analyze evidence` or `dataset` is unknown | Run `command -v rosalind` and `rosalind --version`. Put the current source binary on PATH using [installation](installation.md). |
| Cargo fails while compiling native dependencies | Check Rust 1.83+, C/C++, CMake, pkg-config, and compression development headers in [build prerequisites](analyzer-sdk.md#build-prerequisites). |
| Python cannot import `rosalind` | Install the wheel in the interpreter's virtual environment and run outside the checkout. The distribution name is `rosalind-bio`; see [Python installation](../python/README.md). |
| Python and CLI versions differ | Use the native executable bundled with the installed wheel, or reinstall the matching candidate. Different RC numbers are different producers. |
| Output already exists | Choose a fresh output name/directory. Use `--force` only when intentionally replacing an existing artifact; it is not needed for a first run. |
| Multiple or ambiguous samples | Inspect alignment read groups. Select `--sample NAME`, or use `--pool-samples` only when pooling is the intended scientific operation. |
| Missing index or reference mismatch | Provide BAI/CSI for BAM, CRAI for CRAM, and matching reference contig names/lengths. Local CRAM FASTA needs FAI. Keep files immutable during the run. |
| Unsupported CRAM request | Check [supported versions, layouts, codecs, and limits](SEMANTICS.md#cram-decoder-admission). A refused layout can be valid CRAM outside this preview profile; use an appropriate BAM input. |
| A sparse CRAM request spends time before producing rows | The complete validation pass precedes indexed extraction, including during `--plan` and cache resume. Offline saved-dataset queries avoid the original CRAM. |
| Candidate REF mismatch or unsupported allele | Confirm the reference build and VCF alleles. The current evidence path accepts single-base A/C/G/T SNVs, not indels or symbolic alleles. |
| Counts differ from another pileup tool | Align quality/flag filters, missing-quality handling, read versus fragment units, BAQ, overlapping mates, and coordinate selection using [scientific semantics](SEMANTICS.md). |
| Resource plan refuses the request | Read the reported requirements; use a sufficient native budget or request fewer stored field groups. Reducing fields must still meet the analysis requirements. |
| RSS exceeds the declared budget or a record envelope fails | Treat this as a failed run, retain the error and measurements, and revisit the envelope/budget. A cooperative model is not a universal allocation cap. |
| Saved-dataset query refuses a locus or field | The stored selection and fields must cover the query. Missing data is not zero. Use the [reuse reference](reusable-evidence.md) to fill missing loci from compatible verified inputs. |
| Verification fails after copying a dataset | Copy the entire directory containing `evidence-dataset.manifest.json`, including its descriptor and partition directories; do not copy only the receipt. |
| Replay cannot locate inputs or reproduce the output | Supply the recorded immutable dependencies and matching producer. For an external analyzer use an explicit `--binary`; see [SDK replay](analyzer-sdk.md). |
| Python memory grows despite bounded native batches | Process each batch and release it. Lists, dataframes, and arrays retained by the consumer are outside the native budget. |

The [reuse quickstart](reuse-quickstart.md) is an executable check of dataset
relocation and queries after deleting only generated tutorial input files.
The [researcher tutorial](../examples/research-filter/README.md) is a small
reproducible extraction check with pinned public sources.

If a problem persists, include the source/candidate identity, platform, exact
command, error text, and a minimal shareable reproduction in a
[GitHub issue](https://github.com/logannye/rosalind/issues). Keep credentials and
private genomic inputs out of public reports. For a vulnerability, follow the
[security policy](../SECURITY.md).
