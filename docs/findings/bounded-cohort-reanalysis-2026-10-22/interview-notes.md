# Engineering mechanisms and evidence boundaries

Read these notes with the [measured report](README.md), which binds the completed
checks to the actual hosted source and binary. The G2 decision is to retain the
correctness implementation and stop optional performance expansion.

## What the representation preserves

A saved row contains exact read observations under a particular extraction
contract. Retained A/C/G/T counts can answer a new ALT question at an already
stored position. The same callable depth can be screened at a different positive
threshold. Lowering the original base-quality threshold cannot recover reads that
were excluded before the counts were stored. Per-position counts also cannot
recover haplotype linkage, fragment identity or UMI consensus.

Explain the boundary as `Q(D) = g(E_phi(D))`: a downstream question is recoverable
only when the retained representation under extraction settings `phi` contains
what that question needs. The project demonstrates a narrow recoverable family,
not a universally sufficient statistic or a new variant caller.

## Why comparison identity differs from source reuse identity

Same-source reuse must stay bound to the original alignment content and selected
sample. Cross-specimen comparison expects different alignment sources, while
requiring compatible measurement semantics. A single reuse key cannot express
both relationships correctly.

The native `compatibility_blake3` stays unchanged. The research comparison checks
explicit filters, exact read counting, ordered contigs, required fields, and the
effective analysis-reference role and byte identity. The separate duplicate guard
rejects the same alignment content and complete named scope under two specimen
aliases. Equal sample labels across different files do not establish biological
identity or independence.

Equal contig lengths or assembly names are insufficient reference evidence. The
implementation conservatively requires identical original reference bytes and
source role. Descriptor v1 lacks an explicit format enum, so the code does not
guess FASTA versus packed reference from a pathname or declare sequence-level
equivalence between different representations.

## Why missingness belongs in the calculation

An observed zero-depth row and an absent saved position have different meanings.
The former increments the observation denominator; the latter does not. Both can
add zero to a summary, but those additive zeros describe different sets of
observations. The output therefore retains `n_requested` and `n_observed` alongside
eligible/support counts and exact observation totals.

A candidate absent from every member has zero additive totals over an empty
observed set, undefined fractions, and `reference_validation=unavailable`.
A stored row with zero callable depth still permits REF validation. Fractions
retain numerator and denominator as integers and become null at denominator zero.
Multiple ALT rows repeat the locus denominator; adding those repeated denominators
would double count. Inputs and sums remain checked uint64 values, including values
beyond binary64's exact-integer range, and overflow refuses publication.

## What makes execution bounded

Execution admits a bounded candidate window, then visits one specimen at a time
through one native reader. Window counters and the current specimen's observation
mask replace a dense specimen-by-locus matrix. This is not simply an `O(W)` claim:
member metadata, normalized candidates, coverage, descriptors, native decoders,
receipts and output buffers also occupy space and have admission limits.

The parent model and native planner are cooperative controls, not OS allocation
caps. Process-tree RSS is sampled and may miss short peaks. The internal analyzer
sample excludes final publication; the benchmark's outer worker sampler covers
the full measured process lifetime. Keep these measurement scopes separate.

Small windows can repeatedly touch one physical partition. Native saved-data
reads verify and decode complete touched partitions, even when Python requests
fewer columns. Each fresh native child also has startup and metadata costs. These
mechanisms predict a memory–runtime tradeoff; the pilot measures whether it is
useful before adding a cache or query coalescing.

## What the receipts establish

Native inspection and coverage validate metadata and stored-coordinate membership.
They do not certify partition payloads or candidate REF. Execution performs those
checks when it consumes rows. A successful native stream must be exhausted and
have a completed receipt before its observations support publication.

The research input ledger accumulates actual consumed dependencies and rechecks
their union before publication, including dependencies used by earlier children.
Generated coverage masks and descriptor copies are execution inputs too; they
must remain bound to the versions accepted during planning. This prevents a
later modified mask from silently changing an observed denominator. It does not
provide snapshot isolation against adversarial concurrent mutation.

The outer manifest binds the research outputs and recorded provenance, separately
from native receipts. Hashes prove byte identity under a stated comparison; they
do not authenticate authorship, validate biology, or replace an independent
numerical oracle. Published native child receipts retain original staging-path
provenance. Use the wrapper's verify/replay flow for the outer artifact; do not
claim generic native replay of the aggregate or direct native verification of
relocated child paths.

## Why the baseline must be competent

The synthetic diagnostic compares existing verified saved queries, the bounded
research analyzer, and Parquet with projection, predicate/file pruning, bounded
batches and controlled readahead. The columnar counts come from an independent
simple-read oracle on the same generated observations and saved masks. The cache
does not obtain missing observations that the analyzer was forbidden to read.

Content-checked and performance-only columnar measurements have separate labels.
Neither has native receipt semantics. Each prefix retains one Python worker, while
native processes remain fresh. Setup is included in total-cost estimates, but is
measured once per representation in this diagnostic; query repetitions do not
create a repeated setup distribution. Cache conditions are recorded without a
cold-cache claim. Logical retained file bytes are not a measured peak temporary
storage trace.

An observed loss is informative when its mechanism is clear. If process launches,
repeated decoding, hashing or richer integrity checks explain the result, report
their cost. Do not remove verification solely to produce a favorable speedup.
Stop optional optimization when the comparison supplies no useful advantage or
would consume the application time budget.

## Distinguish the three evidence sources

1. The unchanged three-specimen authored fixture establishes a small known-answer
   target, including measured zero and unmeasured positions. Current execution
   claims require the matching hosted report.
2. Archived HG002 reuse checks current saved-data behavior against retained
   historical independent oracle outputs for one real specimen. They are not
   newly acquired data, a fresh raw-input oracle, or a throughput/amortization study.
3. The prior 108-run representative HG002 engine study belongs to its historical
   revision and binaries, with the recorded CRAM failure and corrected rerun kept
   intact. It does not certify this extension's current revision or create new
   cohort evidence.

Genomic evidence is a different modality from neural recordings. The relevant ML
engineering lesson concerns preserving observation identity, preprocessing,
missingness and exact evaluation inputs. This work establishes no BCI decoding
accuracy, online adaptation or participant benefit.

## Evidence-backed application draft

> Extended a Rust/Python biological-evidence pipeline with exact, windowed
> reanalysis that preserves specimen identity and distinguishes missing data from
> measured zero. Passed 58 research tests and reproduced 3,000 HG002 summary rows
> across two window configurations against retained historical independent oracle
> outputs, without reopening raw reads. A 45-trial synthetic comparison showed
> that the added integrity workflow cost more than existing saved queries and a
> bounded Parquet cache, so I retained the correctness guarantees and stopped
> optional optimization.

The HG002 statement describes saved-only current reanalysis; it does not claim a
new raw extraction or newly executed oracle. The 45 trials are synthetic engineering
comparisons, not biological samples. This is a correctness-and-judgment example,
not a performance breakthrough. If the application needs a stronger throughput
result, use the separately identified historical engine study with its own tested
revision and limitations. Logan should explain the implementation choices and
tool-assisted workflow in his own words.
