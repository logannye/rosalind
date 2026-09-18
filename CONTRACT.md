# The resource and evidence contract

Rosalind separates the evidence a scientific question requires from the resources
used to extract it. A smaller admitted budget can mean smaller computation tiles
and more indexed I/O; it cannot silently remove eligible observations from a
successful exact-evidence result.

This guide describes the current source and unpublished evidence candidate. Public
stable v0.1.0 predates these APIs. See [implementation status](docs/implementation-status.md)
for validation and publication evidence.

## Declare and plan the actual analysis

New analyses should use `analyze evidence`, `analyze panel-qc`, or the
[managed evidence SDK](docs/analyzer-sdk.md). With local indexed BAM and FASTA:

```sh
rosalind analyze evidence --reference genome.fa --alignments sample.sorted.bam \
  --sites candidates.vcf --fields depths,alleles --memory-budget-mb 256 --plan
rosalind analyze evidence --reference genome.fa --alignments sample.sorted.bam \
  --sites candidates.vcf --fields depths,alleles --memory-budget-mb 256 \
  --format arrow-ipc --output evidence.arrow
```

The evidence CLI uses a supplied budget for admission and cooperative governance.
Its plan includes the actual fields, decoder, reference windows, metadata, consumer
and encoder. A budget below the required working set refuses the run. Plans do real
input work: CRAM planning includes a complete sequential record-validation pass.
See [CRAM decoder admission](docs/SEMANTICS.md#cram-decoder-admission).

The legacy top-level `rosalind plan` models legacy pileup using reference metadata
and declared read/depth capacities. It does not substitute for an evidence or
custom-analyzer plan. Legacy commands and the SDK's `RecordOnly` mode distinguish
recorded budgets from enforcement; the generated evidence analyzer uses `--enforce`
to request cooperative enforcement.

## Preserve scientific meaning

The evidence engine aggregates indexed reads into bounded tiles instead of keeping
all active observations. It never downsamples to meet a budget. Successful evidence
bytes are invariant to admitted budgets, microtiles and worker counts for the same
scientific request and producer. Canonical ownership, fixed Arrow batches, checked
integer arithmetic and ordered reduction preserve that invariant.

Fields, reference, sample scope, filters and loci are scientific choices. Tile width,
workers and repeated indexed-record visits are execution choices. Omitted fields
remain absent; an unmeasured locus differs from a stored zero-depth row. Saved
evidence cannot recover excluded observations or change sample scope. See
[SEMANTICS.md](docs/SEMANTICS.md).

The unit is an aligned read base, not a molecule. Overlapping mates count separately.
The profile does not provide UMI consensus, phasing, methylation, local realignment
or clinical interpretation. The built-in SNV caller remains experimental.

## Admission and runtime assurance

Prediction is a declared working-set model, not a proof about every allocation.
Unknown custom-analyzer bounds are permitted in observation mode and refused under
enforcement. A bound must cover retained reducer state and flush transients; the
runner separately reserves its own buffers and receipt state.

- Admission refusal precedes creation of the requested successful output.
- Cooperative checks and sampled RSS can detect a runtime breach. Managed artifact
  runners identify retained resource-failure output as partial; failed work never
  becomes a successful artifact/receipt pair.
- Native decoders can allocate before the next checkpoint. CRAM's whole-file
  validation also has this limitation before indexed extraction. A declared
  envelope and sampled peak do not prove a universal hard RAM cap or guarantee
  recovery before an OS termination.
- `--require-os-limit` additionally checks an existing finite Linux cgroup-v2
  `memory.max` no greater than the declared budget. Rosalind does not create the
  cgroup. A scheduler memory request alone is not this assurance.
- Python-retained batches and downstream R/SQL allocations are outside the native
  child-process budget. Metadata grows with sites, targets and partition inventory.

The governor and managed cancellation scope are process-wide: serialize independent
managed runs in a process. Evidence workers produce first-party partitions; custom
reducers consume them in canonical order. Separate processes are another boundary.

## Legacy compatibility

`features`, `analyze coverage`, experimental `variants` and `ColumnAnalyzer` retain
legacy scientific defaults. Current pileup execution is **exact-or-fail**: exceeding
declared read/depth capacity fails instead of selecting a subset. Historical
sampling receipts keep their original verification capabilities; replaying old
behavior requires the matching producer.

Legacy regions/reference-span shards have deterministic locus ownership and
first-party merge. Reference-pack construction uses bounded buffers; legacy
FM-index construction still uses memory proportional to the reference. Neither
index building nor region-retaining somatic routines inherit the evidence tile
model. See [architecture](ARCHITECTURE.md).

## Verify, replay and publish

```sh
rosalind verify --manifest evidence.arrow.manifest.json
rosalind reproduce --manifest evidence.arrow.manifest.json --inputs input-directory
```

Verification checks receipt integrity and recorded artifact bytes; replay reruns
an operation and compares output bytes. Scientific identity is narrower than the
whole receipt claim: resource declarations and replay settings may change a claim
while successful scientific output remains identical. Host-local measurements
have a separate integrity layer. Paths are relocatable metadata in supported
receipt schemas.

Materialized artifacts use staged output and atomic publication. Create-new is the
default; `--force` explicitly requests replacement where supported. Failure and
cancellation do not publish success. Retain the dependencies required by each
receipt: a JSON sidecar cannot verify missing artifact bytes. Portable dataset
queries instead rely on saved evidence and its lineage without reopening original
alignments.

A receipt is tamper-evident, not a signature, proof of authorship or biological
accuracy. Read [receipt trust](docs/receipts-and-trust.md),
[determinism](docs/determinism.md) and [reusable evidence](docs/reusable-evidence.md).

Start a new reducer with the [builder quickstart](docs/builder-quickstart.md).
Lower-level `EvidenceEngine` and `run_bounded_whole_genome` calls do not themselves
supply managed publication, receipts or a bound for arbitrary custom state.
