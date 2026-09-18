# Determinism contract

For the same scientific inputs, reducer and producer, successful exact-evidence
artifacts have deterministic bytes. Admitted budgets, microtiles and workers do
not change those bytes. This describes the unpublished evidence source/candidate;
it does not retrofit its semantics into public stable v0.1.0 or old receipts.

## Scientific and execution identity

Scientific inputs include alignment/reference content, normalized loci, requested
ALT annotations, sample scope, quality/flag profile, fields, schema and scientific
analyzer parameters. Original variant-record identity also matters for
record-preserving annotation. Offline queries retain their stored extraction profile.

Execution settings include budget, microtiles, workers and cache state. Smaller
tiles can repeat indexed reads, without sampling. Canonical ownership and fixed
encoding batches separate scheduling from successful evidence bytes. Verification
and CRAM whole-file validation still perform real work, including on reused paths.

Keep three identities distinct:

- **Scientific identity** describes what is counted and reduced.
- **Receipt claim** also binds outputs, producer/build identity and replay
  declarations. A resource declaration can change this full claim while the
  scientific output stays identical.
- **Measurements** describe host-local time, memory and work. Their separately
  protected values can differ. Whole receipt JSON is not promised to match across
  machines.

Recorded paths are relocatable metadata in supported schemas. Verification/replay
still require the relevant file bytes. Claims are tamper-evident, not signatures.

## Covered artifacts and limits

Maintained equality checks cover evidence Arrow/TSV, panel summaries,
record-preserving annotation and managed deterministic reducer outputs. Legacy
feature/caller outputs retain their own defaults and fixtures. The evidence path
consumes BAM/CRAM; it does not claim deterministic alignment from FASTQ.

Replay uses the matching producer. Historical sampling, changed scientific profiles
and different encoders are not implicitly equivalent. Cross-platform checks are
evidence for tested builds, not every compiler, library or hardware combination.
Parquet exports preserve integer values and lineage; native byte replay of a
Parquet directory is unsupported. Use materialized Arrow/TSV for byte replay.

Failed/cancelled runs are not successful scientific artifacts. Breach timing may
vary and partials have no success-byte guarantee. Native allocation can precede a
checkpoint, including during CRAM validation. See [the contract](../CONTRACT.md).

## Implementation rules

1. Use dictionary, coordinate and documented allele/target order. Sort unordered
   keys before output and define tie-breakers.
2. Reduce checked integers in canonical order. For floating-point calculations,
   specify reduction order, missingness and serialization.
3. Keep timestamps, random IDs and host measurements out of scientific output.
   Absent loci/fields never become zeros.
4. Separate worker scheduling from output ownership and reducer order. Custom
   reducers receive canonical batches independent of computation tiles.
5. Record all scientific parameters in metadata and replay. Declare retained
   state, including copies kept beyond a borrowed callback.

## Verification

Maintain scientific oracles and repeat-run checks, plus budget/tile/worker,
projection and native/persisted equality for affected paths. Refusal/corruption
must leave no successful artifact. External analyzer conformance tests the runner;
the analyzer author supplies a scientific oracle for a new statistic.

Replay passes tokenized `command_argv` directly without a shell. External receipts
require an explicit `--binary`; receipts cannot implicitly choose executables.
See [the SDK](analyzer-sdk.md), [receipt schema](receipt-schema.md) and
[receipt trust](receipts-and-trust.md).
