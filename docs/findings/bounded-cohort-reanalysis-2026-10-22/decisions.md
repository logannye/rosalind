# Decisions and falsifiable mechanisms

## Preserve the observation mask

The fixture contains raw reads at specimen B:40, but its portable leaf does not
store that position. Reading those raw observations during partial reanalysis
would change the question. The reducer therefore visits only native-normalized
covered intersections, counts stored zero-depth rows as observations, and retains
the requested and observed denominators. The unchanged authored fixture and
independent SAM counting are the small falsifying experiment: a filled missing
position or a dropped measured-zero row must change the expected answer and fail.

## Use a separate comparison contract

Different specimens normally have different source reuse identities. Requiring
the same reuse key would reject useful comparisons; weakening that key would
damage existing saved-data guarantees. Keep native reuse unchanged and compare
measurement semantics, required fields, ordered contigs and effective reference
role/bytes separately. Reject an aliased duplicate alignment and named scope.
Contract tests independently vary source, scope, profile, mask and reference.

## Pin generated execution inputs

Review found that a generated coverage report was reread after preflight without
being pinned to its original bytes. Removing a future member's covered position
could otherwise alter the observation denominator while the final inventory
merely hashed the modified report. Coverage and descriptor copies now enter the
input ledger immediately after creation. The hosted mutation test changes a
future coverage report after the first successful child and requires refusal
before publication. Earlier consumed partition mutation is a separate test.

The same boundary rejects FIFOs, devices and directories before hashing. A file's
reported zero size does not imply a finite readable regular file. Pure tests
ensure rejection occurs before opening, without a potentially blocking FIFO read.

## Correct the timing comparison before interpreting it

Early development runs gave C an extra Python interpreter for every request,
while B/D reused their worker within a prefix. That asymmetry is removed in the
frozen repeat: all paths retain one Python worker, and each native child remains
fresh. D projects only position, callable depth and the requested ALT. Its bounded
scanner uses predicate/file pruning, row-group statistics, one batch/fragment of
readahead, and no scan threads. Performance-only D performs no query file hashing;
content-checked D includes explicitly labeled before/after checks.

The early trials remain in `development/` with their own identities. They are not
combined with the frozen measurements or used as application speed results.

## Preserve native telemetry representation

The first frozen rerun completed native build/tests and the existing Python suite,
but the retained demo refused because a new memory-metric guard expected a JSON
integer. Native receipt measurements are a string map, and `peak_rss_bytes` is an
unsigned decimal string. The parser now accepts bounded ASCII decimal strings
without floating-point conversion and preserves absent values as null. A pure
regression covers values above 2^53, invalid forms and absence; the empty-query
integration checks that no child launch or invented peak appears. The failed run
39815498177 and its stderr remain in `development/`. No scientific count was
published by that failed demo.

## Test the window mechanism within the declared pilot

Smaller windows retain fewer reducer counters but can revisit one physical
partition across multiple children. The predicted effect is more process launches,
hashing and decoding, potentially in exchange for lower parent memory. The
smallest comparison is already frozen in the pilot: W=64 versus W=256 over the
same three members, masks, queries and thresholds, with three rotating repetitions.
Inspect the full setup-plus-query totals, sampled whole-worker peak, launch counts
and repeated partition decodes together. Do not infer a decoder-only timing from
these combined costs or add overlapping internal timers.

The final report records the result and keep/stop decision. Query coalescing is
not silently credited to this sequentially revealed workload. It would require a
separate known-in-advance family and an equivalently coalesced columnar baseline.
No integrity checks are removed to manufacture a favorable comparison.

## Preserve the resource and evidence boundaries

No ready compatible raw HG002 cache was present in the inspected local locations.
A tracked historical archive contains a small portable leaf and independent
oracle outputs, so bounded extraction from that archive provides current
saved-data correctness coverage without bulk acquisition. It cannot supply a
new fresh-extraction baseline, a new raw oracle, or raw-versus-saved amortization.

The pilot is diagnostic and uses synthetic observations. Setup is measured once,
cache state is uncontrolled, three rotations are not a fully balanced order
design, and 10 ms RSS sampling can miss short peaks. Retained file sizes describe
logical storage; accumulated outputs are not peak temporary occupancy. A failed
performance hypothesis therefore supports a narrow stop decision, not a claim
that every larger workload must behave the same way.
