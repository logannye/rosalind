# GIAB metadata preflight and explicit source-only admission

Authorized by the user on 2026-10-23 after the original bounded study stopped.
The user requested metadata preflight, resolution of GIAB input compatibility and
a hosted rerun with frozen models, thresholds and scientific gates unchanged.
The [prior stop](RESULTS.md) and its three runs remain preserved. This continuation
does not change `protocol.json`, model files, threshold policies, sources, intervals,
thinning fractions, seed or scientific decision rules.

## Cheap admission before expensive acquisition

The workflow runs a hosted preflight after tests and before the native build or
full reference/truth download. It acquires the pinned 12,006,264-byte BAI and closed
BAM ranges, capped at 64 MiB combined and 120 seconds. The fixed probe inspects at
most 64 indexed alignments from the first 1 kb of the declared test interval.
It reads no variant labels or model scores. Native record buffers contain read
data, but the probe inspects and retains only identity/flag/coordinate metadata.

The actual SAM header must match the previously retained SHA256
`2a507e6a549c0e6b62bdccb14da7c4b10dec3c7882fb04e357782672ceabc856` and its complete
ordered sequence dictionary. The exact URL, object size, Last-Modified, published
BAI MD5 and byte-range receipts remain enforced. Indexed records must have complete
seven-field Illumina QNAMEs, no original RG tags, monotone coordinates, and valid
paired read-end flags on primary alignments. No probe-region relocation is allowed.

A success-only manifest binds every preflight receipt, the actual header and BAI.
The main adapter validates those bytes, copies the verified current BAI path,
charges preflight bytes to the unchanged combined acquisition budget, and rechecks
the original receipt snapshot after analysis. The full reference checksum and
ordered dictionary must then pass before truth files are acquired. Cheap preflight
alone does not verify reference sequence or biological specimen identity.

## Missing original read groups

The source has no header read groups; the retained novoalign command line cannot
establish the identity of every library in the merged BAM. GIAB's official source
metadata binds this exact file to HG001/NA12878. The adapter admits it using that
file-level attribution, with original library/read-group provenance explicitly
unavailable. It does not infer a library from the program header.

Only derived BAMs receive a technical RG with ID `rosalind_source_` followed by the
SHA256 of the pinned source metadata, `SM:HG001`, and a description identifying it
as synthetic source grouping. No biological library ID is invented. Derived `SO`
becomes coordinate only with full selected-stream monotonicity checks and successful
index construction. Original files and headers remain unmodified and retained by
identity. Read sequence, qualities, CIGAR, QNAME, flags, coordinates and original
auxiliary tags remain unchanged; only the technical RG tag is added to derivatives.

The template key is the pinned source identity, fixed technical group and complete
unmodified QNAME. This explicitly defines the previously unsupported missing-RG
namespace; it is an intake-contract extension, not a claim that original grouping
was recovered. The existing SHA256 thinning function, seed and nested fractions
are unchanged. The namespace is fixed before candidate labels or scores are used.

[SAM](https://github.com/samtools/hts-specs/blob/master/SAMv1.tex) defines QNAME as
the template name. The [Illumina read-name convention](https://support.illumina.com/content/dam/illumina-support/help/BaseSpace_Sequence_Hub_OLH_009008_2/Source/Informatics/BS/UploadFastqReq_swBS.htm)
encodes instrument, run, flowcell, lane, tile and cluster coordinates. These names
support source-local consistency checks; they do not authenticate the donor or
prove global uniqueness across all sequencing data.

## Complete bounded-region audit

Every admitted alignment is audited before thinning. Repeated primary-end core
signatures retain physical multiplicity. Conflicting primary-end identities fail;
both observed primary mates must have compatible flags and reciprocal reference,
start and strand metadata. Secondary/supplementary alignments may repeat QNAMEs.
Unobserved mates are reported as absent from the bounded scope, without asserting
where they are. No read is deduplicated or QNAME rewritten by this audit.

Compact signatures retain no sequences or qualities. The audit stops above eight
million template entries; no entries are silently pruned. The existing process-tree
RSS and scratch ceilings remain in force. Its count must reconcile with the number
of admitted physical records, separately from traversed and excluded records.
The capped preflight is not a substitute for this full-region audit.

## Execution and interpretation

The new [execution authorization](execution-window.json) permits at most two new
hosted attempts, including at most one execution/provenance repair, until 21:20 UTC
(2:20 p.m. Pacific) on 2026-10-23. Provider run inventory and job durations are checked
before work: prior 618 runner seconds remain charged, new usage is bounded by
18,000 seconds, and overall usage by 21,600 seconds. One workflow runs at a time.
No automation is recreated. No local native builds, tests or scientific execution.

The existing policy and protocol SHA256 assertions must pass before fresh
evaluation. A new marker records the start of candidate truth-label access before
the truth index is constructed; model scoring retains its separate marker. No
scientific failure authorizes another candidate pool, threshold, model or gate.
All original absolute recall, noninferiority, practical-volume and cost requirements
remain unchanged. Even a passing screen only nominates further workflow validation.
