# Analysis reference pack (`.rref`) format

Version 1 is a deterministic little-endian mmap format with a fixed 128-byte header,
an ordered contig table, alignment padding, and 64-base sequence blocks. Each block
contains two little-endian 2-bit words (A/C/G/T) followed by one 64-bit ambiguity
mask (`N`). The header records format version, section extents, contig/base counts,
the BLAKE3 of normalized concatenated source bases, and a BLAKE3 over all bytes after
the header.

Builders normalize case, convert U to T, reject unsupported bases, duplicate names,
empty records, oversized contigs, and lines beyond the fixed 1 MiB parser ceiling.
They scan FASTA twice and reject a source that changes between passes. Output uses a
same-directory temporary file and atomic publication.

`reference build` and `reference convert` also create a schema-5 sidecar receipt by
default. The receipt hashes the source and `.rref`, records normalized source
identity and total bases, verifies offline, and can replay the deterministic pack
byte-for-byte.

Readers reject bad magic, unsupported version/endian, corrupt offsets, truncation,
checksum mismatch, invalid UTF-8 names, inconsistent global offsets, or a contig
table whose lengths disagree with the declared reference length before exposing any
base view.

`.rref` intentionally has no suffix array, BWT, or FM-index. Use `.idx` only for
search/alignment workloads. `reference convert` extracts normalized reference bytes
and source identity from a validated legacy index.
