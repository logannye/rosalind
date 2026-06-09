//! # Rosalind — a deterministic, low-memory genomics engine
//!
//! Call variants across a whole genome on a laptop, with memory you can **predict
//! and verify**, and results that are **byte-for-byte reproducible**. Rosalind
//! treats memory as a *contract*: you declare a RAM budget, `rosalind plan` tells
//! you up front whether the job fits, the run honors it (fits-or-refuses cleanly —
//! never a silent OOM-kill), and `rosalind verify` re-checks a BLAKE3 receipt
//! proving the realized peak landed inside your budget.
//!
//! The kernel is a streaming, CIGAR-aware **pileup column stream** bounded by local
//! coverage, not input size — a substrate you can compute arbitrary per-locus
//! analytics on. Variant calling is the first consumer, not the whole product.
//!
//! ```
//! use std::sync::Arc;
//! use rosalind::{PileupEngine, PileupParams, SliceSource};
//! use rosalind::core::{AlignedRead, CigarOp, CigarOpKind, Position, SamFlags};
//!
//! // One 4bp read "ACGT" aligned at chr0:0 over the reference "ACGT".
//! let read = AlignedRead {
//!     contig: 0,
//!     pos: Position(0),
//!     mapq: 60,
//!     flags: SamFlags(0),
//!     cigar: vec![CigarOp::new(CigarOpKind::Match, 4)],
//!     seq: Arc::from(b"ACGT".to_vec().into_boxed_slice()),
//!     qual: Arc::from(vec![40u8; 4].into_boxed_slice()),
//! };
//! let reference: Arc<[u8]> = Arc::from(b"ACGT".to_vec().into_boxed_slice());
//!
//! // The bounded pileup substrate: one PileupColumn per covered position.
//! let mut engine =
//!     PileupEngine::new(SliceSource::new(vec![read]), reference, 0, 0..4, PileupParams::default());
//! let first = engine.next().unwrap().unwrap();
//! assert_eq!(first.depth(), 1);
//! ```
//!
//! ## Research direction (Phase D)
//!
//! Rosalind is also a research vehicle for **space-bounded genomics** — sublinear-space
//! index *construction* along a `~√t` space/time curve, extending the memory contract to the
//! index build itself (today's build is `O(reference)`). That is a direction, not yet shipped;
//! it is tracked in `docs/OPEN_PROBLEMS.md`.

#![warn(missing_docs, missing_debug_implementations)]
#![allow(clippy::new_without_default)]

// Each module is a layer of the genomics engine.
/// The calling layer: probabilistically-grounded, abstention-aware variant calls from pileup columns.
pub mod call;
/// Core types: the lingua franca shared by every layer (io, index, align, pileup, call).
pub mod core;
/// Genomics primitives: the FM-index, persisted memory-mapped index, alignment, sort, eval.
pub mod genomics;
/// IO layer: spec-valid VCF writer + streaming FASTA/FASTQ/BAM readers.
pub mod io;
/// The streaming pileup kernel: one CIGAR-aware, filtered, bounded-memory engine.
pub mod pileup;
/// Reproducibility receipts: canonical-JSON BLAKE3 manifests for every run.
/// Extracted to the `rosalind-receipt` leaf crate (no htslib — wasm-friendly) and
/// re-exported here, so `rosalind::provenance::*` is unchanged.
pub use rosalind_receipt as provenance;
/// Third-party byte re-derivation from a receipt (the `reproduce` verb).
pub mod reproduce;
/// Helper utilities: read-only mmap + peak-RSS measurement.
pub mod util;

// ── Genomics product surface — what builders compose on ───────────────────────
// The bounded streaming substrate:
pub use io::bam::StreamingBamSource;
pub use pileup::{Obs, PileupColumn, PileupEngine, PileupParams, ReadSource, SliceSource};
// The bounded whole-genome germline drive + calls:
pub use call::{
    call_germline_region_streaming, call_germline_whole_genome, GermlineCall, GermlineParams,
};
// ColumnKit: implement one trait, inherit the bounded contract (SDK front door).
pub use call::{run_bounded_whole_genome, ColumnAnalyzer, CoverageTrack, FeatureAnalyzer};
// The memory contract (declare → plan → honor → verify), incl. fleet packing:
pub use call::{
    estimate_variants_working_set, first_fit_decreasing, predicted_peak_rss_bytes, PackJob,
    PackOutcome,
};
pub use core::{MemoryBudget, WorkingSet};
// Build-once → mmap index + the reproducibility receipt:
pub use genomics::{GenomeIndex, IndexReader, ReferenceView};
pub use provenance::{verify_receipt, CommandCapture, RunManifest, VerifyOpts, VerifyReport};
