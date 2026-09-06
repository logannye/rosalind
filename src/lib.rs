//! # Rosalind — a deterministic, low-memory genomics engine
//!
//! Extract reusable short-read evidence with explicit filtering, bounded execution
//! state, and deterministic output. The exact indexed evidence path uses integer
//! summaries over planned tiles; the legacy pileup path retains observations and
//! fails at its declared capacity. Neither may silently sample a successful run.
//!
//! Planning, cooperative memory monitoring, and an existing OS memory limit are
//! distinct assurances. Receipts record claims and measurements for verification;
//! an unsigned receipt does not prove authorship, biological accuracy, or a hard
//! allocation bound. Custom analyzers must declare their own retained memory.
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
//! Search-index construction is deliberately separate from per-locus analysis.
//! External-memory index research is conditional on user evidence; analyzers use
//! lightweight analysis references and do not require an FM-index.

#![warn(missing_docs, missing_debug_implementations)]
#![allow(clippy::new_without_default)]

// Each module is a layer of the genomics engine.
/// The calling layer: probabilistically-grounded, abstention-aware variant calls from pileup columns.
pub mod call;
/// Offline conformance harness for external analyzer binaries.
pub mod conformance;
/// Public orchestration for inheriting planning, enforcement, and receipts.
pub mod contract;
/// Core types: the lingua franca shared by every layer (io, index, align, pileup, call).
pub mod core;
/// Canonical exact-evidence cache and bounded first-party workers.
pub mod dataset;
/// Verified streaming differences between exact evidence datasets.
pub mod dataset_diff;
/// Read-only preflight diagnostics with actionable remediation.
pub mod doctor;
/// Exact indexed read evidence and panel summaries with bounded batch execution.
pub mod evidence;
/// Genomics primitives: the FM-index, persisted memory-mapped index, alignment, sort, eval.
pub mod genomics;
/// IO layer: spec-valid VCF writer + streaming FASTA/FASTQ/BAM readers.
pub mod io;
/// Receipt-driven canonical first-party shard merge.
pub mod merge;
/// The streaming pileup kernel: one CIGAR-aware, filtered, bounded-memory engine.
pub mod pileup;
/// Reproducibility receipts: canonical-JSON BLAKE3 manifests for every run.
/// Extracted to the `rosalind-receipt` leaf crate (no htslib — wasm-friendly) and
/// re-exported here, so `rosalind::provenance::*` is unchanged.
pub use rosalind_receipt as provenance;
/// Receipt inspection, sanitization, and standards export.
pub mod receipt_tools;
/// Third-party byte re-derivation from a receipt (the `reproduce` verb).
pub mod reproduce;
/// Generate standalone downstream analyzer projects.
pub mod scaffold;
/// Canonical interval and deterministic shard selection.
pub mod selection;
/// Loopback-only embedded Receipt Studio server.
pub mod studio;
/// Helper utilities: read-only mmap + peak-RSS measurement.
pub mod util;
/// Record-preserving annotation from verified SNV evidence.
pub mod variant_annotation;
/// Checked VCF/BCF input and output.
pub mod variant_io;

// ── Genomics product surface — what builders compose on ───────────────────────
// The bounded streaming substrate:
pub use io::bam::StreamingBamSource;
pub use pileup::{Obs, PileupColumn, PileupEngine, PileupParams, ReadSource, SliceSource};
// The bounded whole-genome germline drive + calls:
pub use call::{
    call_germline_region_streaming, call_germline_selected_bam, call_germline_whole_genome,
    GermlineCall, GermlineParams,
};
// ColumnKit: implement one trait, inherit the bounded contract (SDK front door).
pub use call::{
    run_bounded_selected_bam, run_bounded_whole_genome, ColumnAnalyzer, CoverageTrack,
    FeatureAnalyzer, FeatureArrowAnalyzer, FEATURE_ARROW_BATCH_ROWS, FEATURE_ARROW_SCHEMA_VERSION,
};
// The memory contract (declare → plan → honor → verify), incl. fleet packing:
pub use call::{
    estimate_variants_working_set, first_fit_decreasing, predicted_peak_rss_bytes, PackJob,
    PackOutcome,
};
pub use core::{MemoryBudget, WorkingSet};
pub use doctor::{run_doctor, run_doctor_selected, DoctorReport, DoctorSpec};
pub use merge::{merge_shards, MergeError, MergeOutcome};
pub use receipt_tools::{export_intoto, inspect_receipt, sanitize_receipt, ReceiptInspection};
pub use selection::{AnalysisSelection, GenomicInterval, IntervalSet, SelectionError};
pub use studio::{serve_studio, StudioSpec};
// Build-once → mmap index + the reproducibility receipt:
pub use conformance::{conform_analyzer, ConformanceReport};
pub use contract::{
    detected_os_memory_limit_bytes, run_column_analysis, run_column_analysis_selected,
    AnalyzerIdentity, AnalyzerMemoryModel, ContractRunError, ContractRunOutcome, ContractRunSpec,
    ContractVerdict, EnforcementAssurance, EnforcementMode, GovernorState, OutputPolicy,
    OutputTarget, ProducerIdentity, RefusalReport, ReplayInvocation,
};
pub use genomics::{
    AnalysisReference, GenomeIndex, IndexReader, ReferencePackBuilder, ReferencePackReader,
    ReferenceProvider, ReferenceSequence, ReferenceView,
};
pub use provenance::{verify_receipt, CommandCapture, RunManifest, VerifyOpts, VerifyReport};
