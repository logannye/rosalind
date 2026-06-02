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
//! Rosalind is also a research vehicle for **space-bounded genomics**: a `~√t`
//! (square-root-space) evaluation framework (Williams 2025; Cook–Mertz 2024) as a
//! continuous space/time knob, aimed at **sublinear-space index construction**.
//! That layer is future work — not yet load-bearing — tracked in
//! `docs/OPEN_PROBLEMS.md`.

#![warn(missing_docs, missing_debug_implementations)]
#![allow(clippy::new_without_default)]

// Core modules - each implements a key component of the algorithm
pub mod algebra; // Algebraic replay engine
pub mod blocking; // Block-respecting simulation
/// The calling layer: calibrated, abstention-aware variant calls from pileup columns.
pub mod call;
/// Core types: the lingua franca shared by every layer (io, index, align, pileup, call).
pub mod core;
pub mod framework; // Generic compressed evaluation
pub mod genomics; // Genomics primitives and algorithms
/// IO layer: spec-valid VCF writer (FASTA/FASTQ/BAM readers arrive in later phases).
pub mod io;
pub mod ledger; // Streaming progress tracking
pub mod machine; // Turing machine representation
/// The streaming pileup kernel: one CIGAR-aware, filtered, bounded-memory engine.
pub mod pileup;
pub mod plugin; // Plugin system
/// Reproducibility receipts: canonical-JSON BLAKE3 manifests for every run.
pub mod provenance;
/// Python bindings for exposing Rosalind components to external runtimes.
#[cfg(feature = "python-bindings")]
pub mod python_bindings;
pub mod space; // Space accounting utilities
pub mod tree; // Height-compressed evaluation tree
pub mod util; // Helper functions

// ── Genomics product surface — what builders compose on ───────────────────────
// The bounded streaming substrate:
pub use io::bam::StreamingBamSource;
pub use pileup::{Obs, PileupColumn, PileupEngine, PileupParams, ReadSource, SliceSource};
// The bounded whole-genome germline drive + calls:
pub use call::{
    call_germline_region_streaming, call_germline_whole_genome, GermlineCall, GermlineParams,
};
// The memory contract (declare → plan → honor → verify):
pub use call::{estimate_variants_working_set, predicted_peak_rss_bytes};
pub use core::{MemoryBudget, WorkingSet};
// Build-once → mmap index + the reproducibility receipt:
pub use genomics::{GenomeIndex, IndexReader, ReferenceView};
pub use provenance::RunManifest;

// ── Research layer (√t space-bounded simulation; Phase D — see OPEN_PROBLEMS) ──
pub use algebra::{AlgebraicEngine, FiniteField};
pub use blocking::{BlockSummary, MovementLog};
pub use ledger::StreamingLedger;
pub use machine::{Configuration, Move, State, Symbol, Transition, TuringMachine};
pub use tree::{CompressedTree, TreeNode};

use thiserror::Error;

/// Main simulation orchestrator
///
/// Coordinates all components to achieve O(√t) space bound
#[derive(Debug)]
pub struct Simulator {
    machine: TuringMachine,
    config: SimulationConfig,
    space_tracker: space::SpaceTracker,
}

/// Configuration parameters for simulation
#[derive(Debug, Clone)]
pub struct SimulationConfig {
    /// Block size parameter b (optimal: √t)
    pub block_size: usize,

    /// Time bound t
    pub time_bound: usize,

    /// Number of blocks T = ⌈t/b⌉
    pub num_blocks: usize,

    /// Field size for algebraic operations (constant)
    pub field_characteristic: u8,

    /// Enable space profiling
    pub profile_space: bool,

    /// Enable detailed logging
    pub verbose: bool,
}

impl SimulationConfig {
    /// Create optimal configuration: b = ⌈√t⌉
    pub fn optimal_for_time(time_bound: usize) -> Self {
        let block_size = (time_bound as f64).sqrt().ceil() as usize;
        let num_blocks = (time_bound + block_size - 1) / block_size;

        Self {
            block_size,
            time_bound,
            num_blocks,
            field_characteristic: 8, // 𝔽_{2^8} constant
            profile_space: false,
            verbose: false,
        }
    }

    /// Theoretical space bound: O(b + t/b + log(t/b))
    pub fn space_bound(&self) -> usize {
        self.block_size + self.num_blocks + (self.num_blocks as f64).log2().ceil() as usize
    }

    /// Simplified O(√t) bound
    pub fn sqrt_t_bound(&self) -> usize {
        (self.time_bound as f64).sqrt().ceil() as usize * 2
    }
}

/// Result of simulation
#[derive(Debug, Clone)]
pub struct SimulationResult {
    /// Whether machine accepted
    pub accepted: bool,

    /// Final configuration
    pub final_config: Configuration,

    /// Space used (in cells)
    pub space_used: usize,

    /// Time steps simulated
    pub time_steps: usize,

    /// Space profile (if enabled)
    pub space_profile: Option<space::SpaceProfile>,
}

impl SimulationResult {
    /// Verify result satisfies theoretical bound
    pub fn satisfies_bound(&self, config: &SimulationConfig) -> bool {
        self.space_used <= config.sqrt_t_bound()
    }
}

/// Errors that can occur during simulation
#[derive(Error, Debug)]
pub enum SimulationError {
    /// Invalid machine configuration
    #[error("Invalid machine configuration: {0}")]
    InvalidMachine(String),

    /// Invalid block size for given time bound
    #[error("Invalid block size {0} for time bound {1}")]
    InvalidBlockSize(usize, usize),

    /// Interface verification failed at specified block
    #[error("Interface verification failed at block {0}")]
    InterfaceCheckFailed(usize),

    /// Algebraic operation failed
    #[error("Algebraic operation failed: {0}")]
    AlgebraError(String),

    /// Space bound exceeded
    #[error("Space bound exceeded: used {used} > bound {bound}")]
    SpaceBoundExceeded {
        /// Space actually used
        used: usize,
        /// Space bound that was exceeded
        bound: usize,
    },
}

impl Simulator {
    /// Create new simulator
    pub fn new(machine: TuringMachine, config: SimulationConfig) -> Self {
        Self {
            machine,
            config: config.clone(),
            space_tracker: space::SpaceTracker::new(config.profile_space),
        }
    }

    /// Run simulation on input
    ///
    /// This is the main entry point that orchestrates:
    /// 1. Height-compressed tree construction
    /// 2. Pointerless DFS evaluation
    /// 3. Space bound verification
    pub fn run(&mut self, input: &[Symbol]) -> Result<SimulationResult, SimulationError> {
        // Validate configuration
        if self.config.block_size == 0 {
            return Err(SimulationError::InvalidBlockSize(0, self.config.time_bound));
        }

        // Initialize algebraic replay engine
        let mut engine = algebra::AlgebraicEngine::new(&self.config, &self.machine);

        // Allocate ledger space
        let ledger_size = self.config.num_blocks * 2 / 8; // 2 bits per block
        self.space_tracker.allocate_ledger(ledger_size);

        // Create implicit compressed tree root
        let root = tree::TreeNode::root(1, self.config.num_blocks);

        // Evaluate root via DFS
        let final_summary = engine.evaluate_dfs(root, input, &mut self.space_tracker)?;

        // Extract result
        let accepted = final_summary.exit_state() == self.machine.accept_state();
        let space_used = self.space_tracker.max_space_used();

        // Verify ledger completion (all merges should be complete)
        if !engine.ledger().all_merges_complete() {
            return Err(SimulationError::AlgebraError(
                "Not all merges completed in streaming ledger".to_string(),
            ));
        }

        // Verify bound
        if space_used > self.config.sqrt_t_bound() {
            return Err(SimulationError::SpaceBoundExceeded {
                used: space_used,
                bound: self.config.sqrt_t_bound(),
            });
        }

        Ok(SimulationResult {
            accepted,
            final_config: final_summary.into_configuration(input, self.machine.blank()),
            space_used,
            time_steps: self.config.time_bound,
            space_profile: self.space_tracker.take_profile(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_optimal_block_size() {
        let config = SimulationConfig::optimal_for_time(10_000);
        assert_eq!(config.block_size, 100);
        assert_eq!(config.num_blocks, 100);
    }

    #[test]
    fn test_space_bound_formula() {
        let config = SimulationConfig::optimal_for_time(10_000);
        let bound = config.space_bound();
        // O(b + T + log T) = 100 + 100 + 7 = 207
        assert!(bound <= 250); // Allow some constant factor
    }
}
