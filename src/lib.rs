//! # O(√t) Space Simulation via Height Compression
//!
//! This library implements the breakthrough algorithm for simulating
//! time-t Turing machine computations in O(√t) space.
//!
//! ## Core Algorithm
//!
//! 1. **Block-respecting simulation**: Partition computation into T = ⌈t/b⌉ blocks
//! 2. **Height compression**: Transform tree from height Θ(T) → O(log T)
//! 3. **Pointerless evaluation**: O(1) bits per level instead of O(log b)
//! 4. **Streaming ledger**: Track T merges with constant tokens
//!
//! Result: Space = O(b + T + log T) = O(b + t/b), optimal at b = √t
//!
//! ## Usage Example
//!
//! ```ignore
//! use sqrt_space_sim::{TuringMachine, Simulator, SimulationConfig};
//!
//! let config = SimulationConfig::optimal_for_time(10_000);
//! let mut sim = Simulator::new(machine, config);
//! let result = sim.run(&input)?;
//! assert!(result.space_used <= O(√10_000));
//! ```

#![warn(missing_docs, missing_debug_implementations)]
#![allow(clippy::new_without_default)]

// Core modules - each implements a key component of the algorithm
pub mod machine;    // Turing machine representation
pub mod blocking;   // Block-respecting simulation
pub mod tree;       // Height-compressed evaluation tree
pub mod algebra;    // Algebraic replay engine
pub mod ledger;     // Streaming progress tracking
pub mod space;      // Space accounting utilities
pub mod util;       // Helper functions
pub mod framework;  // Generic compressed evaluation
pub mod genomics;   // Genomics primitives and algorithms
pub mod plugin;     // Plugin system
/// Python bindings for exposing Rosalind components to external runtimes.
pub mod python_bindings;

// Re-exports for convenience
pub use machine::{TuringMachine, Configuration, Transition, Move, Symbol, State};
pub use blocking::{BlockSummary, MovementLog};
pub use tree::{CompressedTree, TreeNode};
pub use algebra::{AlgebraicEngine, FiniteField};
pub use ledger::StreamingLedger;

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
        self.block_size + self.num_blocks + 
            (self.num_blocks as f64).log2().ceil() as usize
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
        let final_summary = engine.evaluate_dfs(
            root,
            input,
            &mut self.space_tracker,
        )?;
        
        // Extract result
        let accepted = final_summary.exit_state() == self.machine.accept_state();
        let space_used = self.space_tracker.max_space_used();
        
        // Verify ledger completion (all merges should be complete)
        if !engine.ledger().all_merges_complete() {
            return Err(SimulationError::AlgebraError(
                "Not all merges completed in streaming ledger".to_string()
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