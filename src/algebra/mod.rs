//! Algebraic Replay Engine (ARE)
//!
//! Evaluates compressed tree with O(1) per-level storage:
//! - Constant-size field 𝔽_{2^c}
//! - Low-degree polynomial extensions
//! - Constant evaluation grid

mod field;
mod polynomial;
mod combiner;

pub use field::FiniteField;
pub use polynomial::{PolynomialEncoding, EvaluationGrid};
pub use combiner::Combiner;

use crate::{
    blocking::{BlockSummary, simulate_block, BlockId},
    tree::TreeNode,
    ledger::StreamingLedger,
    space::SpaceTracker,
    machine::{Symbol, TuringMachine, Configuration},
    SimulationConfig, SimulationError,
};
use std::collections::HashMap;

/// Algebraic Replay Engine
///
/// Core evaluation with space tracking
#[derive(Debug)]
pub struct AlgebraicEngine {
    config: SimulationConfig,
    #[allow(dead_code)]
    field: FiniteField,
    grid: EvaluationGrid,
    ledger: StreamingLedger,
    combiner: Combiner,
    machine: TuringMachine,
    /// Cache of block summaries for boundary reconstruction
    /// Maps block_id -> BlockSummary
    boundary_cache: HashMap<BlockId, BlockSummary>,
}

impl AlgebraicEngine {
    /// Create new engine
    pub fn new(config: &SimulationConfig, machine: &TuringMachine) -> Self {
        let field = FiniteField::new(config.field_characteristic);
        let grid = EvaluationGrid::new(&field);
        let ledger = StreamingLedger::new(config.num_blocks);
        let combiner = Combiner::new(&field);
        
        Self {
            config: config.clone(),
            field,
            grid,
            ledger,
            combiner,
            machine: machine.clone(),
            boundary_cache: HashMap::new(),
        }
    }
    
    /// Get reference to ledger for verification
    pub fn ledger(&self) -> &StreamingLedger {
        &self.ledger
    }
    
    /// Evaluate tree via pointerless DFS
    ///
    /// Space breakdown:
    /// - Path stack: O(log T) cells
    /// - Streaming ledger: O(T) cells  
    /// - Leaf buffer: O(b) cells (one at a time)
    /// - Field accumulators: O(1) cells
    /// Total: O(b + T + log T) = O(b + t/b)
    pub fn evaluate_dfs(
        &mut self,
        node: TreeNode,
        input: &[Symbol],
        tracker: &mut SpaceTracker,
    ) -> Result<BlockSummary, SimulationError> {
        // Track space for path token (2 bits per level = 1 byte)
        tracker.push_stack_frame(1);
        
        if node.is_leaf() {
            // Base case: leaf → evaluate_leaf()
            let summary = self.evaluate_leaf(node, input, tracker)?;
            
            // Cache the summary for boundary reconstruction of subsequent blocks
            let block_id = node.leaf_block_id();
            self.boundary_cache.insert(block_id, summary.clone());
            
            tracker.pop_stack_frame();
            return Ok(summary);
        }
        
        // Recursive case: internal node
        let (left_child, right_child) = node.children();
        
        // Evaluate left subtree recursively
        let left_summary = self.evaluate_dfs(left_child, input, tracker)?;
        
        // Mark progress in streaming ledger
        self.ledger.mark_left_complete(node);
        
        // Evaluate right subtree recursively
        let right_summary = self.evaluate_dfs(right_child, input, tracker)?;
        
        // Mark progress
        self.ledger.mark_right_complete(node);
        
        // Merge summaries via combiner
        // CRITICAL: This uses exact interface checking + algebraic encoding
        let merged = self.combiner.merge(
            &left_summary,
            &right_summary,
            &self.grid,
            tracker,
        )?;
        
        tracker.pop_stack_frame();
        Ok(merged)
    }
    
    /// Evaluate leaf block
    ///
    /// Space: O(b) for replay buffer + O(1) for accumulators
    fn evaluate_leaf(
        &mut self,
        node: TreeNode,
        input: &[Symbol],
        tracker: &mut SpaceTracker,
    ) -> Result<BlockSummary, SimulationError> {
        let block_id = node.leaf_block_id();
        
        // Allocate O(b) buffer
        tracker.allocate_leaf_buffer(self.config.block_size);
        
        // Reconstruct block summary by simulation
        // Get initial configuration
        let config = if block_id == 1 {
            // First block: use initial configuration
            Configuration::initial(input, self.machine.num_tapes() + 1) // +1 for input tape
        } else {
            // Get previous block's summary from cache
            let prev_block_id = block_id - 1;
            let prev_summary = self.boundary_cache.get(&prev_block_id)
                .ok_or_else(|| SimulationError::InvalidMachine(
                    format!("Missing boundary for block {}", prev_block_id)
                ))?;
            
            // Reconstruct configuration from previous block
            self.reconstruct_boundary_config(prev_summary, input, self.machine.num_tapes() + 1)
        };
        
        // Simulate block
        let summary = simulate_block(
            &self.machine,
            &config,
            block_id,
            self.config.block_size,
        )?;
        
        // Free buffer
        tracker.free_leaf_buffer(self.config.block_size);
        
        Ok(summary)
    }
    
    /// Reconstruct boundary configuration from previous block's summary
    ///
    /// This is used to chain blocks together - each block after the first
    /// starts from the previous block's exit state.
    fn reconstruct_boundary_config(
        &self,
        previous_summary: &BlockSummary,
        input: &[Symbol],
        _num_tapes: usize,
    ) -> Configuration {
        // Reconstruct configuration from previous block's summary
        // Clone summary since into_configuration takes ownership
        let config = previous_summary.clone().into_configuration(input, self.machine.blank());
        
        // Verify state matches (sanity check)
        assert_eq!(config.state(), previous_summary.exit_state());
        
        // Verify head positions match (sanity check)
        let head_positions = config.head_positions();
        assert_eq!(head_positions, previous_summary.exit_heads());
        
        config
    }
}