//! Block summary representation
//!
//! Each σ_k encodes:
//! - Entry/exit states and head positions
//! - Movement log (deterministic micro-ops)
//! - Window bounds per tape
//! - Advisory fingerprints (not used for correctness)

use super::{BlockId, WindowBounds};
use crate::machine::{State, Symbol, Move};

/// Micro-operation: (tape_id, position_offset, symbol_written, move_direction)
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MicroOp {
    pub tape_id: usize,
    pub offset: i64,
    pub symbol: Symbol,
    pub move_dir: Move,
}

/// Summary of one time block
///
/// Space: O(b) cells total
#[derive(Debug, Clone)]
pub struct BlockSummary {
    /// Block identifier
    pub block_id: BlockId,
    
    /// Entry state at start of block
    entry_state: State,
    
    /// Exit state at end of block
    exit_state: State,
    
    /// Head positions at start of block (all tapes including input)
    entry_heads: Vec<i64>,
    
    /// Head positions at end of block (all tapes including input)
    exit_heads: Vec<i64>,
    
    /// Movement log: sequence of micro-operations (≤ b entries)
    movement_log: MovementLog,
    
    /// Window bounds per tape (one per tape)
    windows: Vec<WindowBounds>,
    
    /// Advisory fingerprints (optional, not used for correctness)
    #[allow(dead_code)]
    fingerprints: Vec<u64>,
}

impl BlockSummary {
    /// Create new block summary
    pub fn new(
        block_id: BlockId,
        entry_state: State,
        exit_state: State,
        entry_heads: Vec<i64>,
        exit_heads: Vec<i64>,
        movement_log: MovementLog,
        windows: Vec<WindowBounds>,
    ) -> Self {
        Self {
            block_id,
            entry_state,
            exit_state,
            entry_heads,
            exit_heads,
            movement_log,
            windows,
            fingerprints: Vec::new(),
        }
    }
    
    /// Create default (for testing)
    pub fn default() -> Self {
        Self {
            block_id: 1,
            entry_state: 0,
            exit_state: 0,
            entry_heads: vec![0],
            exit_heads: vec![0],
            movement_log: MovementLog::new(),
            windows: vec![WindowBounds { left: 0, right: 0 }],
            fingerprints: Vec::new(),
        }
    }
    
    /// Get entry state
    pub fn entry_state(&self) -> State {
        self.entry_state
    }
    
    /// Get exit state
    pub fn exit_state(&self) -> State {
        self.exit_state
    }
    
    /// Get entry head positions
    pub fn entry_heads(&self) -> &[i64] {
        &self.entry_heads
    }
    
    /// Get exit head positions
    pub fn exit_heads(&self) -> &[i64] {
        &self.exit_heads
    }
    
    /// Get movement log
    pub fn movement_log(&self) -> &MovementLog {
        &self.movement_log
    }
    
    /// Get window bounds
    pub fn windows(&self) -> &[WindowBounds] {
        &self.windows
    }
    
    /// Convert to full configuration (for final result)
    ///
    /// Reconstructs complete configuration from block summary by replaying movement log.
    /// This is used for boundary reconstruction when chaining blocks.
    pub fn into_configuration(self, input: &[Symbol], blank_symbol: Symbol) -> crate::machine::Configuration {
        use crate::machine::{Configuration, Tape};
        
        let num_tapes = self.entry_heads.len();
        
        // Create input tape (read-only)
        let mut input_tape = Tape::blank(blank_symbol);
        for (i, &symbol) in input.iter().enumerate() {
            input_tape.write(i as i64, symbol);
        }
        
        // Create work tapes and reconstruct contents by replaying movement log
        let mut tapes = vec![input_tape];
        
        // For each work tape (tape_id 1..num_tapes-1)
        for tape_id in 1..num_tapes {
            // Get entry head position for this tape
            let entry_head_pos = self.entry_heads[tape_id];
            
            // Replay movement log to reconstruct tape contents
            let tape_contents = self.movement_log.replay_full_tape(tape_id, entry_head_pos, blank_symbol);
            
            // Create tape with blank symbol
            let mut tape = Tape::blank(blank_symbol);
            
            // Fill in non-blank cells from replay
            for (pos, symbol) in tape_contents {
                tape.write(pos, symbol);
            }
            
            tapes.push(tape);
        }
        
        // Create configuration with reconstructed tapes
        let mut config = Configuration::initial(input, num_tapes);
        config.set_state(self.exit_state);
        config.set_head_positions(self.exit_heads[0], self.exit_heads[1..].to_vec());
        config.set_tapes(tapes);
        
        config
    }
    
    /// Space usage in cells
    pub fn space_usage(&self) -> usize {
        // Movement log is the dominant term (≤ b entries)
        self.movement_log.space_usage()
            // Plus O(1) for states, head positions, windows
            + 2 * std::mem::size_of::<State>() / std::mem::size_of::<Symbol>()
            + self.entry_heads.len() * std::mem::size_of::<i64>() / std::mem::size_of::<Symbol>()
            + self.exit_heads.len() * std::mem::size_of::<i64>() / std::mem::size_of::<Symbol>()
            + self.windows.len() * std::mem::size_of::<WindowBounds>() / std::mem::size_of::<Symbol>()
    }
}

/// Movement log: sequence of micro-operations
///
/// Each micro-op = (tape_id, position_offset, symbol, move)
/// Total: ≤ b operations, O(b) cells
#[derive(Debug, Clone)]
pub struct MovementLog {
    /// Sequence of micro-operations
    ops: Vec<MicroOp>,
}

impl MovementLog {
    /// Create empty log
    pub fn new() -> Self {
        Self {
            ops: Vec::new(),
        }
    }
    
    /// Record one micro-operation
    pub fn record(&mut self, tape_id: usize, offset: i64, symbol: Symbol, move_dir: Move) {
        self.ops.push(MicroOp {
            tape_id,
            offset,
            symbol,
            move_dir,
        });
    }
    
    /// Replay log to reconstruct tape contents within a window
    ///
    /// This is used for interface checking (exact bounded-window replay).
    /// Returns a map of position -> symbol for the window region.
    pub fn replay_on_window(
        &self,
        window: WindowBounds,
        tape_id: usize,
        initial_head_pos: i64,
    ) -> std::collections::HashMap<i64, Symbol> {
        let mut contents = std::collections::HashMap::new();
        let mut current_pos = initial_head_pos;
        
        // Replay operations that affect this tape and window
        for op in &self.ops {
            if op.tape_id == tape_id {
                // Apply the operation
                if window.left <= current_pos && current_pos <= window.right {
                    contents.insert(current_pos, op.symbol);
                }
                
                // Move head
                current_pos = op.move_dir.apply(current_pos);
            }
        }
        
        contents
    }
    
    /// Replay log to get final head position for a tape
    pub fn replay_head_position(&self, tape_id: usize, initial_head_pos: i64) -> i64 {
        let mut current_pos = initial_head_pos;
        
        for op in &self.ops {
            if op.tape_id == tape_id {
                current_pos = op.move_dir.apply(current_pos);
            }
        }
        
        current_pos
    }
    
    /// Replay full tape to reconstruct complete tape contents
    ///
    /// Returns a sparse HashMap of position -> symbol for all non-blank cells
    /// This is used for boundary reconstruction when chaining blocks.
    pub fn replay_full_tape(
        &self,
        tape_id: usize,
        initial_head_pos: i64,
        blank_symbol: Symbol,
    ) -> std::collections::HashMap<i64, Symbol> {
        let mut contents = std::collections::HashMap::new();
        let mut current_pos = initial_head_pos;
        
        // Replay all operations that affect this tape
        for op in &self.ops {
            if op.tape_id == tape_id {
                // Record the symbol written at current position
                // Only store non-blank symbols (sparse representation)
                if op.symbol != blank_symbol {
                    contents.insert(current_pos, op.symbol);
                } else {
                    // If writing blank, remove from sparse representation
                    contents.remove(&current_pos);
                }
                
                // Move head according to operation
                current_pos = op.move_dir.apply(current_pos);
            }
        }
        
        contents
    }
    
    /// Get all operations
    pub fn operations(&self) -> &[MicroOp] {
        &self.ops
    }
    
    /// Space usage
    pub fn space_usage(&self) -> usize {
        // Each MicroOp is O(1) cells: tape_id (usize), offset (i64), symbol (char), move (enum)
        // Approximate: each op is roughly 16-20 bytes = 16-20 cells (assuming 1 byte per cell)
        // For simplicity, we count as 1 cell per op (the log is the dominant term)
        self.ops.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::machine::Move;
    
    #[test]
    fn test_block_summary_size() {
        // Create a movement log with b operations
        let b = 100;
        let mut log = MovementLog::new();
        for i in 0..b {
            log.record(0, i as i64, '1', Move::Right);
        }
        
        let summary = BlockSummary::new(
            1,
            0,
            1,
            vec![0],
            vec![b as i64],
            log,
            vec![WindowBounds { left: 0, right: b as i64 }],
        );
        
        // Verify space usage is O(b)
        let space = summary.space_usage();
        assert!(space >= b, "Space usage should be at least b={} operations", b);
        assert!(space <= b * 2, "Space usage should be O(b), got {}", space);
    }
    
    #[test]
    fn test_movement_log_replay() {
        let mut log = MovementLog::new();
        log.record(0, 0, '1', Move::Right);
        log.record(0, 1, '0', Move::Right);
        
        let window = WindowBounds { left: 0, right: 2 };
        let contents = log.replay_on_window(window, 0, 0);
        
        assert_eq!(contents.get(&0), Some(&'1'));
        assert_eq!(contents.get(&1), Some(&'0'));
        assert_eq!(log.replay_head_position(0, 0), 2);
    }
}