//! Block-respecting simulation framework
//!
//! Partitions time-t computation into T = ⌈t/b⌉ blocks of size b:
//! - Each block touches ≤ b cells per tape
//! - Block summary σ_k: O(b) cells encoding
//! - Movement log: deterministic micro-ops (write, move)
//! - Interface checking: exact bounded-window replay

mod interface;
mod summary;

pub use interface::InterfaceChecker;
pub use summary::{BlockSummary, MovementLog};

use crate::{
    machine::{Configuration, TuringMachine},
    SimulationError,
};

/// Block ID (index in [1, T])
pub type BlockId = usize;

/// Window bounds on a tape (interval of cells)
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WindowBounds {
    /// Left boundary of window
    pub left: i64,
    /// Right boundary of window
    pub right: i64,
}

impl WindowBounds {
    /// Create window for visited cells
    pub fn from_positions(positions: &[i64]) -> Self {
        if positions.is_empty() {
            return Self { left: 0, right: 0 };
        }

        let min_pos = *positions.iter().min().unwrap();
        let max_pos = *positions.iter().max().unwrap();
        Self {
            left: min_pos,
            right: max_pos,
        }
    }

    /// Window length (number of cells)
    pub fn length(&self) -> usize {
        if self.right < self.left {
            0
        } else {
            (self.right - self.left + 1) as usize
        }
    }

    /// Check if windows overlap
    pub fn overlaps(&self, other: &WindowBounds) -> bool {
        self.left <= other.right && other.left <= self.right
    }
}

/// Simulate one block and produce summary
pub fn simulate_block(
    machine: &TuringMachine,
    config_in: &Configuration,
    block_id: BlockId,
    block_size: usize,
) -> Result<BlockSummary, SimulationError> {
    // Clone configuration to avoid mutating input
    let mut config = config_in.clone();

    // Record entry state and head positions
    let entry_state = config.state();
    let entry_heads = config.head_positions();

    // Track all head positions visited during block
    let mut all_head_positions: Vec<Vec<i64>> = vec![Vec::new(); config.num_tapes()];
    for (i, &head_pos) in entry_heads.iter().enumerate() {
        if i < all_head_positions.len() {
            all_head_positions[i].push(head_pos);
        }
    }

    // Create movement log
    let mut movement_log = MovementLog::new();

    // 1. Execute b steps
    for _step in 0..block_size {
        // Check if machine has halted
        if machine.is_halted(config.state()) {
            break;
        }

        // Record current head positions and state before step
        let current_heads = config.head_positions();
        let current_state = config.state();
        let symbols_read = config.read_symbols();

        // Look up transition to know what will be written
        let transition = machine.transition(current_state, &symbols_read);

        // Execute one step
        machine.step(&mut config)?;

        // Record head positions after step
        let new_heads = config.head_positions();

        // Record micro-operations
        // For work tapes: record writes and moves
        // For input tape: only record moves (read-only)
        for tape_id in 0..config.num_tapes() {
            let old_pos = current_heads[tape_id];
            let new_pos = new_heads[tape_id];

            // Determine move direction
            let move_dir = if new_pos < old_pos {
                crate::machine::Move::Left
            } else if new_pos > old_pos {
                crate::machine::Move::Right
            } else {
                crate::machine::Move::Stay
            };

            if tape_id == 0 {
                // Input tape: only track position for window bounds (read-only)
            } else {
                // Work tape: record what was written
                let work_tape_idx = tape_id - 1; // Work tapes are 1-indexed in config
                let symbol_written = if let Some(trans) = transition {
                    if work_tape_idx < trans.writes.len() {
                        trans.writes[work_tape_idx]
                    } else {
                        // No write for this tape, read current
                        config.tape(tape_id).read(old_pos)
                    }
                } else {
                    // No transition found, read current
                    config.tape(tape_id).read(old_pos)
                };

                // Record operation (offset is relative to entry position)
                let offset = old_pos - entry_heads[tape_id];
                movement_log.record(tape_id, offset, symbol_written, move_dir);
            }

            // Track head position
            all_head_positions[tape_id].push(new_pos);
        }
    }

    // Record exit state and head positions
    let exit_state = config.state();
    let exit_heads = config.head_positions();

    // 3. Compute window bounds for each tape
    let windows: Vec<WindowBounds> = all_head_positions
        .iter()
        .map(|positions| WindowBounds::from_positions(positions))
        .collect();

    // 4. Return BlockSummary
    Ok(BlockSummary::new(
        block_id,
        entry_state,
        exit_state,
        entry_heads,
        exit_heads,
        movement_log,
        windows,
    ))
}
