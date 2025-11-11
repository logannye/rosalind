//! Configuration (instantaneous description) of Turing machine
//!
//! Represents complete state at a point in time:
//! - Control state q ∈ Q
//! - Head positions for all tapes
//! - Tape contents (sparse representation)

use super::{Symbol, State};
use std::collections::HashMap;

/// Complete instantaneous description of machine
#[derive(Debug, Clone, PartialEq)]
pub struct Configuration {
    /// Current control state
    state: State,
    
    /// Head position on read-only input tape (tape 0)
    input_head: i64,
    
    /// Head positions on work tapes (tapes 1..num_tapes)
    work_heads: Vec<i64>,
    
    /// Sparse tape contents (tape 0 = input, 1.. = work tapes)
    tapes: Vec<Tape>,
    
    /// Input tape contents (read-only, stored separately)
    input: Vec<Symbol>,
}

impl Configuration {
    /// Create initial configuration for input
    pub fn initial(input: &[Symbol], num_tapes: usize) -> Self {
        let blank = '_';
        
        // Create input tape (read-only, but we'll store it for reading)
        let mut input_tape = Tape::blank(blank);
        for (i, &symbol) in input.iter().enumerate() {
            input_tape.write(i as i64, symbol);
        }
        
        // Create work tapes (all blank initially)
        let mut tapes = vec![input_tape];
        for _ in 1..num_tapes {
            tapes.push(Tape::blank(blank));
        }
        
        Self {
            state: 0, // q_0 = initial state
            input_head: 0,
            work_heads: vec![0; num_tapes - 1], // num_tapes - 1 work tapes
            tapes,
            input: input.to_vec(),
        }
    }
    
    /// Get current state
    pub fn state(&self) -> State {
        self.state
    }
    
    /// Set state (for transitions)
    pub fn set_state(&mut self, state: State) {
        self.state = state;
    }
    
    /// Set head positions (for boundary reconstruction)
    pub fn set_head_positions(&mut self, input_head: i64, work_heads: Vec<i64>) {
        self.input_head = input_head;
        self.work_heads = work_heads;
    }
    
    /// Set tape contents (for boundary reconstruction)
    pub fn set_tapes(&mut self, tapes: Vec<Tape>) {
        self.tapes = tapes;
    }
    
    /// Read symbols under all tape heads
    pub fn read_symbols(&self) -> Vec<Symbol> {
        let mut symbols = Vec::with_capacity(self.tapes.len());
        
        // Read from input tape (tape 0)
        let input_symbol = if (self.input_head as usize) < self.input.len() {
            self.input[self.input_head as usize]
        } else {
            '_' // Blank beyond input
        };
        symbols.push(input_symbol);
        
        // Read from work tapes (tapes 1..)
        for (i, &head_pos) in self.work_heads.iter().enumerate() {
            let tape_idx = i + 1; // Work tape indices start at 1
            symbols.push(self.tapes[tape_idx].read(head_pos));
        }
        
        symbols
    }
    
    /// Apply transition to this configuration
    pub fn apply_transition(&mut self, transition: &super::Transition) {
        // Update state
        self.state = transition.next_state;
        
        // Write symbols to work tapes (input tape is read-only)
        for (i, &symbol) in transition.writes.iter().enumerate() {
            let tape_idx = i + 1; // Work tapes start at index 1
            let head_pos = self.work_heads[i];
            self.tapes[tape_idx].write(head_pos, symbol);
        }
        
        // Move heads
        // Input tape head
        self.input_head = transition.moves[0].apply(self.input_head);
        
        // Work tape heads
        for (i, &move_dir) in transition.moves[1..].iter().enumerate() {
            self.work_heads[i] = move_dir.apply(self.work_heads[i]);
        }
    }
    
    /// Get head positions (input + work tapes)
    pub fn head_positions(&self) -> Vec<i64> {
        let mut positions = vec![self.input_head];
        positions.extend_from_slice(&self.work_heads);
        positions
    }
    
    /// Get tape by index (0 = input, 1.. = work tapes)
    pub fn tape(&self, idx: usize) -> &Tape {
        &self.tapes[idx]
    }
    
    /// Get mutable tape by index
    pub fn tape_mut(&mut self, idx: usize) -> &mut Tape {
        &mut self.tapes[idx]
    }
    
    /// Number of tapes (including input)
    pub fn num_tapes(&self) -> usize {
        self.tapes.len()
    }
}

/// Sparse representation of single tape
#[derive(Debug, Clone, PartialEq)]
pub struct Tape {
    /// Sparse storage: only non-blank cells
    cells: HashMap<i64, Symbol>,
    
    /// Blank symbol (default for unvisited cells)
    blank: Symbol,
}

impl Tape {
    /// Create blank tape
    pub fn blank(blank_symbol: Symbol) -> Self {
        Self {
            cells: HashMap::new(),
            blank: blank_symbol,
        }
    }
    
    /// Read symbol at position
    pub fn read(&self, position: i64) -> Symbol {
        self.cells.get(&position).copied().unwrap_or(self.blank)
    }
    
    /// Write symbol at position
    pub fn write(&mut self, position: i64, symbol: Symbol) {
        if symbol == self.blank {
            // Remove from map if writing blank (space optimization)
            self.cells.remove(&position);
        } else {
            self.cells.insert(position, symbol);
        }
    }
    
    /// Get bounds of written region
    pub fn bounds(&self) -> (i64, i64) {
        if self.cells.is_empty() {
            return (0, 0);
        }
        
        let min_pos = *self.cells.keys().min().unwrap();
        let max_pos = *self.cells.keys().max().unwrap();
        (min_pos, max_pos)
    }
    
    /// Space usage: number of non-blank cells
    pub fn space_usage(&self) -> usize {
        self.cells.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_sparse_tape() {
        // Verify that blank cells don't consume space
        let mut tape = Tape::blank('_');
        
        // Write to a far position
        tape.write(1000000, '1');
        
        // Only one non-blank cell should be stored
        assert_eq!(tape.space_usage(), 1);
        
        // Verify read works
        assert_eq!(tape.read(1000000), '1');
        assert_eq!(tape.read(500000), '_'); // Blank cell
        
        // Writing blank should remove from map
        tape.write(1000000, '_');
        assert_eq!(tape.space_usage(), 0);
    }
    
    #[test]
    fn test_configuration_initial() {
        let input = vec!['1', '0', '1'];
        let config = Configuration::initial(&input, 2); // 1 input + 1 work tape
        
        assert_eq!(config.state(), 0);
        assert_eq!(config.head_positions().len(), 2); // input + 1 work tape
        assert_eq!(config.num_tapes(), 2);
        
        // Verify input reading
        let symbols = config.read_symbols();
        assert_eq!(symbols[0], '1'); // First input symbol
        assert_eq!(symbols[1], '_'); // Work tape is blank
    }
}