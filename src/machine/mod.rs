//! Turing machine representation and execution
//!
//! Provides abstractions for deterministic multitape Turing machines:
//! - Constant number of tapes τ
//! - Finite alphabet Γ (constant size)
//! - Finite state set Q
//! - Deterministic transition function δ

mod config;
mod transition;

pub use config::{Configuration, Tape};
pub use transition::{Move, Transition};

use crate::SimulationError;
use std::collections::HashMap;

/// Tape symbol (element of alphabet Γ)
pub type Symbol = char;

/// Machine state (element of Q)
pub type State = u32;

/// Tape identifier (0 = input tape, 1..τ = work tapes)
pub type TapeId = usize;

/// Deterministic multitape Turing machine
///
/// Model: τ work tapes + 1 read-only input tape over alphabet Γ
#[derive(Debug, Clone)]
pub struct TuringMachine {
    /// Number of work tapes (constant τ)
    num_tapes: usize,

    /// Alphabet Γ (constant size)
    #[allow(dead_code)]
    alphabet: Vec<Symbol>,

    /// Blank symbol
    blank: Symbol,

    /// Set of states Q (used for validation)
    #[allow(dead_code)]
    states: Vec<State>,

    /// Initial state q_0
    initial_state: State,

    /// Accept state
    accept_state: State,

    /// Reject state
    reject_state: State,

    /// Transition function δ: Q × Γ^(τ+1) → Q × Γ^τ × {L,R,S}^(τ+1)
    transitions: HashMap<(State, Vec<Symbol>), Transition>,
}

impl TuringMachine {
    /// Create fluent builder
    pub fn builder() -> TuringMachineBuilder {
        TuringMachineBuilder::new()
    }

    /// Get number of work tapes
    pub fn num_tapes(&self) -> usize {
        self.num_tapes
    }

    /// Get blank symbol
    pub fn blank(&self) -> Symbol {
        self.blank
    }

    /// Get accept state
    pub fn accept_state(&self) -> State {
        self.accept_state
    }

    /// Get reject state
    pub fn reject_state(&self) -> State {
        self.reject_state
    }

    /// Get initial state
    pub fn initial_state(&self) -> State {
        self.initial_state
    }

    /// Look up transition for (state, symbols)
    pub fn transition(&self, state: State, symbols: &[Symbol]) -> Option<&Transition> {
        self.transitions.get(&(state, symbols.to_vec()))
    }

    /// Execute one step from configuration
    pub fn step(&self, config: &mut Configuration) -> Result<(), SimulationError> {
        // 1. Read symbols under all tape heads
        let symbols = config.read_symbols();

        // 2. Look up transition
        let transition = self.transition(config.state(), &symbols).ok_or_else(|| {
            SimulationError::InvalidMachine(format!(
                "No transition for state {} and symbols {:?}",
                config.state(),
                symbols
            ))
        })?;

        // 3. Apply transition to configuration
        config.apply_transition(transition);

        Ok(())
    }

    /// Check if in halting state
    pub fn is_halted(&self, state: State) -> bool {
        state == self.accept_state || state == self.reject_state
    }
}

/// Builder for Turing machines (fluent API)
#[derive(Debug)]
pub struct TuringMachineBuilder {
    num_tapes: Option<usize>,
    alphabet: Option<Vec<Symbol>>,
    blank: Symbol,
    #[allow(dead_code)]
    states: Vec<State>,
    initial_state: Option<State>,
    accept_state: Option<State>,
    reject_state: Option<State>,
    transitions: Vec<(State, Vec<Symbol>, Transition)>,
}

impl TuringMachineBuilder {
    /// Create new builder
    pub fn new() -> Self {
        Self {
            num_tapes: None,
            alphabet: None,
            blank: '_',
            states: Vec::new(),
            initial_state: None,
            accept_state: None,
            reject_state: None,
            transitions: Vec::new(),
        }
    }

    /// Set number of work tapes
    pub fn num_tapes(mut self, n: usize) -> Self {
        self.num_tapes = Some(n);
        self
    }

    /// Set alphabet
    pub fn alphabet(mut self, alphabet: Vec<Symbol>) -> Self {
        self.alphabet = Some(alphabet);
        self
    }

    /// Set blank symbol
    pub fn blank(mut self, blank: Symbol) -> Self {
        self.blank = blank;
        self
    }

    /// Set initial state
    pub fn initial_state(mut self, state: State) -> Self {
        self.initial_state = Some(state);
        self
    }

    /// Set accept state
    pub fn accept_state(mut self, state: State) -> Self {
        self.accept_state = Some(state);
        self
    }

    /// Set reject state
    pub fn reject_state(mut self, state: State) -> Self {
        self.reject_state = Some(state);
        self
    }

    /// Add a transition rule
    pub fn add_transition(
        mut self,
        from_state: State,
        read_symbols: Vec<Symbol>,
        to_state: State,
        write_symbols: Vec<Symbol>,
        moves: Vec<Move>,
    ) -> Self {
        let transition = Transition {
            next_state: to_state,
            writes: write_symbols,
            moves,
        };
        self.transitions
            .push((from_state, read_symbols, transition));
        self
    }

    /// Build the Turing machine
    pub fn build(self) -> Result<TuringMachine, SimulationError> {
        let num_tapes = self.num_tapes.ok_or_else(|| {
            SimulationError::InvalidMachine("Number of tapes not set".to_string())
        })?;

        let alphabet = self
            .alphabet
            .ok_or_else(|| SimulationError::InvalidMachine("Alphabet not set".to_string()))?;

        let initial_state = self.initial_state.unwrap_or(0);
        let accept_state = self
            .accept_state
            .ok_or_else(|| SimulationError::InvalidMachine("Accept state not set".to_string()))?;
        let reject_state = self
            .reject_state
            .ok_or_else(|| SimulationError::InvalidMachine("Reject state not set".to_string()))?;

        // Collect all states
        let mut states = std::collections::HashSet::new();
        states.insert(initial_state);
        states.insert(accept_state);
        states.insert(reject_state);
        for (from, _, ref trans) in &self.transitions {
            states.insert(*from);
            states.insert(trans.next_state);
        }
        let states: Vec<State> = states.into_iter().collect();

        // Build transition map
        let mut transitions = HashMap::new();
        for (from_state, read_symbols, transition) in self.transitions {
            // Validate: read_symbols length should match num_tapes + 1 (input + work tapes)
            if read_symbols.len() != num_tapes + 1 {
                return Err(SimulationError::InvalidMachine(format!(
                    "Transition read symbols length {} doesn't match num_tapes+1 {}",
                    read_symbols.len(),
                    num_tapes + 1
                )));
            }

            // Validate: write_symbols length should match num_tapes (work tapes only)
            if transition.writes.len() != num_tapes {
                return Err(SimulationError::InvalidMachine(format!(
                    "Transition write symbols length {} doesn't match num_tapes {}",
                    transition.writes.len(),
                    num_tapes
                )));
            }

            // Validate: moves length should match num_tapes + 1 (input + work tapes)
            if transition.moves.len() != num_tapes + 1 {
                return Err(SimulationError::InvalidMachine(format!(
                    "Transition moves length {} doesn't match num_tapes+1 {}",
                    transition.moves.len(),
                    num_tapes + 1
                )));
            }

            transitions.insert((from_state, read_symbols), transition);
        }

        Ok(TuringMachine {
            num_tapes,
            alphabet,
            blank: self.blank,
            states,
            initial_state,
            accept_state,
            reject_state,
            transitions,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_single_step() {
        // Create a simple machine: reads '1' and writes '0', moves right
        let machine = TuringMachine::builder()
            .num_tapes(1) // 1 work tape
            .alphabet(vec!['_', '0', '1'])
            .initial_state(0)
            .accept_state(1)
            .reject_state(2)
            .add_transition(
                0,                              // from state
                vec!['1', '_'],                 // read: input='1', work='_'
                1,                              // to state
                vec!['0'],                      // write '0' to work tape
                vec![Move::Right, Move::Right], // move both heads right
            )
            .build()
            .unwrap();

        // Create initial configuration
        let input = vec!['1'];
        let mut config = Configuration::initial(&input, 2); // input + 1 work tape

        // Execute one step
        machine.step(&mut config).unwrap();

        // Verify state changed
        assert_eq!(config.state(), 1);

        // Verify work tape was written
        assert_eq!(config.tape(1).read(0), '0');

        // Verify heads moved right
        assert_eq!(config.head_positions()[0], 1); // input head
        assert_eq!(config.head_positions()[1], 1); // work head
    }
}
