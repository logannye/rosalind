//! Transition function δ: Q × Γ^(τ+1) → Q × Γ^τ × {L,R,S}^(τ+1)

use super::{State, Symbol};

/// Single transition rule
#[derive(Debug, Clone, PartialEq)]
pub struct Transition {
    /// Next state
    pub next_state: State,
    
    /// Symbols to write (one per work tape; input tape is read-only)
    pub writes: Vec<Symbol>,
    
    /// Head movements for all tapes (including input)
    pub moves: Vec<Move>,
}

/// Head movement direction
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Move {
    /// Move left (decrement position)
    Left,
    
    /// Move right (increment position)
    Right,
    
    /// Stay (no movement)
    Stay,
}

impl Move {
    /// Apply move to position
    pub fn apply(&self, position: i64) -> i64 {
        match self {
            Move::Left => position - 1,
            Move::Right => position + 1,
            Move::Stay => position,
        }
    }
    
    /// Encode as integer for movement log
    pub fn to_i8(&self) -> i8 {
        match self {
            Move::Left => -1,
            Move::Right => 1,
            Move::Stay => 0,
        }
    }
    
    /// Decode from integer
    pub fn from_i8(val: i8) -> Self {
        match val {
            -1 => Move::Left,
            1 => Move::Right,
            _ => Move::Stay,
        }
    }
}