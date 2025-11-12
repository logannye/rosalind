//! Test helper functions for creating test machines

#![allow(dead_code)]
#![allow(dead_code)]
use rosalind::{Move, TuringMachine};

/// Create a simple accepting machine
/// Reads '1' and accepts, otherwise rejects
pub fn create_accept_machine() -> TuringMachine {
    TuringMachine::builder()
        .num_tapes(1) // 1 work tape
        .alphabet(vec!['_', '0', '1'])
        .initial_state(0)
        .accept_state(1)
        .reject_state(2)
        .add_transition(
            0,
            vec!['1', '_'],               // read: input='1', work='_'
            1,                            // accept
            vec!['_'],                    // no write
            vec![Move::Stay, Move::Stay], // no move
        )
        .add_transition(
            0,
            vec!['0', '_'], // read: input='0', work='_'
            2,              // reject
            vec!['_'],
            vec![Move::Stay, Move::Stay],
        )
        .build()
        .unwrap()
}

/// Create a simple machine that moves right
pub fn create_right_move_machine() -> TuringMachine {
    TuringMachine::builder()
        .num_tapes(1)
        .alphabet(vec!['_', '0', '1'])
        .initial_state(0)
        .accept_state(1)
        .reject_state(2)
        .add_transition(
            0,
            vec!['1', '_'],
            0,                              // stay in state 0
            vec!['1'],                      // write '1' to work tape
            vec![Move::Right, Move::Right], // move both heads right
        )
        .add_transition(
            0,
            vec!['_', '_'], // end of input
            1,              // accept
            vec!['_'],
            vec![Move::Stay, Move::Stay],
        )
        .build()
        .unwrap()
}
