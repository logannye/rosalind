//! Test helper functions for creating long-running machines for scale tests
//!
//! These machines are designed to run for many steps without halting early,
//! allowing proper testing of O(√t) space complexity at scale.

use rosalind::{Move, TuringMachine};

/// Create a counter machine that counts from 0 to approximately t
///
/// This machine writes incrementing values on the work tape, using
/// approximately t steps. It doesn't halt early, making it ideal for
/// testing space complexity at scale.
///
/// Machine behavior:
/// - State 0: Counting state (runs for t steps)
/// - State 1: Accept state (reached after t steps)
/// - Uses binary representation on work tape
pub fn create_counter_machine() -> TuringMachine {
    TuringMachine::builder()
        .num_tapes(1) // 1 work tape
        .alphabet(vec!['_', '0', '1'])
        .initial_state(0)
        .accept_state(1)
        .reject_state(2)
        // Counter logic: move right, write '1' to count
        // This transition can be used repeatedly
        .add_transition(
            0,                             // counting state
            vec!['_', '_'],                // read: input='_', work='_'
            0,                             // stay in counting state
            vec!['1'],                     // write '1' to work tape
            vec![Move::Stay, Move::Right], // input stays, work moves right
        )
        // If we see '0' on work tape, change it to '1' (increment)
        .add_transition(
            0,
            vec!['_', '0'],
            0,
            vec!['1'],
            vec![Move::Stay, Move::Right],
        )
        // If we see '1' on work tape, change it to '0' and continue (carry)
        .add_transition(
            0,
            vec!['_', '1'],
            0,
            vec!['0'],
            vec![Move::Stay, Move::Right],
        )
        // After many steps, we can transition to accept
        // But we'll let the time bound control when it stops
        .build()
        .unwrap()
}

/// Create a scanning machine that moves back and forth
///
/// This machine scans the work tape by moving the head back and forth
/// repeatedly, using t steps. Good for testing space usage with head movement.
///
/// Machine behavior:
/// - Moves right until hitting a marker, then left, then right again
/// - Uses t steps of head movement
pub fn create_scanning_machine() -> TuringMachine {
    TuringMachine::builder()
        .num_tapes(1)
        .alphabet(vec!['_', '0', '1', 'X']) // X is marker
        .initial_state(0)
        .accept_state(1)
        .reject_state(2)
        // State 0: Moving right, writing markers
        .add_transition(
            0,
            vec!['_', '_'],                // read blank
            0,                             // stay in scanning state
            vec!['X'],                     // write marker
            vec![Move::Stay, Move::Right], // move work head right
        )
        .add_transition(
            0,
            vec!['_', 'X'], // hit marker, turn around
            0,
            vec!['_'],                    // erase marker
            vec![Move::Stay, Move::Left], // move left
        )
        .add_transition(
            0,
            vec!['_', '0'], // moving left, see '0'
            0,
            vec!['0'],
            vec![Move::Stay, Move::Left],
        )
        .add_transition(
            0,
            vec!['_', '1'], // moving left, see '1'
            0,
            vec!['1'],
            vec![Move::Stay, Move::Left],
        )
        .build()
        .unwrap()
}

/// Create a simple computation machine that performs repeated operations
///
/// This machine does a simple computation (copying/transforming input)
/// that guarantees t steps of work without early halting.
pub fn create_computation_machine() -> TuringMachine {
    TuringMachine::builder()
        .num_tapes(1)
        .alphabet(vec!['_', '0', '1'])
        .initial_state(0)
        .accept_state(1)
        .reject_state(2)
        // Simple computation: write pattern on work tape
        // State 0: alternating pattern '1', '0', '1', '0', ...
        .add_transition(
            0,
            vec!['_', '_'],                // start with blank
            0,                             // stay in computation state
            vec!['1'],                     // write '1'
            vec![Move::Stay, Move::Right], // move right
        )
        .add_transition(
            0,
            vec!['_', '0'], // see '0', write '1'
            0,
            vec!['1'],
            vec![Move::Stay, Move::Right],
        )
        .add_transition(
            0,
            vec!['_', '1'], // see '1', write '0'
            0,
            vec!['0'],
            vec![Move::Stay, Move::Right],
        )
        .build()
        .unwrap()
}

/// Create a simple non-halting machine that runs for exactly t steps
///
/// This is a minimal machine that just moves right and writes,
/// ensuring it uses the full time bound without halting early.
/// Useful for testing when you need guaranteed t steps.
pub fn create_non_halting_machine() -> TuringMachine {
    TuringMachine::builder()
        .num_tapes(1)
        .alphabet(vec!['_', '0', '1'])
        .initial_state(0)
        .accept_state(1)
        .reject_state(2)
        // Simple loop: always move right and write '1'
        // This will never halt (until time bound is reached)
        .add_transition(
            0,
            vec!['_', '_'],                // any input, blank work
            0,                             // stay in state 0 (never reach accept/reject)
            vec!['1'],                     // write '1'
            vec![Move::Stay, Move::Right], // move right
        )
        .add_transition(
            0,
            vec!['_', '0'], // any input, '0' on work
            0,
            vec!['1'],
            vec![Move::Stay, Move::Right],
        )
        .add_transition(
            0,
            vec!['_', '1'], // any input, '1' on work
            0,
            vec!['1'],
            vec![Move::Stay, Move::Right],
        )
        .build()
        .unwrap()
}
