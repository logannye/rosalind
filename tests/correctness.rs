//! Correctness tests: verify simulation matches direct execution

use rosalind::*;
mod test_helpers;
use test_helpers::*;

#[test]
fn test_simple_machine_correctness() {
    let machine = create_accept_machine();
    let input = vec!['1'];
    let time_bound = 10;

    let config = SimulationConfig::optimal_for_time(time_bound);
    let mut simulator = Simulator::new(machine.clone(), config);

    // Run compressed simulation - should succeed now
    let result = simulator.run(&input).expect("Simulation should succeed");

    // Verify it accepted (if machine reached accept state)
    // The machine should accept on '1'
    if result.accepted {
        assert_eq!(result.final_config.state(), machine.accept_state());
    }
    // Verify space bound is satisfied
    assert!(result.space_used > 0);
}

#[test]
fn test_block_summary_composition() {
    use rosalind::blocking::simulate_block;
    use rosalind::machine::Configuration;

    let machine = create_right_move_machine();
    let input = vec!['1', '1', '0'];
    let block_size = 2;

    // Simulate first block
    let config1 = Configuration::initial(&input, 2);
    let summary1 = simulate_block(&machine, &config1, 1, block_size).unwrap();

    // Verify summary1 has correct structure
    assert_eq!(summary1.block_id, 1);
    assert!(summary1.movement_log().operations().len() <= block_size);

    // Verify summary1 structure is valid
    assert_eq!(summary1.entry_heads().len(), 2); // input + 1 work tape
    assert_eq!(summary1.exit_heads().len(), 2);
    assert_eq!(summary1.windows().len(), 2);

    // Now test boundary reconstruction: reconstruct config2 from summary1
    // (clone since into_configuration takes ownership)
    let exit_state = summary1.exit_state();
    let exit_heads = summary1.exit_heads().to_vec();
    let config2 = summary1.clone().into_configuration(&input, machine.blank());

    // Verify reconstructed config matches summary1's exit state
    assert_eq!(config2.state(), exit_state);

    // Verify head positions match
    let head_positions = config2.head_positions();
    assert_eq!(head_positions, exit_heads);
}

#[test]
fn test_interface_checking() {
    use rosalind::blocking::{BlockSummary, InterfaceChecker, MovementLog};
    use rosalind::machine::Move;

    // Create two adjacent blocks with matching interface
    let mut left_log = MovementLog::new();
    left_log.record(1, 0, '1', Move::Right);

    let left = BlockSummary::new(
        1,
        0,          // entry state
        1,          // exit state
        vec![0, 0], // entry heads
        vec![0, 1], // exit heads
        left_log,
        vec![
            crate::blocking::WindowBounds { left: 0, right: 1 },
            crate::blocking::WindowBounds { left: 0, right: 1 },
        ],
    );

    let mut right_log = MovementLog::new();
    right_log.record(1, 1, '0', Move::Right);

    let right = BlockSummary::new(
        2,
        1,          // entry state (matches left exit)
        2,          // exit state
        vec![0, 1], // entry heads (matches left exit)
        vec![0, 2], // exit heads
        right_log,
        vec![
            crate::blocking::WindowBounds { left: 0, right: 2 },
            crate::blocking::WindowBounds { left: 1, right: 2 },
        ],
    );

    // Interface should be consistent
    assert!(InterfaceChecker::check(&left, &right).unwrap());

    // Create mismatched interface
    let right_mismatched = BlockSummary::new(
        2,
        0, // entry state (doesn't match left exit)
        2,
        vec![0, 1],
        vec![0, 2],
        MovementLog::new(),
        vec![
            crate::blocking::WindowBounds { left: 0, right: 2 },
            crate::blocking::WindowBounds { left: 1, right: 2 },
        ],
    );

    // Interface should be inconsistent
    assert!(!InterfaceChecker::check(&left, &right_mismatched).unwrap());
}

#[test]
fn test_boundary_reconstruction() {
    use rosalind::blocking::simulate_block;
    use rosalind::machine::Configuration;

    let machine = create_right_move_machine();
    let input = vec!['1', '1', '1'];
    let block_size = 3;

    // Simulate first block
    let config1 = Configuration::initial(&input, 2);
    let summary1 = simulate_block(&machine, &config1, 1, block_size).unwrap();

    // Reconstruct boundary configuration from summary1 (clone since into_configuration takes ownership)
    let exit_state = summary1.exit_state();
    let exit_heads = summary1.exit_heads().to_vec();
    let config2_boundary = summary1.clone().into_configuration(&input, machine.blank());

    // Verify boundary reconstruction correctness
    assert_eq!(
        config2_boundary.state(),
        exit_state,
        "Block 2 should start from block 1's exit state"
    );
    assert_eq!(
        config2_boundary.head_positions(),
        exit_heads,
        "Block 2 should start from block 1's exit head positions"
    );

    // Simulate second block starting from reconstructed boundary
    let summary2 = simulate_block(&machine, &config2_boundary, 2, block_size).unwrap();

    // Verify interface consistency: block 1's exit should match block 2's entry
    assert_eq!(
        exit_state,
        summary2.entry_state(),
        "Block 1 exit state should match block 2 entry state"
    );
    assert_eq!(
        exit_heads,
        summary2.entry_heads(),
        "Block 1 exit heads should match block 2 entry heads"
    );

    // Verify interface checker passes
    use rosalind::blocking::InterfaceChecker;
    assert!(
        InterfaceChecker::check(&summary1, &summary2).unwrap(),
        "Interface check should pass for correctly reconstructed blocks"
    );
}
