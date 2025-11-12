//! Full integration tests

use sqrt_space_sim::*;
mod test_helpers;
use test_helpers::*;

#[test]
fn test_end_to_end_simulation() {
    // Create a simple machine
    let machine = create_accept_machine();
    let input = vec!['1'];
    let time_bound = 10;

    // Create simulator
    let config = SimulationConfig::optimal_for_time(time_bound);
    let mut simulator = Simulator::new(machine.clone(), config.clone());

    // Run simulation - should succeed now with boundary reconstruction
    let result = simulator.run(&input).expect("Simulation should succeed");

    // Verify correctness
    assert!(result.accepted, "Machine should accept on '1'");
    assert_eq!(result.time_steps, time_bound);
    assert!(
        result.satisfies_bound(&config),
        "Space bound should be satisfied"
    );
    assert!(result.space_used > 0, "Space should be used");

    // Verify final state
    assert_eq!(result.final_config.state(), machine.accept_state());
}

#[test]
fn test_large_computation() {
    let machine = create_accept_machine();
    let input = vec!['1'];
    let time_bound = 1000; // Larger time bound to test multi-block

    let config = SimulationConfig {
        profile_space: true,
        ..SimulationConfig::optimal_for_time(time_bound)
    };

    let mut simulator = Simulator::new(machine, config.clone());
    let result = simulator
        .run(&input)
        .expect("Large simulation should succeed");

    // Verify correctness
    assert!(result.accepted, "Machine should accept");

    // Verify space bound
    assert!(
        result.space_used <= config.sqrt_t_bound(),
        "Space used {} exceeds bound {}",
        result.space_used,
        config.sqrt_t_bound()
    );

    // Verify space is indeed O(√t) - should be much less than t
    let sqrt_t = (time_bound as f64).sqrt() as usize;
    assert!(
        result.space_used <= sqrt_t * 20, // Allow constant factor
        "Space used {} should be O(√t) = O({}), not O(t)",
        result.space_used,
        sqrt_t
    );

    // Verify space is significantly less than t (efficiency claim)
    assert!(
        result.space_used < time_bound / 10,
        "Space {} should be much less than t={} for efficiency",
        result.space_used,
        time_bound
    );
}

#[test]
fn test_multiple_blocks() {
    let machine = create_right_move_machine();
    let input = vec!['1', '1', '1'];
    let time_bound = 20;

    let config = SimulationConfig::optimal_for_time(time_bound);
    let mut simulator = Simulator::new(machine, config.clone());

    // Run simulation - should succeed with boundary reconstruction
    let result = simulator
        .run(&input)
        .expect("Multi-block simulation should succeed");

    // Verify bounds
    assert!(
        result.satisfies_bound(&config),
        "Space bound should be satisfied"
    );
    assert!(result.space_used > 0, "Space should be used");

    // Verify we used multiple blocks
    assert!(
        config.num_blocks > 1,
        "Should have multiple blocks for t={}",
        time_bound
    );
}

#[test]
fn test_large_multi_block_simulation() {
    let machine = create_accept_machine();
    let input = vec!['1'];
    let time_bound = 10_000; // Large time bound to test many blocks

    let config = SimulationConfig {
        profile_space: true,
        ..SimulationConfig::optimal_for_time(time_bound)
    };

    let mut simulator = Simulator::new(machine, config.clone());
    let result = simulator
        .run(&input)
        .expect("Large multi-block simulation should succeed");

    // Verify correctness
    assert!(result.accepted, "Machine should accept");

    // Verify space bound
    assert!(
        result.space_used <= config.sqrt_t_bound(),
        "Space used {} exceeds bound {} for t={}",
        result.space_used,
        config.sqrt_t_bound(),
        time_bound
    );

    // Verify O(√t) scaling: space should scale with √t, not t
    let sqrt_t = (time_bound as f64).sqrt() as usize;
    assert!(
        result.space_used <= sqrt_t * 50, // Allow constant factor
        "Space used {} should be O(√t) = O({}), not O(t) = O({})",
        result.space_used,
        sqrt_t,
        time_bound
    );

    // Verify we have multiple blocks
    assert!(
        config.num_blocks > 1,
        "Should have multiple blocks for t={}",
        time_bound
    );

    // Verify space efficiency: space << t
    assert!(
        result.space_used < time_bound / 20,
        "Space {} should be much less than t={} to demonstrate efficiency",
        result.space_used,
        time_bound
    );
}

#[test]
fn test_space_scaling() {
    // Test that space scales as O(√t), not O(t)
    let machine = create_accept_machine();
    let input = vec!['1'];

    let mut space_measurements = Vec::new();

    // Test with increasing time bounds
    for t in [100, 400, 1600, 6400] {
        let config = SimulationConfig::optimal_for_time(t);
        let mut simulator = Simulator::new(machine.clone(), config.clone());
        let result = simulator.run(&input).expect("Simulation should succeed");

        space_measurements.push((t, result.space_used));

        // Verify bound for each
        assert!(
            result.space_used <= config.sqrt_t_bound(),
            "Space {} exceeds bound {} for t={}",
            result.space_used,
            config.sqrt_t_bound(),
            t
        );
    }

    // Verify scaling: when t increases by 4x, space should increase by ~2x (√4 = 2)
    // Not by 4x (which would be linear)
    for i in 1..space_measurements.len() {
        let (t1, s1) = space_measurements[i - 1];
        let (t2, s2) = space_measurements[i];

        let t_ratio = t2 as f64 / t1 as f64;
        let s_ratio = s2 as f64 / s1 as f64;
        let sqrt_t_ratio = t_ratio.sqrt();

        // Space should scale closer to √t than to t
        // Allow some variance but should be much closer to √t ratio
        assert!(
            s_ratio <= sqrt_t_ratio * 2.0, // Allow 2x factor
            "Space scaling {} should be closer to √t scaling {} than linear {}",
            s_ratio,
            sqrt_t_ratio,
            t_ratio
        );

        // Space should not scale linearly
        assert!(
            s_ratio < t_ratio / 2.0,
            "Space scaling {} should be sublinear (much less than t scaling {})",
            s_ratio,
            t_ratio
        );
    }
}
