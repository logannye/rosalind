//! Space bound verification tests

use sqrt_space_sim::*;
mod test_helpers;
use test_helpers::*;

#[test]
fn test_space_bound_satisfaction() {
    let machine = create_accept_machine();
    let input = vec!['1'];
    
    // Test with various time bounds - should all succeed now
    for t in [100, 1_000, 10_000] {
        let config = SimulationConfig {
            profile_space: true,
            ..SimulationConfig::optimal_for_time(t)
        };
        
        let mut simulator = Simulator::new(machine.clone(), config.clone());
        let result = simulator.run(&input).expect(&format!("Simulation should succeed for t={}", t));
        
        // Verify space bound: used ≤ O(√t)
        let sqrt_bound = config.sqrt_t_bound();
        assert!(
            result.space_used <= sqrt_bound,
            "t={}: space_used {} > bound {}",
            t, result.space_used, sqrt_bound
        );
        
        // Verify space is sublinear in t (should be much less than t)
        // For t=100, space should be << 100 (e.g., < 50)
        // For t=1000, space should be << 1000 (e.g., < 200)
        let max_space = t / 2; // Allow up to 50% of t (still sublinear)
        assert!(
            result.space_used < max_space,
            "t={}: space {} should be much less than t (max {}) for efficiency",
            t,
            result.space_used,
            max_space
        );
    }
}

#[test]
fn test_component_space_breakdown() {
    let machine = create_accept_machine();
    let input = vec!['1'];
    let time_bound = 1000; // Larger time bound to test components
    
    let config = SimulationConfig {
        profile_space: true,
        ..SimulationConfig::optimal_for_time(time_bound)
    };
    
    let mut simulator = Simulator::new(machine, config.clone());
    let result = simulator.run(&input).expect("Simulation should succeed");
    
    if let Some(profile) = result.space_profile {
        // Leaf buffer: O(b) = O(√t)
        assert!(
            profile.leaf_buffer_max <= config.block_size * 2,
            "Leaf buffer {} > 2*b = {}",
            profile.leaf_buffer_max,
            config.block_size * 2
        );
        
        // Stack depth: O(log T) where T = num_blocks = O(√t)
        let log_t = (config.num_blocks as f64).log2().ceil() as usize;
        assert!(
            profile.stack_depth_max <= log_t * 3,
            "Stack depth {} > 3*log(T) = {}",
            profile.stack_depth_max,
            log_t * 3
        );
        
        // Ledger: O(T) = O(√t) since T = t/b = t/√t = √t
        assert!(
            profile.ledger_size <= config.num_blocks * 2,
            "Ledger size {} > 2*T = {}",
            profile.ledger_size,
            config.num_blocks * 2
        );
        
        // Total space should be O(√t)
        let total_space = profile.leaf_buffer_max + 
                         (profile.stack_depth_max * std::mem::size_of::<usize>()) +
                         profile.ledger_size;
        
        let sqrt_t = (time_bound as f64).sqrt() as usize;
        assert!(
            total_space <= sqrt_t * 30, // Allow constant factor
            "Total component space {} should be O(√t) = O({})",
            total_space,
            sqrt_t
        );
    }
}

#[test]
fn test_optimal_block_size() {
    let _machine = create_accept_machine();
    let _input = vec!['1'];
    let time_bound = 10_000;
    
    // Test optimal block size
    let optimal_config = SimulationConfig::optimal_for_time(time_bound);
    let expected_block_size = (time_bound as f64).sqrt().ceil() as usize;
    assert_eq!(optimal_config.block_size, expected_block_size,
               "Block size should be ⌈√t⌉");
    assert_eq!(optimal_config.num_blocks, 
               (time_bound + optimal_config.block_size - 1) / optimal_config.block_size,
               "Number of blocks should be ⌈t/b⌉");
    
    // Verify space bound formula
    let bound = optimal_config.space_bound();
    let sqrt_t = (time_bound as f64).sqrt() as usize;
    assert!(
        bound <= sqrt_t * 5,
        "Space bound {} should be O(√t) = O({})",
        bound,
        sqrt_t
    );
}

#[test]
fn test_space_efficiency_claim() {
    // Verify the core efficiency claim: O(√t) space vs O(t) naive
    let machine = create_accept_machine();
    let input = vec!['1'];
    
    // Test with large time bound
    let time_bound = 10_000;
    let config = SimulationConfig::optimal_for_time(time_bound);
    let mut simulator = Simulator::new(machine, config.clone());
    let result = simulator.run(&input).expect("Simulation should succeed");
    
    // Verify space is O(√t), not O(t)
    let sqrt_t = (time_bound as f64).sqrt() as usize; // ~100
    let naive_space = time_bound; // O(t) would be ~10,000
    
    assert!(
        result.space_used <= sqrt_t * 50,
        "Space {} should be O(√t) ≈ {}, not O(t) = {}",
        result.space_used,
        sqrt_t,
        naive_space
    );
    
    // Efficiency ratio: should use much less than naive
    let efficiency_ratio = result.space_used as f64 / naive_space as f64;
    assert!(
        efficiency_ratio < 0.1, // Less than 10% of naive
        "Space efficiency ratio {} should be << 1 (ideally < 0.1)",
        efficiency_ratio
    );
}
