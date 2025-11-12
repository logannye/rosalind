//! Basic end-to-end simulation example

use sqrt_space_sim::{SimulationConfig, Simulator, TuringMachine};

fn main() -> anyhow::Result<()> {
    // Define a simple binary increment machine
    // Input: binary number (LSB first)
    // Output: incremented binary number

    let machine = TuringMachine::builder()
        .num_tapes(1)
        .alphabet(vec!['_', '0', '1'])
        // TODO: Add transition rules for binary increment
        .build()?;

    let input = vec!['1', '1', '0', '1']; // 1011 = 13
    let time_bound = 100;

    let config = SimulationConfig::optimal_for_time(time_bound);
    println!("Block size b = {}", config.block_size);
    println!("Number of blocks T = {}", config.num_blocks);
    println!("Expected space ≤ O(√t) = O({})", config.sqrt_t_bound());

    let mut simulator = Simulator::new(machine, config.clone());
    let result = simulator.run(&input)?;

    println!(
        "\nResult: {}",
        if result.accepted { "ACCEPT" } else { "REJECT" }
    );
    println!("Space used: {} cells", result.space_used);
    println!("Space bound: {} cells", config.sqrt_t_bound());
    println!("Bound satisfied: {}", result.satisfies_bound(&config));

    Ok(())
}
