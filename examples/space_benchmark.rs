//! Benchmark space usage for various time bounds

use rosalind::{SimulationConfig, Simulator, TuringMachine};

fn main() -> anyhow::Result<()> {
    println!("Space Benchmark: O(√t) Verification");
    println!("=====================================\n");

    // TODO: Create a standard test machine
    let machine = TuringMachine::builder().build()?;

    let time_bounds = [100, 1_000, 10_000, 100_000];

    println!(
        "{:>10} | {:>10} | {:>10} | {:>10} | {:>10}",
        "t", "b", "Space", "√t bound", "Ratio"
    );
    println!("{}", "-".repeat(65));

    for &t in &time_bounds {
        let config = SimulationConfig::optimal_for_time(t);
        let _sim = Simulator::new(machine.clone(), config.clone());

        // TODO: Run with appropriate input
        // let result = sim.run(&input)?;

        // Placeholder for now
        println!(
            "{:>10} | {:>10} | {:>10} | {:>10} | {:>10.2}",
            t,
            config.block_size,
            "TODO", // result.space_used,
            config.sqrt_t_bound(),
            0.0 // ratio
        );
    }

    Ok(())
}
