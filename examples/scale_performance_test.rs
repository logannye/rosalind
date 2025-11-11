//! Scale Performance Test: O(√t) Space Complexity Verification
//!
//! This program runs simulations with progressively larger time bounds
//! and measures space usage, execution time, and scaling behavior to verify
//! that the implementation achieves the claimed O(√t) space complexity.
#![allow(dead_code)]
#![allow(unused_variables)]
use sqrt_space_sim::{Simulator, SimulationConfig};
use std::time::Instant;

/// Format number with commas for readability
fn format_number(n: usize) -> String {
    let s = n.to_string();
    let mut result = String::new();
    for (i, c) in s.chars().rev().enumerate() {
        if i > 0 && i % 3 == 0 {
            result.push(',');
        }
        result.push(c);
    }
    result.chars().rev().collect()
}

/// Format ratio with scientific notation for small values
/// Uses scientific notation when ratio < 0.001, otherwise uses decimal notation
fn format_ratio(ratio: f64) -> String {
    if ratio.abs() < 0.001 && ratio != 0.0 {
        // Use scientific notation for very small values
        format!("{:>10.4e}", ratio)
    } else {
        // Use decimal notation with sufficient precision for larger values
        format!("{:>10.8}", ratio)
    }
}

// Import scale test helpers
#[path = "../tests/scale_test_helpers.rs"]
mod scale_test_helpers;

use scale_test_helpers::create_non_halting_machine;

/// Metrics collected for a single test run
#[derive(Debug, Clone)]
struct TestMetrics {
    time_bound: usize,
    block_size: usize,
    num_blocks: usize,
    space_used: usize,
    sqrt_t_bound: usize,
    execution_time_secs: f64,
    space_efficiency_ratio: f64, // space / t
    scaling_ratio: Option<f64>, // space growth / time growth
}

impl TestMetrics {
    fn new(
        time_bound: usize,
        config: &SimulationConfig,
        space_used: usize,
        execution_time_secs: f64,
    ) -> Self {
        let space_efficiency_ratio = if time_bound > 0 {
            space_used as f64 / time_bound as f64
        } else {
            0.0
        };

        Self {
            time_bound,
            block_size: config.block_size,
            num_blocks: config.num_blocks,
            space_used,
            sqrt_t_bound: config.sqrt_t_bound(),
            execution_time_secs,
            space_efficiency_ratio,
            scaling_ratio: None,
        }
    }

    fn calculate_scaling_ratio(&mut self, previous: &TestMetrics) {
        if previous.time_bound > 0 && self.time_bound > 0 {
            let time_ratio = self.time_bound as f64 / previous.time_bound as f64;
            let space_ratio = self.space_used as f64 / previous.space_used as f64;
            
            if time_ratio > 0.0 {
                self.scaling_ratio = Some(space_ratio);
            }
        }
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("Scale Performance Test: O(√t) Space Complexity Verification");
    println!("==========================================================\n");

    // Time bound progression: each increases by 4x (to test √t scaling)
    // Extended to much larger bounds to test memory efficiency at scale
    let time_bounds = vec![
        100,        // 10^2
        400,        // 4×10^2
        1_600,      // 16×10^2
        6_400,      // 64×10^2
        25_600,     // 256×10^2
        100_000,    // 10^5
        400_000,    // 4×10^5
        1_600_000,  // 16×10^5 = 1.6×10^6
        6_400_000,  // 64×10^5 = 6.4×10^6
        25_600_000, // 256×10^5 = 2.56×10^7
        100_000_000, // 10^8
        400_000_000, // 4×10^8
        1_600_000_000, // 16×10^8
        6_400_000_000, // 64×10^8
        25_600_000_000, // 256×10^8
        // 100_000_000_000, // 10^9
    ];
    
    // Create a machine that won't halt early
    let machine = create_non_halting_machine();
    let input = vec!['_']; // Minimal input since machine doesn't use it
    
    // Collect metrics for each time bound
    let mut all_metrics = Vec::new();
    
    println!("Running scale tests...");
    println!("Note: Large time bounds (t > 1,000,000) may take several minutes to complete.\n");
    println!("{:<15} | {:>10} | {:>10} | {:>12} | {:>12} | {:>10} | {:>10} | {:>10}",
             "Time Bound", "Block Size", "Blocks", "Space Used", "√t Bound", "Space/t", "Time(s)", "Scaling");
    println!("{}", "-".repeat(120));
    
    for (i, &time_bound) in time_bounds.iter().enumerate() {
        // Print progress for large tests
        if time_bound >= 1_000_000 {
            println!("Running test for t={} (this may take a while)...", time_bound);
        }
        
        let config = SimulationConfig {
            profile_space: false, // Set to true for detailed component breakdown
            ..SimulationConfig::optimal_for_time(time_bound)
        };
        
        let mut simulator = Simulator::new(machine.clone(), config.clone());
        
        // Measure execution time
        let start = Instant::now();
        let result = match simulator.run(&input) {
            Ok(result) => result,
            Err(e) => {
                eprintln!("Error at t={}: {:?}", time_bound, e);
                eprintln!("Skipping remaining larger time bounds due to error.");
                break; // Stop on error to avoid wasting time on larger bounds
            }
        };
        let elapsed = start.elapsed().as_secs_f64();
        
        // Create metrics
        let mut metrics = TestMetrics::new(
            time_bound,
            &config,
            result.space_used,
            elapsed,
        );
        
        // Calculate scaling ratio if we have previous metrics
        if i > 0 {
            metrics.calculate_scaling_ratio(&all_metrics[i - 1]);
        }
        
        all_metrics.push(metrics.clone());
        
        // Print formatted row
        let scaling_str = if let Some(ratio) = metrics.scaling_ratio {
            format!("{:.2}x", ratio)
        } else {
            "-".to_string()
        };
        
        // Format large numbers with commas for readability
        let time_str = format_number(time_bound);
        let space_str = format_number(metrics.space_used);
        let bound_str = format_number(metrics.sqrt_t_bound);
        let ratio_str = format_ratio(metrics.space_efficiency_ratio);
        
        println!("{:<15} | {:>10} | {:>10} | {:>12} | {:>12} | {:>10} | {:>10.3} | {:>10}",
                 time_str,
                 metrics.block_size,
                 format_number(metrics.num_blocks),
                 space_str,
                 bound_str,
                 ratio_str,
                 metrics.execution_time_secs,
                 scaling_str);
    }
    
    println!();
    
    // Summary statistics
    println!("Summary:");
    println!("--------");
    
    // Verify all tests satisfy space bound
    let all_within_bound = all_metrics.iter()
        .all(|m| m.space_used <= m.sqrt_t_bound);
    
    if all_within_bound {
        println!("✓ All tests satisfy space bound (space ≤ O(√t))");
    } else {
        println!("✗ Some tests exceed space bound!");
        for m in &all_metrics {
            if m.space_used > m.sqrt_t_bound {
                println!("  t={}: space {} > bound {}", 
                         m.time_bound, m.space_used, m.sqrt_t_bound);
            }
        }
    }
    
    // Calculate average scaling ratio
    let scaling_ratios: Vec<f64> = all_metrics.iter()
        .filter_map(|m| m.scaling_ratio)
        .collect();
    
    if !scaling_ratios.is_empty() {
        let avg_scaling = scaling_ratios.iter().sum::<f64>() / scaling_ratios.len() as f64;
        println!("✓ Space scaling: {:.2}x average (expected: ~2.00x for 4x time increase)", avg_scaling);
        
        // Verify scaling is closer to √t (2x) than linear (4x)
        let expected_linear = 4.0;
        let expected_sqrt = 2.0;
        let distance_from_linear = (avg_scaling - expected_linear).abs();
        let distance_from_sqrt = (avg_scaling - expected_sqrt).abs();
        
        if distance_from_sqrt < distance_from_linear {
            println!("✓ Scaling is sublinear (closer to √t than linear)");
        } else {
            println!("⚠ Warning: Scaling may not be sublinear");
        }
    }
    
    // Efficiency analysis
    if let Some(last) = all_metrics.last() {
        let naive_space = last.time_bound;
        let efficiency_pct = (1.0 - last.space_efficiency_ratio) * 100.0;
        println!("✓ Efficiency: {:.2}% space savings vs naive O(t) at t={}",
                 efficiency_pct, last.time_bound);
        println!("✓ Largest test: t={} completed in {:.2}s",
                 last.time_bound, last.execution_time_secs);
    }
    
    // Verify O(√t) behavior
    println!();
    println!("O(√t) Verification:");
    println!("-------------------");
    
    if all_metrics.len() >= 2 {
        for i in 1..all_metrics.len() {
            let prev = &all_metrics[i - 1];
            let curr = &all_metrics[i];
            
            let time_ratio = curr.time_bound as f64 / prev.time_bound as f64;
            let space_ratio = curr.space_used as f64 / prev.space_used as f64;
            let sqrt_ratio = time_ratio.sqrt();
            
            println!("t: {} → {} ({}x), space: {} → {} ({:.2}x), √t ratio: {:.2}x",
                     prev.time_bound, curr.time_bound, time_ratio as usize,
                     prev.space_used, curr.space_used, space_ratio,
                     sqrt_ratio);
            
            // Verify space ratio is closer to √t ratio than linear
            let distance_from_linear = (space_ratio - time_ratio).abs();
            let distance_from_sqrt = (space_ratio - sqrt_ratio).abs();
            
            if distance_from_sqrt < distance_from_linear {
                println!("  ✓ Space scales as O(√t), not O(t)");
            } else {
                println!("  ⚠ Space scaling may not be O(√t)");
            }
        }
    }
    
    // Optional: CSV export
    let csv_output = std::env::var("SCALE_TEST_CSV").is_ok();
    if csv_output {
        println!("\nCSV Export:");
        println!("time_bound,block_size,num_blocks,space_used,sqrt_t_bound,space_efficiency,execution_time,scaling_ratio");
        for m in &all_metrics {
            let scaling = m.scaling_ratio.map(|r| r.to_string()).unwrap_or_else(|| "-".to_string());
            println!("{},{},{},{},{},{:.4},{:.3},{}",
                     m.time_bound, m.block_size, m.num_blocks,
                     m.space_used, m.sqrt_t_bound, m.space_efficiency_ratio,
                     m.execution_time_secs, scaling);
        }
    }
    
    Ok(())
}

