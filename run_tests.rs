//! Simple test runner binary for quick testing
//! Run with: cargo run --bin run_tests

use std::process::Command;

fn main() {
    println!("==========================================");
    println!("sqrt-space-sim Quick Test Runner");
    println!("==========================================");
    println!();
    
    // Run cargo test
    let status = Command::new("cargo")
        .args(&["test", "--", "--nocapture"])
        .status()
        .expect("Failed to run cargo test");
    
    if status.success() {
        println!("\n✓ All tests passed!");
        std::process::exit(0);
    } else {
        println!("\n✗ Some tests failed!");
        std::process::exit(1);
    }
}

