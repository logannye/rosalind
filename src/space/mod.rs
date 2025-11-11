//! Space accounting and profiling
//!
//! Tracks space usage to verify O(√t) bound

mod allocator;

pub use allocator::SpaceTracker;

/// Detailed space profile (if enabled)
#[derive(Debug, Clone)]
pub struct SpaceProfile {
    /// Maximum space used
    pub max_space: usize,
    
    /// Space over time (snapshots)
    pub timeline: Vec<(usize, usize)>, // (time_point, space_used)
    
    /// Breakdown by component
    pub leaf_buffer_max: usize,
    /// Maximum stack depth reached
    pub stack_depth_max: usize,
    /// Size of streaming ledger
    pub ledger_size: usize,
}

impl SpaceProfile {
    /// Verify bound is satisfied
    pub fn satisfies_bound(&self, bound: usize) -> bool {
        self.max_space <= bound
    }
    
    /// Generate report
    pub fn report(&self) -> String {
        // TODO: Format detailed report
        format!(
            "Max space: {} cells\nComponents:\n  Leaf: {}\n  Stack: {}\n  Ledger: {}",
            self.max_space,
            self.leaf_buffer_max,
            self.stack_depth_max,
            self.ledger_size
        )
    }
}