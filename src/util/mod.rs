//! Utility functions

mod index_free;

pub use index_free::IndexFreeIterator;

/// Marker-based scanning utilities
#[allow(dead_code)]
#[derive(Debug)]
pub struct MarkerScanner {
    // TODO: Implement marker-based iteration
    // Used for O(1) loop control instead of O(log b) counters
}