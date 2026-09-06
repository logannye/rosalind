//! Utility functions

/// Transactional creation of user-facing artifacts.
pub mod atomic;
/// Create-new transactional publication of an artifact directory.
pub mod atomic_directory;
pub mod mmap;
pub mod rss;
