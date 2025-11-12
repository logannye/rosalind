//! Plugin framework enabling custom genomic analyses to reuse the O(√t) engine.

mod api;
mod examples;
mod registry;

pub use api::{GenomicPlugin, PluginExecutor};
pub use examples::RNASeqQuantification;
pub use registry::PluginRegistry;
