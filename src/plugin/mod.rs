//! Plugin framework enabling custom genomic analyses to reuse the O(√t) engine.

mod api;
mod registry;
mod examples;

pub use api::{GenomicPlugin, PluginExecutor};
pub use registry::PluginRegistry;
pub use examples::RNASeqQuantification;


