use std::any::Any;
use std::collections::HashMap;
use std::sync::Arc;

use crate::plugin::GenomicPlugin;

/// Metadata describing a registered plugin.
#[derive(Debug, Clone)]
pub struct PluginInfo {
    pub name: String,
    pub description: String,
}

#[derive(Debug)]
struct PluginEntry {
    plugin: Box<dyn Any + Send + Sync>,
    description: String,
}

/// Registry storing available genomic plugins.
#[derive(Debug, Default)]
pub struct PluginRegistry {
    entries: HashMap<String, PluginEntry>,
}

impl PluginRegistry {
    /// Create an empty registry.
    pub fn new() -> Self {
        Self {
            entries: HashMap::new(),
        }
    }

    /// Register a plugin and return an `Arc` handle to it.
    pub fn register<P>(&mut self, plugin: P) -> Arc<P>
    where
        P: GenomicPlugin,
    {
        let arc = Arc::new(plugin);
        let name = arc.name().to_string();
        let description = arc.description().to_string();
        self.entries.insert(
            name.clone(),
            PluginEntry {
                plugin: Box::new(arc.clone()),
                description,
            },
        );
        arc
    }

    /// Retrieve a plugin by name, downcasting to the requested type.
    pub fn get<P>(&self, name: &str) -> Option<Arc<P>>
    where
        P: GenomicPlugin,
    {
        self.entries.get(name).and_then(|entry| {
            entry
                .plugin
                .downcast_ref::<Arc<P>>()
                .map(|arc| Arc::clone(arc))
        })
    }

    /// List all registered plugins.
    pub fn list(&self) -> Vec<PluginInfo> {
        self.entries
            .iter()
            .map(|(name, entry)| PluginInfo {
                name: name.clone(),
                description: entry.description.clone(),
            })
            .collect()
    }
}
