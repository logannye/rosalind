//! Python bindings that expose selected Rosalind components via PyO3.
use std::ops::Range;
use std::sync::Arc;

use pyo3::{exceptions::PyRuntimeError, prelude::*, types::PyModule};

use crate::genomics::{AlignedRead, CigarOp, CigarOpKind, PileupWorkload};
use crate::plugin::{PluginExecutor, PluginRegistry, RNASeqQuantification};

/// Python-facing entry point for running Rosalind plugins.
#[pyclass]
#[derive(Debug)]
pub struct PyGenomicEngine {
    registry: PluginRegistry,
}

#[pymethods]
impl PyGenomicEngine {
    #[new]
    /// Create the engine with built-in plugins registered.
    pub fn new() -> Self {
        let mut registry = PluginRegistry::new();
        registry.register(RNASeqQuantification::default());
        Self { registry }
    }

    /// List the names of all registered plugins.
    pub fn list_plugins(&self) -> Vec<String> {
        self.registry
            .list()
            .into_iter()
            .map(|info| info.name)
            .collect()
    }

    /// Run the RNA-seq quantification example plugin.
    ///
    /// Args:
    ///     region_start: Start position (inclusive).
    ///     region_end: End position (exclusive).
    ///     reads: List of `(position, sequence)` tuples.
    ///     block_size: Number of bases per evaluation block.
    ///
    /// Returns:
    ///     List of `(position, depth)` tuples.
    pub fn run_rna_seq_plugin(
        &self,
        region_start: u32,
        region_end: u32,
        reads: Vec<(u32, String)>,
        block_size: usize,
    ) -> PyResult<Vec<(u32, u32)>> {
        let plugin = self
            .registry
            .get::<RNASeqQuantification>("rna_seq_quant")
            .ok_or_else(|| PyRuntimeError::new_err("RNA-seq plugin not registered"))?;

        let workload = build_workload(region_start..region_end, reads, block_size)?;
        let executor = PluginExecutor::new(plugin);
        let evaluation = executor
            .run(&workload)
            .map_err(|err| PyRuntimeError::new_err(err.to_string()))?;

        Ok(evaluation.output)
    }
}

fn build_workload(
    region: Range<u32>,
    reads: Vec<(u32, String)>,
    block_size: usize,
) -> PyResult<PileupWorkload> {
    let mut aligned_reads = Vec::with_capacity(reads.len());
    for (pos, seq) in reads {
        let sequence = seq.trim().to_ascii_uppercase();
        let seq_bytes = sequence.as_bytes().to_vec();
        if seq_bytes.is_empty() {
            continue;
        }
        aligned_reads.push(AlignedRead::new(
            Arc::from("chr1"),
            pos,
            60,
            vec![CigarOp::new(CigarOpKind::Match, seq_bytes.len() as u32)],
            seq_bytes,
            vec![30u8; sequence.len()],
            false,
        ));
    }

    Ok(PileupWorkload::new(aligned_reads, region, block_size))
}

/// Create Python module.
#[pymodule]
pub fn rosalind_py(_py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyGenomicEngine>()?;
    Ok(())
}
