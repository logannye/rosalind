use std::sync::Arc;

use crate::framework::{BlockContext, BlockProcessor, CompressedEvaluator, EvaluatorConfig, EvaluationResult, FrameworkError};

/// Trait implemented by plugins that run on top of the compressed evaluation engine.
pub trait GenomicPlugin: Send + Sync + 'static {
    /// Type describing the full plugin input (e.g. reads, reference window).
    type Input: Send + Sync;
    /// Type emitted per processed block.
    type BlockSummary: Send + Sync;
    /// Final result type returned after evaluation.
    type Output: Send + Sync;

    /// Unique plugin name.
    fn name(&self) -> &'static str;

    /// Human-readable description.
    fn description(&self) -> &'static str;

    /// Compute total units (logical length) for the given input.
    fn total_units(&self, input: &Self::Input) -> usize;

    /// Select block size for the given input (defaults to √n in many cases).
    fn block_size(&self, input: &Self::Input) -> usize;

    /// Allocate workspace size in bytes (defaults to block size).
    fn workspace_bytes(&self, _input: &Self::Input) -> usize {
        self.block_size(_input)
    }

    /// Process a single block with O(b) workspace.
    fn process_block(
        &self,
        input: &Self::Input,
        context: &BlockContext,
        workspace: &mut [u8],
    ) -> Result<Self::BlockSummary, FrameworkError>;

    /// Merge summaries from adjacent blocks.
    fn merge_summaries(
        &self,
        left: &Self::BlockSummary,
        right: &Self::BlockSummary,
    ) -> Result<Self::BlockSummary, FrameworkError>;

    /// Finalize the result at the root.
    fn finalize(
        &self,
        root: &Self::BlockSummary,
        input: &Self::Input,
    ) -> Result<Self::Output, FrameworkError>;
}

/// Executes a plugin by adapting it to the block-processing interface.
#[derive(Debug)]
pub struct PluginExecutor<P: GenomicPlugin> {
    plugin: Arc<P>,
}

impl<P: GenomicPlugin> PluginExecutor<P> {
    /// Create a new executor for the supplied plugin.
    pub fn new(plugin: Arc<P>) -> Self {
        Self { plugin }
    }

    fn build_config(&self, input: &P::Input) -> Result<EvaluatorConfig, FrameworkError> {
        let total_units = self.plugin.total_units(input);
        let block_size = self.plugin.block_size(input);
        EvaluatorConfig::with_block_size(total_units, block_size).map(|config| {
            config
                .with_workspace_bytes(self.plugin.workspace_bytes(input))
                .with_space_profiling(false)
        })
    }

    /// Run plugin evaluation and return the full evaluation result.
    pub fn run(
        &self,
        input: &P::Input,
    ) -> Result<EvaluationResult<P::Output, P::BlockSummary>, FrameworkError> {
        let config = self.build_config(input)?;
        let adapter = PluginAdapter {
            plugin: Arc::clone(&self.plugin),
        };
        let mut evaluator = CompressedEvaluator::new(adapter, config);
        evaluator.evaluate(input)
    }
}

struct PluginAdapter<P: GenomicPlugin> {
    plugin: Arc<P>,
}

impl<P: GenomicPlugin> BlockProcessor for PluginAdapter<P> {
    type Input = P::Input;
    type BlockSummary = P::BlockSummary;
    type Output = P::Output;

    fn process_block(
        &mut self,
        input: &Self::Input,
        context: &BlockContext,
        workspace: &mut [u8],
    ) -> Result<Self::BlockSummary, FrameworkError> {
        self.plugin.process_block(input, context, workspace)
    }

    fn merge(
        &mut self,
        left: &Self::BlockSummary,
        right: &Self::BlockSummary,
    ) -> Result<Self::BlockSummary, FrameworkError> {
        self.plugin.merge_summaries(left, right)
    }

    fn finalize(
        &mut self,
        root: &Self::BlockSummary,
        input: &Self::Input,
    ) -> Result<Self::Output, FrameworkError> {
        self.plugin.finalize(root, input)
    }
}

