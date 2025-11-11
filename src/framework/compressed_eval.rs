use std::ops::Range;

use crate::{
    ledger::StreamingLedger,
    space::{SpaceProfile, SpaceTracker},
    tree::TreeNode,
};
use thiserror::Error;

/// Errors that can occur while using the generic compressed evaluator.
#[derive(Debug, Error)]
pub enum FrameworkError {
    /// Configuration invalid (e.g., zero block size).
    #[error("invalid evaluator configuration: {0}")]
    InvalidConfiguration(String),

    /// Requested block index is out of range for the configured number of blocks.
    #[error("block id {block_id} out of range (max {max_blocks})")]
    BlockOutOfRange {
        /// Block identifier (1-indexed) that was requested.
        block_id: usize,
        /// Maximum valid block identifier.
        max_blocks: usize,
    },

    /// User-supplied processor reported an error.
    #[error("processor error: {0}")]
    Processor(String),

    /// Encountered more space usage than the theoretical bound allows.
    #[error("space bound exceeded: used {used} > bound {bound}")]
    SpaceBoundExceeded {
        /// Actual space used.
        used: usize,
        /// Theoretical space bound.
        bound: usize,
    },

    /// Streaming ledger did not observe all merges.
    #[error("incomplete streaming ledger for {num_blocks} blocks")]
    LedgerIncomplete {
        /// Number of blocks expected.
        num_blocks: usize,
    },
}

impl FrameworkError {
    /// Helper for constructing processor-originated errors.
    pub fn processor_failure(msg: impl Into<String>) -> Self {
        FrameworkError::Processor(msg.into())
    }
}

/// Configuration parameters for generic compressed evaluation.
#[derive(Debug, Clone)]
pub struct EvaluatorConfig {
    /// Block size parameter `b`.
    pub block_size: usize,
    /// Number of blocks `T = ⌈n / b⌉`.
    pub num_blocks: usize,
    /// Total number of logical units (e.g., symbols, bases).
    pub total_units: usize,
    /// Size of reusable workspace buffer for block processing.
    pub workspace_bytes: usize,
    /// Enable detailed space profiling.
    pub profile_space: bool,
    /// Enable verbose logging/hooks for debugging.
    pub verbose: bool,
}

impl EvaluatorConfig {
    /// Construct configuration given total units; chooses `b = ⌈√n⌉`.
    pub fn optimal_for_units(total_units: usize) -> Result<Self, FrameworkError> {
        if total_units == 0 {
            return Err(FrameworkError::InvalidConfiguration(
                "total units must be > 0".to_string(),
            ));
        }
        let block_size = (total_units as f64).sqrt().ceil() as usize;
        Self::with_block_size(total_units, block_size)
    }

    /// Construct configuration with explicit block size.
    pub fn with_block_size(
        total_units: usize,
        block_size: usize,
    ) -> Result<Self, FrameworkError> {
        if block_size == 0 {
            return Err(FrameworkError::InvalidConfiguration(
                "block size must be > 0".to_string(),
            ));
        }
        if total_units == 0 {
            return Err(FrameworkError::InvalidConfiguration(
                "total units must be > 0".to_string(),
            ));
        }
        let num_blocks = (total_units + block_size - 1) / block_size;
        Ok(Self {
            block_size,
            num_blocks,
            total_units,
            workspace_bytes: block_size,
            profile_space: false,
            verbose: false,
        })
    }

    /// Set workspace bytes (defaults to block size).
    pub fn with_workspace_bytes(mut self, workspace_bytes: usize) -> Self {
        self.workspace_bytes = workspace_bytes.max(1);
        self
    }

    /// Enable space profiling.
    pub fn with_space_profiling(mut self, enabled: bool) -> Self {
        self.profile_space = enabled;
        self
    }

    /// Enable verbose mode.
    pub fn with_verbose(mut self, enabled: bool) -> Self {
        self.verbose = enabled;
        self
    }

    /// Theoretical space bound `O(b + T + log T)`.
    pub fn space_bound(&self) -> usize {
        if self.num_blocks == 0 {
            return self.block_size;
        }
        let log_term = (self.num_blocks as f64).log2().ceil() as usize;
        self.block_size + self.num_blocks + log_term
    }

    /// Compute block context (start/end indices) for a given block id (1-indexed).
    pub fn block_context(&self, block_id: usize) -> Result<BlockContext, FrameworkError> {
        if block_id == 0 || block_id > self.num_blocks {
            return Err(FrameworkError::BlockOutOfRange {
                block_id,
                max_blocks: self.num_blocks,
            });
        }
        let start = (block_id - 1) * self.block_size;
        let end = (start + self.block_size).min(self.total_units);
        Ok(BlockContext {
            block_id,
            range: start..end,
        })
    }
}

/// Per-block metadata supplied to processors.
#[derive(Debug, Clone)]
pub struct BlockContext {
    /// 1-indexed block identifier.
    pub block_id: usize,
    /// Range within the logical input covered by the block.
    pub range: Range<usize>,
}

impl BlockContext {
    /// Number of logical units covered by this block.
    pub fn len(&self) -> usize {
        self.range.end.saturating_sub(self.range.start)
    }

    /// Whether the block is empty (should not normally happen unless padding is used).
    pub fn is_empty(&self) -> bool {
        self.range.end <= self.range.start
    }
}

/// Trait implemented by domain-specific processors that operate on blocks.
pub trait BlockProcessor {
    /// Type describing the full input workload.
    type Input;
    /// Summary emitted per block (must be mergeable).
    type BlockSummary;
    /// Final output type produced after evaluation.
    type Output;

    /// Process a single block using `O(b)` workspace.
    fn process_block(
        &mut self,
        input: &Self::Input,
        context: &BlockContext,
        workspace: &mut [u8],
    ) -> Result<Self::BlockSummary, FrameworkError>;

    /// Merge summaries from adjacent blocks (associative).
    fn merge(
        &mut self,
        left: &Self::BlockSummary,
        right: &Self::BlockSummary,
    ) -> Result<Self::BlockSummary, FrameworkError>;

    /// Finalize result at the root.
    fn finalize(
        &mut self,
        root: &Self::BlockSummary,
        input: &Self::Input,
    ) -> Result<Self::Output, FrameworkError>;
}

/// Result returned by the compressed evaluator.
#[derive(Debug)]
pub struct EvaluationResult<O, S> {
    /// Final output synthesized from the root summary.
    pub output: O,
    /// Summary stored at the root of the compressed tree.
    pub root_summary: S,
    /// Maximum space observed during evaluation.
    pub space_used: usize,
    /// Theoretical space bound for this configuration.
    pub space_bound: usize,
    /// Optional detailed space profile.
    pub space_profile: Option<SpaceProfile>,
}

/// Generic evaluator that reuses the height-compressed tree and streaming ledger.
#[derive(Debug)]
pub struct CompressedEvaluator<P: BlockProcessor> {
    processor: P,
    config: EvaluatorConfig,
    workspace: Vec<u8>,
}

impl<P: BlockProcessor> CompressedEvaluator<P> {
    /// Create a new evaluator with the provided processor and configuration.
    pub fn new(processor: P, config: EvaluatorConfig) -> Self {
        Self {
            processor,
            config,
            workspace: Vec::new(),
        }
    }

    /// Access configuration.
    pub fn config(&self) -> &EvaluatorConfig {
        &self.config
    }

    /// Mutable access to configuration (for tuning).
    pub fn config_mut(&mut self) -> &mut EvaluatorConfig {
        &mut self.config
    }

    /// Execute evaluation using the block processor and return the final output and metadata.
    pub fn evaluate(
        &mut self,
        input: &P::Input,
    ) -> Result<EvaluationResult<P::Output, P::BlockSummary>, FrameworkError> {
        if self.config.num_blocks == 0 {
            return Err(FrameworkError::InvalidConfiguration(
                "number of blocks must be > 0".to_string(),
            ));
        }

        let mut tracker = SpaceTracker::new(self.config.profile_space);
        let mut ledger = StreamingLedger::new(self.config.num_blocks);

        let mut workspace = std::mem::take(&mut self.workspace);
        if workspace.len() < self.config.workspace_bytes {
            workspace.resize(self.config.workspace_bytes.max(1), 0);
        }

        let root = TreeNode::root(1, self.config.num_blocks);
        let root_summary = self.evaluate_dfs(
            root,
            input,
            &mut tracker,
            &mut ledger,
            workspace.as_mut_slice(),
        )?;

        self.workspace = workspace;

        if !ledger.all_merges_complete() {
            return Err(FrameworkError::LedgerIncomplete {
                num_blocks: self.config.num_blocks,
            });
        }

        let space_used = tracker.max_space_used();
        let space_bound = self.config.space_bound();
        if space_used > space_bound {
            return Err(FrameworkError::SpaceBoundExceeded {
                used: space_used,
                bound: space_bound,
            });
        }

        let output = self.processor.finalize(&root_summary, input)?;
        let space_profile = tracker.take_profile();

        Ok(EvaluationResult {
            output,
            root_summary,
            space_used,
            space_bound,
            space_profile,
        })
    }

    fn evaluate_dfs(
        &mut self,
        node: TreeNode,
        input: &P::Input,
        tracker: &mut SpaceTracker,
        ledger: &mut StreamingLedger,
        workspace: &mut [u8],
    ) -> Result<P::BlockSummary, FrameworkError> {
        tracker.push_stack_frame(1);

        let result = if node.is_leaf() {
            let leaf_block_id = node.leaf_block_id();
            let context = self.config.block_context(leaf_block_id)?;

            tracker.allocate_leaf_buffer(self.config.workspace_bytes);
            let block_result = self
                .processor
                .process_block(input, &context, workspace);
            tracker.free_leaf_buffer(self.config.workspace_bytes);

            block_result
        } else {
            let (left_child, right_child) = node.children();

            let left_summary =
                self.evaluate_dfs(left_child, input, tracker, ledger, workspace)?;
            ledger.mark_left_complete(node);

            let right_summary =
                self.evaluate_dfs(right_child, input, tracker, ledger, workspace)?;
            ledger.mark_right_complete(node);

            self.processor.merge(&left_summary, &right_summary)
        };

        tracker.pop_stack_frame();
        result
    }
}

