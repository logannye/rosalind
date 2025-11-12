use std::ops::Range;
use std::sync::Arc;

use crate::framework::{
    BlockContext, BlockProcessor, CompressedEvaluator, EvaluatorConfig, FrameworkError,
};
use crate::genomics::{align_within_block, BlockAlignmentSummary, BlockedFMIndex, FMInterval};
use thiserror::Error;

/// Result of aligning a single read.
#[derive(Debug, Clone)]
pub struct AlignmentResult {
    /// FM interval representing candidate suffixes matching the read.
    pub interval: FMInterval,
    /// Number of bases from the read that participated in the search.
    pub processed_bases: usize,
    /// Count of mismatches or aborted iterations.
    pub mismatches: usize,
    /// Span of block identifiers that contributed to the result.
    pub span: (usize, usize),
    /// Heuristic score used for ranking alignments.
    pub score: f32,
}

impl AlignmentResult {
    /// Returns `true` when there are remaining candidate suffixes.
    pub fn has_candidates(&self) -> bool {
        !self.interval.is_empty()
    }
}

/// Alignment workload passed to the block processor.
#[derive(Debug, Clone)]
struct AlignmentWorkload {
    read: Arc<[u8]>,
}

impl AlignmentWorkload {
    fn new(read: &[u8]) -> Self {
        Self {
            read: Arc::from(read.to_vec()),
        }
    }

    fn as_slice(&self) -> &[u8] {
        &self.read
    }
}

/// Errors surfaced by the aligner.
#[derive(Debug, Error)]
pub enum AlignerError {
    /// Failure while constructing or querying the FM index.
    #[error("fm-index error: {0}")]
    FMIndex(#[from] crate::genomics::FMIndexError),

    /// Failure reported by the compressed evaluation framework.
    #[error("framework error: {0}")]
    Framework(#[from] FrameworkError),
}

/// FM-index backed aligner that executes alignment using the O(√t) framework.
#[derive(Debug)]
pub struct BWTAligner {
    fm_index: Arc<BlockedFMIndex>,
    evaluator: CompressedEvaluator<BWTAlignmentProcessor>,
}

impl BWTAligner {
    /// Construct a new aligner from the reference genome.
    pub fn new(reference: &[u8]) -> Result<Self, AlignerError> {
        let suggested_block_size = ((reference.len() as f64).sqrt().ceil() as usize).max(64);
        let fm_index = Arc::new(BlockedFMIndex::build(reference, suggested_block_size)?);

        let config = EvaluatorConfig::with_block_size(fm_index.len(), fm_index.block_size())?
            .with_workspace_bytes(fm_index.block_size())
            .with_space_profiling(false);

        let processor = BWTAlignmentProcessor::new(Arc::clone(&fm_index));
        let evaluator = CompressedEvaluator::new(processor, config);

        Ok(Self {
            fm_index,
            evaluator,
        })
    }

    /// Align a read and return its FM-interval summary.
    pub fn align_read(&mut self, read: &[u8]) -> Result<AlignmentResult, AlignerError> {
        let workload = AlignmentWorkload::new(read);
        let evaluation = self.evaluator.evaluate(&workload)?;
        Ok(evaluation.output)
    }

    /// Align multiple reads, returning a vector of results.
    pub fn align_batch<I>(&mut self, reads: I) -> Result<Vec<AlignmentResult>, AlignerError>
    where
        I: IntoIterator,
        I::Item: AsRef<[u8]>,
    {
        let mut results = Vec::new();
        for read in reads {
            results.push(self.align_read(read.as_ref())?);
        }
        Ok(results)
    }

    /// Reference to the underlying FM-index.
    pub fn fm_index(&self) -> &BlockedFMIndex {
        &self.fm_index
    }
}

#[derive(Debug, Clone)]
struct BWTAlignmentProcessor {
    fm_index: Arc<BlockedFMIndex>,
}

impl BWTAlignmentProcessor {
    fn new(fm_index: Arc<BlockedFMIndex>) -> Self {
        Self { fm_index }
    }

    fn segment_for_block(&self, read_len: usize, block_id: usize) -> Range<usize> {
        let num_blocks = self.fm_index.num_blocks().max(1);
        let block_idx = block_id.saturating_sub(1);
        let chunk = (read_len + num_blocks - 1) / num_blocks;
        let start = (block_idx).saturating_mul(chunk);
        let end = start.saturating_add(chunk).min(read_len);
        start..end
    }
}

impl BlockProcessor for BWTAlignmentProcessor {
    type Input = AlignmentWorkload;
    type BlockSummary = BlockAlignmentSummary;
    type Output = AlignmentResult;

    fn process_block(
        &mut self,
        input: &Self::Input,
        context: &BlockContext,
        workspace: &mut [u8],
    ) -> Result<Self::BlockSummary, FrameworkError> {
        let read = input.as_slice();
        let segment = self.segment_for_block(read.len(), context.block_id);

        if segment.is_empty() {
            return Ok(BlockAlignmentSummary::empty(
                context.block_id - 1,
                FMInterval::full(self.fm_index.len()),
            ));
        }

        align_within_block(
            &read[segment.clone()],
            &self.fm_index,
            context.block_id - 1,
            workspace,
        )
        .map_err(|err| FrameworkError::processor_failure(err.to_string()))
    }

    fn merge(
        &mut self,
        left: &Self::BlockSummary,
        right: &Self::BlockSummary,
    ) -> Result<Self::BlockSummary, FrameworkError> {
        Ok(left.merge_with(right))
    }

    fn finalize(
        &mut self,
        root: &Self::BlockSummary,
        _input: &Self::Input,
    ) -> Result<Self::Output, FrameworkError> {
        Ok(AlignmentResult {
            interval: root.interval,
            processed_bases: root.processed_bases,
            mismatches: root.mismatches,
            span: root.span,
            score: root.score,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn aligner_returns_interval_for_read() {
        let reference = b"ACGTCGTAACGT";
        let mut aligner = BWTAligner::new(reference).expect("builder should succeed");
        let result = aligner
            .align_read(b"CGTA")
            .expect("alignment should succeed");
        assert!(result.has_candidates());
        assert!(result.processed_bases > 0);
    }
}
