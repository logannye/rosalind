use crate::framework::{FrameworkError, BlockProcessor};
use crate::genomics::{PileupProcessor, PileupSummary, PileupWorkload};
use crate::plugin::GenomicPlugin;

/// Example plugin that produces expression depth per genomic position.
#[derive(Debug, Default)]
pub struct RNASeqQuantification;

impl GenomicPlugin for RNASeqQuantification {
    type Input = PileupWorkload;
    type BlockSummary = PileupSummary;
    type Output = Vec<(u32, u32)>; // (position, depth)

    fn name(&self) -> &'static str {
        "rna_seq_quant"
    }

    fn description(&self) -> &'static str {
        "Compute per-base coverage for RNA-seq reads."
    }

    fn total_units(&self, input: &Self::Input) -> usize {
        (input.region.end - input.region.start) as usize
    }

    fn block_size(&self, input: &Self::Input) -> usize {
        input.bases_per_block
    }

    fn process_block(
        &self,
        input: &Self::Input,
        context: &crate::framework::BlockContext,
        workspace: &mut [u8],
    ) -> Result<Self::BlockSummary, FrameworkError> {
        let mut processor = PileupProcessor::new();
        processor.process_block(input, context, workspace)
    }

    fn merge_summaries(
        &self,
        left: &Self::BlockSummary,
        right: &Self::BlockSummary,
    ) -> Result<Self::BlockSummary, FrameworkError> {
        Ok(left.merge(right))
    }

    fn finalize(
        &self,
        root: &Self::BlockSummary,
        _input: &Self::Input,
    ) -> Result<Self::Output, FrameworkError> {
        Ok(root
            .nodes
            .iter()
            .map(|node| (node.position, node.depth))
            .collect())
    }

    fn workspace_bytes(&self, input: &Self::Input) -> usize {
        input.bases_per_block
    }
}

