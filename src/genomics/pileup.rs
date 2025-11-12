use std::ops::Range;
use std::sync::Arc;

use crate::framework::{BlockContext, BlockProcessor, FrameworkError};
use crate::genomics::AlignedRead;

const NUM_BASES: usize = 4; // A, C, G, T

fn base_index(base: u8) -> Option<usize> {
    match base {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' | b'U' | b'u' => Some(3),
        _ => None,
    }
}

/// Aggregated pileup statistics for a genomic position.
#[derive(Debug, Clone, PartialEq)]
pub struct PileupNode {
    /// Genomic coordinate (0-based) of the pileup position.
    pub position: u32,
    /// Per-base observation counts [A, C, G, T].
    pub base_counts: [u32; NUM_BASES],
    /// Sum of normalised quality scores per base.
    pub quality_sums: [f32; NUM_BASES],
    /// Total number of reads covering this position.
    pub depth: u32,
}

impl PileupNode {
    pub(crate) fn new(position: u32) -> Self {
        Self {
            position,
            base_counts: [0; NUM_BASES],
            quality_sums: [0.0; NUM_BASES],
            depth: 0,
        }
    }

    pub(crate) fn observe(&mut self, base_idx: usize, quality: u8) {
        self.base_counts[base_idx] += 1;
        self.quality_sums[base_idx] += (quality as f32) / 93.0; // Normalize to [0,1]
        self.depth += 1;
    }

    pub(crate) fn merge(&self, other: &Self) -> Self {
        debug_assert_eq!(self.position, other.position);
        let mut merged = Self::new(self.position);
        for i in 0..NUM_BASES {
            merged.base_counts[i] = self.base_counts[i] + other.base_counts[i];
            merged.quality_sums[i] = self.quality_sums[i] + other.quality_sums[i];
        }
        merged.depth = self.depth + other.depth;
        merged
    }
}

/// Summary of pileup statistics for a block.
#[derive(Debug, Clone)]
pub struct PileupSummary {
    /// Identifier of the block that generated the summary.
    pub block_id: usize,
    /// Genomic range covered by this summary.
    pub region: Range<u32>,
    /// Ordered pileup nodes for the covered region.
    pub nodes: Vec<PileupNode>,
}

impl PileupSummary {
    /// Construct an empty summary for a given block.
    pub fn empty(block_id: usize, region: Range<u32>) -> Self {
        Self {
            block_id,
            region,
            nodes: Vec::new(),
        }
    }

    /// Merge two adjacent summaries, preserving positional ordering.
    pub fn merge(&self, other: &Self) -> Self {
        let mut merged_nodes = Vec::with_capacity(self.nodes.len() + other.nodes.len());
        let mut i = 0;
        let mut j = 0;

        while i < self.nodes.len() && j < other.nodes.len() {
            let left = &self.nodes[i];
            let right = &other.nodes[j];
            match left.position.cmp(&right.position) {
                std::cmp::Ordering::Less => {
                    merged_nodes.push(left.clone());
                    i += 1;
                }
                std::cmp::Ordering::Greater => {
                    merged_nodes.push(right.clone());
                    j += 1;
                }
                std::cmp::Ordering::Equal => {
                    merged_nodes.push(left.merge(right));
                    i += 1;
                    j += 1;
                }
            }
        }

        merged_nodes.extend_from_slice(&self.nodes[i..]);
        merged_nodes.extend_from_slice(&other.nodes[j..]);

        let start = self.region.start.min(other.region.start);
        let end = self.region.end.max(other.region.end);

        Self {
            block_id: other.block_id,
            region: start..end,
            nodes: merged_nodes,
        }
    }
}

/// Workload consumed by the pileup processor.
#[derive(Debug, Clone)]
pub struct PileupWorkload {
    /// Collection of aligned reads to include in the pileup.
    pub reads: Arc<[AlignedRead]>,
    /// Genomic window to evaluate.
    pub region: Range<u32>,
    /// Desired number of bases per evaluation block.
    pub bases_per_block: usize,
}

impl PileupWorkload {
    /// Construct a new workload given reads and a target region.
    pub fn new(reads: Vec<AlignedRead>, region: Range<u32>, bases_per_block: usize) -> Self {
        Self {
            reads: Arc::from(reads.into_boxed_slice()),
            region,
            bases_per_block,
        }
    }

    fn block_region(&self, block_id: usize) -> Range<u32> {
        let block_idx = block_id.saturating_sub(1) as u32;
        let start = self.region.start + block_idx * self.bases_per_block as u32;
        let end = (start + self.bases_per_block as u32).min(self.region.end);
        start..end
    }
}

/// Block processor constructing pileup summaries.
#[derive(Debug, Clone)]
pub struct PileupProcessor;

impl PileupProcessor {
    /// Create a new pileup processor.
    pub fn new() -> Self {
        Self
    }

    fn build_summary(&self, reads: &[AlignedRead], region: &Range<u32>) -> PileupSummary {
        if region.start >= region.end {
            return PileupSummary::empty(0, region.clone());
        }

        let window_len = (region.end - region.start) as usize;
        let mut nodes: Vec<PileupNode> = (0..window_len)
            .map(|offset| PileupNode::new(region.start + offset as u32))
            .collect();

        for read in reads {
            let read_start = read.pos;
            let read_end = read.end();

            if read_end <= region.start || read_start >= region.end {
                continue;
            }

            let overlap_start = region.start.max(read_start);
            let overlap_end = region.end.min(read_end);
            let read_offset_start = (overlap_start - read_start) as usize;
            let read_offset_end = (overlap_end - read_start) as usize;

            for (offset, node_idx) in
                (read_offset_start..read_offset_end).zip((overlap_start - region.start) as usize..)
            {
                if let Some(base) = read.base_at(offset) {
                    if let Some(idx) = base_index(base) {
                        let qual = read.quality_at(offset).unwrap_or(30);
                        nodes[node_idx].observe(idx, qual);
                    }
                }
            }
        }

        nodes.retain(|node| node.depth > 0);

        PileupSummary {
            block_id: 0,
            region: region.clone(),
            nodes,
        }
    }
}

impl BlockProcessor for PileupProcessor {
    type Input = PileupWorkload;
    type BlockSummary = PileupSummary;
    type Output = PileupSummary;

    fn process_block(
        &mut self,
        input: &Self::Input,
        context: &BlockContext,
        _workspace: &mut [u8],
    ) -> Result<Self::BlockSummary, FrameworkError> {
        let region = input.block_region(context.block_id);
        if region.start >= region.end {
            return Ok(PileupSummary::empty(context.block_id, region));
        }

        let summary = self.build_summary(input.reads.as_ref(), &region);
        Ok(PileupSummary {
            block_id: context.block_id,
            ..summary
        })
    }

    fn merge(
        &mut self,
        left: &Self::BlockSummary,
        right: &Self::BlockSummary,
    ) -> Result<Self::BlockSummary, FrameworkError> {
        Ok(left.merge(right))
    }

    fn finalize(
        &mut self,
        root: &Self::BlockSummary,
        _input: &Self::Input,
    ) -> Result<Self::Output, FrameworkError> {
        Ok(root.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genomics::{CigarOp, CigarOpKind};

    #[test]
    fn pileup_processor_aggregates_counts() {
        let reads = vec![
            AlignedRead::new(
                "chr1",
                100,
                vec![CigarOp::new(CigarOpKind::Match, 4)],
                b"ACGT".to_vec(),
                vec![30; 4],
                false,
            ),
            AlignedRead::new(
                "chr1",
                101,
                vec![CigarOp::new(CigarOpKind::Match, 4)],
                b"CGTA".to_vec(),
                vec![25; 4],
                false,
            ),
        ];

        let workload = PileupWorkload::new(reads, 100..110, 5);
        let mut processor = PileupProcessor::new();
        let context = BlockContext {
            block_id: 1,
            range: 0..5,
        };
        let mut workspace = vec![0u8; 5];

        let summary = processor
            .process_block(&workload, &context, workspace.as_mut_slice())
            .expect("pileup should succeed");

        assert!(!summary.nodes.is_empty());
        let first = &summary.nodes[0];
        assert_eq!(first.position, 100);
        assert_eq!(first.depth, 1);
    }
}
