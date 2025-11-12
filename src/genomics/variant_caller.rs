use std::sync::Arc;

use crate::framework::{CompressedEvaluator, EvaluatorConfig, FrameworkError};
use crate::genomics::{
    bayesian_variant_caller, PileupProcessor, PileupSummary, PileupWorkload, VariantCall,
};
use crate::genomics::{AlignedRead, PileupNode};
use thiserror::Error;

/// Variant identified from streaming pileup analysis.
#[derive(Debug, Clone, PartialEq)]
pub struct Variant {
    /// Chromosome/contig name.
    pub chrom: Arc<str>,
    /// Genomic coordinate (0-based).
    pub position: u32,
    /// Reference base observed in the input genome.
    pub reference: u8,
    /// Alternate base inferred from reads.
    pub alternate: u8,
    /// Total read depth supporting reference or alternate.
    pub depth: u32,
    /// Heuristic quality score for the variant.
    pub quality: f32,
    /// Fraction of reads supporting the alternate allele.
    pub allele_fraction: f32,
}

/// Errors originating from variant calling.
#[derive(Debug, Error)]
pub enum VariantCallerError {
    /// Error propagated from the compressed evaluation framework.
    #[error("framework error: {0}")]
    Framework(#[from] FrameworkError),
}

/// Streaming variant caller built on the compressed evaluation framework.
#[derive(Debug)]
pub struct StreamingVariantCaller {
    evaluator: CompressedEvaluator<PileupProcessor>,
    reference: Arc<[u8]>,
    chrom: Arc<str>,
    region_start: u32,
    bases_per_block: usize,
    quality_threshold: f32,
    prior: f32,
}

impl StreamingVariantCaller {
    /// Create a new caller from the reference window and configuration.
    pub fn new(
        chrom: impl Into<Arc<str>>,
        reference: Arc<[u8]>,
        region_start: u32,
        bases_per_block: usize,
        quality_threshold: f32,
        prior: f32,
    ) -> Result<Self, FrameworkError> {
        let total_units = reference.len();
        let config = EvaluatorConfig::with_block_size(total_units, bases_per_block)?
            .with_workspace_bytes(bases_per_block)
            .with_space_profiling(false);

        let processor = PileupProcessor::new();
        let evaluator = CompressedEvaluator::new(processor, config);

        Ok(Self {
            evaluator,
            reference,
            chrom: chrom.into(),
            region_start,
            bases_per_block,
            quality_threshold,
            prior,
        })
    }

    /// Call variants from a batch of aligned reads.
    pub fn call_variants(
        &mut self,
        reads: Vec<AlignedRead>,
    ) -> Result<Vec<Variant>, VariantCallerError> {
        let region = self.region_start..(self.region_start + self.reference.len() as u32);
        let workload = PileupWorkload::new(reads, region, self.bases_per_block);
        let evaluation = self.evaluator.evaluate(&workload)?;
        let summary = evaluation.output;
        Ok(self.extract_variants(&summary))
    }

    fn extract_variants(&self, summary: &PileupSummary) -> Vec<Variant> {
        summary
            .nodes
            .iter()
            .filter_map(|node| self.call_variant_for_node(node))
            .collect()
    }

    fn call_variant_for_node(&self, node: &PileupNode) -> Option<Variant> {
        let offset = (node.position - self.region_start) as usize;
        if offset >= self.reference.len() {
            return None;
        }
        let reference_base = self.reference[offset];
        let VariantCall {
            alt_base,
            quality,
            allele_fraction,
        } = bayesian_variant_caller(node, reference_base, self.prior)?;

        if quality < self.quality_threshold {
            return None;
        }

        Some(Variant {
            chrom: Arc::clone(&self.chrom),
            position: node.position,
            reference: reference_base,
            alternate: alt_base,
            depth: node.depth,
            quality,
            allele_fraction,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genomics::{CigarOp, CigarOpKind};

    #[test]
    fn streaming_variant_caller_produces_variants() {
        let reference = Arc::from(b"ACGTACGT".to_vec().into_boxed_slice());
        let chrom = Arc::from("chr1");

        let reads = vec![
            AlignedRead::new(
                Arc::clone(&chrom),
                0,
                60,
                vec![CigarOp::new(CigarOpKind::Match, 4)],
                b"ACGT".to_vec(),
                vec![30; 4],
                false,
            ),
            AlignedRead::new(
                Arc::clone(&chrom),
                2,
                55,
                vec![CigarOp::new(CigarOpKind::Match, 4)],
                b"GTAA".to_vec(),
                vec![25; 4],
                false,
            ),
        ];

        let mut caller = StreamingVariantCaller::new(
            Arc::clone(&chrom),
            Arc::clone(&reference),
            0,
            4,
            5.0,
            1e-6,
        )
        .unwrap();

        let variants = caller.call_variants(reads).unwrap();
        assert!(!variants.is_empty());
    }
}
