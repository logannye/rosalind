//! Standalone, bounded research reducer using only public Rosalind APIs.
use rosalind::evidence::{
    EvidenceAnalyzer, EvidenceBatch, EvidenceEngine, EvidenceError, EvidenceFields,
    EvidenceRequest, EvidenceRequirements, EvidenceSelection,
};

#[derive(Debug, Default)]
struct CandidateSummary {
    loci: u64,
    callable_reads: u64,
    candidate_alt_reads: u64,
}

fn add(counter: &mut u64, value: u64) -> Result<(), EvidenceError> {
    *counter = counter
        .checked_add(value)
        .ok_or(EvidenceError::CounterOverflow)?;
    Ok(())
}

impl EvidenceAnalyzer for CandidateSummary {
    fn requirements(&self) -> EvidenceRequirements {
        EvidenceRequirements {
            fields: EvidenceFields::DEPTHS.union(EvidenceFields::ALLELES),
            requires_reference: true,
            context_bases: 0,
            retained_bytes: Some(std::mem::size_of::<Self>() as u64),
        }
    }

    fn on_batch(&mut self, batch: &EvidenceBatch) -> Result<(), EvidenceError> {
        for row in &batch.rows {
            add(&mut self.loci, 1)?;
            add(&mut self.callable_reads, row.callable_depth)?;
            for alternate in &row.requested_alts {
                let index = b"ACGT"
                    .iter()
                    .position(|base| base == alternate)
                    .ok_or_else(|| EvidenceError::Analyzer("invalid SNV allele".into()))?;
                add(&mut self.candidate_alt_reads, row.allele_counts[index])?;
            }
        }
        Ok(())
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args_os().skip(1).collect();
    if args.len() != 3 {
        return Err(
            "usage: rosalind-example-evidence-analyzer ALIGNMENTS REFERENCE CANDIDATES.vcf".into(),
        );
    }
    let mut request = EvidenceRequest::new(&args[0], &args[1]);
    request.execution.memory_budget_bytes = Some(128 << 20);
    let mut engine = EvidenceEngine::open(request)?;
    let selection = EvidenceSelection::from_vcf(&args[2], engine.contigs())?;
    engine.set_selection(selection)?;
    let mut analyzer = CandidateSummary::default();
    eprintln!("plan: {:?}", engine.plan_for_analyzer(&analyzer)?);
    engine.run(&mut analyzer)?;
    // Publish this tiny result only after the complete traversal succeeded.
    println!("selected_loci\tcallable_reads\tcandidate_alt_reads");
    println!(
        "{}\t{}\t{}",
        analyzer.loci, analyzer.callable_reads, analyzer.candidate_alt_reads
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use rosalind::evidence::EvidenceRow;

    #[test]
    fn counts_requested_alts_once_and_keeps_zero_loci() {
        let mut analyzer = CandidateSummary::default();
        let batch = EvidenceBatch {
            contig_id: 0,
            contig: "chr1".into(),
            canonical_tile_start: 0,
            rows: vec![
                EvidenceRow {
                    callable_depth: 10,
                    allele_counts: [6, 3, 1, 0],
                    requested_alts: vec![b'C', b'G'],
                    ..EvidenceRow::default()
                },
                EvidenceRow::default(),
            ],
        };
        analyzer.on_batch(&batch).unwrap();
        assert_eq!(
            (
                analyzer.loci,
                analyzer.callable_reads,
                analyzer.candidate_alt_reads
            ),
            (2, 10, 4)
        );
        assert_eq!(analyzer.requirements().retained_bytes, Some(24));
    }
}
