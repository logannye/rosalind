//! Extract a panel's evidence once and aggregate every original BED target.
//! Usage: cargo run --example exact_panel -- sample.bam genome.fa panel.bed
//! The FASTA requires .fai and BAM requires BAI/CSI. Use `rosalind evidence` /
//! `rosalind panel-qc` for transactional CLI artifacts and portable receipts.
use rosalind::evidence::{EvidenceEngine, EvidenceRequest, PanelQcAnalyzer, PanelTarget};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() != 3 {
        return Err("usage: exact_panel ALIGNMENTS REFERENCE TARGETS.bed".into());
    }
    let request = EvidenceRequest::new(&args[0], &args[1]);
    let mut engine = EvidenceEngine::open(request)?;
    let targets = PanelTarget::from_bed(&args[2], engine.contigs())?;
    let mut analyzer = PanelQcAnalyzer::new(targets)?.with_min_callable_depth(20);
    engine.set_selection(analyzer.selection())?;
    let plan = engine.plan_for_analyzer(&analyzer)?;
    eprintln!(
        "{} requested loci; {}-base execution windows",
        plan.selected_loci, plan.microtile_bases
    );
    engine.run(&mut analyzer)?;
    analyzer.write_tsv(std::io::stdout().lock(), engine.contigs())?;
    Ok(())
}
