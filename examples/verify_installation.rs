use std::fs;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use anyhow::{Context, Result};
use rosalind::call::{call_germline_region, GermlineCall, GermlineParams};
use rosalind::core::{AlignedRead, CigarOp, CigarOpKind, Locus, Position, SamFlags};
use rosalind::genomics::BWTAligner;
use rosalind::pileup::{PileupParams, SliceSource};

fn main() -> Result<()> {
    let data_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("examples/data");
    ensure_data_exists(&data_dir)?;

    let reference = load_reference(&data_dir.join("ref.fa"))?;
    let mut aligner = BWTAligner::new(&reference).context("failed to build aligner")?;

    let read = b"ACGTACGT";
    let summary = aligner
        .align_read(read)
        .context("alignment of sanity check read failed")?;

    println!(
        "Aligner ready: interval=[{}, {}) width={} mismatches={}",
        summary.interval.lower,
        summary.interval.upper,
        summary.interval.width(),
        summary.mismatches
    );

    let variants = call_variants(&reference).context("variant caller failed")?;
    println!("Variant caller ready: {} variants emitted", variants.len());

    println!("Rosalind installation looks good ✅");
    Ok(())
}

fn ensure_data_exists(dir: &Path) -> Result<()> {
    for name in ["ref.fa", "reads.fastq"] {
        let path = dir.join(name);
        if !path.exists() {
            anyhow::bail!("missing required file: {}", path.display());
        }
    }
    Ok(())
}

fn load_reference(path: &PathBuf) -> Result<Vec<u8>> {
    let contents = fs::read_to_string(path)
        .with_context(|| format!("failed to read reference at {}", path.display()))?;
    let sequence = contents
        .lines()
        .filter(|line| !line.starts_with('>') && !line.trim().is_empty())
        .collect::<String>();
    Ok(sequence.to_ascii_uppercase().into_bytes())
}

fn call_variants(reference: &[u8]) -> Result<Vec<(Locus, u8, GermlineCall)>> {
    let reference_arc: Arc<[u8]> = Arc::from(reference.to_vec().into_boxed_slice());

    let reads = vec![AlignedRead {
        contig: 0,
        pos: Position(0),
        mapq: 60,
        flags: SamFlags::default(),
        cigar: vec![CigarOp::new(CigarOpKind::Match, 8)],
        seq: Arc::from(b"ACGTACGT".to_vec().into_boxed_slice()),
        qual: Arc::from(vec![30u8; 8].into_boxed_slice()),
    }];

    let sites = call_germline_region(
        SliceSource::new(reads),
        Arc::clone(&reference_arc),
        0,
        0..reference.len() as u32,
        PileupParams::default(),
        &GermlineParams::default(),
    )?;
    Ok(sites)
}
