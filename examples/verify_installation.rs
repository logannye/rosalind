use std::fs;
use std::path::PathBuf;
use std::sync::Arc;

use anyhow::{Context, Result};
use rosalind::genomics::{
    AlignedRead, BWTAligner, CigarOp, CigarOpKind, StreamingVariantCaller, Variant,
};

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

fn ensure_data_exists(dir: &PathBuf) -> Result<()> {
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

fn call_variants(reference: &[u8]) -> Result<Vec<Variant>> {
    let chrom = Arc::from("chr1");
    let reference_arc = Arc::from(reference.to_vec().into_boxed_slice());

    let reads = vec![AlignedRead::new(
        Arc::clone(&chrom),
        0,
        60,
        vec![CigarOp::new(CigarOpKind::Match, 8)],
        b"ACGTACGT".to_vec(),
        vec![30u8; 8],
        false,
    )];

    let mut caller =
        StreamingVariantCaller::new(Arc::clone(&chrom), reference_arc, 0, 32, 5.0, 1e-6)?;

    Ok(caller.call_variants(reads)?)
}
