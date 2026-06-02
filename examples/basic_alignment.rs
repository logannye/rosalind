//! Basic genomics example demonstrating the BWT aligner.

use rosalind::genomics::BWTAligner;

fn main() -> anyhow::Result<()> {
    // Minimal reference contig (re-used in README quick start).
    let reference = include_bytes!("data/ref.fa")
        .split(|&b| b == b'\n')
        .filter(|line| !line.starts_with(&[b'>']))
        .flatten()
        .copied()
        .collect::<Vec<u8>>();

    // Reads that match positions in the reference sequence.
    let reads = ["ACGTACGT", "TTTACGT", "ACGTACGTACGT", "GGGACGT"];

    let mut aligner = BWTAligner::new(&reference)?;

    for (idx, read) in reads.iter().enumerate() {
        let result = aligner.align_read(read.as_bytes())?;
        println!(
            "read {idx}: interval=[{}, {}) width={} mismatches={}",
            result.interval.lower,
            result.interval.upper,
            result.interval.width(),
            result.mismatches
        );
    }

    Ok(())
}
