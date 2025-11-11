use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;
use std::sync::Arc;

use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use sqrt_space_sim::genomics::{
    AlignedRead, BWTAligner, CigarOp, CigarOpKind, StreamingVariantCaller, Variant,
};

#[derive(Parser, Debug)]
#[command(name = "rosalind", about = "Genomic analysis engine using O(√t) space")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Align reads against a reference genome using the BWT aligner.
    Align {
        /// Reference genome (plain FASTA without headers or raw sequence file).
        reference: PathBuf,
        /// Reads file (one sequence per line).
        reads: PathBuf,
    },
    /// Call variants from aligned reads using the streaming variant caller.
    Variants {
        /// Reference genome (plain FASTA or raw sequence).
        reference: PathBuf,
        /// Alignments file (`<position>\t<sequence>` per line).
        alignments: PathBuf,
        /// Chromosome name (default: chr1).
        #[arg(long, default_value = "chr1")]
        chrom: String,
        /// Bases per block for streaming evaluation.
        #[arg(long, default_value_t = 1024)]
        block_size: usize,
        /// Minimum quality threshold for reporting variants.
        #[arg(long, default_value_t = 10.0)]
        quality_threshold: f32,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Align { reference, reads } => run_align(reference, reads)?,
        Commands::Variants {
            reference,
            alignments,
            chrom,
            block_size,
            quality_threshold,
        } => run_variants(reference, alignments, chrom, block_size, quality_threshold)?,
    }

    Ok(())
}

fn run_align(reference_path: PathBuf, reads_path: PathBuf) -> Result<()> {
    let reference = read_sequence_file(&reference_path)
        .with_context(|| format!("failed to read reference from {}", reference_path.display()))?;
    let mut aligner =
        BWTAligner::new(&reference).context("failed to initialize BWT aligner")?;

    let reader = BufReader::new(File::open(&reads_path).with_context(|| {
        format!("failed to open reads file {}", reads_path.display())
    })?);

    for (idx, line) in reader.lines().enumerate() {
        let read = line?.trim().to_string();
        if read.is_empty() {
            continue;
        }
        let result = aligner
            .align_read(read.as_bytes())
            .with_context(|| format!("alignment failed for read {}", idx + 1))?;

        println!(
            "read {}\tinterval=[{}, {})\twidth={}\tscore={:.2}\tmismatches={}",
            idx + 1,
            result.interval.lower,
            result.interval.upper,
            result.interval.width(),
            result.score,
            result.mismatches
        );
    }

    Ok(())
}

fn run_variants(
    reference_path: PathBuf,
    alignments_path: PathBuf,
    chrom: String,
    block_size: usize,
    quality_threshold: f32,
) -> Result<()> {
    let reference_vec = read_sequence_file(&reference_path).with_context(|| {
        format!(
            "failed to read reference from {}",
            reference_path.display()
        )
    })?;
    let reference = Arc::from(reference_vec.into_boxed_slice());
    let region_start = 0u32;

    let chrom_arc = Arc::from(chrom);
    let reads = read_alignment_file(&alignments_path, &chrom_arc)?;

    let mut caller = StreamingVariantCaller::new(
        Arc::clone(&chrom_arc),
        Arc::clone(&reference),
        region_start,
        block_size,
        quality_threshold,
        1e-6,
    )
    .context("failed to initialize variant caller")?;

    let variants = caller
        .call_variants(reads)
        .context("variant calling failed")?;

    if variants.is_empty() {
        println!("No variants detected above threshold.");
    } else {
        for variant in variants {
            print_variant(&variant);
        }
    }

    Ok(())
}

fn read_sequence_file(path: &PathBuf) -> Result<Vec<u8>> {
    let contents = std::fs::read_to_string(path)?;
    let sequence: String = contents
        .lines()
        .filter(|line| !line.starts_with('>') && !line.trim().is_empty())
        .collect::<Vec<_>>()
        .join("");
    Ok(sequence.trim().to_ascii_uppercase().into_bytes())
}

fn read_alignment_file(path: &PathBuf, chrom: &Arc<str>) -> Result<Vec<AlignedRead>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut reads = Vec::new();

    for (line_no, line) in reader.lines().enumerate() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }

        let mut fields = line.split_whitespace();
        let pos_str = fields
            .next()
            .ok_or_else(|| anyhow::anyhow!("missing position on line {}", line_no + 1))?;
        let seq = fields
            .next()
            .ok_or_else(|| anyhow::anyhow!("missing sequence on line {}", line_no + 1))?;

        let pos: u32 = pos_str.parse().with_context(|| {
            format!(
                "invalid position '{}' on line {}",
                pos_str,
                line_no + 1
            )
        })?;

        let sequence = seq.to_ascii_uppercase().into_bytes();
        let qualities = vec![30u8; sequence.len()];

        reads.push(AlignedRead::new(
            Arc::clone(chrom),
            pos,
            vec![CigarOp::new(CigarOpKind::Match, sequence.len() as u32)],
            sequence,
            qualities,
            false,
        ));
    }

    Ok(reads)
}

fn print_variant(variant: &Variant) {
    println!(
        "{}\t{}\t{}\t{}\tdepth={}\tAF={:.3}\tQUAL={:.2}",
        variant.chrom,
        variant.position,
        variant.reference as char,
        variant.alternate as char,
        variant.depth,
        variant.allele_fraction,
        variant.quality
    );
}

