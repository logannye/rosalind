use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::PathBuf;
use std::sync::Arc;

use anyhow::{anyhow, bail, Context, Result};
use clap::{Parser, Subcommand, ValueEnum};
use rosalind::genomics::{
    create_bam_writer, write_vcf, AlignedRead, BWTAligner, CigarOp, CigarOpKind,
    StreamingVariantCaller,
};
use rust_htslib::bam::{
    self, record::Aux, record::Cigar as BamCigar, record::CigarString, record::Record,
};

#[derive(Parser, Debug)]
#[command(name = "rosalind", about = "Genomic analysis engine using O(√t) space")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Align reads against a reference genome and emit SAM records.
    Align {
        /// Reference genome in FASTA format (only the first record is used).
        #[arg(long)]
        reference: PathBuf,
        /// Reads file in FASTQ format.
        #[arg(long)]
        reads: PathBuf,
        /// Maximum mismatches permitted when seeding alignments.
        #[arg(long, default_value_t = 2)]
        max_mismatches: usize,
        /// Offset applied to reported reference positions (1-based in SAM).
        #[arg(long, default_value_t = 0)]
        reference_offset: u32,
        /// Output format for the alignment.
        #[arg(long, value_enum, default_value_t = OutputFormat::Sam)]
        format: OutputFormat,
        /// Optional path to write the output (stdout if omitted for SAM).
        #[arg(short, long)]
        output: Option<PathBuf>,
    },
    /// Call variants from aligned reads using the streaming variant caller.
    Variants {
        /// Reference genome (FASTA).
        #[arg(long)]
        reference: PathBuf,
        /// Alignments in SAM format (primary alignments only).
        #[arg(long)]
        alignments: PathBuf,
        /// Chromosome name (defaults to the first FASTA record if omitted).
        #[arg(long)]
        chrom: Option<String>,
        /// Starting offset (0-based) for the reference region.
        #[arg(long, default_value_t = 0)]
        region_start: u32,
        /// Minimum MAPQ required for a read to be considered.
        #[arg(long, default_value_t = 0)]
        mapq_threshold: u8,
        /// Optional VCF output path (stdout if omitted).
        #[arg(short, long)]
        output: Option<PathBuf>,
        /// Bases per block for streaming evaluation.
        #[arg(long, default_value_t = 1024)]
        block_size: usize,
        /// Minimum quality threshold for reporting variants.
        #[arg(long, default_value_t = 10.0)]
        quality_threshold: f32,
    },
}

#[derive(Copy, Clone, Debug, ValueEnum, Eq, PartialEq)]
enum OutputFormat {
    Sam,
    Bam,
}

struct FastaRecord {
    name: String,
    sequence: Vec<u8>,
}

struct FastqRecord {
    name: String,
    sequence: Vec<u8>,
    qualities: Vec<u8>,
}

struct AlignmentCandidate {
    position: usize,
    mismatches: usize,
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Align {
            reference,
            reads,
            max_mismatches,
            reference_offset,
            format,
            output,
        } => run_align(
            reference,
            reads,
            max_mismatches,
            reference_offset,
            format,
            output,
        )?,
        Commands::Variants {
            reference,
            alignments,
            chrom,
            region_start,
            mapq_threshold,
            output,
            block_size,
            quality_threshold,
        } => run_variants(
            reference,
            alignments,
            chrom,
            region_start,
            mapq_threshold,
            output,
            block_size,
            quality_threshold,
        )?,
    }

    Ok(())
}

fn run_align(
    reference_path: PathBuf,
    reads_path: PathBuf,
    max_mismatches: usize,
    reference_offset: u32,
    format: OutputFormat,
    output: Option<PathBuf>,
) -> Result<()> {
    let fasta = read_fasta(&reference_path)
        .with_context(|| format!("failed to read reference from {}", reference_path.display()))?;
    let records = read_fastq(&reads_path)
        .with_context(|| format!("failed to read reads from {}", reads_path.display()))?;
    let mut aligner = BWTAligner::new(&fasta.sequence)
        .context("failed to initialize the FM-index aligner for the reference")?;
    let alignments =
        align_reads(&mut aligner, &records, max_mismatches).context("aligning reads failed")?;

    match format {
        OutputFormat::Sam => {
            if let Some(path) = output {
                let file = File::create(&path)
                    .with_context(|| format!("failed to create SAM file {}", path.display()))?;
                let mut writer = io::BufWriter::new(file);
                write_sam_alignments(
                    &mut writer,
                    &fasta.name,
                    fasta.sequence.len(),
                    reference_offset,
                    &records,
                    &alignments,
                )?;
            } else {
                let stdout = io::stdout();
                let mut handle = stdout.lock();
                write_sam_alignments(
                    &mut handle,
                    &fasta.name,
                    fasta.sequence.len(),
                    reference_offset,
                    &records,
                    &alignments,
                )?;
            }
        }
        OutputFormat::Bam => {
            let path = output.ok_or_else(|| {
                anyhow!("--output <FILE> must be provided when writing BAM output")
            })?;
            let mut writer = create_bam_writer(&path, &fasta.name, fasta.sequence.len())
                .with_context(|| format!("failed to create BAM writer for {}", path.display()))?;
            write_bam_alignments(&mut writer, reference_offset, &records, &alignments)?;
        }
    }

    Ok(())
}

fn align_reads(
    aligner: &mut BWTAligner,
    records: &[FastqRecord],
    max_mismatches: usize,
) -> Result<Vec<Option<AlignmentCandidate>>> {
    let mut results = Vec::with_capacity(records.len());
    for record in records {
        let alignment = aligner
            .align_read(&record.sequence)
            .with_context(|| format!("failed to align read {}", record.name))?;

        if alignment.has_candidates() && alignment.mismatches <= max_mismatches {
            let position = aligner.fm_index().sa_at(alignment.interval.lower as usize);
            results.push(Some(AlignmentCandidate {
                position,
                mismatches: alignment.mismatches,
            }));
        } else {
            results.push(None);
        }
    }

    Ok(results)
}

fn run_variants(
    reference_path: PathBuf,
    alignments_path: PathBuf,
    chrom: Option<String>,
    region_start: u32,
    mapq_threshold: u8,
    output: Option<PathBuf>,
    block_size: usize,
    quality_threshold: f32,
) -> Result<()> {
    let fasta = read_fasta(&reference_path)
        .with_context(|| format!("failed to read reference from {}", reference_path.display()))?;
    let chrom_name = chrom.unwrap_or_else(|| fasta.name.clone());
    let chrom_arc: Arc<str> = chrom_name.into();
    let reference = Arc::from(fasta.sequence.into_boxed_slice());
    let reads = read_alignment_file(&alignments_path, Some(&chrom_arc)).with_context(|| {
        format!(
            "failed to read SAM alignments from {}",
            alignments_path.display()
        )
    })?;

    let mut caller = StreamingVariantCaller::new(
        Arc::clone(&chrom_arc),
        Arc::clone(&reference),
        region_start,
        block_size,
        quality_threshold,
        1e-6,
    )
    .context("failed to initialize variant caller")?;

    let filtered_reads: Vec<AlignedRead> = reads
        .into_iter()
        .filter(|read| read.mapq() >= mapq_threshold)
        .collect();

    let variants = caller
        .call_variants(filtered_reads)
        .context("variant calling failed")?;

    if let Some(path) = output {
        let file = File::create(&path)
            .with_context(|| format!("failed to create VCF file {}", path.display()))?;
        let mut writer = io::BufWriter::new(file);
        write_vcf(&mut writer, &variants)?;
    } else {
        let stdout = io::stdout();
        let mut handle = stdout.lock();
        write_vcf(&mut handle, &variants)?;
    }

    Ok(())
}

fn read_fasta(path: &PathBuf) -> Result<FastaRecord> {
    let contents = std::fs::read_to_string(path)
        .with_context(|| format!("failed to open {}", path.display()))?;

    let mut name = None;
    let mut sequence = String::new();

    for line in contents.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }

        if trimmed.starts_with('>') {
            if name.is_some() {
                eprintln!(
                    "warning: only the first FASTA record is currently used; ignoring '{}'",
                    trimmed
                );
                break;
            }
            let header = trimmed.trim_start_matches('>').trim();
            let primary_name = header
                .split_whitespace()
                .next()
                .ok_or_else(|| anyhow!("FASTA header '{}' has no name", header))?
                .to_string();
            name = Some(primary_name);
        } else {
            sequence.push_str(trimmed);
        }
    }

    let name = name.ok_or_else(|| anyhow!("FASTA file {} is missing a header", path.display()))?;
    if sequence.is_empty() {
        bail!("FASTA record {} has no sequence data", name);
    }

    Ok(FastaRecord {
        name,
        sequence: sequence.to_ascii_uppercase().into_bytes(),
    })
}

fn read_fastq(path: &PathBuf) -> Result<Vec<FastqRecord>> {
    let file = File::open(path)
        .with_context(|| format!("failed to open FASTQ file {}", path.display()))?;
    let mut reader = BufReader::new(file).lines();
    let mut records = Vec::new();

    while let Some(header) = reader.next() {
        let header = header?;
        if header.trim().is_empty() {
            continue;
        }
        if !header.starts_with('@') {
            bail!(
                "expected FASTQ header starting with '@', got '{}' in {}",
                header,
                path.display()
            );
        }
        let header_rest = header[1..].trim();
        let name = header_rest
            .split_whitespace()
            .next()
            .ok_or_else(|| anyhow!("FASTQ header '{}' has no read name", header))?
            .to_string();

        let seq_line = reader
            .next()
            .ok_or_else(|| anyhow!("unexpected end of FASTQ file while reading {}", name))??;
        let plus_line = reader
            .next()
            .ok_or_else(|| anyhow!("unexpected end of FASTQ file after sequence {}", name))??;
        if !plus_line.trim().starts_with('+') {
            bail!(
                "expected '+' separator after sequence for read {}, found '{}'",
                name,
                plus_line
            );
        }
        let qual_line = reader
            .next()
            .ok_or_else(|| anyhow!("unexpected end of FASTQ file after '+' for {}", name))??;

        let sequence = seq_line.trim().to_ascii_uppercase().into_bytes();
        let qualities = qual_line.trim().as_bytes().to_vec();

        if sequence.len() != qualities.len() {
            bail!(
                "sequence/quality length mismatch for read {} ({} vs {})",
                name,
                sequence.len(),
                qualities.len()
            );
        }

        records.push(FastqRecord {
            name,
            sequence,
            qualities,
        });
    }

    Ok(records)
}

fn mapq_from_mismatches(mismatches: usize) -> u8 {
    match mismatches {
        0 => 60,
        1 => 40,
        2 => 20,
        _ => 0,
    }
}

fn write_sam_alignments<W: Write>(
    writer: &mut W,
    reference_name: &str,
    reference_len: usize,
    reference_offset: u32,
    reads: &[FastqRecord],
    alignments: &[Option<AlignmentCandidate>],
) -> Result<()> {
    writeln!(writer, "@HD\tVN:1.6\tSO:unknown")?;
    writeln!(writer, "@SQ\tSN:{reference_name}\tLN:{reference_len}")?;

    if reads.len() != alignments.len() {
        bail!("internal error: read-alignment count mismatch");
    }

    for (record, alignment) in reads.iter().zip(alignments.iter()) {
        let seq_str = String::from_utf8(record.sequence.clone())
            .map_err(|_| anyhow!("FASTQ sequence for {} is not valid ASCII", record.name))?;
        let qual_str = String::from_utf8(record.qualities.clone())
            .map_err(|_| anyhow!("FASTQ qualities for {} are not valid ASCII", record.name))?;

        if let Some(hit) = alignment {
            let mapq = mapq_from_mismatches(hit.mismatches);
            let pos = hit.position + reference_offset as usize + 1;
            writeln!(
                writer,
                "{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{seq}\t{qual}\tNM:i:{nm}",
                qname = record.name,
                flag = 0,
                rname = reference_name,
                pos = pos,
                cigar = format!("{}M", record.sequence.len()),
                seq = seq_str,
                qual = qual_str,
                nm = hit.mismatches
            )?;
        } else {
            writeln!(
                writer,
                "{qname}\t{flag}\t*\t0\t0\t*\t*\t0\t0\t{seq}\t{qual}",
                qname = record.name,
                flag = 4,
                seq = seq_str,
                qual = qual_str
            )?;
        }
    }

    writer.flush()?;
    Ok(())
}

fn write_bam_alignments(
    writer: &mut bam::Writer,
    reference_offset: u32,
    reads: &[FastqRecord],
    alignments: &[Option<AlignmentCandidate>],
) -> Result<()> {
    if reads.len() != alignments.len() {
        bail!("internal error: read-alignment count mismatch");
    }

    for (record, alignment) in reads.iter().zip(alignments.iter()) {
        if let Some(hit) = alignment {
            let mut bam_record = Record::new();
            let cigar = CigarString::from(vec![BamCigar::Match(record.sequence.len() as u32)]);
            bam_record.set(
                record.name.as_bytes(),
                Some(&cigar),
                &record.sequence,
                &record.qualities,
            );
            bam_record.set_tid(0);
            bam_record.set_pos((hit.position + reference_offset as usize) as i64);
            bam_record.set_flags(0);
            bam_record.set_mapq(mapq_from_mismatches(hit.mismatches));
            bam_record.set_mtid(-1);
            bam_record.set_mpos(-1);
            bam_record.set_insert_size(0);
            bam_record.push_aux(b"NM", Aux::I32(hit.mismatches as i32))?;
            writer.write(&bam_record)?;
        } else {
            let mut bam_record = Record::new();
            bam_record.set(
                record.name.as_bytes(),
                None,
                &record.sequence,
                &record.qualities,
            );
            bam_record.set_tid(-1);
            bam_record.set_pos(-1);
            bam_record.set_flags(0x4);
            bam_record.set_mapq(0);
            bam_record.set_mtid(-1);
            bam_record.set_mpos(-1);
            bam_record.set_insert_size(0);
            writer.write(&bam_record)?;
        }
    }

    Ok(())
}

fn read_alignment_file(
    path: &PathBuf,
    target_chrom: Option<&Arc<str>>,
) -> Result<Vec<AlignedRead>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut reads = Vec::new();

    for (line_no, line) in reader.lines().enumerate() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('@') {
            continue;
        }

        let fields: Vec<&str> = trimmed.split('\t').collect();
        if fields.len() < 11 {
            bail!(
                "SAM record on line {} has {} fields (expected ≥ 11)",
                line_no + 1,
                fields.len()
            );
        }

        let flag: u16 = fields[1]
            .parse()
            .with_context(|| format!("invalid FLAG '{}' on line {}", fields[1], line_no + 1))?;
        if flag & 0x4 != 0 {
            continue;
        }

        let rname = fields[2];
        if rname == "*" {
            continue;
        }
        if let Some(target) = target_chrom {
            if rname != target.as_ref() {
                continue;
            }
        }

        let pos: u32 = fields[3]
            .parse()
            .with_context(|| format!("invalid POS '{}' on line {}", fields[3], line_no + 1))?;
        let mapq: u8 = fields[4]
            .parse()
            .with_context(|| format!("invalid MAPQ '{}' on line {}", fields[4], line_no + 1))?;
        let cigar = parse_cigar(fields[5], fields[9].len())
            .with_context(|| format!("invalid CIGAR '{}' on line {}", fields[5], line_no + 1))?;

        let sequence = fields[9].to_ascii_uppercase().into_bytes();
        let qual_field = fields[10].as_bytes();
        if qual_field.len() != sequence.len() {
            bail!(
                "sequence/quality length mismatch on line {} ({} vs {})",
                line_no + 1,
                sequence.len(),
                qual_field.len()
            );
        }
        let qualities: Vec<u8> = qual_field.iter().map(|q| q.saturating_sub(33)).collect();

        let is_reverse = flag & 0x10 != 0;

        reads.push(AlignedRead::new(
            rname.to_string(),
            pos.saturating_sub(1),
            mapq,
            cigar,
            sequence,
            qualities,
            is_reverse,
        ));
    }

    Ok(reads)
}

fn parse_cigar(cigar: &str, read_len: usize) -> Result<Vec<CigarOp>> {
    if cigar == "*" {
        bail!("variant calling requires mapped reads (CIGAR cannot be '*')");
    }
    if !cigar.ends_with('M') {
        bail!("only match operations (M) are supported in CIGAR strings");
    }
    let len: u32 = cigar[..cigar.len() - 1]
        .parse()
        .with_context(|| format!("invalid CIGAR length in '{}'", cigar))?;
    if len as usize != read_len {
        bail!(
            "CIGAR length {} does not match read length {}",
            len,
            read_len
        );
    }
    Ok(vec![CigarOp::new(CigarOpKind::Match, len)])
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::env;
    use std::fs;
    use std::path::PathBuf;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn temp_file_path(suffix: &str) -> PathBuf {
        let timestamp = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("time went backwards")
            .as_nanos();
        let mut path = env::temp_dir();
        path.push(format!("rosalind-test-{suffix}-{timestamp}.tmp"));
        path
    }

    fn write_temp_file(contents: &str, suffix: &str) -> PathBuf {
        let path = temp_file_path(suffix);
        fs::write(&path, contents).expect("failed to write temp file");
        path
    }

    #[test]
    fn fasta_parser_extracts_primary_name() {
        let path = write_temp_file(">chr1 CP068277.2 description\nACGT\n", "fasta");
        let record = read_fasta(&path).expect("FASTA read should succeed");
        assert_eq!(record.name, "chr1");
        let _ = fs::remove_file(path);
    }

    #[test]
    fn fastq_parser_trims_after_space() {
        let contents = "@read1 1:N:0:CG\nAC\n+\n!!\n";
        let path = write_temp_file(contents, "fastq");
        let records = read_fastq(&path).expect("FASTQ read should succeed");
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].name, "read1");
        let _ = fs::remove_file(path);
    }

    #[test]
    fn align_reads_with_fm_index() {
        let reference = b"ACGTACGT";
        let mut aligner =
            BWTAligner::new(reference).expect("should build FM-index aligner for test");
        let records = vec![FastqRecord {
            name: "read1".to_string(),
            sequence: b"ACGT".to_vec(),
            qualities: b"IIII".to_vec(),
        }];
        let alignments = align_reads(&mut aligner, &records, 2).expect("alignment should run");
        assert!(alignments[0].is_some());
    }
}
