use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::PathBuf;
use std::sync::Arc;
use std::time::Instant;

use anyhow::{anyhow, bail, Context, Result};
use clap::{Parser, Subcommand, ValueEnum};
use rosalind::genomics::{
    compare_callsets, create_bam_writer, read_vcf_variants, sort_bam_deterministic, AlignedRead,
    BWTAligner, BedIndex, CigarOp, CigarOpKind,
};
use rosalind::util::rss::peak_rss_bytes;
use rust_htslib::bam::Read as BamRead;
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
        /// Reads file in FASTQ format (single-end).
        #[arg(long)]
        reads: Option<PathBuf>,
        /// Reads R1 file in FASTQ format (paired-end).
        #[arg(long)]
        reads_r1: Option<PathBuf>,
        /// Reads R2 file in FASTQ format (paired-end).
        #[arg(long)]
        reads_r2: Option<PathBuf>,
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
    /// Deterministically coordinate-sort a BAM file using bounded memory.
    Sort {
        /// Input BAM path.
        #[arg(long)]
        input: PathBuf,
        /// Output BAM path.
        #[arg(short, long)]
        output: PathBuf,
        /// Memory budget (MiB) for in-memory sorting chunks.
        #[arg(long, default_value_t = 1024)]
        memory_mb: usize,
    },
    /// End-to-end tumor/normal somatic calling (align + sort + call).
    Somatic {
        /// Reference FASTA (single-contig for now).
        #[arg(long)]
        reference: PathBuf,
        /// Tumor FASTQ reads (single-end).
        #[arg(long)]
        tumor: Option<PathBuf>,
        /// Tumor FASTQ reads R1 (paired-end).
        #[arg(long)]
        tumor_r1: Option<PathBuf>,
        /// Tumor FASTQ reads R2 (paired-end).
        #[arg(long)]
        tumor_r2: Option<PathBuf>,
        /// Normal FASTQ reads (single-end).
        #[arg(long)]
        normal: Option<PathBuf>,
        /// Normal FASTQ reads R1 (paired-end).
        #[arg(long)]
        normal_r1: Option<PathBuf>,
        /// Normal FASTQ reads R2 (paired-end).
        #[arg(long)]
        normal_r2: Option<PathBuf>,
        /// Output VCF path.
        #[arg(short, long)]
        output: PathBuf,
        /// Working directory for intermediate BAMs.
        #[arg(long)]
        workdir: Option<PathBuf>,
        /// Sorting memory budget (MiB).
        #[arg(long, default_value_t = 1024)]
        memory_mb: usize,
    },
    /// Compare a called VCF against a truth VCF (optionally masked by BED).
    EvalSomatic {
        /// Reference FASTA (single contig slice used by the VCFs).
        #[arg(long)]
        reference: PathBuf,
        /// Called VCF path.
        #[arg(long)]
        calls: PathBuf,
        /// Truth VCF path.
        #[arg(long)]
        truth: PathBuf,
        /// Optional BED mask (0-based half-open).
        #[arg(long)]
        regions: Option<PathBuf>,
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

#[derive(Debug, Clone)]
struct FastqRecord {
    name: String,
    sequence: Vec<u8>,
    qualities: Vec<u8>,
}

#[derive(Debug, Clone)]
struct FastqPair {
    name: String,
    r1: FastqRecord,
    r2: FastqRecord,
}

#[derive(Debug)]
enum ResolvedReads {
    Single(Vec<FastqRecord>),
    Paired(Vec<FastqPair>),
}

struct AlignmentCandidate {
    position: usize,
    mismatches: usize,
    mapq: u8,
    is_reverse: bool,
    cigar: Vec<CigarOp>,
    as_score: i32,
    md: String,
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Align {
            reference,
            reads,
            reads_r1,
            reads_r2,
            max_mismatches,
            reference_offset,
            format,
            output,
        } => run_align(
            reference,
            reads,
            reads_r1,
            reads_r2,
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
        Commands::Sort {
            input,
            output,
            memory_mb,
        } => {
            let bytes = memory_mb.saturating_mul(1024 * 1024).max(1024 * 1024);
            sort_bam_deterministic(input, output, bytes).context("sorting BAM failed")?;
        }
        Commands::Somatic {
            reference,
            tumor,
            tumor_r1,
            tumor_r2,
            normal,
            normal_r1,
            normal_r2,
            output,
            workdir,
            memory_mb,
        } => {
            run_somatic(
                reference, tumor, tumor_r1, tumor_r2, normal, normal_r1, normal_r2, output,
                workdir, memory_mb,
            )?;
        }
        Commands::EvalSomatic {
            reference,
            calls,
            truth,
            regions,
        } => {
            run_eval_somatic(reference, calls, truth, regions)?;
        }
    }

    Ok(())
}

fn run_eval_somatic(
    reference_path: PathBuf,
    calls_path: PathBuf,
    truth_path: PathBuf,
    regions_path: Option<PathBuf>,
) -> Result<()> {
    let fasta = read_fasta(&reference_path)
        .with_context(|| format!("failed to read reference from {}", reference_path.display()))?;
    let reference = fasta.sequence;

    let calls_txt = std::fs::read_to_string(&calls_path)
        .with_context(|| format!("failed to read calls VCF {}", calls_path.display()))?;
    let truth_txt = std::fs::read_to_string(&truth_path)
        .with_context(|| format!("failed to read truth VCF {}", truth_path.display()))?;

    let calls = read_vcf_variants(&calls_txt)
        .with_context(|| format!("failed to parse calls VCF {}", calls_path.display()))?;
    let truth = read_vcf_variants(&truth_txt)
        .with_context(|| format!("failed to parse truth VCF {}", truth_path.display()))?;

    let bed = if let Some(path) = regions_path {
        let bed_txt = std::fs::read_to_string(&path)
            .with_context(|| format!("failed to read BED {}", path.display()))?;
        Some(BedIndex::from_str(&bed_txt).with_context(|| "failed to parse BED")?)
    } else {
        None
    };

    let report = compare_callsets(&reference, &calls, &truth, bed.as_ref())?;
    println!("truth_total={}", report.total_truth);
    println!("calls_total={}", report.total_calls);
    println!("tp={}", report.true_positive);
    println!("fp={}", report.false_positive);
    println!("fn={}", report.false_negative);
    println!("precision={:.6}", report.precision());
    println!("recall={:.6}", report.recall());
    for (ty, (tp, fp, fn_)) in report.by_type.iter() {
        println!("type={:?} tp={} fp={} fn={}", ty, tp, fp, fn_);
    }
    Ok(())
}

fn run_somatic(
    reference_path: PathBuf,
    tumor_fastq: Option<PathBuf>,
    tumor_r1: Option<PathBuf>,
    tumor_r2: Option<PathBuf>,
    normal_fastq: Option<PathBuf>,
    normal_r1: Option<PathBuf>,
    normal_r2: Option<PathBuf>,
    output_vcf: PathBuf,
    workdir: Option<PathBuf>,
    memory_mb: usize,
) -> Result<()> {
    let start_total = Instant::now();
    let fasta = read_fasta(&reference_path)
        .with_context(|| format!("failed to read reference from {}", reference_path.display()))?;
    let reference: Arc<[u8]> = Arc::from(fasta.sequence.clone().into_boxed_slice());

    let workdir = workdir.unwrap_or_else(|| {
        output_vcf
            .parent()
            .map(|p| p.to_path_buf())
            .unwrap_or_else(|| PathBuf::from("."))
    });
    std::fs::create_dir_all(&workdir)
        .with_context(|| format!("failed to create workdir {}", workdir.display()))?;

    let tumor_bam = workdir.join("tumor.bam");
    let tumor_sorted = workdir.join("tumor.sorted.bam");
    let normal_bam = workdir.join("normal.bam");
    let normal_sorted = workdir.join("normal.sorted.bam");

    // Align tumor reads.
    let start_align_tumor = Instant::now();
    {
        let mut writer =
            create_bam_writer(&tumor_bam, fasta.name.as_str(), fasta.sequence.len())
                .with_context(|| format!("failed to create tumor BAM {}", tumor_bam.display()))?;
        let mut aligner = BWTAligner::new(&reference)?;

        match resolve_reads(
            tumor_fastq.clone(),
            tumor_r1.clone(),
            tumor_r2.clone(),
            "tumor",
        )? {
            ResolvedReads::Single(reads) => {
                let alignments = align_reads(&mut aligner, &reads, 10)?;
                write_bam_alignments(&mut writer, 0, &reads, &alignments)?;
            }
            ResolvedReads::Paired(pairs) => {
                let alignments = align_pairs(&mut aligner, &pairs, 10)?;
                write_bam_alignments_paired(&mut writer, 0, &pairs, &alignments)?;
            }
        }
    }
    let dur_align_tumor = start_align_tumor.elapsed();

    // Align normal reads.
    let start_align_normal = Instant::now();
    {
        let mut writer = create_bam_writer(&normal_bam, fasta.name.as_str(), fasta.sequence.len())
            .with_context(|| format!("failed to create normal BAM {}", normal_bam.display()))?;
        let mut aligner = BWTAligner::new(&reference)?;

        match resolve_reads(
            normal_fastq.clone(),
            normal_r1.clone(),
            normal_r2.clone(),
            "normal",
        )? {
            ResolvedReads::Single(reads) => {
                let alignments = align_reads(&mut aligner, &reads, 10)?;
                write_bam_alignments(&mut writer, 0, &reads, &alignments)?;
            }
            ResolvedReads::Paired(pairs) => {
                let alignments = align_pairs(&mut aligner, &pairs, 10)?;
                write_bam_alignments_paired(&mut writer, 0, &pairs, &alignments)?;
            }
        }
    }
    let dur_align_normal = start_align_normal.elapsed();

    // Sort BAMs deterministically.
    let start_sort = Instant::now();
    let bytes = memory_mb.saturating_mul(1024 * 1024).max(1024 * 1024);
    sort_bam_deterministic(&tumor_bam, &tumor_sorted, bytes)?;
    sort_bam_deterministic(&normal_bam, &normal_sorted, bytes)?;
    let dur_sort = start_sort.elapsed();

    // Call somatic SNVs on the new engine (indels deferred to Phase C).
    use rosalind::call::{call_somatic_region, SomaticParams};
    use rosalind::core::ContigSet;
    use rosalind::io::bam::BamSource;
    use rosalind::io::vcf::write_somatic_vcf;
    use rosalind::pileup::PileupParams;
    use rosalind::provenance::{blake3_file, write_manifest, FileHash, RunManifest};

    let start_call = Instant::now();
    let mut contigs = ContigSet::new();
    let contig_id = contigs.push(fasta.name.clone(), reference.len() as u32);
    let region = 0..(reference.len() as u32);

    let tumor_src = BamSource::new(&tumor_sorted, &contigs)
        .map_err(|e| anyhow!("failed to read tumor BAM {}: {e}", tumor_sorted.display()))?;
    let normal_src = BamSource::new(&normal_sorted, &contigs)
        .map_err(|e| anyhow!("failed to read normal BAM {}: {e}", normal_sorted.display()))?;
    let calls = call_somatic_region(
        tumor_src,
        normal_src,
        Arc::clone(&reference),
        contig_id,
        region,
        PileupParams::default(),
        &SomaticParams::default(),
    )
    .map_err(|e| anyhow!("somatic calling failed: {e}"))?;
    let dur_call = start_call.elapsed();

    // Write spec-valid somatic VCF (TUMOR/NORMAL).
    {
        let file = File::create(&output_vcf)
            .with_context(|| format!("failed to create somatic VCF {}", output_vcf.display()))?;
        let mut writer = io::BufWriter::new(file);
        write_somatic_vcf(&mut writer, &contigs, &calls)?;
        writer.flush()?;
    }

    // Reproducibility receipt (BLAKE3, canonical JSON).
    let tumor_inputs: Vec<PathBuf> = match (&tumor_fastq, &tumor_r1, &tumor_r2) {
        (Some(p), None, None) => vec![p.clone()],
        (None, Some(r1), Some(r2)) => vec![r1.clone(), r2.clone()],
        _ => Vec::new(),
    };
    let normal_inputs: Vec<PathBuf> = match (&normal_fastq, &normal_r1, &normal_r2) {
        (Some(p), None, None) => vec![p.clone()],
        (None, Some(r1), Some(r2)) => vec![r1.clone(), r2.clone()],
        _ => Vec::new(),
    };
    let mut manifest = RunManifest::new("somatic");
    manifest.inputs.push(FileHash {
        path: reference_path.display().to_string(),
        blake3: blake3_file(&reference_path)?,
    });
    for p in tumor_inputs.iter().chain(normal_inputs.iter()) {
        if p.exists() {
            manifest.inputs.push(FileHash {
                path: p.display().to_string(),
                blake3: blake3_file(p)?,
            });
        }
    }
    manifest.outputs.push(FileHash {
        path: output_vcf.display().to_string(),
        blake3: blake3_file(&output_vcf)?,
    });
    manifest
        .params
        .insert("somatic_snv_only".to_string(), "true".to_string());
    let manifest_path = write_manifest(&output_vcf, &manifest)?;
    eprintln!("wrote reproducibility receipt: {}", manifest_path.display());

    // Performance/RSS report (kept separate from determinism checks).
    let perf_path = workdir.join("somatic.perf.txt");
    let mut f = io::BufWriter::new(
        File::create(&perf_path)
            .with_context(|| format!("failed to create perf report {}", perf_path.display()))?,
    );
    writeln!(f, "align_tumor_ms={}", dur_align_tumor.as_millis())?;
    writeln!(f, "align_normal_ms={}", dur_align_normal.as_millis())?;
    writeln!(f, "sort_ms={}", dur_sort.as_millis())?;
    writeln!(f, "call_ms={}", dur_call.as_millis())?;
    writeln!(f, "total_ms={}", start_total.elapsed().as_millis())?;
    writeln!(f, "peak_rss_bytes={}", peak_rss_bytes())?;
    f.flush()?;

    Ok(())
}

fn run_align(
    reference_path: PathBuf,
    reads_path: Option<PathBuf>,
    reads_r1: Option<PathBuf>,
    reads_r2: Option<PathBuf>,
    max_mismatches: usize,
    reference_offset: u32,
    format: OutputFormat,
    output: Option<PathBuf>,
) -> Result<()> {
    let fasta = read_fasta(&reference_path)
        .with_context(|| format!("failed to read reference from {}", reference_path.display()))?;
    let reads = resolve_reads(reads_path, reads_r1, reads_r2, "reads")
        .context("failed to resolve reads inputs")?;
    let mut aligner = BWTAligner::new(&fasta.sequence)
        .context("failed to initialize the FM-index aligner for the reference")?;

    match reads {
        ResolvedReads::Single(records) => {
            let alignments = align_reads(&mut aligner, &records, max_mismatches)
                .context("aligning reads failed")?;
            match format {
                OutputFormat::Sam => {
                    if let Some(path) = output {
                        let file = File::create(&path).with_context(|| {
                            format!("failed to create SAM file {}", path.display())
                        })?;
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
                        .with_context(|| {
                            format!("failed to create BAM writer for {}", path.display())
                        })?;
                    write_bam_alignments(&mut writer, reference_offset, &records, &alignments)?;
                }
            }
        }
        ResolvedReads::Paired(pairs) => {
            let alignments = align_pairs(&mut aligner, &pairs, max_mismatches)
                .context("aligning reads failed")?;
            match format {
                OutputFormat::Sam => {
                    if let Some(path) = output {
                        let file = File::create(&path).with_context(|| {
                            format!("failed to create SAM file {}", path.display())
                        })?;
                        let mut writer = io::BufWriter::new(file);
                        write_sam_alignments_paired(
                            &mut writer,
                            &fasta.name,
                            fasta.sequence.len(),
                            reference_offset,
                            &pairs,
                            &alignments,
                        )?;
                    } else {
                        let stdout = io::stdout();
                        let mut handle = stdout.lock();
                        write_sam_alignments_paired(
                            &mut handle,
                            &fasta.name,
                            fasta.sequence.len(),
                            reference_offset,
                            &pairs,
                            &alignments,
                        )?;
                    }
                }
                OutputFormat::Bam => {
                    let path = output.ok_or_else(|| {
                        anyhow!("--output <FILE> must be provided when writing BAM output")
                    })?;
                    let mut writer = create_bam_writer(&path, &fasta.name, fasta.sequence.len())
                        .with_context(|| {
                            format!("failed to create BAM writer for {}", path.display())
                        })?;
                    write_bam_alignments_paired(
                        &mut writer,
                        reference_offset,
                        &pairs,
                        &alignments,
                    )?;
                }
            }
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
            let position = alignment
                .primary_position
                .map(|p| p as usize)
                .or_else(|| {
                    if alignment.interval.is_empty() {
                        None
                    } else {
                        Some(aligner.fm_index().sa_at(alignment.interval.lower as usize))
                    }
                })
                .ok_or_else(|| {
                    anyhow!("alignment reported candidates but no position could be located")
                })?;
            results.push(Some(AlignmentCandidate {
                position,
                mismatches: alignment.mismatches,
                mapq: alignment.mapq,
                is_reverse: alignment.is_reverse,
                cigar: alignment.cigar.clone(),
                as_score: alignment.as_score,
                md: alignment.md.clone(),
            }));
        } else {
            results.push(None);
        }
    }

    Ok(results)
}

type PairAlignments = (Option<AlignmentCandidate>, Option<AlignmentCandidate>);

fn align_pairs(
    aligner: &mut BWTAligner,
    pairs: &[FastqPair],
    max_mismatches: usize,
) -> Result<Vec<PairAlignments>> {
    let mut out = Vec::with_capacity(pairs.len());
    for pair in pairs {
        let a1 = aligner
            .align_read(&pair.r1.sequence)
            .with_context(|| format!("failed to align read {}", pair.name))?;
        let a2 = aligner
            .align_read(&pair.r2.sequence)
            .with_context(|| format!("failed to align read {}", pair.name))?;

        let c1 = if a1.has_candidates() && a1.mismatches <= max_mismatches {
            Some(AlignmentCandidate {
                position: a1
                    .primary_position
                    .map(|p| p as usize)
                    .or_else(|| {
                        if a1.interval.is_empty() {
                            None
                        } else {
                            Some(aligner.fm_index().sa_at(a1.interval.lower as usize))
                        }
                    })
                    .ok_or_else(|| {
                        anyhow!("alignment reported candidates but no position could be located")
                    })?,
                mismatches: a1.mismatches,
                mapq: a1.mapq,
                is_reverse: a1.is_reverse,
                cigar: a1.cigar.clone(),
                as_score: a1.as_score,
                md: a1.md.clone(),
            })
        } else {
            None
        };

        let c2 = if a2.has_candidates() && a2.mismatches <= max_mismatches {
            Some(AlignmentCandidate {
                position: a2
                    .primary_position
                    .map(|p| p as usize)
                    .or_else(|| {
                        if a2.interval.is_empty() {
                            None
                        } else {
                            Some(aligner.fm_index().sa_at(a2.interval.lower as usize))
                        }
                    })
                    .ok_or_else(|| {
                        anyhow!("alignment reported candidates but no position could be located")
                    })?,
                mismatches: a2.mismatches,
                mapq: a2.mapq,
                is_reverse: a2.is_reverse,
                cigar: a2.cigar.clone(),
                as_score: a2.as_score,
                md: a2.md.clone(),
            })
        } else {
            None
        };

        out.push((c1, c2));
    }
    Ok(out)
}

fn run_variants(
    reference_path: PathBuf,
    alignments_path: PathBuf,
    chrom: Option<String>,
    region_start: u32,
    mapq_threshold: u8,
    output: Option<PathBuf>,
    _block_size: usize,
    quality_threshold: f32,
) -> Result<()> {
    use rosalind::call::{call_germline_region, GermlineParams};
    use rosalind::core::ContigSet;
    use rosalind::io::bam::BamSource;
    use rosalind::io::vcf::{write_germline_vcf, GermlineRow};
    use rosalind::pileup::{PileupParams, SliceSource};
    use rosalind::provenance::{blake3_file, write_manifest, FileHash, RunManifest};

    let fasta = read_fasta(&reference_path)
        .with_context(|| format!("failed to read reference from {}", reference_path.display()))?;
    let chrom_name = chrom.unwrap_or_else(|| fasta.name.clone());
    let chrom_arc: Arc<str> = chrom_name.clone().into();
    let reference: Arc<[u8]> = Arc::from(fasta.sequence.into_boxed_slice());

    // Single-contig run: one contig spanning the reference window.
    let mut contigs = ContigSet::new();
    let contig_id = contigs.push(chrom_name.clone(), region_start + reference.len() as u32);

    let region = region_start..(region_start + reference.len() as u32);
    let pileup_params = PileupParams {
        min_mapq: mapq_threshold, // now actually honored by the engine
        ..PileupParams::default()
    };
    let germline_params = GermlineParams {
        min_qual: quality_threshold as f64,
        ..GermlineParams::default()
    };

    let is_bam = alignments_path
        .extension()
        .and_then(|ext| ext.to_str())
        .map(|ext| ext.eq_ignore_ascii_case("bam"))
        .unwrap_or(false);

    let sites = if is_bam {
        let source = BamSource::new(&alignments_path, &contigs)
            .map_err(|e| anyhow!("failed to read BAM {}: {e}", alignments_path.display()))?;
        call_germline_region(
            source,
            Arc::clone(&reference),
            contig_id,
            region,
            pileup_params,
            &germline_params,
        )
        .map_err(|e| anyhow!("variant calling failed (BAM): {e}"))?
    } else {
        let legacy =
            read_alignment_file(&alignments_path, Some(&chrom_arc)).with_context(|| {
                format!(
                    "failed to read SAM alignments from {}",
                    alignments_path.display()
                )
            })?;
        let core_reads: Vec<rosalind::core::AlignedRead> = legacy
            .iter()
            .map(|r| legacy_read_to_core(r, contig_id))
            .collect();
        let source = SliceSource::new(core_reads);
        call_germline_region(
            source,
            Arc::clone(&reference),
            contig_id,
            region,
            pileup_params,
            &germline_params,
        )
        .map_err(|e| anyhow!("variant calling failed (SAM): {e}"))?
    };

    let rows: Vec<GermlineRow> = sites
        .into_iter()
        .map(|(locus, ref_base, call)| GermlineRow {
            locus,
            ref_base,
            call,
        })
        .collect();

    match output {
        Some(path) => {
            let file = File::create(&path)
                .with_context(|| format!("failed to create VCF file {}", path.display()))?;
            let mut writer = io::BufWriter::new(file);
            write_germline_vcf(&mut writer, &contigs, &chrom_name, &rows)?;
            writer.flush()?;
            drop(writer);

            // Reproducibility receipt next to the VCF.
            let mut manifest = RunManifest::new("variants");
            manifest.inputs.push(FileHash {
                path: reference_path.display().to_string(),
                blake3: blake3_file(&reference_path)?,
            });
            manifest.inputs.push(FileHash {
                path: alignments_path.display().to_string(),
                blake3: blake3_file(&alignments_path)?,
            });
            manifest.outputs.push(FileHash {
                path: path.display().to_string(),
                blake3: blake3_file(&path)?,
            });
            manifest
                .params
                .insert("mapq_threshold".to_string(), mapq_threshold.to_string());
            manifest.params.insert(
                "min_qual".to_string(),
                (quality_threshold as f64).to_string(),
            );
            manifest
                .params
                .insert("region_start".to_string(), region_start.to_string());
            let manifest_path = write_manifest(&path, &manifest)?;
            eprintln!("wrote reproducibility receipt: {}", manifest_path.display());
        }
        None => {
            let stdout = io::stdout();
            let mut handle = stdout.lock();
            write_germline_vcf(&mut handle, &contigs, &chrom_name, &rows)?;
        }
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

fn normalize_read_name(raw: &str) -> String {
    let name = raw.trim();
    let name = name.split_whitespace().next().unwrap_or(name);
    if name.ends_with("/1") || name.ends_with("/2") {
        name[..name.len() - 2].to_string()
    } else {
        name.to_string()
    }
}

fn read_fastq_pairs(r1: &PathBuf, r2: &PathBuf) -> Result<Vec<FastqPair>> {
    let mut r1_records =
        read_fastq(r1).with_context(|| format!("failed to read FASTQ R1 from {}", r1.display()))?;
    let mut r2_records =
        read_fastq(r2).with_context(|| format!("failed to read FASTQ R2 from {}", r2.display()))?;

    if r1_records.len() != r2_records.len() {
        bail!(
            "paired FASTQ length mismatch: {} has {} records, {} has {} records",
            r1.display(),
            r1_records.len(),
            r2.display(),
            r2_records.len()
        );
    }

    let mut out = Vec::with_capacity(r1_records.len());
    for (mut rec1, mut rec2) in r1_records.drain(..).zip(r2_records.drain(..)) {
        let name1 = normalize_read_name(&rec1.name);
        let name2 = normalize_read_name(&rec2.name);
        if name1 != name2 {
            bail!("mate name mismatch: '{}' vs '{}'", name1, name2);
        }
        rec1.name = name1.clone();
        rec2.name = name1.clone();
        out.push(FastqPair {
            name: name1,
            r1: rec1,
            r2: rec2,
        });
    }
    Ok(out)
}

fn resolve_reads(
    single: Option<PathBuf>,
    r1: Option<PathBuf>,
    r2: Option<PathBuf>,
    label: &str,
) -> Result<ResolvedReads> {
    match (single, r1, r2) {
        (Some(path), None, None) => Ok(ResolvedReads::Single(read_fastq(&path)?)),
        (None, Some(r1), Some(r2)) => Ok(ResolvedReads::Paired(read_fastq_pairs(&r1, &r2)?)),
        (None, None, None) => bail!(
            "missing reads for {}: provide --{} <FASTQ> or --{}-r1/--{}-r2",
            label,
            label,
            label,
            label
        ),
        _ => bail!(
            "invalid reads args for {}: provide either single-end or paired-end (not both/incomplete)",
            label
        ),
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
    writeln!(
        writer,
        "@PG\tID:rosalind\tPN:rosalind\tVN:{}",
        env!("CARGO_PKG_VERSION")
    )?;

    if reads.len() != alignments.len() {
        bail!("internal error: read-alignment count mismatch");
    }

    for (record, alignment) in reads.iter().zip(alignments.iter()) {
        let seq_str = String::from_utf8(record.sequence.clone())
            .map_err(|_| anyhow!("FASTQ sequence for {} is not valid ASCII", record.name))?;
        let qual_str = String::from_utf8(record.qualities.clone())
            .map_err(|_| anyhow!("FASTQ qualities for {} are not valid ASCII", record.name))?;

        if let Some(hit) = alignment {
            let mapq = hit.mapq;
            let pos = hit.position + reference_offset as usize + 1;
            let flag = if hit.is_reverse { 0x10 } else { 0 };
            let cigar = if hit.cigar.is_empty() {
                format!("{}M", record.sequence.len())
            } else {
                format_cigar(&hit.cigar)
            };
            writeln!(
                writer,
                "{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{seq}\t{qual}\tNM:i:{nm}\tAS:i:{as_score}\tMD:Z:{md}",
                qname = record.name,
                flag = flag,
                rname = reference_name,
                pos = pos,
                cigar = cigar,
                seq = seq_str,
                qual = qual_str,
                nm = hit.mismatches,
                as_score = hit.as_score,
                md = hit.md
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

fn write_sam_alignments_paired<W: Write>(
    writer: &mut W,
    reference_name: &str,
    reference_len: usize,
    reference_offset: u32,
    reads: &[FastqPair],
    alignments: &[PairAlignments],
) -> Result<()> {
    writeln!(writer, "@HD\tVN:1.6\tSO:unknown")?;
    writeln!(writer, "@SQ\tSN:{reference_name}\tLN:{reference_len}")?;
    writeln!(
        writer,
        "@PG\tID:rosalind\tPN:rosalind\tVN:{}",
        env!("CARGO_PKG_VERSION")
    )?;

    if reads.len() != alignments.len() {
        bail!("internal error: read-alignment count mismatch");
    }

    for (pair, (a1, a2)) in reads.iter().zip(alignments.iter()) {
        write_one_sam_mate(
            writer,
            reference_name,
            reference_offset,
            &pair.r1,
            a1,
            a2,
            true,
        )?;
        write_one_sam_mate(
            writer,
            reference_name,
            reference_offset,
            &pair.r2,
            a2,
            a1,
            false,
        )?;
    }

    writer.flush()?;
    Ok(())
}

fn write_one_sam_mate<W: Write>(
    writer: &mut W,
    reference_name: &str,
    reference_offset: u32,
    record: &FastqRecord,
    this: &Option<AlignmentCandidate>,
    mate: &Option<AlignmentCandidate>,
    is_first: bool,
) -> Result<()> {
    let seq_str = String::from_utf8(record.sequence.clone())
        .map_err(|_| anyhow!("FASTQ sequence for {} is not valid ASCII", record.name))?;
    let qual_str = String::from_utf8(record.qualities.clone())
        .map_err(|_| anyhow!("FASTQ qualities for {} are not valid ASCII", record.name))?;

    let mut flag: u16 = 0x1; // paired
    flag |= if is_first { 0x40 } else { 0x80 };

    let (rname, pos, mapq, cigar) = if let Some(hit) = this {
        if hit.is_reverse {
            flag |= 0x10;
        }
        let pos = hit.position + reference_offset as usize + 1;
        let cigar = if hit.cigar.is_empty() {
            format!("{}M", record.sequence.len())
        } else {
            format_cigar(&hit.cigar)
        };
        (reference_name.to_string(), pos, hit.mapq, cigar)
    } else {
        flag |= 0x4;
        ("*".to_string(), 0usize, 0u8, "*".to_string())
    };

    let (rnext, pnext, tlen, mate_reverse_flag) = if let Some(m) = mate {
        let pnext = m.position + reference_offset as usize + 1;
        let tlen = compute_tlen(this.as_ref(), mate.as_ref(), record.sequence.len());
        let mate_is_rev = m.is_reverse;
        ("=".to_string(), pnext, tlen, mate_is_rev)
    } else {
        flag |= 0x8;
        ("*".to_string(), 0usize, 0i64, false)
    };
    if mate_reverse_flag {
        flag |= 0x20;
    }

    let nm = this.as_ref().map(|h| h.mismatches).unwrap_or(0);
    let as_score = this.as_ref().map(|h| h.as_score).unwrap_or(0);
    let md = this.as_ref().map(|h| h.md.as_str()).unwrap_or("");

    writeln!(
        writer,
        "{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}\tNM:i:{nm}\tAS:i:{as_score}\tMD:Z:{md}",
        qname = record.name,
        flag = flag,
        rname = rname,
        pos = pos,
        mapq = mapq,
        cigar = cigar,
        rnext = rnext,
        pnext = pnext,
        tlen = tlen,
        seq = seq_str,
        qual = qual_str,
        nm = nm,
        as_score = as_score,
        md = md
    )?;
    Ok(())
}

fn compute_tlen(
    this: Option<&AlignmentCandidate>,
    mate: Option<&AlignmentCandidate>,
    read_len: usize,
) -> i64 {
    let (Some(a), Some(b)) = (this, mate) else {
        return 0;
    };
    let a_start = a.position as i64;
    let b_start = b.position as i64;
    let a_end = a_start + read_len as i64;
    let b_end = b_start + read_len as i64;
    let left = a_start.min(b_start);
    let right = a_end.max(b_end);
    let tlen = right - left;
    if a_start <= b_start {
        tlen
    } else {
        -tlen
    }
}

/// Convert FASTQ-ASCII quality bytes (Phred+33, e.g. `b'I'` = 73 for Q40) to
/// raw Phred values (e.g. 40) as required by the BAM QUAL field.
///
/// The FASTQ parser stores quality bytes verbatim from the quality line, so
/// in-memory `FastqRecord.qualities` are ASCII-encoded (Phred+33). BAM/htslib
/// expects raw Phred (0–93) in `Record::set`; this function decodes them.
fn fastq_quals_to_phred(ascii: &[u8]) -> Vec<u8> {
    ascii.iter().map(|q| q.saturating_sub(33)).collect()
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
        let phred_quals = fastq_quals_to_phred(&record.qualities);
        if let Some(hit) = alignment {
            let mut bam_record = Record::new();
            let cigar_ops = if hit.cigar.is_empty() {
                vec![BamCigar::Match(record.sequence.len() as u32)]
            } else {
                hit.cigar
                    .iter()
                    .map(|op| match op.kind {
                        CigarOpKind::Match => BamCigar::Match(op.len),
                        CigarOpKind::Insertion => BamCigar::Ins(op.len),
                        CigarOpKind::Deletion => BamCigar::Del(op.len),
                        CigarOpKind::SoftClip => BamCigar::SoftClip(op.len),
                        CigarOpKind::HardClip => BamCigar::HardClip(op.len),
                    })
                    .collect()
            };
            let cigar = CigarString::from(cigar_ops);
            bam_record.set(
                record.name.as_bytes(),
                Some(&cigar),
                &record.sequence,
                &phred_quals,
            );
            bam_record.set_tid(0);
            bam_record.set_pos((hit.position + reference_offset as usize) as i64);
            let flags = if hit.is_reverse { 0x10 } else { 0 };
            bam_record.set_flags(flags);
            bam_record.set_mapq(hit.mapq);
            bam_record.set_mtid(-1);
            bam_record.set_mpos(-1);
            bam_record.set_insert_size(0);
            bam_record.push_aux(b"NM", Aux::I32(hit.mismatches as i32))?;
            bam_record.push_aux(b"AS", Aux::I32(hit.as_score))?;
            if !hit.md.is_empty() {
                bam_record.push_aux(b"MD", Aux::String(hit.md.as_str()))?;
            }
            writer.write(&bam_record)?;
        } else {
            let mut bam_record = Record::new();
            bam_record.set(record.name.as_bytes(), None, &record.sequence, &phred_quals);
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

fn write_bam_alignments_paired(
    writer: &mut bam::Writer,
    reference_offset: u32,
    reads: &[FastqPair],
    alignments: &[PairAlignments],
) -> Result<()> {
    if reads.len() != alignments.len() {
        bail!("internal error: read-alignment count mismatch");
    }

    for (pair, (a1, a2)) in reads.iter().zip(alignments.iter()) {
        write_one_bam_mate(writer, reference_offset, &pair.r1, a1, a2, true)?;
        write_one_bam_mate(writer, reference_offset, &pair.r2, a2, a1, false)?;
    }
    Ok(())
}

fn write_one_bam_mate(
    writer: &mut bam::Writer,
    reference_offset: u32,
    record: &FastqRecord,
    this: &Option<AlignmentCandidate>,
    mate: &Option<AlignmentCandidate>,
    is_first: bool,
) -> Result<()> {
    let mut bam_record = Record::new();

    let mut flags: u16 = 0x1; // paired
    flags |= if is_first { 0x40 } else { 0x80 };

    let (tid, pos) = if let Some(hit) = this {
        if hit.is_reverse {
            flags |= 0x10;
        }
        (0i32, (hit.position + reference_offset as usize) as i64)
    } else {
        flags |= 0x4;
        (-1i32, -1i64)
    };

    if mate.is_none() {
        flags |= 0x8;
    } else if mate.as_ref().is_some_and(|m| m.is_reverse) {
        flags |= 0x20;
    }

    let cigar_ops = if let Some(hit) = this {
        if hit.cigar.is_empty() {
            vec![BamCigar::Match(record.sequence.len() as u32)]
        } else {
            hit.cigar
                .iter()
                .map(|op| match op.kind {
                    CigarOpKind::Match => BamCigar::Match(op.len),
                    CigarOpKind::Insertion => BamCigar::Ins(op.len),
                    CigarOpKind::Deletion => BamCigar::Del(op.len),
                    CigarOpKind::SoftClip => BamCigar::SoftClip(op.len),
                    CigarOpKind::HardClip => BamCigar::HardClip(op.len),
                })
                .collect()
        }
    } else {
        Vec::new()
    };
    let cigar = if cigar_ops.is_empty() {
        None
    } else {
        Some(CigarString::from(cigar_ops))
    };

    let phred_quals = fastq_quals_to_phred(&record.qualities);
    bam_record.set(
        record.name.as_bytes(),
        cigar.as_ref(),
        &record.sequence,
        &phred_quals,
    );
    bam_record.set_flags(flags);
    bam_record.set_tid(tid);
    bam_record.set_pos(pos);
    bam_record.set_mapq(this.as_ref().map(|h| h.mapq).unwrap_or(0));

    if let Some(m) = mate {
        bam_record.set_mtid(0);
        bam_record.set_mpos((m.position + reference_offset as usize) as i64);
        bam_record.set_insert_size(compute_tlen(
            this.as_ref(),
            mate.as_ref(),
            record.sequence.len(),
        ));
    } else {
        bam_record.set_mtid(-1);
        bam_record.set_mpos(-1);
        bam_record.set_insert_size(0);
    }

    if let Some(hit) = this {
        bam_record.push_aux(b"NM", Aux::I32(hit.mismatches as i32))?;
        bam_record.push_aux(b"AS", Aux::I32(hit.as_score))?;
        if !hit.md.is_empty() {
            bam_record.push_aux(b"MD", Aux::String(hit.md.as_str()))?;
        }
    }

    writer.write(&bam_record)?;
    Ok(())
}

fn format_cigar(cigar: &[CigarOp]) -> String {
    let mut out = String::new();
    for op in cigar {
        out.push_str(&op.len.to_string());
        let ch = match op.kind {
            CigarOpKind::Match => 'M',
            CigarOpKind::Insertion => 'I',
            CigarOpKind::Deletion => 'D',
            CigarOpKind::SoftClip => 'S',
            CigarOpKind::HardClip => 'H',
        };
        out.push(ch);
    }
    out
}

fn read_alignment_file(
    path: &PathBuf,
    target_chrom: Option<&Arc<str>>,
) -> Result<Vec<AlignedRead>> {
    if path
        .extension()
        .and_then(|ext| ext.to_str())
        .map(|ext| ext.eq_ignore_ascii_case("bam"))
        .unwrap_or(false)
    {
        return read_bam_alignment_file(path, target_chrom);
    }

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

/// Convert a legacy `genomics::AlignedRead` into a canonical `core::AlignedRead`
/// for the new calling vertical. The legacy type has no RefSkip/Pad CIGAR ops.
fn legacy_read_to_core(r: &AlignedRead, contig: u32) -> rosalind::core::AlignedRead {
    use rosalind::core::{CigarOp as CoreOp, CigarOpKind as CoreKind};
    let cigar = r
        .cigar
        .iter()
        .map(|op| {
            let kind = match op.kind {
                CigarOpKind::Match => CoreKind::Match,
                CigarOpKind::Insertion => CoreKind::Insertion,
                CigarOpKind::Deletion => CoreKind::Deletion,
                CigarOpKind::SoftClip => CoreKind::SoftClip,
                CigarOpKind::HardClip => CoreKind::HardClip,
            };
            CoreOp::new(kind, op.len)
        })
        .collect();
    let flags = if r.is_reverse {
        rosalind::core::SamFlags(rosalind::core::SamFlags::REVERSE)
    } else {
        rosalind::core::SamFlags::default()
    };
    rosalind::core::AlignedRead {
        contig,
        pos: rosalind::core::Position(r.pos),
        mapq: r.mapq,
        flags,
        cigar,
        seq: std::sync::Arc::clone(&r.sequence),
        qual: std::sync::Arc::clone(&r.qualities),
    }
}

fn read_bam_alignment_file(
    path: &PathBuf,
    target_chrom: Option<&Arc<str>>,
) -> Result<Vec<AlignedRead>> {
    let mut reader = bam::Reader::from_path(path)
        .with_context(|| format!("failed to open BAM {}", path.display()))?;
    let header = reader.header().to_owned();

    let mut reads = Vec::new();
    for record_result in reader.records() {
        let record = record_result?;

        if record.is_unmapped() {
            continue;
        }

        let tid = record.tid();
        if tid < 0 {
            continue;
        }
        let rname_bytes = header.tid2name(tid as u32);
        let rname = std::str::from_utf8(rname_bytes)
            .map_err(|_| anyhow!("BAM reference name is not valid UTF-8"))?;
        if let Some(target) = target_chrom {
            if rname != target.as_ref() {
                continue;
            }
        }

        let pos0 = record.pos();
        if pos0 < 0 {
            continue;
        }

        let mut cigar = Vec::new();
        for c in record.cigar().iter() {
            let (kind, len) = match *c {
                BamCigar::Match(l) | BamCigar::Equal(l) | BamCigar::Diff(l) => {
                    (CigarOpKind::Match, l)
                }
                BamCigar::Ins(l) => (CigarOpKind::Insertion, l),
                BamCigar::Del(l) => (CigarOpKind::Deletion, l),
                BamCigar::SoftClip(l) => (CigarOpKind::SoftClip, l),
                BamCigar::HardClip(l) => (CigarOpKind::HardClip, l),
                _ => continue,
            };
            cigar.push(CigarOp::new(kind, len));
        }

        let seq = record.seq().as_bytes();
        let sequence: Vec<u8> = seq.iter().map(|b| b.to_ascii_uppercase()).collect();
        let qualities: Vec<u8> = record.qual().to_vec();

        reads.push(AlignedRead::new(
            rname.to_string(),
            pos0 as u32,
            record.mapq(),
            cigar,
            sequence,
            qualities,
            record.is_reverse(),
        ));
    }

    Ok(reads)
}

fn parse_cigar(cigar: &str, read_len: usize) -> Result<Vec<CigarOp>> {
    if cigar == "*" {
        bail!("variant calling requires mapped reads (CIGAR cannot be '*')");
    }

    let mut ops = Vec::new();
    let mut run = String::new();
    let mut read_consuming: usize = 0;

    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            run.push(ch);
            continue;
        }
        let len: u32 = run
            .parse()
            .with_context(|| format!("invalid CIGAR run length in '{cigar}'"))?;
        run.clear();
        // Map operations the same way `read_bam_alignment_file` does; skip ops
        // that `CigarOpKind` does not represent (e.g. N/P), mirroring that reader.
        let kind = match ch {
            'M' | '=' | 'X' => CigarOpKind::Match,
            'I' => CigarOpKind::Insertion,
            'D' => CigarOpKind::Deletion,
            'S' => CigarOpKind::SoftClip,
            'H' => CigarOpKind::HardClip,
            _ => continue,
        };
        // M/=/X, I, and S consume read bases; D and H do not.
        if matches!(
            kind,
            CigarOpKind::Match | CigarOpKind::Insertion | CigarOpKind::SoftClip
        ) {
            read_consuming += len as usize;
        }
        ops.push(CigarOp::new(kind, len));
    }

    if !run.is_empty() {
        bail!("CIGAR '{cigar}' ends with a run length but no operation");
    }
    if ops.is_empty() {
        bail!("CIGAR '{cigar}' contains no supported operations");
    }
    if read_consuming != read_len {
        bail!("CIGAR '{cigar}' consumes {read_consuming} read bases but the read has {read_len}");
    }

    Ok(ops)
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
    fn parse_cigar_handles_multi_op_indels() {
        // The exact CIGAR the SAM reader previously choked on (regression).
        // 146M + 2I + 2M consumes 150 read bases.
        let ops = parse_cigar("146M2I2M", 150).expect("multi-op CIGAR should parse");
        assert_eq!(
            ops,
            vec![
                CigarOp::new(CigarOpKind::Match, 146),
                CigarOp::new(CigarOpKind::Insertion, 2),
                CigarOp::new(CigarOpKind::Match, 2),
            ]
        );
    }

    #[test]
    fn parse_cigar_handles_match_softclip_and_deletion() {
        assert_eq!(
            parse_cigar("150M", 150).unwrap(),
            vec![CigarOp::new(CigarOpKind::Match, 150)]
        );
        // Soft-clips consume read bases (5 + 140 + 5 = 150).
        let sc = parse_cigar("5S140M5S", 150).unwrap();
        assert_eq!(sc[0], CigarOp::new(CigarOpKind::SoftClip, 5));
        assert_eq!(sc[2], CigarOp::new(CigarOpKind::SoftClip, 5));
        // Deletions consume reference only, so read-consuming length is 20.
        let del = parse_cigar("10M2D10M", 20).unwrap();
        assert_eq!(del[1], CigarOp::new(CigarOpKind::Deletion, 2));
    }

    #[test]
    fn parse_cigar_rejects_length_mismatch_and_star() {
        assert!(parse_cigar("100M", 150).is_err());
        assert!(parse_cigar("*", 0).is_err());
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

    // ── fastq_quals_to_phred unit tests ────────────────────────────────────────
    // BAM/htslib expect raw Phred (0–93) in Record::set; the FASTQ parser keeps
    // ASCII-encoded (Phred+33) bytes. These tests guard the conversion helper.

    #[test]
    fn fastq_quals_to_phred_decodes_single_byte() {
        // ASCII b'I' = 73; Phred = 73 - 33 = 40 (Q40, a typical high-quality base).
        assert_eq!(fastq_quals_to_phred(b"I"), vec![40]);
    }

    #[test]
    fn fastq_quals_to_phred_decodes_lowest_quality() {
        // ASCII b'!' = 33; Phred = 33 - 33 = 0.
        assert_eq!(fastq_quals_to_phred(b"!"), vec![0]);
    }

    #[test]
    fn fastq_quals_to_phred_decodes_run() {
        // b"!I~" → [0, 40, 93]  (b'~' = 126; 126 - 33 = 93, the maximum valid Phred score).
        assert_eq!(fastq_quals_to_phred(b"!I~"), vec![0, 40, 93]);
    }

    #[test]
    fn fastq_quals_to_phred_empty_input_returns_empty() {
        assert_eq!(fastq_quals_to_phred(b""), Vec::<u8>::new());
    }

    #[test]
    fn fastq_quals_to_phred_saturating_does_not_underflow() {
        // A byte below 33 should saturate to 0 rather than wrapping (defensive).
        assert_eq!(fastq_quals_to_phred(&[0u8]), vec![0]);
        assert_eq!(fastq_quals_to_phred(&[32u8]), vec![0]);
    }
}
