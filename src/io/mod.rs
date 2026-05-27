//! IO layer: standards-compliant readers and writers. Phase A4 lands the
//! spec-valid VCF writer; readers (FASTA/FASTQ/BAM) migrate here in later
//! phases.

pub mod vcf;

pub mod bam;
