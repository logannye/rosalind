use anyhow::{anyhow, Result};
use std::io::Write;

use super::SomaticVariant;

const HEADER: &str = "##fileformat=VCFv4.3\n##source=Rosalind\n##INFO=<ID=DP_T,Number=1,Type=Integer,Description=\"Tumor depth\">\n##INFO=<ID=DP_N,Number=1,Type=Integer,Description=\"Normal depth\">\n##INFO=<ID=AF_T,Number=1,Type=Float,Description=\"Tumor alt allele fraction\">\n##INFO=<ID=AF_N,Number=1,Type=Float,Description=\"Normal alt allele fraction\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// Write somatic SNVs to VCF (INFO-only, no sample columns yet).
pub fn write_somatic_vcf<W: Write>(writer: &mut W, variants: &[SomaticVariant]) -> Result<()> {
    writer.write_all(HEADER.as_bytes())?;

    // Deterministic canonical ordering.
    let mut ordered: Vec<&SomaticVariant> = variants.iter().collect();
    ordered.sort_by(|a, b| {
        a.chrom
            .as_ref()
            .cmp(b.chrom.as_ref())
            .then_with(|| a.position.cmp(&b.position))
            .then_with(|| a.reference.cmp(&b.reference))
            .then_with(|| a.alternate.cmp(&b.alternate))
    });

    for v in ordered {
        let line = format!(
            "{chrom}\t{pos}\t.\t{ref_base}\t{alt_base}\t{qual:.2}\t{filter}\tDP_T={dp_t};DP_N={dp_n};AF_T={af_t:.3};AF_N={af_n:.3}\n",
            chrom = v.chrom,
            pos = v.position + 1,
            ref_base = v.reference as char,
            alt_base = v.alternate as char,
            qual = v.quality,
            filter = v.filter,
            dp_t = v.tumor_depth,
            dp_n = v.normal_depth,
            af_t = v.tumor_af,
            af_n = v.normal_af
        );
        writer.write_all(line.as_bytes())?;
    }

    writer.flush()?;
    Ok(())
}

/// Render somatic SNVs into a VCF string (useful for tests and snapshots).
pub fn render_somatic_vcf(variants: &[SomaticVariant]) -> Result<String> {
    let mut buffer = Vec::new();
    write_somatic_vcf(&mut buffer, variants)?;
    String::from_utf8(buffer).map_err(|_| anyhow!("rendered VCF is not valid UTF-8"))
}
