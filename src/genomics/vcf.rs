use anyhow::{anyhow, Result};
use std::io::Write;

use super::Variant;

const HEADER: &str =
    "##fileformat=VCFv4.3\n##source=Rosalind\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// Write variants in a normalized VCF form.
pub fn write_vcf<W: Write>(writer: &mut W, variants: &[Variant]) -> Result<()> {
    writer.write_all(HEADER.as_bytes())?;

    // Determinism requirement: callers may provide variants in any order.
    // Emit records in a canonical, stable order.
    let mut ordered: Vec<&Variant> = variants.iter().collect();
    ordered.sort_by(|a, b| {
        a.chrom
            .as_ref()
            .cmp(b.chrom.as_ref())
            .then_with(|| a.position.cmp(&b.position))
            .then_with(|| a.reference.cmp(&b.reference))
            .then_with(|| a.alternate.cmp(&b.alternate))
            .then_with(|| a.depth.cmp(&b.depth))
    });

    for variant in ordered {
        let line = format!(
            "{chrom}\t{pos}\t.\t{ref_base}\t{alt_base}\t{qual:.2}\tPASS\tDP={depth};AF={af:.3}\n",
            chrom = variant.chrom,
            pos = variant.position + 1,
            ref_base = variant.reference as char,
            alt_base = variant.alternate as char,
            depth = variant.depth,
            af = variant.allele_fraction,
            qual = variant.quality
        );
        writer.write_all(line.as_bytes())?;
    }

    writer.flush()?;
    Ok(())
}

/// Render variants into a VCF string (useful for tests and snapshots).
pub fn render_vcf(variants: &[Variant]) -> Result<String> {
    let mut buffer = Vec::new();
    write_vcf(&mut buffer, variants)?;
    String::from_utf8(buffer).map_err(|_| anyhow!("rendered VCF is not valid UTF-8"))
}
