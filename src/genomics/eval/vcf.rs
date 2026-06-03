use std::collections::BTreeMap;

use thiserror::Error;

/// Diploid zygosity of a biallelic call, parsed from a VCF `GT`. Phase-insensitive
/// (`0|1` and `0/1` are both `Het`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Zygosity {
    /// One alt copy (`0/1`).
    Het,
    /// Two alt copies (`1/1`).
    Hom,
}

/// A minimal parsed VCF variant record (one biallelic record per ALT allele;
/// multi-allelic sites are decomposed on read).
#[derive(Debug, Clone, PartialEq)]
pub struct VcfVariant {
    /// Contig name.
    pub chrom: String,
    /// 0-based genomic position (VCF POS-1).
    pub pos0: u32,
    /// Reference allele as bytes.
    pub reference: Vec<u8>,
    /// Alternate allele as bytes.
    pub alternate: Vec<u8>,
    /// QUAL, if present.
    pub qual: Option<f32>,
    /// FILTER column (raw).
    pub filter: String,
    /// INFO key/value pairs (best-effort parsing).
    pub info: BTreeMap<String, String>,
    /// Diploid zygosity from the sample `GT`, if a FORMAT/sample column is present.
    pub genotype: Option<Zygosity>,
}

#[derive(Debug, Error)]
/// Errors that can occur while parsing VCF.
pub enum VcfParseError {
    /// A VCF data line could not be parsed.
    #[error("invalid VCF line {line}: {msg}")]
    InvalidLine {
        /// 1-based line number in the input.
        line: usize,
        /// Human-readable parse error.
        msg: String,
    },
}

/// Parse a (possibly headered) VCF into a vector of single-allelic variants.
pub fn read_vcf_variants(contents: &str) -> Result<Vec<VcfVariant>, VcfParseError> {
    let mut out = Vec::new();

    for (i, line) in contents.lines().enumerate() {
        let line_no = i + 1;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = trimmed.split('\t').collect();
        if cols.len() < 8 {
            return Err(VcfParseError::InvalidLine {
                line: line_no,
                msg: "expected at least 8 tab-separated columns".to_string(),
            });
        }
        let chrom = cols[0].to_string();
        let pos1: u32 = cols[1].parse().map_err(|_| VcfParseError::InvalidLine {
            line: line_no,
            msg: "POS is not a u32".to_string(),
        })?;
        if pos1 == 0 {
            return Err(VcfParseError::InvalidLine {
                line: line_no,
                msg: "POS must be 1-based".to_string(),
            });
        }
        let pos0 = pos1 - 1;
        let reference = cols[3].as_bytes().to_vec();
        let qual = if cols[5] == "." {
            None
        } else {
            Some(cols[5].parse().map_err(|_| VcfParseError::InvalidLine {
                line: line_no,
                msg: "QUAL is not a float".to_string(),
            })?)
        };
        let filter = cols[6].to_string();
        let info = parse_info(cols[7]);

        // GT allele indices from the FORMAT (cols[8]) + first sample (cols[9]), if present.
        let gt_indices = parse_gt_indices(cols.get(8).copied(), cols.get(9).copied());

        // Decompose multi-allelic ALT into one biallelic record per ALT allele.
        for (alt_idx, alt) in cols[4].split(',').enumerate() {
            let allele_index = (alt_idx + 1) as u8; // GT indices: ref=0, first alt=1, ...
            let genotype = match &gt_indices {
                // No GT column: emit the allele for detection, unscored for genotype.
                None => None,
                // GT present: zygosity = how many times this allele index appears in it.
                Some(idx) => match idx.iter().filter(|&&a| a == allele_index).count() {
                    0 => continue, // allele absent from the genotype -> drop this record
                    1 => Some(Zygosity::Het),
                    _ => Some(Zygosity::Hom),
                },
            };
            out.push(VcfVariant {
                chrom: chrom.clone(),
                pos0,
                reference: reference.clone(),
                alternate: alt.as_bytes().to_vec(),
                qual,
                filter: filter.clone(),
                info: info.clone(),
                genotype,
            });
        }
    }

    Ok(out)
}

/// Parse the GT allele indices (e.g. `0/1` -> [0,1], `1|2` -> [1,2]) from the
/// FORMAT + first sample columns. Returns `None` if no GT field is present.
fn parse_gt_indices(format: Option<&str>, sample: Option<&str>) -> Option<Vec<u8>> {
    let (format, sample) = (format?, sample?);
    let gt_pos = format.split(':').position(|k| k == "GT")?;
    let gt = sample.split(':').nth(gt_pos)?;
    let indices: Vec<u8> = gt
        .split(['/', '|'])
        .filter_map(|a| a.parse::<u8>().ok()) // '.' (no-call) is filtered out
        .collect();
    if indices.is_empty() {
        None
    } else {
        Some(indices)
    }
}

fn parse_info(info: &str) -> BTreeMap<String, String> {
    let mut map = BTreeMap::new();
    if info == "." {
        return map;
    }
    for item in info.split(';') {
        if item.is_empty() {
            continue;
        }
        if let Some((k, v)) = item.split_once('=') {
            map.insert(k.to_string(), v.to_string());
        } else {
            map.insert(item.to_string(), "true".to_string());
        }
    }
    map
}

#[cfg(test)]
mod gt_tests {
    use super::*;

    #[test]
    fn parses_gt_zygosity_from_format_sample() {
        let vcf = "chr1\t10\t.\tA\tC\t.\tPASS\t.\tGT:DP\t0/1:30\n\
                   chr1\t20\t.\tG\tT\t.\tPASS\t.\tGT\t1/1\n\
                   chr1\t30\t.\tT\tA\t.\tPASS\t.\tGT\t0|1\n";
        let v = read_vcf_variants(vcf).unwrap();
        assert_eq!(v[0].genotype, Some(Zygosity::Het));
        assert_eq!(v[1].genotype, Some(Zygosity::Hom));
        assert_eq!(v[2].genotype, Some(Zygosity::Het)); // phase-insensitive
    }

    #[test]
    fn missing_sample_is_genotype_none_and_still_parses() {
        // An 8-column VCF (no FORMAT/sample) still parses, genotype unknown.
        let vcf = "chr1\t10\t.\tA\tC\t.\tPASS\t.\n";
        let v = read_vcf_variants(vcf).unwrap();
        assert_eq!(v.len(), 1);
        assert_eq!(v[0].genotype, None);
    }

    #[test]
    fn multiallelic_decomposes_per_alt_with_decomposed_genotype() {
        // ALT A,C with GT 1/2 -> two biallelic records, each Het.
        let vcf = "chr1\t10\t.\tG\tA,C\t.\tPASS\t.\tGT\t1/2\n";
        let v = read_vcf_variants(vcf).unwrap();
        assert_eq!(v.len(), 2);
        assert_eq!(v[0].alternate, b"A");
        assert_eq!(v[0].genotype, Some(Zygosity::Het));
        assert_eq!(v[1].alternate, b"C");
        assert_eq!(v[1].genotype, Some(Zygosity::Het));
    }

    #[test]
    fn multiallelic_hom_second_allele_drops_the_absent_first() {
        // ALT A,C with GT 2/2 -> only the C record, Hom.
        let vcf = "chr1\t10\t.\tG\tA,C\t.\tPASS\t.\tGT\t2/2\n";
        let v = read_vcf_variants(vcf).unwrap();
        assert_eq!(v.len(), 1);
        assert_eq!(v[0].alternate, b"C");
        assert_eq!(v[0].genotype, Some(Zygosity::Hom));
    }

    #[test]
    fn multiallelic_without_gt_emits_all_alleles_unscored() {
        let vcf = "chr1\t10\t.\tG\tA,C\t.\tPASS\t.\n";
        let v = read_vcf_variants(vcf).unwrap();
        assert_eq!(v.len(), 2);
        assert_eq!(v[0].genotype, None);
        assert_eq!(v[1].genotype, None);
    }
}
