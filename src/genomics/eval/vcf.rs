use std::collections::BTreeMap;

use thiserror::Error;

/// A minimal parsed VCF variant record (single-allelic only).
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
        let alt_field = cols[4];
        // Single-allele only for v1 harness.
        if alt_field.contains(',') {
            continue;
        }
        let alternate = alt_field.as_bytes().to_vec();
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

        out.push(VcfVariant {
            chrom,
            pos0,
            reference,
            alternate,
            qual,
            filter,
            info,
        });
    }

    Ok(out)
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


