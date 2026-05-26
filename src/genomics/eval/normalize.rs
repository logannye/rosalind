use thiserror::Error;

use super::VcfVariant;

#[derive(Debug, Error)]
/// Errors that can occur during variant normalization.
pub enum NormalizeError {
    /// The variant POS is outside the provided reference slice.
    #[error("variant position {pos0} out of bounds for reference length {reference_len}")]
    OutOfBounds {
        /// 0-based position.
        pos0: u32,
        /// Reference length.
        reference_len: usize,
    },
    /// An allele became empty during normalization.
    #[error("empty allele after normalization")]
    EmptyAllele,
}

/// A normalized variant key suitable for deterministic comparisons.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct NormalizedVariant {
    /// Contig name.
    pub chrom: String,
    /// 0-based coordinate.
    pub pos0: u32,
    /// Normalized reference allele.
    pub reference: Vec<u8>,
    /// Normalized alternate allele.
    pub alternate: Vec<u8>,
}

/// Deterministically normalize a VCF variant against a reference contig slice.
///
/// This performs:
/// - uppercase alleles
/// - trim common suffix then prefix (minimal representation)
/// - best-effort left-alignment for indels within the provided reference slice
pub fn normalize_variant(reference: &[u8], v: &VcfVariant) -> Result<NormalizedVariant, NormalizeError> {
    let mut pos0 = v.pos0;
    let mut r = v.reference.iter().map(|b| b.to_ascii_uppercase()).collect::<Vec<u8>>();
    let mut a = v.alternate.iter().map(|b| b.to_ascii_uppercase()).collect::<Vec<u8>>();

    if pos0 as usize >= reference.len() {
        return Err(NormalizeError::OutOfBounds {
            pos0,
            reference_len: reference.len(),
        });
    }

    // Minimal representation: trim common suffix.
    while r.len() > 1 && a.len() > 1 && r.last() == a.last() {
        r.pop();
        a.pop();
    }
    // Trim common prefix, adjusting position.
    while r.len() > 1 && a.len() > 1 && r.first() == a.first() {
        r.remove(0);
        a.remove(0);
        pos0 += 1;
    }

    if r.is_empty() || a.is_empty() {
        return Err(NormalizeError::EmptyAllele);
    }

    // Left-align indels (best-effort).
    if r.len() != a.len() {
        left_align(reference, &mut pos0, &mut r, &mut a)?;
    }

    Ok(NormalizedVariant {
        chrom: v.chrom.clone(),
        pos0,
        reference: r,
        alternate: a,
    })
}

fn left_align(
    reference: &[u8],
    pos0: &mut u32,
    r: &mut Vec<u8>,
    a: &mut Vec<u8>,
) -> Result<(), NormalizeError> {
    // Classic left alignment rotates the differing suffix as long as the base to the left matches.
    // This is correct for simple representations and sufficient for a deterministic harness.
    while *pos0 > 0 {
        let left_base = reference[(*pos0 as usize) - 1].to_ascii_uppercase();
        let r_last = *r.last().ok_or(NormalizeError::EmptyAllele)?;
        let a_last = *a.last().ok_or(NormalizeError::EmptyAllele)?;

        // Only rotate if the left base matches the last base of either allele (cyclic shift).
        // This works for the common insert/delete normalization patterns.
        if left_base != r_last && left_base != a_last {
            break;
        }

        // Shift left by 1: prepend left_base and drop last base.
        *pos0 -= 1;

        r.insert(0, left_base);
        r.pop();

        a.insert(0, left_base);
        a.pop();

        // Re-trim suffix/prefix again (keep minimal).
        while r.len() > 1 && a.len() > 1 && r.last() == a.last() {
            r.pop();
            a.pop();
        }
        while r.len() > 1 && a.len() > 1 && r.first() == a.first() {
            r.remove(0);
            a.remove(0);
            *pos0 += 1;
        }

        if r.is_empty() || a.is_empty() {
            return Err(NormalizeError::EmptyAllele);
        }
    }
    Ok(())
}


