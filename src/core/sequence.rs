//! Canonical DNA base model. Exact bases are A/C/G/T (U folds to T); every
//! other byte — including IUPAC ambiguity codes — folds to `N` (lossy but
//! explicit). `N` is not a callable allele.

/// A canonical DNA base. `N` represents any non-ACGT / ambiguous base.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BaseCode {
    /// Adenine.
    A,
    /// Cytosine.
    C,
    /// Guanine.
    G,
    /// Thymine (also the fold target for uracil `U`).
    T,
    /// Any non-ACGT / ambiguous base.
    N,
}

impl BaseCode {
    /// Allele index 0..=3 for A/C/G/T; `None` for `N` (not a callable allele).
    pub fn allele_index(self) -> Option<usize> {
        match self {
            BaseCode::A => Some(0),
            BaseCode::C => Some(1),
            BaseCode::G => Some(2),
            BaseCode::T => Some(3),
            BaseCode::N => None,
        }
    }

    /// Canonical uppercase ASCII byte for this base.
    pub fn to_ascii(self) -> u8 {
        match self {
            BaseCode::A => b'A',
            BaseCode::C => b'C',
            BaseCode::G => b'G',
            BaseCode::T => b'T',
            BaseCode::N => b'N',
        }
    }

    /// Map an ASCII byte to a `BaseCode`. A/C/G/T (any case) map directly,
    /// `U`/`u` fold to `T`, and every other byte folds to `N`.
    pub fn from_ascii_lossy(b: u8) -> BaseCode {
        match b.to_ascii_uppercase() {
            b'A' => BaseCode::A,
            b'C' => BaseCode::C,
            b'G' => BaseCode::G,
            b'T' | b'U' => BaseCode::T,
            _ => BaseCode::N,
        }
    }

    /// Whether `b` is a recognized exact base (A/C/G/T/U). Bytes that fold to
    /// `N` return `false`; callers use this to count ambiguity warnings.
    pub fn is_exact(b: u8) -> bool {
        matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'U')
    }
}

/// Convenience: allele index 0..=3 for an ASCII base, or `None` for N/other.
pub fn allele_index(base: u8) -> Option<usize> {
    BaseCode::from_ascii_lossy(base).allele_index()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_bases_map_directly_and_u_folds_to_t() {
        assert_eq!(BaseCode::from_ascii_lossy(b'a'), BaseCode::A);
        assert_eq!(BaseCode::from_ascii_lossy(b'C'), BaseCode::C);
        assert_eq!(BaseCode::from_ascii_lossy(b'u'), BaseCode::T);
        assert_eq!(BaseCode::from_ascii_lossy(b'T'), BaseCode::T);
    }

    #[test]
    fn ambiguous_and_unknown_fold_to_n_and_are_not_exact() {
        assert_eq!(BaseCode::from_ascii_lossy(b'R'), BaseCode::N); // IUPAC purine
        assert_eq!(BaseCode::from_ascii_lossy(b'.'), BaseCode::N);
        assert!(!BaseCode::is_exact(b'R'));
        assert!(BaseCode::is_exact(b'g'));
    }

    #[test]
    fn allele_index_is_some_for_acgt_and_none_for_n() {
        assert_eq!(allele_index(b'A'), Some(0));
        assert_eq!(allele_index(b'T'), Some(3));
        assert_eq!(allele_index(b'N'), None);
        assert_eq!(allele_index(b'R'), None);
    }
}
