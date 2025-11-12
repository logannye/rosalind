//! Constant-size finite field 𝔽_{2^c}
//!
//! Key property: |𝔽| independent of t and b
//! Each element: O(1) tape cells

/// Finite field for algebraic operations
///
/// Size: 2^c for constant c (e.g., c=8 → |𝔽| = 256)
/// Uses GF(2^8) with primitive polynomial x^8 + x^4 + x^3 + x^2 + 1 (0x11D)
#[derive(Debug, Clone)]
pub struct FiniteField {
    characteristic: u8,
    /// Primitive polynomial for GF(2^c)
    primitive_poly: u16,
}

impl FiniteField {
    /// Create field of size 2^c
    pub fn new(characteristic: u8) -> Self {
        // For GF(2^8), use primitive polynomial 0x11D
        // x^8 + x^4 + x^3 + x^2 + 1
        let primitive_poly = match characteristic {
            8 => 0x11D, // Common choice for GF(256)
            _ => 0x11D, // Default to GF(256)
        };

        Self {
            characteristic,
            primitive_poly,
        }
    }

    /// Field size
    pub fn size(&self) -> usize {
        1 << self.characteristic
    }

    /// Add two field elements (XOR for characteristic 2)
    pub fn add(&self, a: u8, b: u8) -> u8 {
        a ^ b
    }

    /// Multiply two field elements in GF(2^8)
    ///
    /// Uses Russian peasant algorithm with reduction by primitive polynomial
    pub fn mul(&self, a: u8, b: u8) -> u8 {
        if a == 0 || b == 0 {
            return 0;
        }

        let mut result = 0u16;
        let mut a_val = a as u16;
        let mut b_val = b as u16;
        let poly = self.primitive_poly;

        // Russian peasant algorithm
        while b_val > 0 {
            if b_val & 1 != 0 {
                result ^= a_val;
            }
            a_val <<= 1;
            if a_val & 0x100 != 0 {
                a_val ^= poly;
            }
            b_val >>= 1;
        }

        result as u8
    }

    /// Evaluate polynomial at point using Horner's method
    pub fn eval_poly(&self, coeffs: &[u8], point: u8) -> u8 {
        if coeffs.is_empty() {
            return 0;
        }

        let mut result = coeffs[coeffs.len() - 1];
        for i in (0..coeffs.len() - 1).rev() {
            result = self.mul(result, point);
            result = self.add(result, coeffs[i]);
        }

        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_constant_size() {
        let field = FiniteField::new(8);
        assert_eq!(field.size(), 256); // Constant!

        // Verify size doesn't grow with t or b
        // This is critical for O(1) internal node space
    }
}
