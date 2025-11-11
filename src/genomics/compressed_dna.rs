use std::fmt;

use thiserror::Error;

/// Number of bases encoded per `u64` chunk.
const BASES_PER_WORD: usize = 32;
/// Bits used to encode a single DNA base (A/C/G/T).
const BITS_PER_BASE: usize = 2;

/// Bitmask tracking ambiguous bases (currently only `N`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AmbiguityMask {
    bits: Vec<u64>,
    len: usize,
}

impl AmbiguityMask {
    fn new(len: usize) -> Self {
        let words = words_for_len(len);
        Self {
            bits: vec![0; words],
            len,
        }
    }

    fn set(&mut self, idx: usize) {
        let (word_idx, bit_idx) = bit_position(idx);
        if word_idx >= self.bits.len() {
            self.bits.resize(word_idx + 1, 0);
        }
        self.bits[word_idx] |= 1u64 << bit_idx;
        self.len = self.len.max(idx + 1);
    }

    fn test(&self, idx: usize) -> bool {
        let (word_idx, bit_idx) = bit_position(idx);
        self.bits
            .get(word_idx)
            .map(|word| (word >> bit_idx) & 1 == 1)
            .unwrap_or(false)
    }

    /// Returns the number of positions tracked by this mask.
    pub fn len(&self) -> usize {
        self.len
    }

    /// Access the underlying bit words (useful for serialization).
    pub fn bits(&self) -> &[u64] {
        &self.bits
    }
}

/// Errors that can occur while working with compressed DNA sequences.
#[derive(Debug, Error)]
pub enum CompressedDNAError {
    /// Encountered a base that cannot be represented in the 2-bit alphabet.
    #[error("unsupported nucleotide '{0}' at position {1}")]
    UnsupportedBase(char, usize),
}

/// DNA sequence compressed using 2-bit encoding per base.
///
/// Only the canonical bases (A, C, G, T) are stored directly; ambiguous bases
/// (currently only `N`) are tracked separately via an ambiguity mask so that
/// decoding restores the original character.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CompressedDNA {
    data: Vec<u64>,
    len: usize,
    ambiguity: AmbiguityMask,
}

impl CompressedDNA {
    /// Compress a DNA string (ASCII bases) into 2-bit representation.
    pub fn compress(sequence: &[u8]) -> Result<Self, CompressedDNAError> {
        let len = sequence.len();
        let words = words_for_len(len);
        let mut data = vec![0u64; words];
        let mut ambiguity = AmbiguityMask::new(len);

        for (idx, &base) in sequence.iter().enumerate() {
            let (code, is_ambiguous) =
                encode_base(base).ok_or_else(|| CompressedDNAError::UnsupportedBase(
                    base as char,
                    idx,
                ))?;

            if is_ambiguous {
                ambiguity.set(idx);
            }

            let (word_idx, bit_shift) = word_position(idx);
            data[word_idx] |= (code as u64) << bit_shift;
        }

        Ok(Self { data, len, ambiguity })
    }

    /// Create an owned compressed DNA sequence from a vector of packed words.
    ///
    /// # Panics
    /// Panics if `len` exceeds the capacity implied by `data`.
    pub fn from_parts(data: Vec<u64>, len: usize, ambiguity: AmbiguityMask) -> Self {
        let capacity = data.len() * BASES_PER_WORD;
        assert!(
            len <= capacity,
            "length {} exceeds backing capacity {} ({} words)",
            len,
            capacity,
            data.len()
        );
        Self { data, len, ambiguity }
    }

    /// Number of bases in the sequence.
    pub fn len(&self) -> usize {
        self.len
    }

    /// Returns `true` when the sequence is empty.
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Access the underlying packed words (little-endian base order).
    pub fn words(&self) -> &[u64] {
        &self.data
    }

    /// Access the ambiguity mask.
    pub fn ambiguity(&self) -> &AmbiguityMask {
        &self.ambiguity
    }

    /// Retrieve the base at `idx` as an uppercase ASCII byte.
    pub fn base_at(&self, idx: usize) -> Option<u8> {
        if idx >= self.len {
            return None;
        }
        if self.ambiguity.test(idx) {
            return Some(b'N');
        }
        let (word_idx, bit_shift) = word_position(idx);
        let word = self.data[word_idx];
        let code = ((word >> bit_shift) & 0b11) as u8;
        Some(decode_base(code))
    }

    /// Decode into a newly allocated vector of uppercase ASCII bases.
    pub fn to_vec(&self) -> Vec<u8> {
        let mut out = vec![b'N'; self.len];
        self.decode_into(&mut out);
        out
    }

    /// Decode into the provided buffer (must be at least `self.len()`).
    pub fn decode_into(&self, out: &mut [u8]) {
        assert!(
            out.len() >= self.len,
            "output buffer too small: {} < {}",
            out.len(),
            self.len
        );

        for idx in 0..self.len {
            out[idx] = if self.ambiguity.test(idx) {
                b'N'
            } else {
                let (word_idx, bit_shift) = word_position(idx);
                let code = ((self.data[word_idx] >> bit_shift) & 0b11) as u8;
                decode_base(code)
            };
        }
    }

    /// Append a single base to the compressed sequence.
    pub fn push(&mut self, base: u8) -> Result<(), CompressedDNAError> {
        let idx = self.len;
        let (code, is_ambiguous) =
            encode_base(base).ok_or_else(|| CompressedDNAError::UnsupportedBase(
                base as char,
                idx,
            ))?;

        let (word_idx, bit_shift) = word_position(idx);
        if word_idx >= self.data.len() {
            self.data.push(0);
        }
        self.data[word_idx] |= (code as u64) << bit_shift;
        if is_ambiguous {
            self.ambiguity.set(idx);
        }
        self.len += 1;
        Ok(())
    }

    /// Extend the sequence by appending the provided bases.
    pub fn extend_from_slice(&mut self, sequence: &[u8]) -> Result<(), CompressedDNAError> {
        for &base in sequence {
            self.push(base)?;
        }
        Ok(())
    }

    /// Iterate over decoded bases (uppercase ASCII).
    pub fn iter(&self) -> CompressedDNAIter<'_> {
        CompressedDNAIter {
            dna: self,
            index: 0,
        }
    }
}

impl fmt::Display for CompressedDNA {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let decoded = self.to_vec();
        let as_str = String::from_utf8_lossy(&decoded);
        write!(f, "{as_str}")
    }
}

/// Iterator over decoded bases in a `CompressedDNA` sequence.
#[derive(Debug)]
pub struct CompressedDNAIter<'a> {
    dna: &'a CompressedDNA,
    index: usize,
}

impl<'a> Iterator for CompressedDNAIter<'a> {
    type Item = u8;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index >= self.dna.len {
            return None;
        }
        let base = self.dna.base_at(self.index);
        self.index += 1;
        base
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = self.dna.len.saturating_sub(self.index);
        (remaining, Some(remaining))
    }
}

impl ExactSizeIterator for CompressedDNAIter<'_> {}

fn encode_base(base: u8) -> Option<(u8, bool)> {
    match base {
        b'A' | b'a' => Some((0b00, false)),
        b'C' | b'c' => Some((0b01, false)),
        b'G' | b'g' => Some((0b10, false)),
        b'T' | b't' | b'U' | b'u' => Some((0b11, false)),
        b'N' | b'n' => Some((0b00, true)),
        _ => None,
    }
}

fn decode_base(code: u8) -> u8 {
    match code & 0b11 {
        0b00 => b'A',
        0b01 => b'C',
        0b10 => b'G',
        0b11 => b'T',
        _ => b'N',
    }
}

fn words_for_len(len: usize) -> usize {
    if len == 0 {
        0
    } else {
        (len + BASES_PER_WORD - 1) / BASES_PER_WORD
    }
}

fn word_position(idx: usize) -> (usize, usize) {
    let word_idx = idx / BASES_PER_WORD;
    let offset = idx % BASES_PER_WORD;
    let bit_shift = offset * BITS_PER_BASE;
    (word_idx, bit_shift)
}

fn bit_position(idx: usize) -> (usize, usize) {
    let word_idx = idx / 64;
    let bit_idx = idx % 64;
    (word_idx, bit_idx)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn compress_and_decode_roundtrip() {
        let seq = b"ACGTACGTNNACGT";
        let compressed = CompressedDNA::compress(seq).expect("compression should succeed");
        assert_eq!(compressed.len(), seq.len());

        let decoded = compressed.to_vec();
        assert_eq!(decoded, seq);
    }

    #[test]
    fn push_and_extend() {
        let mut dna =
            CompressedDNA::compress(b"ACG").expect("initial compression should succeed");
        dna.push(b'T').unwrap();
        dna.extend_from_slice(b"NN").unwrap();
        assert_eq!(dna.to_vec(), b"ACGTNN");
    }

    #[test]
    fn iterator_yields_correct_bases() {
        let seq = b"AACCGGTTNN";
        let compressed = CompressedDNA::compress(seq).unwrap();
        let collected: Vec<u8> = compressed.iter().collect();
        assert_eq!(collected, seq);
    }

    #[test]
    fn unsupported_base_returns_error() {
        let result = CompressedDNA::compress(b"ABCD");
        assert!(matches!(
            result,
            Err(CompressedDNAError::UnsupportedBase('B', 1))
        ));
    }
}

