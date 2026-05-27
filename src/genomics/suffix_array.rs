//! Suffix array construction.
//!
//! This module provides a production-grade suffix array builder for large
//! references. We use SA-IS induced sorting over an integer alphabet.
//!
//! Notes:
//! - Input is an integer slice where 0 is reserved for the unique sentinel and
//!   all other symbols are in `1..=K`.
//! - Output is a suffix array of length `text.len()` with indices in `0..n`.
//! - For human-scale references, indices are stored as `u32` (n must fit).

use thiserror::Error;

/// Errors returned by suffix array construction.
#[derive(Debug, Error)]
pub enum SuffixArrayError {
    #[error("text too large for u32 indices: length={0}")]
    TooLarge(usize),
    #[error("invalid alphabet: symbol {symbol} exceeds max {max_symbol}")]
    InvalidAlphabet { symbol: u32, max_symbol: u32 },
    #[error("text must contain at least one symbol")]
    Empty,
}

/// Build a suffix array using SA-IS.
///
/// Requirements:
/// - `text.len() > 0`
/// - `text` contains a unique sentinel `0` at the final position.
/// - All symbols are in `[0..=max_symbol]`.
pub fn sais_u32(text: &[u32], max_symbol: u32) -> Result<Vec<u32>, SuffixArrayError> {
    let n = text.len();
    if n == 0 {
        return Err(SuffixArrayError::Empty);
    }
    if n > (u32::MAX as usize) {
        return Err(SuffixArrayError::TooLarge(n));
    }
    for &sym in text {
        if sym > max_symbol {
            return Err(SuffixArrayError::InvalidAlphabet {
                symbol: sym,
                max_symbol,
            });
        }
    }
    // Common SA-IS contract: last symbol is sentinel and must be minimal.
    debug_assert_eq!(*text.last().unwrap(), 0);

    Ok(sais_impl(text, max_symbol))
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum SuffixType {
    S,
    L,
}

fn classify_types(text: &[u32]) -> Vec<SuffixType> {
    let n = text.len();
    let mut t = vec![SuffixType::S; n];
    t[n - 1] = SuffixType::S;
    for i in (0..n - 1).rev() {
        let a = text[i];
        let b = text[i + 1];
        if a < b {
            t[i] = SuffixType::S;
        } else if a > b {
            t[i] = SuffixType::L;
        } else {
            t[i] = t[i + 1];
        }
    }
    t
}

#[inline]
fn is_lms(types: &[SuffixType], i: usize) -> bool {
    i > 0 && types[i] == SuffixType::S && types[i - 1] == SuffixType::L
}

fn bucket_sizes(text: &[u32], max_symbol: u32) -> Vec<usize> {
    let mut sizes = vec![0usize; (max_symbol as usize) + 1];
    for &sym in text {
        sizes[sym as usize] += 1;
    }
    sizes
}

fn bucket_heads(sizes: &[usize]) -> Vec<usize> {
    let mut heads = vec![0usize; sizes.len()];
    let mut sum = 0usize;
    for (i, &sz) in sizes.iter().enumerate() {
        heads[i] = sum;
        sum += sz;
    }
    heads
}

fn bucket_tails(sizes: &[usize]) -> Vec<usize> {
    let mut tails = vec![0usize; sizes.len()];
    let mut sum = 0usize;
    for (i, &sz) in sizes.iter().enumerate() {
        sum += sz;
        tails[i] = sum - 1;
    }
    tails
}

fn induce_sort(text: &[u32], max_symbol: u32, types: &[SuffixType], lms: &[usize]) -> Vec<i32> {
    let n = text.len();
    let sizes = bucket_sizes(text, max_symbol);
    let mut sa = vec![-1i32; n];

    // 1) Place LMS at bucket tails.
    {
        let mut tails = bucket_tails(&sizes);
        for &pos in lms.iter().rev() {
            let sym = text[pos] as usize;
            sa[tails[sym]] = pos as i32;
            if tails[sym] > 0 {
                tails[sym] -= 1;
            }
        }
    }

    // 2) Induce L-type to bucket heads.
    {
        let mut heads = bucket_heads(&sizes);
        for i in 0..n {
            let j = sa[i];
            if j <= 0 {
                continue;
            }
            let p = (j - 1) as usize;
            if types[p] == SuffixType::L {
                let sym = text[p] as usize;
                sa[heads[sym]] = p as i32;
                heads[sym] += 1;
            }
        }
    }

    // 3) Induce S-type to bucket tails.
    {
        let mut tails = bucket_tails(&sizes);
        for i in (0..n).rev() {
            let j = sa[i];
            if j <= 0 {
                continue;
            }
            let p = (j - 1) as usize;
            if types[p] == SuffixType::S {
                let sym = text[p] as usize;
                sa[tails[sym]] = p as i32;
                if tails[sym] > 0 {
                    tails[sym] -= 1;
                }
            }
        }
    }

    sa
}

fn lms_substring_equal(text: &[u32], types: &[SuffixType], a: usize, b: usize) -> bool {
    if a == b {
        return true;
    }
    let n = text.len();
    let mut i = 0usize;
    loop {
        let a_i = a + i;
        let b_i = b + i;
        if a_i >= n || b_i >= n {
            return false;
        }
        if text[a_i] != text[b_i] {
            return false;
        }
        let a_lms = is_lms(types, a_i);
        let b_lms = is_lms(types, b_i);
        if a_lms != b_lms {
            return false;
        }
        if i > 0 && a_lms && b_lms {
            return true;
        }
        i += 1;
    }
}

fn sais_impl(text: &[u32], max_symbol: u32) -> Vec<u32> {
    let n = text.len();
    if n == 1 {
        return vec![0u32];
    }

    let types = classify_types(text);
    let mut lms_positions = Vec::new();
    for i in 1..n {
        if is_lms(&types, i) {
            lms_positions.push(i);
        }
    }

    // Induced sort LMS positions.
    let sa0 = induce_sort(text, max_symbol, &types, &lms_positions);

    // Extract LMS order from SA.
    let mut lms_in_sa_order = Vec::with_capacity(lms_positions.len());
    for &idx in &sa0 {
        if idx >= 0 {
            let p = idx as usize;
            if is_lms(&types, p) {
                lms_in_sa_order.push(p);
            }
        }
    }

    // Name LMS substrings.
    let mut name = 0u32;
    let mut lms_name = vec![u32::MAX; n];
    let mut prev = None;
    for &p in &lms_in_sa_order {
        if let Some(prev_p) = prev {
            if !lms_substring_equal(text, &types, prev_p, p) {
                name += 1;
            }
        }
        lms_name[p] = name;
        prev = Some(p);
    }
    let new_alphabet = name + 1;

    // Build reduced string (in original LMS order).
    let mut reduced = Vec::with_capacity(lms_positions.len() + 1);
    for &p in &lms_positions {
        reduced.push(lms_name[p]);
    }
    // Reduced string must end with sentinel 0.
    // The last LMS is always the sentinel position (n-1) in SA-IS for valid input,
    // but our LMS list excludes n-1 unless it is LMS; ensure reduced ends with 0 by appending.
    // For typical SA-IS, the sentinel at n-1 is LMS by definition (since it is S and previous is L).
    // Still, be defensive:
    if *text.last().unwrap() == 0 {
        // Ensure last LMS exists.
        if !lms_positions.contains(&(n - 1)) {
            reduced.push(0);
        }
    }

    // Recurse or compute SA of reduced directly.
    let reduced_sa: Vec<u32> = if new_alphabet == reduced.len() as u32 {
        // Names are unique: we can construct SA by inverse mapping.
        let mut sa = vec![0u32; reduced.len()];
        for (i, &sym) in reduced.iter().enumerate() {
            sa[sym as usize] = i as u32;
        }
        sa
    } else {
        // Recurse. Max symbol is new_alphabet-1.
        sais_impl(&reduced, new_alphabet - 1)
    };

    // Map reduced SA back to LMS positions (reverse LMS order list).
    let mut ordered_lms = Vec::with_capacity(lms_positions.len());
    for &i in reduced_sa.iter().skip(1) {
        // Skip sentinel in reduced SA.
        let pos = lms_positions[i as usize];
        ordered_lms.push(pos);
    }
    // Append sentinel LMS (n-1).
    ordered_lms.push(n - 1);

    // Final induced sort using ordered LMS.
    let sa_final = induce_sort(text, max_symbol, &types, &ordered_lms);
    sa_final.into_iter().map(|v| v as u32).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn naive_sa(text: &[u32]) -> Vec<u32> {
        let mut sa: Vec<u32> = (0..text.len() as u32).collect();
        sa.sort_by(|&a, &b| text[a as usize..].cmp(&text[b as usize..]));
        sa
    }

    #[test]
    fn sais_matches_naive_small() {
        // Alphabet: {0,1,2,3} with sentinel 0 at end.
        let text = vec![2, 1, 3, 1, 2, 0];
        let sa = sais_u32(&text, 3).expect("sais");
        assert_eq!(sa, naive_sa(&text));
    }

    #[test]
    fn sais_handles_repeated_symbols() {
        let text = vec![1, 1, 1, 1, 0];
        let sa = sais_u32(&text, 1).expect("sais");
        assert_eq!(sa, naive_sa(&text));
    }
}
