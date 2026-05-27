//! Alignment refinement utilities (banded affine-gap DP).
//!
//! This module is intentionally small and deterministic. It is not yet a full
//! BWA-MEM/GATK-quality aligner, but provides the core building block needed to
//! refine seed-derived candidates into a CIGAR with mismatch/indel accounting.

use crate::genomics::{CigarOp, CigarOpKind};

const NEG_INF: i32 = i32::MIN / 4;

/// Scoring parameters for affine-gap alignment.
#[derive(Debug, Clone, Copy)]
pub struct AffineScores {
    pub match_score: i32,
    pub mismatch_score: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
}

impl Default for AffineScores {
    fn default() -> Self {
        Self {
            match_score: 2,
            mismatch_score: -4,
            gap_open: -6,
            gap_extend: -1,
        }
    }
}

/// Result of refining a candidate alignment.
#[derive(Debug, Clone)]
pub struct RefinedAlignment {
    /// 0-based reference start coordinate for the alignment.
    pub ref_start: usize,
    /// CIGAR operations (M/I/D only for now; soft clips will be added later).
    pub cigar: Vec<CigarOp>,
    /// Alignment score under the provided scoring parameters.
    pub score: i32,
    /// Edit distance style mismatch/indel count (NM tag semantics).
    pub nm: u32,
    /// SAM MD string describing mismatches and deletions.
    pub md: String,
}

/// Banded affine-gap alignment that consumes the entire read (semi-global on reference).
///
/// - The alignment consumes all bases of `read` (no clipping yet).
/// - The alignment may start and end anywhere within `reference_window` without penalty.
/// - The DP is restricted to a diagonal band of width `band`.
pub fn banded_affine_align(
    reference_window: &[u8],
    read: &[u8],
    band: usize,
    scores: AffineScores,
) -> Option<RefinedAlignment> {
    if read.is_empty() || reference_window.is_empty() {
        return None;
    }

    let m = read.len();
    let n = reference_window.len();

    // Allocate DP matrices within the full grid but only compute a band around the diagonal.
    // For read-length Illumina (~150) and modest band (~32), this is fine.
    let mut m_mat = vec![vec![NEG_INF; n + 1]; m + 1];
    let mut ix = vec![vec![NEG_INF; n + 1]; m + 1]; // insertion in reference (gap in ref): consumes read
    let mut iy = vec![vec![NEG_INF; n + 1]; m + 1]; // deletion in reference (gap in read): consumes ref

    // Traceback pointers: 0=M,1=Ix,2=Iy for each matrix; store predecessor matrix id.
    let mut tb_m = vec![vec![0u8; n + 1]; m + 1];
    let mut tb_ix = vec![vec![0u8; n + 1]; m + 1];
    let mut tb_iy = vec![vec![0u8; n + 1]; m + 1];

    // Semi-global on reference: allow starting anywhere in reference at no cost.
    for j in 0..=n {
        m_mat[0][j] = 0;
        ix[0][j] = NEG_INF;
        iy[0][j] = NEG_INF;
    }
    // Allow leading insertions in reference (gaps in ref) with penalties.
    m_mat[0][0] = 0;
    for i in 1..=m {
        ix[i][0] = scores.gap_open + (i as i32 - 1) * scores.gap_extend;
        tb_ix[i][0] = 1;
        m_mat[i][0] = NEG_INF;
        iy[i][0] = NEG_INF;
    }

    for i in 1..=m {
        let j_min = i.saturating_sub(band).max(1);
        let j_max = (i + band).min(n);
        for j in j_min..=j_max {
            // Update Ix (gap in ref): from M or Ix at (i-1, j)
            let from_m = m_mat[i - 1][j] + scores.gap_open;
            let from_ix = ix[i - 1][j] + scores.gap_extend;
            if from_m >= from_ix {
                ix[i][j] = from_m;
                tb_ix[i][j] = 0;
            } else {
                ix[i][j] = from_ix;
                tb_ix[i][j] = 1;
            }

            // Update Iy (gap in read): from M or Iy at (i, j-1)
            let from_m = m_mat[i][j - 1] + scores.gap_open;
            let from_iy = iy[i][j - 1] + scores.gap_extend;
            if from_m >= from_iy {
                iy[i][j] = from_m;
                tb_iy[i][j] = 0;
            } else {
                iy[i][j] = from_iy;
                tb_iy[i][j] = 2;
            }

            // Update M: from best of (M, Ix, Iy) at (i-1, j-1) + match/mismatch
            let a = read[i - 1];
            let b = reference_window[j - 1];
            let s = if a == b {
                scores.match_score
            } else {
                scores.mismatch_score
            };

            let (best_prev, prev_id) =
                best3(m_mat[i - 1][j - 1], ix[i - 1][j - 1], iy[i - 1][j - 1]);
            m_mat[i][j] = best_prev + s;
            tb_m[i][j] = prev_id;
        }
    }

    // Choose best end column with free end on reference.
    let mut best_score = NEG_INF;
    let mut best_j = 0usize;
    let mut best_matrix = 0u8;
    for j in 0..=n {
        let (score, which) = best3(m_mat[m][j], ix[m][j], iy[m][j]);
        if score > best_score {
            best_score = score;
            best_j = j;
            best_matrix = which;
        }
    }
    if best_score == NEG_INF {
        return None;
    }

    // Traceback to build ops.
    let mut ops_rev: Vec<(CigarOpKind, u32)> = Vec::new();
    let mut i = m;
    let mut j = best_j;
    let mut matrix = best_matrix;

    while i > 0 {
        match matrix {
            0 => {
                // M matrix consumes both
                push_op(&mut ops_rev, CigarOpKind::Match, 1);
                let prev = tb_m[i][j];
                i -= 1;
                j = j.saturating_sub(1);
                matrix = prev;
            }
            1 => {
                // Ix consumes read
                push_op(&mut ops_rev, CigarOpKind::Insertion, 1);
                let prev = tb_ix[i][j];
                i -= 1;
                matrix = prev;
            }
            2 => {
                // Iy consumes reference
                push_op(&mut ops_rev, CigarOpKind::Deletion, 1);
                let prev = tb_iy[i][j];
                j = j.saturating_sub(1);
                matrix = prev;
            }
            _ => break,
        }
        if j == 0 && i == 0 {
            break;
        }
    }

    ops_rev.reverse();
    let cigar: Vec<CigarOp> = ops_rev
        .into_iter()
        .map(|(kind, len)| CigarOp { kind, len })
        .collect();

    let ref_start = compute_ref_start(reference_window, read, &cigar, best_j)?;
    let (nm, md) = compute_nm_md(reference_window, read, ref_start, &cigar);

    Some(RefinedAlignment {
        ref_start,
        cigar,
        score: best_score,
        nm,
        md,
    })
}

#[inline]
fn best3(a: i32, b: i32, c: i32) -> (i32, u8) {
    if a >= b && a >= c {
        (a, 0)
    } else if b >= c {
        (b, 1)
    } else {
        (c, 2)
    }
}

fn push_op(out: &mut Vec<(CigarOpKind, u32)>, kind: CigarOpKind, len: u32) {
    if let Some((last_kind, last_len)) = out.last_mut() {
        if *last_kind == kind {
            *last_len += len;
            return;
        }
    }
    out.push((kind, len));
}

fn compute_ref_start(
    reference_window: &[u8],
    read: &[u8],
    cigar: &[CigarOp],
    ref_end_j: usize,
) -> Option<usize> {
    // We aligned the full read and ended at reference column ref_end_j (in DP coordinates).
    // Walk backwards through the CIGAR to find how many reference bases were consumed.
    let mut ref_consumed = 0usize;
    let mut read_consumed = 0usize;
    for op in cigar {
        match op.kind {
            CigarOpKind::Match => {
                ref_consumed += op.len as usize;
                read_consumed += op.len as usize;
            }
            CigarOpKind::Insertion => {
                read_consumed += op.len as usize;
            }
            CigarOpKind::Deletion => {
                ref_consumed += op.len as usize;
            }
            _ => {}
        }
    }
    if read_consumed != read.len() {
        return None;
    }
    if ref_consumed > ref_end_j {
        return None;
    }
    let start = ref_end_j - ref_consumed;
    if start <= reference_window.len() {
        Some(start)
    } else {
        None
    }
}

fn compute_nm_md(
    reference_window: &[u8],
    read: &[u8],
    ref_start: usize,
    cigar: &[CigarOp],
) -> (u32, String) {
    let mut nm: u32 = 0;
    let mut md = String::new();
    let mut match_run: u32 = 0;

    let mut r_i = 0usize;
    let mut ref_i = ref_start;

    for op in cigar {
        match op.kind {
            CigarOpKind::Match => {
                for _ in 0..op.len {
                    let rb = read[r_i];
                    let gb = reference_window[ref_i];
                    if rb == gb {
                        match_run += 1;
                    } else {
                        // mismatch
                        nm += 1;
                        md.push_str(&match_run.to_string());
                        match_run = 0;
                        md.push(gb as char);
                    }
                    r_i += 1;
                    ref_i += 1;
                }
            }
            CigarOpKind::Insertion => {
                // insertion contributes to NM but not MD
                nm += op.len;
                r_i += op.len as usize;
            }
            CigarOpKind::Deletion => {
                nm += op.len;
                md.push_str(&match_run.to_string());
                match_run = 0;
                md.push('^');
                for _ in 0..op.len {
                    md.push(reference_window[ref_i] as char);
                    ref_i += 1;
                }
            }
            _ => {}
        }
    }
    md.push_str(&match_run.to_string());
    (nm, md)
}
