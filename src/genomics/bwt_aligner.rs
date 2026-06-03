use std::sync::Arc;

use crate::genomics::alignment::{banded_affine_align, AffineScores};
use crate::genomics::{BlockedFMIndex, CigarOp, FMInterval};
use thiserror::Error;

/// Result of aligning a single read.
#[derive(Debug, Clone)]
pub struct AlignmentResult {
    /// FM interval representing candidate suffixes matching the read.
    pub interval: FMInterval,
    /// Primary mapping position (0-based reference coordinate), if a candidate was found.
    pub primary_position: Option<u32>,
    /// Whether the primary mapping is on the reverse-complement strand.
    pub is_reverse: bool,
    /// CIGAR operations for the refined alignment (if refinement succeeded).
    pub cigar: Vec<CigarOp>,
    /// NM tag value for the refined alignment.
    pub nm: u32,
    /// AS tag value for the refined alignment.
    pub as_score: i32,
    /// MD string for the refined alignment.
    pub md: String,
    /// Mapping quality estimate (placeholder; will be calibrated later).
    pub mapq: u8,
    /// Number of bases from the read that participated in the search.
    pub processed_bases: usize,
    /// Count of mismatches or aborted iterations.
    pub mismatches: usize,
    /// Span of block identifiers that contributed to the result.
    pub span: (usize, usize),
    /// Heuristic score used for ranking alignments.
    pub score: f32,
}

impl AlignmentResult {
    /// Returns `true` when there are remaining candidate suffixes.
    pub fn has_candidates(&self) -> bool {
        self.primary_position.is_some() || !self.interval.is_empty()
    }
}

#[derive(Debug, Clone, Copy)]
struct SeedConfig {
    k: usize,
    stride: usize,
    max_hits_per_seed: usize,
    diagonal_tolerance: i32,
}

/// Errors surfaced by the aligner.
#[derive(Debug, Error)]
pub enum AlignerError {
    /// Failure while constructing or querying the FM index.
    #[error("fm-index error: {0}")]
    FMIndex(#[from] crate::genomics::FMIndexError),
}

/// FM-index backed aligner (seeding + simple deterministic chaining).
#[derive(Debug)]
pub struct BWTAligner {
    fm_index: Arc<BlockedFMIndex>,
    reference: Arc<[u8]>,
    seed_cfg: SeedConfig,
}

impl BWTAligner {
    /// Construct a new aligner from the reference genome.
    pub fn new(reference: &[u8]) -> Result<Self, AlignerError> {
        let suggested_block_size = ((reference.len() as f64).sqrt().ceil() as usize).max(64);
        let fm_index = Arc::new(BlockedFMIndex::build(reference, suggested_block_size)?);
        let reference: Arc<[u8]> = Arc::from(
            reference
                .iter()
                .map(|b| b.to_ascii_uppercase())
                .collect::<Vec<u8>>()
                .into_boxed_slice(),
        );

        Ok(Self {
            fm_index,
            reference,
            seed_cfg: SeedConfig {
                k: 19,
                stride: 5,
                max_hits_per_seed: 64,
                diagonal_tolerance: 3,
            },
        })
    }

    /// Align a read and return a deterministic candidate mapping position (if found).
    ///
    /// Current behavior:
    /// - Uses exact-match FM-index seeding (no mismatches/indels yet).
    /// - Produces candidate mapping by diagonal clustering of seed hits.
    /// - Does not yet compute CIGAR/MAPQ (next milestone).
    pub fn align_read(&mut self, read: &[u8]) -> Result<AlignmentResult, AlignerError> {
        let forward = self.align_read_one_strand(read, false);
        let reverse = {
            let rc = reverse_complement(read);
            self.align_read_one_strand(&rc, true)
        };

        Ok(pick_best(forward, reverse))
    }

    /// Align multiple reads, returning a vector of results.
    pub fn align_batch<I>(&mut self, reads: I) -> Result<Vec<AlignmentResult>, AlignerError>
    where
        I: IntoIterator,
        I::Item: AsRef<[u8]>,
    {
        let mut results = Vec::new();
        for read in reads {
            results.push(self.align_read(read.as_ref())?);
        }
        Ok(results)
    }

    /// Reference to the underlying FM-index.
    pub fn fm_index(&self) -> &BlockedFMIndex {
        &self.fm_index
    }

    fn align_read_one_strand(&self, read: &[u8], is_reverse: bool) -> AlignmentResult {
        let read_upper: Vec<u8> = read.iter().map(|b| b.to_ascii_uppercase()).collect();

        // Compute full-read interval (exact match). This preserves existing test expectations.
        let interval = self.fm_index.backward_search(&read_upper);

        // Seed-based candidate generation for longer reads.
        let mut hits = Vec::new();
        if read_upper.len() >= self.seed_cfg.k {
            for read_pos in (0..=read_upper.len() - self.seed_cfg.k).step_by(self.seed_cfg.stride) {
                let kmer = &read_upper[read_pos..read_pos + self.seed_cfg.k];
                if kmer.contains(&b'N') {
                    continue;
                }
                let k_interval = self.fm_index.backward_search(kmer);
                if k_interval.is_empty() {
                    continue;
                }
                let locs = self
                    .fm_index
                    .locate_interval(k_interval, self.seed_cfg.max_hits_per_seed);
                for ref_pos in locs {
                    hits.push(SeedHit {
                        ref_pos: ref_pos as i64,
                        read_pos: read_pos as i64,
                        len: self.seed_cfg.k as i32,
                    });
                }
            }
        }

        let candidate_start = if hits.is_empty() {
            if interval.is_empty() {
                None
            } else {
                Some(self.fm_index.sa_at(interval.lower as usize))
            }
        } else {
            chain_hits(&hits, self.seed_cfg.diagonal_tolerance).map(|c| c as usize)
        };

        // Refine candidate with banded affine-gap DP.
        let mut primary_position: Option<u32> = None;
        let mut cigar: Vec<CigarOp> = Vec::new();
        let mut nm: u32 = 0;
        let mut as_score: i32 = 0;
        let mut md = String::new();

        if let Some(start) = candidate_start {
            let band = 32usize;
            let window_start = start.saturating_sub(band);
            let window_end = (start + read_upper.len() + band).min(self.reference.len());
            if window_start < window_end {
                let window = &self.reference[window_start..window_end];
                if let Some(refined) =
                    banded_affine_align(window, &read_upper, band, AffineScores::default())
                {
                    primary_position = Some((window_start + refined.ref_start) as u32);
                    cigar = refined.cigar;
                    nm = refined.nm;
                    as_score = refined.score;
                    md = refined.md;
                }
            }
        }

        let score = if primary_position.is_some() {
            as_score as f32
        } else if interval.is_empty() {
            0.0
        } else {
            1.0
        };

        AlignmentResult {
            interval,
            primary_position,
            is_reverse,
            cigar,
            nm,
            as_score,
            md,
            mapq: if primary_position.is_some() { 60 } else { 0 },
            processed_bases: read_upper.len(),
            mismatches: nm as usize,
            span: (0, 0),
            score,
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct SeedHit {
    ref_pos: i64,
    read_pos: i64,
    len: i32,
}

fn chain_hits(hits: &[SeedHit], diag_tol: i32) -> Option<i64> {
    // Group by diagonal d = ref_pos - read_pos. This is a simple deterministic
    // chaining heuristic appropriate for exact-match seeding.
    let mut sorted = hits.to_vec();
    sorted.sort_by(|a, b| {
        (a.ref_pos - a.read_pos)
            .cmp(&(b.ref_pos - b.read_pos))
            .then_with(|| a.read_pos.cmp(&b.read_pos))
            .then_with(|| a.ref_pos.cmp(&b.ref_pos))
    });

    let mut best_score: i64 = -1;
    let mut best_diag: i64 = 0;

    let mut i = 0usize;
    while i < sorted.len() {
        let d0 = sorted[i].ref_pos - sorted[i].read_pos;
        let mut j = i + 1;
        let mut score = sorted[i].len as i64;
        while j < sorted.len() {
            let dj = sorted[j].ref_pos - sorted[j].read_pos;
            if (dj - d0).abs() > diag_tol as i64 {
                break;
            }
            score += sorted[j].len as i64;
            j += 1;
        }
        if score > best_score || (score == best_score && d0 < best_diag) {
            best_score = score;
            best_diag = d0;
        }
        i = j;
    }

    if best_score <= 0 {
        return None;
    }
    if best_diag < 0 {
        return None;
    }
    Some(best_diag)
}

fn pick_best(a: AlignmentResult, b: AlignmentResult) -> AlignmentResult {
    // Deterministic tie-break: higher score, then prefer forward strand, then lower position.
    match (a.primary_position, b.primary_position) {
        (Some(pa), Some(pb)) => {
            if a.score > b.score {
                a
            } else if b.score > a.score {
                b
            } else if !a.is_reverse && b.is_reverse {
                a
            } else if a.is_reverse && !b.is_reverse {
                b
            } else if pa <= pb {
                a
            } else {
                b
            }
        }
        (Some(_), None) => a,
        (None, Some(_)) => b,
        (None, None) => {
            // Fall back to exact interval presence.
            if !a.interval.is_empty() && b.interval.is_empty() {
                a
            } else if a.interval.is_empty() && !b.interval.is_empty() {
                b
            } else {
                a
            }
        }
    }
}

fn reverse_complement(read: &[u8]) -> Vec<u8> {
    read.iter()
        .rev()
        .map(|&b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn aligner_returns_interval_for_read() {
        let reference = b"ACGTCGTAACGT";
        let mut aligner = BWTAligner::new(reference).expect("builder should succeed");
        let result = aligner
            .align_read(b"CGTA")
            .expect("alignment should succeed");
        assert!(result.has_candidates());
        assert!(result.processed_bases > 0);
    }
}
