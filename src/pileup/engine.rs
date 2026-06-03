//! The streaming pileup engine: CIGAR-aware, read-filtered, strand-aware,
//! bounded-memory. Yields one `PileupColumn` per covered reference position.

use std::collections::HashMap;
use std::ops::Range;
use std::sync::Arc;

use crate::core::{
    allele_index, AlignedRead, CoreError, Locus, Position, WorkingSet, PILEUP_ENGINE_OVERHEAD,
    PILEUP_MAP_BYTES_PER_BASE, PILEUP_PER_READ_OVERHEAD,
};
use crate::pileup::column::{Obs, PileupColumn};
use crate::pileup::source::ReadSource;

/// Fixed-seed FNV-1a-64 over a read's identity (position, end, flags, CIGAR, seq,
/// qual). Deterministic across runs and processes (NOT `DefaultHasher`, whose seed
/// is per-process). The hash is uncorrelated with the allele a read carries at any
/// single column — the basis of the unbiased depth-cap sample.
fn read_priority(read: &AlignedRead) -> u64 {
    const FNV_OFFSET: u64 = 0xcbf2_9ce4_8422_2325;
    const FNV_PRIME: u64 = 0x0000_0100_0000_01b3;
    #[inline]
    fn fold(mut h: u64, bytes: &[u8]) -> u64 {
        for &b in bytes {
            h ^= b as u64;
            h = h.wrapping_mul(FNV_PRIME);
        }
        h
    }
    let mut h = FNV_OFFSET;
    h = fold(h, &read.pos.0.to_le_bytes());
    h = fold(h, &read.end().to_le_bytes());
    h = fold(h, &read.flags.0.to_le_bytes());
    for op in &read.cigar {
        h = fold(h, &[op.kind as u8]);
        h = fold(h, &op.len.to_le_bytes());
    }
    h = fold(h, &read.seq);
    h = fold(h, &read.qual);
    h
}

/// Read-level filters applied as reads enter the pileup.
#[derive(Debug, Clone)]
pub struct PileupParams {
    /// Minimum mapping quality; reads below this are skipped.
    pub min_mapq: u8,
    /// Minimum base quality; observations below this are dropped.
    pub min_base_qual: u8,
    /// Skip secondary alignments (SAM flag 0x100).
    pub skip_secondary: bool,
    /// Skip supplementary alignments (SAM flag 0x800).
    pub skip_supplementary: bool,
    /// Skip PCR/optical duplicates (SAM flag 0x400).
    pub skip_duplicate: bool,
    /// Cap on the active read set per position (deterministic downsampling).
    /// `None` = uncapped (default). When `Some(d)`, an unbiased min-hash reservoir
    /// keeps `d` reads covering each position; reads removed are counted
    /// `over_max_depth`.
    pub max_depth: Option<u32>,
    /// When `Some(m)` (set only under `--enforce`), a read whose `seq.len()`
    /// exceeds `m` aborts the run — the predicted memory envelope assumes `<= m`,
    /// so a longer read would silently void it. `None` = no check (default).
    pub max_read_len: Option<u32>,
}

impl Default for PileupParams {
    fn default() -> Self {
        Self {
            min_mapq: 0,
            min_base_qual: 0,
            skip_secondary: true,
            skip_supplementary: true,
            skip_duplicate: true,
            max_depth: None,
            max_read_len: None,
        }
    }
}

/// Counts of reads skipped during pileup, by reason (surfaced to callers/CLI).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SkipCounts {
    /// Reads with the unmapped flag set.
    pub unmapped: u64,
    /// Reads on a contig other than the target.
    pub wrong_contig: u64,
    /// Secondary alignments skipped.
    pub secondary: u64,
    /// Supplementary alignments skipped.
    pub supplementary: u64,
    /// Duplicate reads skipped.
    pub duplicate: u64,
    /// Reads below the MAPQ threshold.
    pub low_mapq: u64,
    /// Reads dropped because the position was already at `max_depth`.
    pub over_max_depth: u64,
}

impl SkipCounts {
    /// Total reads skipped across all reasons.
    pub fn total(&self) -> u64 {
        self.unmapped
            + self.wrong_contig
            + self.secondary
            + self.supplementary
            + self.duplicate
            + self.low_mapq
            + self.over_max_depth
    }

    /// Field-wise add (summing per-contig engines on the whole-genome path).
    pub fn accumulate(&mut self, other: &SkipCounts) {
        self.unmapped += other.unmapped;
        self.wrong_contig += other.wrong_contig;
        self.secondary += other.secondary;
        self.supplementary += other.supplementary;
        self.duplicate += other.duplicate;
        self.low_mapq += other.low_mapq;
        self.over_max_depth += other.over_max_depth;
    }
}

/// A read currently overlapping the cursor, with its CIGAR projection precomputed.
#[derive(Debug)]
struct ActiveRead {
    /// Half-open reference end (CIGAR-derived) — used to expire the read.
    end: u32,
    /// Map from reference position to read offset for this read's Match bases.
    ref_to_read: HashMap<u32, usize>,
    /// Read sequence (forward-reference orientation).
    seq: Arc<[u8]>,
    /// Per-base qualities (parallel to `seq`).
    qual: Arc<[u8]>,
    /// Mapping quality.
    mapq: u8,
    /// Reverse-strand flag (metadata only — never applied to `seq`).
    reverse: bool,
    /// Fixed-seed content hash — the depth-cap reservoir admission/eviction key.
    priority: u64,
}

/// Streaming, CIGAR-aware, bounded-memory pileup over one contig region.
///
/// Yields one [`PileupColumn`] per covered reference position. The working set
/// is bounded by local read coverage, not by input size.
#[derive(Debug)]
pub struct PileupEngine<S: ReadSource> {
    source: S,
    reference: Arc<[u8]>,
    contig: u32,
    region: Range<u32>,
    params: PileupParams,
    active: Vec<ActiveRead>,
    next_read: Option<AlignedRead>,
    pos: u32,
    skips: SkipCounts,
    source_done: bool,
}

impl<S: ReadSource> PileupEngine<S> {
    /// Create an engine over `reference` bytes for `contig`, covering `region`
    /// (0-based half-open; `reference[0]` is the base at `region.start`).
    pub fn new(
        source: S,
        reference: Arc<[u8]>,
        contig: u32,
        region: Range<u32>,
        params: PileupParams,
    ) -> Self {
        let pos = region.start;
        Self {
            source,
            reference,
            contig,
            region,
            params,
            active: Vec::new(),
            next_read: None,
            pos,
            skips: SkipCounts::default(),
            source_done: false,
        }
    }

    /// Reads skipped so far, by reason. Final after iteration completes.
    pub fn skip_counts(&self) -> SkipCounts {
        self.skips
    }

    /// Current working-set estimate: bounded by the active read set (local
    /// coverage), independent of total input size. Foundation for `rosalind plan`.
    pub fn current_working_set(&self) -> WorkingSet {
        // The decoded reference for this contig is resident in the engine.
        let reference_bytes = self.reference.len() as u64;
        // Each active read holds its projection map (16 B/entry) plus its seq and
        // qual byte buffers; count all three (the map alone is a large undercount,
        // especially for long reads).
        let active_bytes: u64 = self
            .active
            .iter()
            .map(|r| {
                (r.ref_to_read.len() as u64) * PILEUP_MAP_BYTES_PER_BASE
                    + r.seq.len() as u64
                    + r.qual.len() as u64
                    + PILEUP_PER_READ_OVERHEAD
            })
            .sum();
        WorkingSet {
            bytes: reference_bytes + active_bytes + PILEUP_ENGINE_OVERHEAD,
        }
    }

    /// Whether the read passes flag/MAPQ filters. (Contig routing and the
    /// unmapped flag are handled in `advance_to`.)
    fn passes_filters(&mut self, read: &AlignedRead) -> bool {
        let f = read.flags;
        if self.params.skip_secondary && f.is_secondary() {
            self.skips.secondary += 1;
            return false;
        }
        if self.params.skip_supplementary && f.is_supplementary() {
            self.skips.supplementary += 1;
            return false;
        }
        if self.params.skip_duplicate && f.is_duplicate() {
            self.skips.duplicate += 1;
            return false;
        }
        if read.mapq < self.params.min_mapq {
            self.skips.low_mapq += 1;
            return false;
        }
        true
    }

    /// Precompute a read's reference→read-offset map and add it to the active set.
    /// `priority` is the precomputed reservoir key (see `read_priority`).
    fn ingest(&mut self, read: AlignedRead, priority: u64) {
        let end = read.end();
        let mut ref_to_read = HashMap::new();
        for rb in read.projected_bases() {
            ref_to_read.insert(rb.ref_pos, rb.read_offset);
        }
        self.active.push(ActiveRead {
            end,
            ref_to_read,
            seq: Arc::clone(&read.seq),
            qual: Arc::clone(&read.qual),
            mapq: read.mapq,
            reverse: read.flags.is_reverse(),
            priority,
        });
    }

    /// Expire reads that no longer cover `pos`, then pull in reads starting at or
    /// before `pos`. Reads are coordinate-sorted, so we stop at the first read
    /// that starts after `pos` on the target contig.
    fn advance_to(&mut self, pos: u32) -> Result<(), CoreError> {
        self.active.retain(|r| r.end > pos);
        loop {
            if self.next_read.is_none() && !self.source_done {
                match self.source.next_read()? {
                    Some(r) => self.next_read = Some(r),
                    None => self.source_done = true,
                }
            }
            // Peek the routing fields without holding a borrow across `take`.
            let (rc, rp, unmapped) = match self.next_read.as_ref() {
                Some(r) => (r.contig, r.pos.0, r.flags.is_unmapped()),
                None => break,
            };
            if unmapped {
                self.next_read.take();
                self.skips.unmapped += 1;
                continue;
            }
            match rc.cmp(&self.contig) {
                std::cmp::Ordering::Greater => break, // sorted: no more target-contig reads
                std::cmp::Ordering::Less => {
                    self.next_read.take();
                    self.skips.wrong_contig += 1;
                }
                std::cmp::Ordering::Equal => {
                    if rp > pos {
                        break; // future read on our contig
                    }
                    let read = self.next_read.take().unwrap();
                    if !self.passes_filters(&read) {
                        continue;
                    }
                    if read.end() <= pos {
                        continue; // does not reach the cursor
                    }
                    // --enforce precondition: a read longer than the declared
                    // --max-read-len voids the predicted envelope. Fail loud here
                    // rather than silently over-allocate (the would-be OOM path).
                    if let Some(maxlen) = self.params.max_read_len {
                        if read.seq.len() as u32 > maxlen {
                            return Err(CoreError::ReadExceedsDeclaredLength {
                                len: read.seq.len() as u32,
                                declared: maxlen,
                            });
                        }
                    }
                    let prio = read_priority(&read);
                    if let Some(max) = self.params.max_depth {
                        if self.active.len() as u32 >= max {
                            // Unbiased min-hash reservoir: keep the `max`
                            // smallest-priority reads covering the cursor. Evict the
                            // greatest-priority resident iff the newcomer ranks below
                            // it; otherwise refuse. Either way the cap removed one
                            // read (counted). Priority is a fixed-seed content hash,
                            // so selection is independent of start position / allele
                            // (the leftmost-arrival bias that silently dropped deep
                            // variants is gone). The active set stays bounded by
                            // `max`, preserving the working-set guarantee.
                            let mut worst = 0usize;
                            for i in 1..self.active.len() {
                                if self.active[i].priority > self.active[worst].priority {
                                    worst = i;
                                }
                            }
                            self.skips.over_max_depth += 1;
                            if prio < self.active[worst].priority {
                                self.active.swap_remove(worst);
                                self.ingest(read, prio);
                            }
                            continue;
                        }
                    }
                    self.ingest(read, prio);
                }
            }
        }
        Ok(())
    }

    /// Build the column at the current cursor from the active set.
    fn build_column(&self) -> PileupColumn {
        // The cursor starts at region.start and only increments, so this
        // subtraction cannot underflow.
        debug_assert!(self.pos >= self.region.start);
        let ref_idx = (self.pos - self.region.start) as usize;
        let ref_base = self.reference.get(ref_idx).copied().unwrap_or(b'N');
        let mut obs = Vec::new();
        let mut raw_depth: u32 = 0;
        for r in &self.active {
            if let Some(&off) = r.ref_to_read.get(&self.pos) {
                raw_depth += 1;
                // Forward orientation — read the stored SEQ byte directly. No
                // complement (this is the reverse-strand fix). `off` is CIGAR-
                // derived and not bounds-checked against seq.len() (per the
                // RefBase contract), so a malformed read folds to N, not a panic.
                let base = r.seq.get(off).copied().unwrap_or(b'N');
                let bq = r.qual.get(off).copied().unwrap_or(0);
                if bq < self.params.min_base_qual {
                    continue;
                }
                if let Some(allele) = allele_index(base) {
                    obs.push(Obs {
                        allele: allele as u8,
                        base_qual: bq,
                        mapq: r.mapq,
                        reverse: r.reverse,
                    });
                }
            }
        }
        // Emit observations in a canonical total order so the downstream diploid
        // log-likelihood accumulation (germline.rs, an order-sensitive f64 sum) is
        // order-independent and the VCF stays byte-identical regardless of the
        // active set's internal order (which the depth-cap reservoir may permute).
        obs.sort_by(|a, b| {
            (a.allele, a.base_qual, a.mapq, a.reverse as u8).cmp(&(
                b.allele,
                b.base_qual,
                b.mapq,
                b.reverse as u8,
            ))
        });
        PileupColumn {
            locus: Locus {
                contig: self.contig,
                pos: Position(self.pos),
            },
            ref_base,
            raw_depth,
            obs,
        }
    }
}

impl<S: ReadSource> Iterator for PileupEngine<S> {
    type Item = Result<PileupColumn, CoreError>;

    fn next(&mut self) -> Option<Self::Item> {
        // Loop over empty positions (no recursion → bounded stack, any gap size).
        while self.pos < self.region.end {
            let pos = self.pos;
            if let Err(e) = self.advance_to(pos) {
                self.pos = self.region.end;
                return Some(Err(e));
            }
            let column = self.build_column();
            self.pos += 1;
            if !column.obs.is_empty() {
                return Some(Ok(column));
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::{CigarOp, CigarOpKind, SamFlags};
    use crate::pileup::source::SliceSource;

    // A fully-matched read at `pos` on contig 0 carrying `seq` (forward orientation).
    fn mread(pos: u32, seq: &[u8], reverse: bool) -> AlignedRead {
        let mut flags = SamFlags::default();
        if reverse {
            flags = SamFlags(SamFlags::REVERSE);
        }
        AlignedRead {
            contig: 0,
            pos: Position(pos),
            mapq: 60,
            flags,
            cigar: vec![CigarOp::new(CigarOpKind::Match, seq.len() as u32)],
            seq: Arc::from(seq.to_vec().into_boxed_slice()),
            qual: Arc::from(vec![30u8; seq.len()].into_boxed_slice()),
        }
    }

    fn engine(reads: Vec<AlignedRead>, reference: &[u8]) -> PileupEngine<SliceSource> {
        PileupEngine::new(
            SliceSource::new(reads),
            Arc::from(reference.to_vec().into_boxed_slice()),
            0,
            0..reference.len() as u32,
            PileupParams::default(),
        )
    }

    fn columns(e: PileupEngine<SliceSource>) -> Vec<PileupColumn> {
        let mut out = Vec::new();
        for c in e {
            out.push(c.expect("pileup column"));
        }
        out
    }

    #[test]
    fn obs_are_emitted_in_canonical_order() {
        // A column with mixed alleles must emit obs sorted by
        // (allele, base_qual, mapq, reverse) — independent of read arrival order,
        // so the downstream f64 likelihood sum is order-stable.
        let reference = b"AAAA";
        // Three reads covering pos 0 with different alleles at offset 0:
        // C (allele 1), A (allele 0, ref), G (allele 2). Arrival order C, A, G.
        let reads = vec![
            mread(0, b"C", false),
            mread(0, b"A", false),
            mread(0, b"G", false),
        ];
        let cols = columns(engine(reads, reference));
        let at0 = cols.iter().find(|c| c.locus.pos.0 == 0).unwrap();
        let alleles: Vec<u8> = at0.obs.iter().map(|o| o.allele).collect();
        // Canonical: ascending by allele (0, 1, 2) regardless of C, A, G arrival.
        assert_eq!(alleles, vec![0, 1, 2]);
    }

    #[test]
    fn read_exceeding_declared_max_read_len_errors() {
        // Under --enforce (max_read_len Some), a read longer than the declared cap
        // aborts the run rather than silently voiding the predicted envelope.
        let reference = b"AAAAAAAA";
        let params = PileupParams {
            max_read_len: Some(4),
            ..PileupParams::default()
        };
        let mut e = PileupEngine::new(
            SliceSource::new(vec![mread(0, b"CCCCCC", false)]), // len 6 > 4
            Arc::from(reference.to_vec().into_boxed_slice()),
            0,
            0..8,
            params,
        );
        let err = loop {
            match e.next() {
                Some(Ok(_)) => continue,
                Some(Err(err)) => break err,
                None => panic!("expected an error, got clean end"),
            }
        };
        assert!(matches!(
            err,
            crate::core::CoreError::ReadExceedsDeclaredLength {
                len: 6,
                declared: 4
            }
        ));
    }

    #[test]
    fn unbiased_cap_keeps_downstream_starting_reads() {
        // The biased (leftmost-arrival) cap dropped reads that START at a deep
        // variant, zeroing the alt allele. The min-hash reservoir keeps a
        // start-position-independent sample, so a variant carried only by reads
        // that begin AT the variant column survives. Reads must be DISTINCT
        // (real reads are) — identical content shares one hash and degenerates.
        let reference = vec![b'A'; 64];
        let v = 32u32; // variant column
        let cap = 30u32;
        let mut reads = Vec::new();
        // 30 ref reads starting upstream (pos 0), distinct lengths -> distinct
        // hashes, all carrying A (ref) at v.
        for i in 0..30u32 {
            let len = 33 + i as usize; // covers v=32 (len>32); 33..62 < 64
            reads.push(mread(0, &vec![b'A'; len], false));
        }
        // 30 alt reads starting AT v, distinct lengths, carrying C at v (offset 0).
        for j in 0..30u32 {
            let len = 1 + j as usize; // start v, covers v
            let mut seq = vec![b'A'; len];
            seq[0] = b'C';
            reads.push(mread(v, &seq, false));
        }
        let params = PileupParams {
            max_depth: Some(cap),
            ..PileupParams::default()
        };
        let e = PileupEngine::new(
            SliceSource::new(reads),
            Arc::from(reference.into_boxed_slice()),
            0,
            0..64,
            params,
        );
        let mut at_v = None;
        for c in e {
            let col = c.unwrap();
            if col.locus.pos.0 == v {
                at_v = Some(col);
            }
        }
        let counts = at_v.expect("a column at v").allele_counts();
        // Both alleles present: ref (A=0) AND alt (C=1) — the biased cap gave alt==0.
        assert!(counts[1] > 0, "alt allele must survive the cap: {counts:?}");
        assert!(counts[0] > 0, "ref allele must remain too: {counts:?}");
        // Capped: total observed at v <= cap.
        assert!(counts.iter().sum::<u32>() <= cap, "depth must be capped");
    }

    #[test]
    fn reverse_strand_reads_are_not_complemented() {
        // Regression: forward and reverse-strand reads carrying the SAME forward
        // SEQ must contribute the SAME allele (legacy code complemented reverse reads).
        let reference = b"AAAAA";
        let fwd = mread(0, b"AAGAA", false);
        let rev = mread(0, b"AAGAA", true);
        let cols = columns(engine(vec![fwd, rev], reference));
        let at2 = cols
            .iter()
            .find(|c| c.locus.pos.0 == 2)
            .expect("column at pos 2");
        // Both observe G (allele 2); none observe C (allele 1, the complement of G).
        assert_eq!(at2.allele_counts(), [0, 0, 2, 0]);
    }

    #[test]
    fn basic_ungapped_pileup_counts() {
        let reference = b"ACGTACGT";
        let reads = vec![mread(0, b"ACGT", false), mread(2, b"GTAC", false)];
        let cols = columns(engine(reads, reference));
        let at2 = cols.iter().find(|c| c.locus.pos.0 == 2).unwrap();
        assert_eq!(at2.depth(), 2);
        assert_eq!(at2.ref_base, b'G');
        assert_eq!(at2.allele_counts(), [0, 0, 2, 0]);
        assert!(cols.iter().all(|c| c.depth() > 0));
    }

    #[test]
    fn sparse_coverage_skips_empty_positions_without_recursion() {
        // A huge gap between two reads must not overflow the stack (loop, not recursion).
        let mut reference = vec![b'A'; 100_000];
        reference[0] = b'C';
        reference[99_999] = b'C';
        let reads = vec![mread(0, b"C", false), mread(99_999, b"C", false)];
        let cols = columns(engine(reads, &reference));
        assert_eq!(cols.len(), 2);
        assert_eq!(cols[0].locus.pos.0, 0);
        assert_eq!(cols[1].locus.pos.0, 99_999);
    }

    #[test]
    fn default_params_skip_noise_and_keep_quality_open() {
        let p = PileupParams::default();
        assert!(p.skip_secondary && p.skip_supplementary && p.skip_duplicate);
        assert_eq!(p.min_mapq, 0);
        assert_eq!(p.min_base_qual, 0);
    }

    #[test]
    fn skip_counts_accumulate_sums_fieldwise() {
        let mut a = SkipCounts {
            over_max_depth: 2,
            low_mapq: 1,
            ..SkipCounts::default()
        };
        let b = SkipCounts {
            over_max_depth: 3,
            duplicate: 4,
            ..SkipCounts::default()
        };
        a.accumulate(&b);
        assert_eq!(a.over_max_depth, 5);
        assert_eq!(a.low_mapq, 1);
        assert_eq!(a.duplicate, 4);
    }

    #[test]
    fn skip_counts_total_sums_all_reasons() {
        let s = SkipCounts {
            unmapped: 1,
            wrong_contig: 2,
            secondary: 3,
            supplementary: 4,
            duplicate: 5,
            low_mapq: 6,
            over_max_depth: 0,
        };
        assert_eq!(s.total(), 21);
    }

    #[test]
    fn insertion_bases_do_not_shift_downstream_reference_positions() {
        // 2M 1I 2M at ref 10: offsets 0,1 -> ref 10,11; offset 2 = inserted (no ref);
        // offsets 3,4 -> ref 12,13.  seq "CCAGT": the inserted base is 'A' (offset 2).
        let reference = b"AAAAAAAAAAAAAAAA"; // 16 'A'
        let read = AlignedRead {
            contig: 0,
            pos: Position(10),
            mapq: 60,
            flags: SamFlags::default(),
            cigar: vec![
                CigarOp::new(CigarOpKind::Match, 2),
                CigarOp::new(CigarOpKind::Insertion, 1),
                CigarOp::new(CigarOpKind::Match, 2),
            ],
            seq: Arc::from(b"CCAGT".to_vec().into_boxed_slice()), // off: C C A(ins) G T
            qual: Arc::from(vec![30u8; 5].into_boxed_slice()),
        };
        let cols = columns(engine(vec![read], reference));
        let get = |p: u32| {
            cols.iter()
                .find(|c| c.locus.pos.0 == p)
                .map(|c| c.allele_counts())
        };
        assert_eq!(get(10), Some([0, 1, 0, 0])); // C (offset 0)
        assert_eq!(get(11), Some([0, 1, 0, 0])); // C (offset 1)
        assert_eq!(get(12), Some([0, 0, 1, 0])); // G (offset 3 — NOT the inserted 'A')
        assert_eq!(get(13), Some([0, 0, 0, 1])); // T (offset 4)
                                                 // The inserted 'A' (offset 2) appears at no reference position.
        assert!(cols.iter().all(|c| c.allele_counts()[0] == 0));
    }

    #[test]
    fn deletion_leaves_a_reference_gap_with_no_observation() {
        // 2M 1D 2M at ref 0: ref 0,1 observed; ref 2 deleted (no obs); ref 3,4 observed.
        let reference = b"AAAAAAAA";
        let read = AlignedRead {
            contig: 0,
            pos: Position(0),
            mapq: 60,
            flags: SamFlags::default(),
            cigar: vec![
                CigarOp::new(CigarOpKind::Match, 2),
                CigarOp::new(CigarOpKind::Deletion, 1),
                CigarOp::new(CigarOpKind::Match, 2),
            ],
            seq: Arc::from(b"GGGG".to_vec().into_boxed_slice()),
            qual: Arc::from(vec![30u8; 4].into_boxed_slice()),
        };
        let cols = columns(engine(vec![read], reference));
        let positions: Vec<u32> = cols.iter().map(|c| c.locus.pos.0).collect();
        assert_eq!(positions, vec![0, 1, 3, 4]); // ref 2 (deleted) emits no column
    }

    #[test]
    fn soft_clipped_bases_are_excluded() {
        // 2S 3M at ref 1: first 2 read bases clipped; only the 3 matched bases pile up.
        let reference = b"AAAAAAA";
        let read = AlignedRead {
            contig: 0,
            pos: Position(1),
            mapq: 60,
            flags: SamFlags::default(),
            cigar: vec![
                CigarOp::new(CigarOpKind::SoftClip, 2),
                CigarOp::new(CigarOpKind::Match, 3),
            ],
            seq: Arc::from(b"TTCGA".to_vec().into_boxed_slice()), // TT clipped; CGA -> ref 1,2,3
            qual: Arc::from(vec![30u8; 5].into_boxed_slice()),
        };
        let cols = columns(engine(vec![read], reference));
        let positions: Vec<u32> = cols.iter().map(|c| c.locus.pos.0).collect();
        assert_eq!(positions, vec![1, 2, 3]);
        let at1 = cols.iter().find(|c| c.locus.pos.0 == 1).unwrap();
        assert_eq!(at1.allele_counts(), [0, 1, 0, 0]); // 'C' (offset 2), not the clipped 'T'
    }

    #[test]
    fn long_read_piles_up_every_matched_base() {
        // Read-length-agnostic: a 5000-base full-match read covers 5000 positions.
        let reference = vec![b'A'; 6000];
        let seq = vec![b'C'; 5000];
        let read = AlignedRead {
            contig: 0,
            pos: Position(1000),
            mapq: 60,
            flags: SamFlags::default(),
            cigar: vec![CigarOp::new(CigarOpKind::Match, 5000)],
            seq: Arc::from(seq.into_boxed_slice()),
            qual: Arc::from(vec![30u8; 5000].into_boxed_slice()),
        };
        let cols = columns(engine(vec![read], &reference));
        assert_eq!(cols.len(), 5000);
        assert_eq!(cols.first().unwrap().locus.pos.0, 1000);
        assert_eq!(cols.last().unwrap().locus.pos.0, 5999);
        assert!(cols.iter().all(|c| c.allele_counts() == [0, 1, 0, 0]));
    }

    fn flagged_read(pos: u32, seq: &[u8], flag_bits: u16) -> AlignedRead {
        let mut r = mread(pos, seq, false);
        r.flags = SamFlags(flag_bits);
        r
    }

    #[test]
    fn filters_skip_secondary_supplementary_and_duplicate_reads() {
        let reference = b"AAAAA";
        let reads = vec![
            mread(0, b"CCCCC", false),                          // kept
            flagged_read(0, b"GGGGG", SamFlags::SECONDARY),     // skipped
            flagged_read(0, b"GGGGG", SamFlags::SUPPLEMENTARY), // skipped
            flagged_read(0, b"GGGGG", SamFlags::DUPLICATE),     // skipped
        ];
        let mut e = engine(reads, reference);
        let mut cols = Vec::new();
        for c in e.by_ref() {
            cols.push(c.unwrap());
        }
        // Only the kept read's 'C' (allele 1) appears at every position.
        assert!(cols.iter().all(|c| c.allele_counts() == [0, 1, 0, 0]));
        let s = e.skip_counts();
        assert_eq!((s.secondary, s.supplementary, s.duplicate), (1, 1, 1));
    }

    #[test]
    fn filters_skip_low_mapq_reads() {
        let reference = b"AAAAA";
        let mut low = mread(0, b"GGGGG", false);
        low.mapq = 3;
        let reads = vec![mread(0, b"CCCCC", false), low];
        let params = PileupParams {
            min_mapq: 10,
            ..PileupParams::default()
        };
        let mut e = PileupEngine::new(
            SliceSource::new(reads),
            Arc::from(reference.to_vec().into_boxed_slice()),
            0,
            0..5,
            params,
        );
        let mut cols = Vec::new();
        for c in e.by_ref() {
            cols.push(c.unwrap());
        }
        assert!(cols.iter().all(|c| c.allele_counts() == [0, 1, 0, 0]));
        assert_eq!(e.skip_counts().low_mapq, 1);
    }

    #[test]
    fn low_base_quality_observations_are_dropped() {
        let reference = b"AAAAA";
        let mut r = mread(0, b"GGGGG", false);
        r.qual = Arc::from(vec![2u8; 5].into_boxed_slice()); // below the floor
        let params = PileupParams {
            min_base_qual: 20,
            ..PileupParams::default()
        };
        let mut e = PileupEngine::new(
            SliceSource::new(vec![r]),
            Arc::from(reference.to_vec().into_boxed_slice()),
            0,
            0..5,
            params,
        );
        // All observations dropped → no columns emitted.
        assert!(e.next().is_none());
    }

    use crate::core::MemoryBudget;

    #[test]
    fn working_set_is_bounded_by_coverage_not_input_size() {
        // 50,000 single-base reads tiled across the reference at depth ~1: the
        // active set (and thus the working set) stays tiny throughout iteration.
        let reference = vec![b'A'; 50_000];
        let reads: Vec<AlignedRead> = (0..50_000u32).map(|p| mread(p, b"C", false)).collect();
        let mut e = PileupEngine::new(
            SliceSource::new(reads),
            Arc::from(reference.into_boxed_slice()),
            0,
            0..50_000,
            PileupParams::default(),
        );
        let budget = MemoryBudget::from_mb(1);
        let mut max_ws = 0u64;
        while let Some(c) = e.next() {
            let _ = c.unwrap();
            let ws = e.current_working_set();
            max_ws = max_ws.max(ws.bytes);
            assert!(ws.fits(budget), "working set {} exceeded 1 MiB", ws.bytes);
        }
        // Sanity: peak working set is far below holding all reads would cost.
        assert!(
            max_ws < 64 * 1024,
            "peak working set unexpectedly large: {max_ws}"
        );
    }

    #[test]
    fn shuffled_source_yields_identical_columns() {
        // SliceSource sorts on construction, so a shuffled input must produce the
        // exact same column sequence (deterministic output).
        let reference = b"ACGTACGTACGT";
        let make = |order: Vec<(u32, &'static [u8])>| {
            let reads: Vec<AlignedRead> =
                order.into_iter().map(|(p, s)| mread(p, s, false)).collect();
            columns(engine(reads, reference))
        };
        let a = make(vec![(0, b"AAAA"), (4, b"CCCC"), (8, b"GGGG")]);
        let b = make(vec![(8, b"GGGG"), (0, b"AAAA"), (4, b"CCCC")]);
        assert_eq!(a, b);
    }

    #[test]
    fn engine_is_an_open_substrate_for_arbitrary_per_locus_analysis() {
        // A non-caller consumer: compute per-position coverage directly from the
        // PileupColumn stream — no variant calling involved.
        let reference = b"AAAAAAAA";
        let reads = vec![mread(0, b"CCCC", false), mread(2, b"CCCC", false)];
        let mut coverage = Vec::new();
        let e = engine(reads, reference);
        for c in e {
            let col = c.unwrap();
            coverage.push((col.locus.pos.0, col.depth()));
        }
        // pos 0,1 depth 1; pos 2,3 depth 2; pos 4,5 depth 1.
        assert_eq!(
            coverage,
            vec![(0, 1), (1, 1), (2, 2), (3, 2), (4, 1), (5, 1)]
        );
    }

    #[test]
    fn raw_depth_counts_covering_reads_even_when_base_is_dropped() {
        // Two reads cover ref 0; one carries 'N' (not a callable allele), so it is
        // excluded from obs/depth() but still counted in raw_depth.
        let reference = b"AA";
        let reads = vec![mread(0, b"C", false), mread(0, b"N", false)];
        let cols = columns(engine(reads, reference));
        let at0 = cols.iter().find(|c| c.locus.pos.0 == 0).unwrap();
        assert_eq!(at0.depth(), 1); // only the callable 'C'
        assert_eq!(at0.allele_counts(), [0, 1, 0, 0]);
        assert_eq!(at0.raw_depth, 2); // both reads cover the position
    }

    #[test]
    fn params_default_max_depth_is_none_and_skipcounts_total_includes_over_max_depth() {
        assert_eq!(PileupParams::default().max_depth, None);
        let s = SkipCounts {
            unmapped: 1,
            wrong_contig: 2,
            secondary: 3,
            supplementary: 4,
            duplicate: 5,
            low_mapq: 6,
            over_max_depth: 7,
        };
        assert_eq!(s.total(), 28);
    }

    #[test]
    fn working_set_counts_reference_and_read_byte_buffers() {
        // One 4-base read fully covering a 10-base reference. After advancing to
        // pos 0 the active set holds that read; the working set must include the
        // reference bytes (10) AND the read's seq+qual buffers (4+4), not just the
        // projection map.
        let reference = b"ACGTACGTAC"; // 10 bytes
        let mut e = engine(vec![mread(0, b"ACGT", false)], reference);
        let first = e.next().expect("a column").expect("ok"); // drives advance_to(0)
        assert_eq!(first.locus.pos.0, 0);
        let ws = e.current_working_set().bytes;
        // reference (10) + map(4*16=64) + seq(4) + qual(4) + per-read(64) + fixed(256)
        // = 10 + 64 + 4 + 4 + 64 + 256 = 402.
        assert_eq!(ws, 402);
    }

    #[test]
    fn max_depth_caps_active_set_deterministically() {
        // 5 reads all covering pos 0..4; cap at 2. Only the first 2 (arrival order)
        // are kept; the other 3 are counted over_max_depth. Capped output is
        // identical regardless of input order (SliceSource sorts on construction).
        let reference = b"AAAA";
        let params = PileupParams {
            max_depth: Some(2),
            ..PileupParams::default()
        };
        let run = |reads: Vec<AlignedRead>| -> (Vec<u32>, u64) {
            let mut e = PileupEngine::new(
                SliceSource::new(reads),
                Arc::from(reference.to_vec().into_boxed_slice()),
                0,
                0..4,
                params.clone(),
            );
            let mut depths = Vec::new();
            for c in e.by_ref() {
                depths.push(c.unwrap().raw_depth);
            }
            (depths, e.skip_counts().over_max_depth)
        };
        let reads_a = vec![
            mread(0, b"CCCC", false),
            mread(0, b"CCCC", false),
            mread(0, b"CCCC", false),
            mread(0, b"CCCC", false),
            mread(0, b"CCCC", false),
        ];
        let mut reads_b = reads_a.clone();
        reads_b.reverse();
        let (depths_a, over_a) = run(reads_a);
        let (depths_b, over_b) = run(reads_b);
        // Capped: every position sees at most 2 reads.
        assert!(
            depths_a.iter().all(|&d| d <= 2),
            "raw depth must be capped at 2"
        );
        assert_eq!(over_a, 3, "3 of 5 reads dropped over max_depth");
        // Deterministic regardless of input order.
        assert_eq!(depths_a, depths_b);
        assert_eq!(over_a, over_b);
    }
}
