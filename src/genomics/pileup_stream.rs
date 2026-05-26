//! Streaming pileup over coordinate-sorted BAM.
//!
//! This provides a bounded-memory pileup engine: it maintains only the active
//! reads overlapping the current coordinate.

use std::path::Path;
use std::sync::Arc;

use anyhow::{anyhow, Context, Result};
use rust_htslib::bam;
use rust_htslib::bam::Read as BamRead;

use crate::genomics::PileupNode;

#[derive(Debug, Clone)]
struct ActiveRead {
    start: u32,
    end: u32,
    seq: Vec<u8>,
    qual: Vec<u8>,
    is_reverse: bool,
}

impl ActiveRead {
    fn base_and_qual_at(&self, position: u32) -> Option<(u8, u8)> {
        if position < self.start || position >= self.end {
            return None;
        }
        let offset = (position - self.start) as usize;
        if offset >= self.seq.len() {
            return None;
        }
        if self.is_reverse {
            let idx = self.seq.len() - 1 - offset;
            let base = complement(self.seq[idx]);
            let qual = *self.qual.get(idx).unwrap_or(&30);
            Some((base, qual))
        } else {
            let base = self.seq[offset];
            let qual = *self.qual.get(offset).unwrap_or(&30);
            Some((base, qual))
        }
    }
}

/// A streaming pileup over a single contig for a given region.
#[derive(Debug)]
pub struct BamPileupStream {
    reader: bam::Reader,
    header: bam::HeaderView,
    chrom: Arc<str>,
    region: std::ops::Range<u32>,
    next_record: Option<bam::Record>,
    active: Vec<ActiveRead>,
    current_pos: u32,
}

impl BamPileupStream {
    /// Create a pileup stream over a sorted BAM.
    pub fn new(
        bam_path: impl AsRef<Path>,
        chrom: Arc<str>,
        region: std::ops::Range<u32>,
    ) -> Result<Self> {
        let reader = bam::Reader::from_path(bam_path.as_ref())
            .with_context(|| format!("failed to open BAM {}", bam_path.as_ref().display()))?;
        let header = reader.header().to_owned();
        Ok(Self {
            reader,
            header,
            chrom,
            region: region.clone(),
            next_record: None,
            active: Vec::new(),
            current_pos: region.start,
        })
    }

    fn tid_matches(&self, tid: i32) -> bool {
        if tid < 0 {
            return false;
        }
        let name = self.header.tid2name(tid as u32);
        std::str::from_utf8(name)
            .map(|s| s == self.chrom.as_ref())
            .unwrap_or(false)
    }

    fn fetch_next_record(&mut self) -> Result<Option<bam::Record>> {
        let mut record = bam::Record::new();
        match self.reader.read(&mut record) {
            None => Ok(None),
            Some(Ok(())) => Ok(Some(record)),
            Some(Err(e)) => Err(anyhow!(e)),
        }
    }

    fn ensure_next_record_loaded(&mut self) -> Result<()> {
        if self.next_record.is_some() {
            return Ok(());
        }
        self.next_record = self.fetch_next_record()?;
        Ok(())
    }

    fn advance_active_to(&mut self, pos: u32) -> Result<()> {
        // Drop expired reads.
        self.active.retain(|r| r.end > pos);

        // Pull in new reads whose start <= pos.
        loop {
            self.ensure_next_record_loaded()?;
            let Some(rec) = self.next_record.take() else { break };

            if rec.is_unmapped() || !self.tid_matches(rec.tid()) {
                // Skip other contigs/unmapped.
                self.next_record = None;
                continue;
            }

            let start = rec.pos();
            if start < 0 {
                self.next_record = None;
                continue;
            }
            let start = start as u32;

            // Stop once we reach reads that start after the current pileup position.
            if start > pos {
                self.next_record = Some(rec);
                break;
            }

            // Only support simple match-style alignments in this streaming pileup for now.
            // (We will generalize to full CIGAR-aware pileup in the somatic caller work.)
            if !is_simple_match_cigar(&rec) {
                self.next_record = None;
                continue;
            }

            let seq = rec.seq().as_bytes();
            let seq: Vec<u8> = seq.iter().map(|b| b.to_ascii_uppercase()).collect();
            let qual: Vec<u8> = rec.qual().to_vec();
            let end = start + seq.len() as u32;

            self.active.push(ActiveRead {
                start,
                end,
                seq,
                qual,
                is_reverse: rec.is_reverse(),
            });

            // Maintain deterministic iteration order in the active set.
            self.active.sort_by(|a, b| {
                a.end.cmp(&b.end).then_with(|| a.start.cmp(&b.start))
            });

            self.next_record = None;
        }

        Ok(())
    }
}

impl Iterator for BamPileupStream {
    type Item = Result<PileupNode>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current_pos >= self.region.end {
            return None;
        }

        let pos = self.current_pos;
        self.current_pos += 1;

        if let Err(e) = self.advance_active_to(pos) {
            return Some(Err(e));
        }

        let mut node = PileupNode::new(pos);
        for read in &self.active {
            if let Some((base, qual)) = read.base_and_qual_at(pos) {
                if let Some(idx) = base_index(base) {
                    node.observe(idx, qual);
                }
            }
        }

        if node.depth == 0 {
            // Skip empty positions (stream remains bounded).
            return self.next();
        }

        Some(Ok(node))
    }
}

fn base_index(base: u8) -> Option<usize> {
    match base {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

fn complement(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => b'N',
    }
}

fn is_simple_match_cigar(rec: &bam::Record) -> bool {
    // Accept exactly one operation consuming the whole read: M/=/X.
    let cigar = rec.cigar();
    if cigar.len() != 1 {
        return false;
    }
    match cigar.iter().next() {
        Some(op) => match *op {
            bam::record::Cigar::Match(l)
            | bam::record::Cigar::Equal(l)
            | bam::record::Cigar::Diff(l) => l as usize == rec.seq_len(),
            _ => false,
        },
        None => false,
    }
}


