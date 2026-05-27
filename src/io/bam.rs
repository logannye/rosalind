//! The BAM → `core::AlignedRead` boundary. `rust_htslib` lives here and never
//! leaks past this module: `BamSource` yields canonical `core::AlignedRead`s and
//! implements `pileup::ReadSource`, so the kernel stays htslib-free.
//!
//! This adapter pre-loads + sorts the records (it is not yet bounded-memory; a
//! streaming `.csi`-fetch BAM source is a later phase). SEQ is taken
//! forward-oriented per the SAM spec — no reverse-complement is applied.

use std::path::Path;
use std::sync::Arc;

use rust_htslib::bam::record::Cigar as BamCigar;
use rust_htslib::bam::{self, Read as BamRead};

use crate::core::{AlignedRead, CigarOp, CigarOpKind, ContigSet, CoreError, Position, SamFlags};
use crate::pileup::ReadSource;

/// Read all mapped records of a BAM into canonical `core::AlignedRead`s, mapping
/// each record's reference name to a contig id via `contigs`. Records that are
/// unmapped, have no tid, or whose reference is absent from `contigs` are
/// skipped. SEQ is uppercased and kept forward-oriented; strand is recorded in
/// `flags` only.
pub fn read_bam_as_core_reads(
    path: &Path,
    contigs: &ContigSet,
) -> Result<Vec<AlignedRead>, CoreError> {
    let mut reader = bam::Reader::from_path(path)
        .map_err(|e| CoreError::MalformedRecord(format!("open BAM {}: {e}", path.display())))?;
    let header = reader.header().to_owned();

    let mut out = Vec::new();
    for rec in reader.records() {
        let rec = rec.map_err(|e| CoreError::MalformedRecord(e.to_string()))?;
        if rec.is_unmapped() {
            continue;
        }
        let tid = rec.tid();
        if tid < 0 {
            continue;
        }
        let name = std::str::from_utf8(header.tid2name(tid as u32))
            .map_err(|_| CoreError::MalformedRecord("BAM reference name is not UTF-8".into()))?;
        let contig = match contigs.by_name(name) {
            Some(c) => c.id,
            None => continue,
        };
        let pos0 = rec.pos();
        if pos0 < 0 {
            continue;
        }

        let mut cigar = Vec::new();
        for c in rec.cigar().iter() {
            let (kind, len) = match *c {
                BamCigar::Match(l) | BamCigar::Equal(l) | BamCigar::Diff(l) => {
                    (CigarOpKind::Match, l)
                }
                BamCigar::Ins(l) => (CigarOpKind::Insertion, l),
                BamCigar::Del(l) => (CigarOpKind::Deletion, l),
                BamCigar::RefSkip(l) => (CigarOpKind::RefSkip, l),
                BamCigar::SoftClip(l) => (CigarOpKind::SoftClip, l),
                BamCigar::HardClip(l) => (CigarOpKind::HardClip, l),
                BamCigar::Pad(l) => (CigarOpKind::Pad, l),
            };
            cigar.push(CigarOp::new(kind, len));
        }

        let seq: Vec<u8> = rec
            .seq()
            .as_bytes()
            .iter()
            .map(|b| b.to_ascii_uppercase())
            .collect();
        let qual: Vec<u8> = rec.qual().to_vec();

        out.push(AlignedRead {
            contig,
            pos: Position(pos0 as u32),
            mapq: rec.mapq(),
            flags: SamFlags(rec.flags()),
            cigar,
            seq: Arc::from(seq.into_boxed_slice()),
            qual: Arc::from(qual.into_boxed_slice()),
        });
    }
    Ok(out)
}

/// A `ReadSource` over a BAM file. Pre-loads and coordinate-sorts the records on
/// construction, then yields them in `(contig, pos)` order.
#[derive(Debug)]
pub struct BamSource {
    reads: std::vec::IntoIter<AlignedRead>,
}

impl BamSource {
    /// Open `path`, convert its mapped records to `core::AlignedRead`s (mapping
    /// reference names via `contigs`), and sort by `(contig, pos)`.
    pub fn new(path: &Path, contigs: &ContigSet) -> Result<Self, CoreError> {
        let mut reads = read_bam_as_core_reads(path, contigs)?;
        reads.sort_by_key(|r| (r.contig, r.pos));
        Ok(Self {
            reads: reads.into_iter(),
        })
    }
}

impl ReadSource for BamSource {
    fn next_read(&mut self) -> Result<Option<AlignedRead>, CoreError> {
        Ok(self.reads.next())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::{Cigar, CigarString, Record};
    use std::time::{SystemTime, UNIX_EPOCH};

    fn tmp(name: &str) -> std::path::PathBuf {
        let ts = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        std::env::temp_dir().join(format!("rosalind-bamsrc-{name}-{ts}.bam"))
    }

    fn one_contig() -> ContigSet {
        let mut c = ContigSet::new();
        c.push("chr1", 1000);
        c
    }

    #[test]
    fn reads_bam_into_core_reads_with_flags_and_cigar() {
        let path = tmp("basic");
        let mut header = bam::Header::new();
        header.push_record(
            bam::header::HeaderRecord::new(b"SQ")
                .push_tag(b"SN", &"chr1")
                .push_tag(b"LN", &1000),
        );
        {
            let mut w = bam::Writer::from_path(&path, &header, bam::Format::Bam).unwrap();
            // Forward read at pos 10, 3M.
            let mut fwd = Record::new();
            fwd.set(
                b"fwd",
                Some(&CigarString::from(vec![Cigar::Match(3)])),
                b"ACG",
                b"III",
            );
            fwd.set_tid(0);
            fwd.set_pos(10);
            fwd.set_flags(0);
            fwd.set_mapq(60);
            w.write(&fwd).unwrap();
            // Reverse read at pos 20, 2M.
            let mut rev = Record::new();
            rev.set(
                b"rev",
                Some(&CigarString::from(vec![Cigar::Match(2)])),
                b"TT",
                b"II",
            );
            rev.set_tid(0);
            rev.set_pos(20);
            rev.set_flags(0x10); // REVERSE
            rev.set_mapq(40);
            w.write(&rev).unwrap();
        }

        let reads = read_bam_as_core_reads(&path, &one_contig()).unwrap();
        assert_eq!(reads.len(), 2);
        let fwd = reads.iter().find(|r| r.pos.0 == 10).unwrap();
        assert_eq!(fwd.contig, 0);
        assert_eq!(fwd.mapq, 60);
        assert!(!fwd.flags.is_reverse());
        assert_eq!(fwd.cigar, vec![CigarOp::new(CigarOpKind::Match, 3)]);
        assert_eq!(&fwd.seq[..], b"ACG");
        let rev = reads.iter().find(|r| r.pos.0 == 20).unwrap();
        assert!(rev.flags.is_reverse());
        // SEQ stays forward-oriented (no reverse-complement applied).
        assert_eq!(&rev.seq[..], b"TT");

        // BamSource yields the same reads, coordinate-sorted.
        let mut src = BamSource::new(&path, &one_contig()).unwrap();
        let first = src.next_read().unwrap().unwrap();
        assert_eq!(first.pos.0, 10);
        let second = src.next_read().unwrap().unwrap();
        assert_eq!(second.pos.0, 20);
        assert!(src.next_read().unwrap().is_none());

        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn unknown_contig_records_are_skipped() {
        let path = tmp("unknown");
        let mut header = bam::Header::new();
        header.push_record(
            bam::header::HeaderRecord::new(b"SQ")
                .push_tag(b"SN", &"chrX")
                .push_tag(b"LN", &1000),
        );
        {
            let mut w = bam::Writer::from_path(&path, &header, bam::Format::Bam).unwrap();
            let mut r = Record::new();
            r.set(
                b"x",
                Some(&CigarString::from(vec![Cigar::Match(2)])),
                b"AC",
                b"II",
            );
            r.set_tid(0);
            r.set_pos(5);
            r.set_flags(0);
            r.set_mapq(60);
            w.write(&r).unwrap();
        }
        // ContigSet has chr1 only → the chrX record is skipped.
        let reads = read_bam_as_core_reads(&path, &one_contig()).unwrap();
        assert!(reads.is_empty());
        std::fs::remove_file(&path).ok();
    }
}
