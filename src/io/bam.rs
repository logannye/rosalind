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

/// Map one BAM `Record` to a canonical `AlignedRead`, or `None` if it is
/// unmapped, has no tid, a negative pos, or a reference name absent from
/// `contigs`. Shared by `read_bam_as_core_reads` and `StreamingBamSource`.
pub(crate) fn record_to_aligned_read(
    rec: &bam::Record,
    header: &bam::HeaderView,
    contigs: &ContigSet,
) -> Result<Option<AlignedRead>, CoreError> {
    if rec.is_unmapped() {
        return Ok(None);
    }
    let tid = rec.tid();
    if tid < 0 {
        return Ok(None);
    }
    let name = std::str::from_utf8(header.tid2name(tid as u32))
        .map_err(|_| CoreError::MalformedRecord("BAM reference name is not UTF-8".into()))?;
    let contig = match contigs.by_name(name) {
        Some(c) => c.id,
        None => return Ok(None),
    };
    let pos0 = rec.pos();
    if pos0 < 0 {
        return Ok(None);
    }

    let mut cigar = Vec::new();
    for c in rec.cigar().iter() {
        let (kind, len) = match *c {
            BamCigar::Match(l) | BamCigar::Equal(l) | BamCigar::Diff(l) => (CigarOpKind::Match, l),
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

    Ok(Some(AlignedRead {
        contig,
        pos: Position(pos0 as u32),
        mapq: rec.mapq(),
        flags: SamFlags(rec.flags()),
        cigar,
        seq: Arc::from(seq.into_boxed_slice()),
        qual: Arc::from(qual.into_boxed_slice()),
    }))
}

/// Cross-check the BAM header's `@SQ` contig lengths against the index `ContigSet`.
/// For every header contig whose name is also in the index, the lengths must match;
/// a mismatch means the alignments were built against a different reference (a
/// different assembly/patch), which would silently produce coordinate-shifted or
/// truncated calls. Names present in the header but absent from the index are left
/// alone (their records are skipped downstream).
pub(crate) fn validate_contig_lengths(
    header: &bam::HeaderView,
    contigs: &ContigSet,
) -> Result<(), CoreError> {
    let mut matched = 0usize;
    for tid in 0..header.target_count() {
        let name = std::str::from_utf8(header.tid2name(tid))
            .map_err(|_| CoreError::MalformedRecord("BAM reference name is not UTF-8".into()))?;
        if let Some(c) = contigs.by_name(name) {
            matched += 1;
            let header_len = header.target_len(tid);
            if header_len != Some(c.length as u64) {
                return Err(CoreError::MalformedRecord(format!(
                    "BAM @SQ length for contig '{name}' ({}) disagrees with the index ({}); \
                     the alignments were built against a different reference",
                    header_len
                        .map(|l| l.to_string())
                        .unwrap_or_else(|| "missing".into()),
                    c.length
                )));
            }
        }
    }
    // If the BAM names a reference but NONE of its @SQ contigs resolve into the
    // index, the two use incompatible naming schemes (the classic UCSC 'chr1' vs
    // Ensembl '1' / GRCh38-vs-hg38 mismatch). Left unchecked, every record is
    // silently dropped and the run writes an empty VCF with exit 0 — the worst
    // kind of failure. Refuse up front with an actionable message instead.
    if header.target_count() > 0 && matched == 0 {
        let bam_names: Vec<String> = (0..header.target_count())
            .filter_map(|tid| {
                std::str::from_utf8(header.tid2name(tid))
                    .ok()
                    .map(str::to_string)
            })
            .take(3)
            .collect();
        return Err(CoreError::MalformedRecord(format!(
            "no BAM @SQ contig name matches the index — the alignments and the index use \
             incompatible naming schemes (e.g. UCSC 'chr1' vs Ensembl '1'). BAM names start \
             with [{}]; check the reference/assembly used for alignment.",
            bam_names.join(", ")
        )));
    }
    Ok(())
}

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
        if let Some(read) = record_to_aligned_read(&rec, &header, contigs)? {
            out.push(read);
        }
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

/// A bounded `ReadSource` over a **coordinate-sorted** BAM: pulls one record at a
/// time (never materialises the file). Enforces sort order via a `(contig_id, pos)`
/// monotonicity guard — `next_read` errors if a read arrives out of order (an
/// unsorted BAM, or one whose `@SQ` order disagrees with `contigs`).
#[derive(Debug)]
pub struct StreamingBamSource<'a> {
    reader: bam::Reader,
    header: bam::HeaderView,
    contigs: &'a ContigSet,
    last: Option<(u32, u32)>,
}

impl<'a> StreamingBamSource<'a> {
    /// Open `path` for streaming, mapping reference names via `contigs`.
    pub fn new(path: &Path, contigs: &'a ContigSet) -> Result<Self, CoreError> {
        let reader = bam::Reader::from_path(path)
            .map_err(|e| CoreError::MalformedRecord(format!("open BAM {}: {e}", path.display())))?;
        let header = reader.header().to_owned();
        // Reject a BAM aligned to a different-length reference up front — names
        // matching the index must agree on length, or every coordinate is suspect.
        validate_contig_lengths(&header, contigs)?;
        Ok(Self {
            reader,
            header,
            contigs,
            last: None,
        })
    }
}

impl ReadSource for StreamingBamSource<'_> {
    fn next_read(&mut self) -> Result<Option<AlignedRead>, CoreError> {
        loop {
            let mut record = bam::Record::new();
            match self.reader.read(&mut record) {
                None => return Ok(None),
                Some(Err(e)) => return Err(CoreError::MalformedRecord(e.to_string())),
                Some(Ok(())) => {}
            }
            let Some(read) = record_to_aligned_read(&record, &self.header, self.contigs)? else {
                continue; // unmapped / absent contig / etc. — skip
            };
            let key = (read.contig, read.pos.0);
            if let Some(prev) = self.last {
                if key < prev {
                    return Err(CoreError::MalformedRecord(format!(
                        "alignments are not coordinate-sorted (or @SQ order disagrees with the \
                         index): read at contig {} pos {} follows contig {} pos {}",
                        key.0, key.1, prev.0, prev.1
                    )));
                }
            }
            self.last = Some(key);
            return Ok(Some(read));
        }
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

    fn test_header(contigs: &[(&str, usize)]) -> bam::Header {
        let mut header = bam::Header::new();
        for (name, len) in contigs {
            let mut rec = bam::header::HeaderRecord::new(b"SQ");
            rec.push_tag(b"SN", name);
            rec.push_tag(b"LN", &(*len as i64));
            header.push_record(&rec);
        }
        header
    }

    fn push_record(
        writer: &mut bam::Writer,
        _header: &bam::HeaderView,
        tid: i32,
        pos: i64,
        seq: &[u8],
    ) {
        let mut rec = Record::new();
        let cigar = CigarString(vec![Cigar::Match(seq.len() as u32)]);
        let qual = vec![30u8; seq.len()];
        rec.set(b"r", Some(&cigar), seq, &qual);
        rec.set_tid(tid);
        rec.set_pos(pos);
        rec.set_mapq(60);
        writer.write(&rec).unwrap();
    }

    fn write_sorted_test_bam(path: &Path) {
        let header = test_header(&[("chr1", 1000), ("chr2", 1000)]);
        let mut writer = bam::Writer::from_path(path, &header, bam::Format::Bam).unwrap();
        let hv = writer.header().clone();
        push_record(&mut writer, &hv, 0, 10, b"ACGT");
        push_record(&mut writer, &hv, 0, 20, b"ACGT");
        push_record(&mut writer, &hv, 1, 5, b"ACGT");
    }

    fn write_unsorted_test_bam(path: &Path) {
        let header = test_header(&[("chr1", 1000)]);
        let mut writer = bam::Writer::from_path(path, &header, bam::Format::Bam).unwrap();
        let hv = writer.header().clone();
        push_record(&mut writer, &hv, 0, 50, b"ACGT");
        push_record(&mut writer, &hv, 0, 10, b"ACGT"); // out of order
    }

    #[test]
    fn streaming_source_yields_same_reads_as_materializing() {
        // Build a tiny coordinate-sorted BAM with 2 contigs, read it both ways.
        let dir = std::env::temp_dir().join(format!(
            "rosalind-stream-{}",
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let bam_path = dir.join("reads.bam");
        write_sorted_test_bam(&bam_path); // helper below

        let mut contigs = ContigSet::new();
        contigs.push("chr1", 1000);
        contigs.push("chr2", 1000);

        let materialized = read_bam_as_core_reads(&bam_path, &contigs).unwrap();
        let mut streamed = Vec::new();
        let mut src = StreamingBamSource::new(&bam_path, &contigs).unwrap();
        while let Some(r) = src.next_read().unwrap() {
            streamed.push(r);
        }
        assert_eq!(streamed.len(), materialized.len());
        for (a, b) in streamed.iter().zip(materialized.iter()) {
            assert_eq!((a.contig, a.pos), (b.contig, b.pos));
        }
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn streaming_source_rejects_disjoint_contig_naming() {
        // BAM @SQ uses Ensembl '1'; the index uses UCSC 'chr1'. No name resolves,
        // so opening must refuse up front rather than silently dropping every read
        // into an empty VCF.
        let bam_path = tmp("naming-mismatch");
        let header = test_header(&[("1", 1000)]); // Ensembl-style
        {
            let mut writer = bam::Writer::from_path(&bam_path, &header, bam::Format::Bam).unwrap();
            let hv = writer.header().clone();
            push_record(&mut writer, &hv, 0, 10, b"ACGT");
        }
        let contigs = one_contig(); // UCSC-style 'chr1'
        let err = StreamingBamSource::new(&bam_path, &contigs)
            .expect_err("disjoint naming must be rejected");
        let msg = err.to_string();
        assert!(
            msg.contains("no BAM @SQ contig name matches the index"),
            "unexpected error: {msg}"
        );
        std::fs::remove_file(&bam_path).ok();
    }

    #[test]
    fn streaming_source_rejects_out_of_order() {
        let dir = std::env::temp_dir().join(format!(
            "rosalind-stream-bad-{}",
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let bam_path = dir.join("unsorted.bam");
        write_unsorted_test_bam(&bam_path); // two reads, descending pos on one contig

        let mut contigs = ContigSet::new();
        contigs.push("chr1", 1000);

        let mut src = StreamingBamSource::new(&bam_path, &contigs).unwrap();
        // First read OK; the second (lower pos) must error.
        let _ = src.next_read().unwrap();
        assert!(
            src.next_read().is_err(),
            "out-of-order read must be rejected"
        );
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn streaming_source_rejects_sq_length_mismatch() {
        // BAM header says chr1 is 2000 bp; the index ContigSet says 1000 bp — the
        // alignments were built against a different reference. Reject at open.
        let dir = std::env::temp_dir().join(format!(
            "rosalind-stream-sqlen-{}",
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let bam_path = dir.join("wrongref.bam");
        let header = test_header(&[("chr1", 2000)]);
        {
            let mut writer = bam::Writer::from_path(&bam_path, &header, bam::Format::Bam).unwrap();
            let hv = writer.header().clone();
            push_record(&mut writer, &hv, 0, 10, b"ACGT");
        }
        let mut contigs = ContigSet::new();
        contigs.push("chr1", 1000);
        let err = StreamingBamSource::new(&bam_path, &contigs);
        assert!(err.is_err(), "length mismatch must be rejected at open");
        let msg = format!("{}", err.err().unwrap());
        assert!(
            msg.contains("disagrees with the index"),
            "clear message: {msg}"
        );
        std::fs::remove_dir_all(&dir).ok();
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
