//! Deterministic coordinate sorting for BAM files (external merge sort).
//!
//! This is the foundation for streaming pileup and calling: downstream stages
//! assume coordinate-sorted alignments.

use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::fs;
use std::path::{Path, PathBuf};

use anyhow::{anyhow, Context, Result};
use rust_htslib::bam;
use rust_htslib::bam::record::Record;
use rust_htslib::bam::Read as BamRead;

/// Deterministically coordinate-sort a BAM file using bounded memory.
///
/// - Uses stable sorting within each chunk.
/// - Uses a deterministic k-way merge across chunk files.
pub fn sort_bam_deterministic(
    input: impl AsRef<Path>,
    output: impl AsRef<Path>,
    memory_bytes: usize,
) -> Result<()> {
    let input = input.as_ref();
    let output = output.as_ref();

    if memory_bytes < 1 << 20 {
        return Err(anyhow!("memory_bytes must be at least 1MiB"));
    }

    let mut reader = bam::Reader::from_path(input)
        .with_context(|| format!("failed to open BAM {}", input.display()))?;
    let header = coordinate_sorted_header(reader.header())?;

    let temp_dir = PathBuf::from(format!("{}.sort_tmp", output.display()));
    fs::create_dir_all(&temp_dir)
        .with_context(|| format!("failed to create temp dir {}", temp_dir.display()))?;

    let mut chunk_paths = Vec::new();
    let mut current_chunk: Vec<Record> = Vec::new();
    let mut current_bytes: usize = 0;

    for rec_result in reader.records() {
        let rec = rec_result?;
        // rust-htslib does not expose raw record byte length; approximate.
        let approx = 128usize + rec.qname().len() + rec.seq_len() + rec.cigar().len() * 8;
        if !current_chunk.is_empty() && current_bytes + approx > memory_bytes {
            let chunk_path =
                spill_chunk(&header, &temp_dir, chunk_paths.len(), &mut current_chunk)?;
            chunk_paths.push(chunk_path);
            current_bytes = 0;
        }
        current_bytes += approx;
        current_chunk.push(rec);
    }

    if !current_chunk.is_empty() {
        let chunk_path = spill_chunk(&header, &temp_dir, chunk_paths.len(), &mut current_chunk)?;
        chunk_paths.push(chunk_path);
    }

    merge_chunks(&header, &chunk_paths, output)?;

    // Best-effort cleanup.
    for p in chunk_paths {
        let _ = fs::remove_file(p);
    }
    let _ = fs::remove_dir(&temp_dir);

    Ok(())
}

/// Preserve the input dictionary and metadata while declaring the property this
/// function has just established. `doctor` intentionally trusts the header for
/// its cheap preflight; callers can request a full record scan separately.
fn coordinate_sorted_header(view: &bam::HeaderView) -> Result<bam::Header> {
    let text = std::str::from_utf8(view.as_bytes()).context("BAM header is not UTF-8")?;
    let lines = text
        .lines()
        .filter(|line| !line.is_empty())
        .collect::<Vec<_>>();
    let mut header = bam::Header::new();

    let mut hd = bam::header::HeaderRecord::new(b"HD");
    let mut has_version = false;
    if let Some(line) = lines
        .iter()
        .find(|line| line.starts_with("@HD\t") || **line == "@HD")
    {
        for field in line.split('\t').skip(1) {
            let (tag, value) = field
                .split_once(':')
                .ok_or_else(|| anyhow!("malformed BAM @HD field {field:?}"))?;
            if tag == "VN" {
                has_version = true;
            }
            hd.push_tag(
                tag.as_bytes(),
                if tag == "SO" { "coordinate" } else { value },
            );
        }
    }
    if !has_version {
        hd.push_tag(b"VN", "1.6");
    }
    if !text.lines().any(|line| {
        line.starts_with("@HD\t") && line.split('\t').any(|field| field.starts_with("SO:"))
    }) {
        hd.push_tag(b"SO", "coordinate");
    }
    header.push_record(&hd);

    for line in lines.into_iter().filter(|line| !line.starts_with("@HD")) {
        if let Some(comment) = line.strip_prefix("@CO") {
            header.push_comment(comment.strip_prefix('\t').unwrap_or(comment).as_bytes());
            continue;
        }
        let mut fields = line.split('\t');
        let kind = fields
            .next()
            .and_then(|value| value.strip_prefix('@'))
            .filter(|value| value.len() == 2)
            .ok_or_else(|| anyhow!("malformed BAM header line {line:?}"))?;
        let mut record = bam::header::HeaderRecord::new(kind.as_bytes());
        for field in fields {
            let (tag, value) = field
                .split_once(':')
                .ok_or_else(|| anyhow!("malformed BAM header field {field:?}"))?;
            record.push_tag(tag.as_bytes(), value);
        }
        header.push_record(&record);
    }
    Ok(header)
}

fn spill_chunk(
    header: &bam::Header,
    dir: &Path,
    idx: usize,
    chunk: &mut Vec<Record>,
) -> Result<PathBuf> {
    chunk.sort_by(sort_key_cmp);

    let path = dir.join(format!("{idx:05}.bam"));
    let mut writer = bam::Writer::from_path(&path, header, bam::Format::Bam)
        .with_context(|| format!("failed to create chunk {}", path.display()))?;
    for rec in chunk.drain(..) {
        writer.write(&rec)?;
    }
    Ok(path)
}

fn merge_chunks(header: &bam::Header, chunk_paths: &[PathBuf], output: &Path) -> Result<()> {
    if chunk_paths.is_empty() {
        return Err(anyhow!("no chunks to merge"));
    }

    let mut readers: Vec<bam::Reader> = Vec::with_capacity(chunk_paths.len());
    for p in chunk_paths {
        readers.push(bam::Reader::from_path(p)?);
    }

    let mut heap: BinaryHeap<HeapItem> = BinaryHeap::new();
    for (idx, reader) in readers.iter_mut().enumerate() {
        if let Some(rec) = next_record(reader)? {
            heap.push(HeapItem::new(idx, rec));
        }
    }

    let mut writer = bam::Writer::from_path(output, header, bam::Format::Bam)
        .with_context(|| format!("failed to create output BAM {}", output.display()))?;

    while let Some(item) = heap.pop() {
        writer.write(&item.record)?;
        let src = item.source_idx;
        if let Some(next) = next_record(&mut readers[src])? {
            heap.push(HeapItem::new(src, next));
        }
    }

    Ok(())
}

fn next_record(reader: &mut bam::Reader) -> Result<Option<Record>> {
    let mut rec = Record::new();
    match reader.read(&mut rec) {
        None => Ok(None),
        Some(Ok(())) => Ok(Some(rec)),
        Some(Err(e)) => Err(anyhow!(e)),
    }
}

#[derive(Debug, Clone, Eq, PartialEq)]
struct SortKey {
    tid: u32,
    pos: i64,
    is_reverse: bool,
    qname: Vec<u8>,
}

impl Ord for SortKey {
    fn cmp(&self, other: &Self) -> Ordering {
        // Reverse ordering for max-heap usage in BinaryHeap (we want min-key first).
        other
            .tid
            .cmp(&self.tid)
            .then_with(|| other.pos.cmp(&self.pos))
            .then_with(|| other.is_reverse.cmp(&self.is_reverse))
            .then_with(|| other.qname.cmp(&self.qname))
    }
}

impl PartialOrd for SortKey {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Debug)]
struct HeapItem {
    key: SortKey,
    source_idx: usize,
    record: Record,
}

impl HeapItem {
    fn new(source_idx: usize, record: Record) -> Self {
        let key = SortKey {
            tid: coordinate_tid(&record),
            pos: record.pos(),
            is_reverse: record.is_reverse(),
            qname: record.qname().to_vec(),
        };
        Self {
            key,
            source_idx,
            record,
        }
    }
}

impl Eq for HeapItem {}

impl PartialEq for HeapItem {
    fn eq(&self, other: &Self) -> bool {
        // Consistent with the now-total `Ord` (key + source_idx).
        self.key == other.key && self.source_idx == other.source_idx
    }
}

impl Ord for HeapItem {
    fn cmp(&self, other: &Self) -> Ordering {
        // Total order: key first (reversed for the min-key-first max-heap, as
        // `SortKey::cmp` already is), then `source_idx` — also reversed so the
        // LOWER source_idx pops first. Records are read in input order and
        // assigned to chunks sequentially, so equal-key records pop in input
        // order regardless of how they were partitioned into chunks — making the
        // merge output independent of the `--memory-mb` budget.
        self.key
            .cmp(&other.key)
            .then_with(|| other.source_idx.cmp(&self.source_idx))
    }
}

impl PartialOrd for HeapItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

fn sort_key_cmp(a: &Record, b: &Record) -> Ordering {
    coordinate_tid(a)
        .cmp(&coordinate_tid(b))
        .then_with(|| a.pos().cmp(&b.pos()))
        .then_with(|| a.is_reverse().cmp(&b.is_reverse()))
        .then_with(|| a.qname().cmp(b.qname()))
}

// SAM coordinate order places records without a reference (tid=-1) after
// every reference-assigned record, including unmapped reads with mate coordinates.
fn coordinate_tid(record: &Record) -> u32 {
    record.tid() as u32
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::{Cigar, CigarString};

    fn rec(qname: &[u8], tid: i32, pos: i64) -> Record {
        let mut r = Record::new();
        let cigar = CigarString(vec![Cigar::Match(1)]);
        r.set(qname, Some(&cigar), b"A", &[30u8]);
        r.set_tid(tid);
        r.set_pos(pos);
        r
    }

    #[test]
    fn merge_tie_break_pops_equal_key_records_in_source_index_order() {
        // Two records with a fully-equal sort key (same tid/pos/strand/qname) must
        // pop in ascending source_idx (= input order), independent of push order —
        // this is what makes the merge output independent of the chunk partition
        // (--memory-mb). A record with a smaller pos pops before both.
        let a = HeapItem::new(0, rec(b"dup", 0, 100));
        let b = HeapItem::new(3, rec(b"dup", 0, 100)); // equal key, higher source_idx
        let early = HeapItem::new(2, rec(b"dup", 0, 50)); // smaller pos -> pops first
        let unplaced = HeapItem::new(1, rec(b"unplaced", -1, -1));

        let mut heap: BinaryHeap<HeapItem> = BinaryHeap::new();
        // Push in an order that does NOT match the desired pop order.
        heap.push(b);
        heap.push(early);
        heap.push(a);
        heap.push(unplaced);

        let p1 = heap.pop().unwrap();
        let p2 = heap.pop().unwrap();
        let p3 = heap.pop().unwrap();
        assert_eq!(p1.record.pos(), 50, "smallest key pops first");
        assert_eq!(
            p2.source_idx, 0,
            "equal-key tie: lower source_idx pops first"
        );
        assert_eq!(p3.source_idx, 3, "then the higher source_idx");
        assert_eq!(heap.pop().unwrap().record.tid(), -1);
    }
}
