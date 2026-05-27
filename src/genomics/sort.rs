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
    let header = bam::Header::from_template(reader.header());

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

fn spill_chunk(
    header: &bam::Header,
    dir: &Path,
    idx: usize,
    chunk: &mut Vec<Record>,
) -> Result<PathBuf> {
    chunk.sort_by(|a, b| sort_key_cmp(a, b));

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
    tid: i32,
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
            tid: record.tid(),
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
        self.key == other.key
    }
}

impl Ord for HeapItem {
    fn cmp(&self, other: &Self) -> Ordering {
        self.key.cmp(&other.key)
    }
}

impl PartialOrd for HeapItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

fn sort_key_cmp(a: &Record, b: &Record) -> Ordering {
    a.tid()
        .cmp(&b.tid())
        .then_with(|| a.pos().cmp(&b.pos()))
        .then_with(|| a.is_reverse().cmp(&b.is_reverse()))
        .then_with(|| a.qname().cmp(b.qname()))
}
