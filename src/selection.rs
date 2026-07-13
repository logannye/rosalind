//! Canonical analysis selections: whole genome, normalized intervals, and
//! deterministic reference-span shards.

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use thiserror::Error;

use crate::core::ContigSet;

/// A zero-based, half-open interval on one reference contig.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct GenomicInterval {
    /// Dense reference contig id.
    pub contig: u32,
    /// Inclusive zero-based start.
    pub start: u32,
    /// Exclusive zero-based end.
    pub end: u32,
}

impl GenomicInterval {
    /// Length in bases.
    pub fn len(self) -> u32 {
        self.end - self.start
    }

    /// Whether this interval contains no bases.
    pub fn is_empty(self) -> bool {
        self.start == self.end
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum SelectionOrigin {
    Programmatic,
    Region(String),
    Bed(PathBuf),
    Shard,
}

/// Ordered, disjoint, normalized genomic intervals.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IntervalSet {
    intervals: Vec<GenomicInterval>,
    origin: SelectionOrigin,
}

impl IntervalSet {
    /// Validate and normalize programmatically supplied intervals.
    pub fn new(
        intervals: impl IntoIterator<Item = GenomicInterval>,
        contigs: &ContigSet,
    ) -> Result<Self, SelectionError> {
        normalize(
            intervals.into_iter().collect(),
            contigs,
            SelectionOrigin::Programmatic,
        )
    }

    /// Parse one samtools-style, 1-based inclusive region.
    pub fn parse_region(region: &str, contigs: &ContigSet) -> Result<Self, SelectionError> {
        let (name, coordinates) = region.rsplit_once(':').ok_or_else(|| {
            SelectionError::InvalidRegion("expected CONTIG:START-END".to_string())
        })?;
        let (start, end) = coordinates.split_once('-').ok_or_else(|| {
            SelectionError::InvalidRegion("expected CONTIG:START-END".to_string())
        })?;
        let parse = |value: &str| {
            value
                .replace(',', "")
                .parse::<u64>()
                .map_err(|_| SelectionError::InvalidRegion(format!("invalid coordinate {value:?}")))
        };
        let start = parse(start)?;
        let end = parse(end)?;
        if start == 0 || end < start {
            return Err(SelectionError::InvalidRegion(
                "region coordinates are 1-based inclusive and require 1 <= START <= END".into(),
            ));
        }
        let contig = contigs
            .by_name(name)
            .ok_or_else(|| SelectionError::UnknownContig(name.to_string()))?;
        if end > contig.length as u64 {
            return Err(SelectionError::OutOfBounds {
                contig: name.to_string(),
                start: start - 1,
                end,
                length: contig.length,
            });
        }
        normalize(
            vec![GenomicInterval {
                contig: contig.id,
                start: u32::try_from(start - 1).expect("validated against u32 contig length"),
                end: u32::try_from(end).expect("validated against u32 contig length"),
            }],
            contigs,
            SelectionOrigin::Region(region.to_string()),
        )
    }

    /// Parse a BED file (zero-based, half-open), normalizing overlap and adjacency.
    pub fn from_bed(path: impl AsRef<Path>, contigs: &ContigSet) -> Result<Self, SelectionError> {
        let path = path.as_ref();
        let file = File::open(path)?;
        let mut intervals = Vec::new();
        for (line_index, line) in BufReader::new(file).lines().enumerate() {
            let line = line?;
            let trimmed = line.trim();
            if trimmed.is_empty()
                || trimmed.starts_with('#')
                || trimmed.starts_with("track ")
                || trimmed.starts_with("browser ")
            {
                continue;
            }
            let fields: Vec<&str> = trimmed.split_whitespace().collect();
            if fields.len() < 3 {
                return Err(SelectionError::InvalidBed {
                    line: line_index + 1,
                    message: "expected at least three columns".into(),
                });
            }
            let contig = contigs
                .by_name(fields[0])
                .ok_or_else(|| SelectionError::UnknownContig(fields[0].to_string()))?;
            let start = fields[1]
                .parse::<u64>()
                .map_err(|_| SelectionError::InvalidBed {
                    line: line_index + 1,
                    message: "start is not an unsigned integer".into(),
                })?;
            let end = fields[2]
                .parse::<u64>()
                .map_err(|_| SelectionError::InvalidBed {
                    line: line_index + 1,
                    message: "end is not an unsigned integer".into(),
                })?;
            if end < start {
                return Err(SelectionError::InvalidBed {
                    line: line_index + 1,
                    message: "end precedes start".into(),
                });
            }
            if end > contig.length as u64 {
                return Err(SelectionError::OutOfBounds {
                    contig: fields[0].to_string(),
                    start,
                    end,
                    length: contig.length,
                });
            }
            if start != end {
                intervals.push(GenomicInterval {
                    contig: contig.id,
                    start: start as u32,
                    end: end as u32,
                });
            }
        }
        normalize(intervals, contigs, SelectionOrigin::Bed(path.to_path_buf()))
    }

    /// Ordered, disjoint intervals.
    pub fn intervals(&self) -> &[GenomicInterval] {
        &self.intervals
    }

    /// Total selected reference bases.
    pub fn total_bases(&self) -> u64 {
        self.intervals
            .iter()
            .map(|interval| interval.len() as u64)
            .sum()
    }

    /// Largest decoded reference interval.
    pub fn largest_interval(&self) -> u64 {
        self.intervals
            .iter()
            .map(|interval| interval.len() as u64)
            .max()
            .unwrap_or(0)
    }

    /// Stable BLAKE3 of normalized little-endian interval triples.
    pub fn blake3(&self) -> String {
        let mut hasher = blake3::Hasher::new();
        for interval in &self.intervals {
            hasher.update(&interval.contig.to_le_bytes());
            hasher.update(&interval.start.to_le_bytes());
            hasher.update(&interval.end.to_le_bytes());
        }
        hasher.finalize().to_hex().to_string()
    }

    /// Original region expression, when constructed from `parse_region`.
    pub fn region_origin(&self) -> Option<&str> {
        match &self.origin {
            SelectionOrigin::Region(region) => Some(region),
            _ => None,
        }
    }

    /// Original BED path, when constructed from `from_bed`.
    pub fn bed_origin(&self) -> Option<&Path> {
        match &self.origin {
            SelectionOrigin::Bed(path) => Some(path),
            _ => None,
        }
    }
}

/// Selection applied to one analysis run.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub enum AnalysisSelection {
    /// Sequentially stream every reference contig; BAM index not required.
    #[default]
    WholeGenome,
    /// Run normalized indexed intervals.
    Intervals(IntervalSet),
    /// Run one deterministic reference-span shard.
    Shard {
        /// Total number of shards.
        count: u32,
        /// Zero-based shard index.
        index: u32,
        /// Intervals owned by this shard.
        intervals: IntervalSet,
    },
}

impl AnalysisSelection {
    /// Construct one deterministic `reference-span-v1` shard.
    pub fn shard(count: u32, index: u32, contigs: &ContigSet) -> Result<Self, SelectionError> {
        let total = contigs.total_length();
        if count == 0 {
            return Err(SelectionError::InvalidShard(
                "shard count must be positive".into(),
            ));
        }
        if index >= count {
            return Err(SelectionError::InvalidShard(format!(
                "shard index {index} must be less than shard count {count}"
            )));
        }
        if count as u64 > total {
            return Err(SelectionError::InvalidShard(format!(
                "shard count {count} exceeds reference length {total}"
            )));
        }
        let start = ((total as u128 * index as u128) / count as u128) as u64;
        let end = ((total as u128 * (index + 1) as u128) / count as u128) as u64;
        let mut intervals = Vec::new();
        for contig in contigs.iter() {
            let contig_start = contig.global_offset;
            let contig_end = contig_start + contig.length as u64;
            let owned_start = start.max(contig_start);
            let owned_end = end.min(contig_end);
            if owned_start < owned_end {
                intervals.push(GenomicInterval {
                    contig: contig.id,
                    start: (owned_start - contig_start) as u32,
                    end: (owned_end - contig_start) as u32,
                });
            }
        }
        Ok(Self::Shard {
            count,
            index,
            intervals: normalize(intervals, contigs, SelectionOrigin::Shard)?,
        })
    }

    /// Selected intervals, or `None` for whole-genome streaming.
    pub fn intervals(&self) -> Option<&IntervalSet> {
        match self {
            Self::WholeGenome => None,
            Self::Intervals(intervals) | Self::Shard { intervals, .. } => Some(intervals),
        }
    }

    /// Largest reference allocation required by this selection.
    pub fn largest_reference_span(&self, contigs: &ContigSet) -> u64 {
        self.intervals().map_or_else(
            || {
                contigs
                    .iter()
                    .map(|contig| contig.length as u64)
                    .max()
                    .unwrap_or(0)
            },
            IntervalSet::largest_interval,
        )
    }

    /// Stable selection kind used in receipts.
    pub fn kind(&self) -> &'static str {
        match self {
            Self::WholeGenome => "whole",
            Self::Intervals(intervals) if intervals.region_origin().is_some() => "region",
            Self::Intervals(_) => "regions",
            Self::Shard { .. } => "shard",
        }
    }
}

/// Selection parsing and validation failures.
#[derive(Debug, Error)]
pub enum SelectionError {
    /// Filesystem or BED read failure.
    #[error("selection I/O error: {0}")]
    Io(#[from] std::io::Error),
    /// Invalid `CONTIG:START-END` syntax or coordinates.
    #[error("invalid region: {0}")]
    InvalidRegion(String),
    /// Invalid BED row.
    #[error("invalid BED at line {line}: {message}")]
    InvalidBed {
        /// One-based BED line number.
        line: usize,
        /// Concrete parsing or validation failure.
        message: String,
    },
    /// Contig absent from the analysis reference.
    #[error("unknown reference contig {0:?}")]
    UnknownContig(String),
    /// Interval exceeds its contig.
    #[error("interval {contig}:{start}-{end} exceeds contig length {length}")]
    OutOfBounds {
        /// Reference contig name.
        contig: String,
        /// Zero-based start.
        start: u64,
        /// Zero-based exclusive end.
        end: u64,
        /// Declared reference contig length.
        length: u32,
    },
    /// Invalid shard count or index.
    #[error("invalid shard: {0}")]
    InvalidShard(String),
}

fn normalize(
    mut intervals: Vec<GenomicInterval>,
    contigs: &ContigSet,
    origin: SelectionOrigin,
) -> Result<IntervalSet, SelectionError> {
    for interval in &intervals {
        let contig = contigs
            .by_id(interval.contig)
            .ok_or_else(|| SelectionError::UnknownContig(interval.contig.to_string()))?;
        if interval.start > interval.end || interval.end > contig.length {
            return Err(SelectionError::OutOfBounds {
                contig: contig.name.to_string(),
                start: interval.start as u64,
                end: interval.end as u64,
                length: contig.length,
            });
        }
    }
    intervals.retain(|interval| !interval.is_empty());
    intervals.sort_unstable();
    let mut merged: Vec<GenomicInterval> = Vec::with_capacity(intervals.len());
    for interval in intervals {
        if let Some(previous) = merged.last_mut() {
            if previous.contig == interval.contig && interval.start <= previous.end {
                previous.end = previous.end.max(interval.end);
                continue;
            }
        }
        merged.push(interval);
    }
    Ok(IntervalSet {
        intervals: merged,
        origin,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn contigs() -> ContigSet {
        let mut contigs = ContigSet::new();
        contigs.push("chr1", 10);
        contigs.push("chr-two", 7);
        contigs
    }

    #[test]
    fn region_is_one_based_inclusive() {
        let set = IntervalSet::parse_region("chr1:2-5", &contigs()).unwrap();
        assert_eq!(
            set.intervals(),
            &[GenomicInterval {
                contig: 0,
                start: 1,
                end: 5
            }]
        );
    }

    #[test]
    fn normalization_merges_overlap_and_adjacency() {
        let set = IntervalSet::new(
            [
                GenomicInterval {
                    contig: 0,
                    start: 5,
                    end: 8,
                },
                GenomicInterval {
                    contig: 0,
                    start: 1,
                    end: 5,
                },
                GenomicInterval {
                    contig: 1,
                    start: 0,
                    end: 2,
                },
            ],
            &contigs(),
        )
        .unwrap();
        assert_eq!(
            set.intervals(),
            &[
                GenomicInterval {
                    contig: 0,
                    start: 1,
                    end: 8
                },
                GenomicInterval {
                    contig: 1,
                    start: 0,
                    end: 2
                },
            ]
        );
    }

    #[test]
    fn shard_union_is_complete_and_disjoint() {
        let contigs = contigs();
        let mut ownership = vec![0u8; contigs.total_length() as usize];
        for index in 0..6 {
            let AnalysisSelection::Shard { intervals, .. } =
                AnalysisSelection::shard(6, index, &contigs).unwrap()
            else {
                unreachable!()
            };
            for interval in intervals.intervals() {
                let contig = contigs.by_id(interval.contig).unwrap();
                for position in interval.start..interval.end {
                    ownership[(contig.global_offset + position as u64) as usize] += 1;
                }
            }
        }
        assert!(ownership.iter().all(|owners| *owners == 1));
    }

    proptest! {
        #[test]
        fn randomized_shards_cover_every_base_once(
            lengths in prop::collection::vec(1u32..200, 1..10),
            requested_count in 1u32..100,
        ) {
            let mut contigs = ContigSet::new();
            for (index, length) in lengths.into_iter().enumerate() {
                contigs.push(format!("contig-{index}:unusual"), length);
            }
            let count = requested_count.min(contigs.total_length() as u32);
            let mut ownership = vec![0u8; contigs.total_length() as usize];
            for index in 0..count {
                let AnalysisSelection::Shard { intervals, .. } =
                    AnalysisSelection::shard(count, index, &contigs).unwrap()
                else { unreachable!() };
                for interval in intervals.intervals() {
                    let contig = contigs.by_id(interval.contig).unwrap();
                    for position in interval.start..interval.end {
                        ownership[(contig.global_offset + position as u64) as usize] += 1;
                    }
                }
            }
            prop_assert!(ownership.iter().all(|owners| *owners == 1));
        }
    }

    #[test]
    fn bed_comments_unusual_names_ordering_and_empty_rows_normalize() {
        let root = std::env::temp_dir().join(format!(
            "rosalind-bed-{}",
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&root).unwrap();
        let bed = root.join("intervals with spaces.bed");
        std::fs::write(
            &bed,
            "# comment\ntrack name=x\nchr-two\t2\t2\nchr-two\t4\t7\nchr1\t5\t8\nchr1\t1\t5\n",
        )
        .unwrap();
        let set = IntervalSet::from_bed(&bed, &contigs()).unwrap();
        assert_eq!(
            set.intervals(),
            &[
                GenomicInterval {
                    contig: 0,
                    start: 1,
                    end: 8
                },
                GenomicInterval {
                    contig: 1,
                    start: 4,
                    end: 7
                },
            ]
        );
        std::fs::write(&bed, "chr1\t0\t18446744073709551615\n").unwrap();
        assert!(matches!(
            IntervalSet::from_bed(&bed, &contigs()),
            Err(SelectionError::OutOfBounds { .. })
        ));
        std::fs::write(&bed, "missing\t0\t1\n").unwrap();
        assert!(matches!(
            IntervalSet::from_bed(&bed, &contigs()),
            Err(SelectionError::UnknownContig(_))
        ));
        std::fs::remove_dir_all(root).unwrap();
    }
}
