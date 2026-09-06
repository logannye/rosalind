use super::*;
use crate::core::ContigSet;
use crate::selection::GenomicInterval;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

/// One original panel target; overlap with other targets is preserved.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PanelTarget {
    /// Original target identifier, preserved in panel output.
    pub id: String,
    /// Reference contig identifier or canonical contig name for a batch.
    pub contig: u32,
    /// Zero-based inclusive target start.
    pub start: u32,
    /// Zero-based exclusive target end.
    pub end: u32,
}
impl PanelTarget {
    /// BED columns1-3 select an interval; column4 supplies the target identifier.
    /// Without column4, `target-{line_number}` is used. Input order is preserved.
    pub fn from_bed(
        path: impl AsRef<Path>,
        contigs: &ContigSet,
    ) -> Result<Vec<Self>, EvidenceError> {
        let mut targets = Vec::new();
        for (line_number, line) in BufReader::new(std::fs::File::open(path)?)
            .lines()
            .enumerate()
        {
            let line = line?;
            let trimmed = line.trim();
            if trimmed.is_empty()
                || trimmed.starts_with('#')
                || trimmed.starts_with("track ")
                || trimmed.starts_with("browser ")
            {
                continue;
            }
            let fields: Vec<_> = trimmed.split_whitespace().collect();
            if fields.len() < 3 {
                return Err(EvidenceError::InvalidRequest(format!(
                    "panel BED line{} needs three columns",
                    line_number + 1
                )));
            }
            let contig = contigs.by_name(fields[0]).ok_or_else(|| {
                EvidenceError::InvalidRequest(format!("unknown panel contig {}", fields[0]))
            })?;
            let parse = |value: &str| {
                value.parse::<u32>().map_err(|_| {
                    EvidenceError::InvalidRequest("invalid panel BED coordinate".into())
                })
            };
            let start = parse(fields[1])?;
            let end = parse(fields[2])?;
            if start >= end || end > contig.length {
                return Err(EvidenceError::InvalidRequest(format!(
                    "panel target line{} is empty or out of bounds",
                    line_number + 1
                )));
            }
            let id = fields.get(3).map_or_else(
                || format!("target-{}", line_number + 1),
                |value| value.to_string(),
            );
            targets.push(Self {
                id,
                contig: contig.id,
                start,
                end,
            });
        }
        Ok(targets)
    }
}
/// Integer panel sufficient statistics; denominators always use full target length.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PanelSummary {
    /// Original, unmerged panel target.
    pub target: PanelTarget,
    /// Number of unique input loci received for this target.
    pub observed_loci: u64,
    /// Sum of depth before profile filters across the entire target.
    pub prefilter_depth_sum: u64,
    /// Sum of depth after read-level filtering across the entire target.
    pub aligned_depth_sum: u64,
    /// Sum of quality-qualified A/C/G/T depth across the entire target.
    pub callable_depth_sum: u64,
    /// Exact sum of callable-observation base qualities.
    pub base_quality_sum: u64,
    /// Exact sum of callable-observation mapping qualities.
    pub mapping_quality_sum: u64,
    /// Target loci reaching the configured minimum callable depth.
    pub callable_positions: u64,
    /// Minimum callable depth, including uncovered target positions.
    pub minimum_callable_depth: u64,
    /// Maximum callable depth across the target.
    pub maximum_callable_depth: u64,
    /// Number of positions reaching1x,10x,20x,30x callable depth.
    pub breadth_positions: [u64; 4],
}
impl PanelSummary {
    /// Return the complete target length, including uncovered positions.
    pub fn target_length(&self) -> u64 {
        (self.target.end - self.target.start) as u64
    }
    /// Mean callable depth over the complete target length.
    pub fn mean_callable_depth(&self) -> f64 {
        self.callable_depth_sum as f64 / self.target_length() as f64
    }
}

/// Bounded panel aggregation over the union of requested target loci. Targets
/// remain distinct even when intervals or identifiers overlap.
#[derive(Debug)]
pub struct PanelQcAnalyzer {
    summaries: Vec<PanelSummary>,
    last_locus: Option<(u32, u32)>,
    /// Configured depth threshold for the callable-position count.
    pub min_callable_depth: u64,
    sorted_targets: Vec<usize>,
    next_target: usize,
    active_targets: Vec<usize>,
}
impl PanelQcAnalyzer {
    /// Construct the configured value without starting evidence extraction.
    pub fn new(targets: Vec<PanelTarget>) -> Result<Self, EvidenceError> {
        if targets
            .iter()
            .any(|target| target.start >= target.end || target.id.contains(['\t', '\n', '\r']))
        {
            return Err(EvidenceError::InvalidRequest(
                "panel targets must be nonempty and IDs cannot contain tabs/newlines".into(),
            ));
        }
        let mut sorted_targets: Vec<_> = (0..targets.len()).collect();
        sorted_targets.sort_by_key(|&index| (targets[index].contig, targets[index].start, index));
        Ok(Self {
            summaries: targets
                .into_iter()
                .map(|target| PanelSummary {
                    target,
                    observed_loci: 0,
                    prefilter_depth_sum: 0,
                    aligned_depth_sum: 0,
                    callable_depth_sum: 0,
                    base_quality_sum: 0,
                    mapping_quality_sum: 0,
                    callable_positions: 0,
                    minimum_callable_depth: u64::MAX,
                    maximum_callable_depth: 0,
                    breadth_positions: [0; 4],
                })
                .collect(),
            last_locus: None,
            min_callable_depth: 20,
            sorted_targets,
            next_target: 0,
            active_targets: Vec::new(),
        })
    }
    /// Set the callability threshold before consuming any batches.
    pub fn with_min_callable_depth(mut self, depth: u64) -> Self {
        self.min_callable_depth = depth;
        self
    }
    /// Return the unmerged targets as an extraction selection; the engine unions overlaps.
    pub fn selection(&self) -> EvidenceSelection {
        EvidenceSelection::Intervals(
            self.summaries
                .iter()
                .map(|summary| GenomicInterval {
                    contig: summary.target.contig,
                    start: summary.target.start,
                    end: summary.target.end,
                })
                .collect(),
        )
    }
    /// Borrow target summaries in the original input order.
    pub fn summaries(&self) -> &[PanelSummary] {
        &self.summaries
    }
    /// Write deterministic panel metrics with full-target-length denominators.
    pub fn write_tsv(
        &self,
        mut output: impl Write,
        contigs: &ContigSet,
    ) -> Result<(), EvidenceError> {
        writeln!(output,"#target\tcontig\tstart\tend\tlength\tprefilter_depth_sum\taligned_depth_sum\tcallable_depth_sum\tbase_quality_sum\tmapping_quality_sum\tcallable_positions\tcallable_threshold\tmin_callable_depth\tmax_callable_depth\tmean_callable_depth\tbreadth_1x\tbreadth_10x\tbreadth_20x\tbreadth_30x")?;
        for summary in &self.summaries {
            let contig = contigs.by_id(summary.target.contig).ok_or_else(|| {
                EvidenceError::InvalidInput("unknown panel summary contig".into())
            })?;
            let length = summary.target_length();
            let minimum = if summary.observed_loci < length {
                0
            } else {
                summary.minimum_callable_depth
            };
            write!(
                output,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t",
                summary.target.id,
                contig.name,
                summary.target.start,
                summary.target.end,
                length,
                summary.prefilter_depth_sum,
                summary.aligned_depth_sum,
                summary.callable_depth_sum,
                summary.base_quality_sum,
                summary.mapping_quality_sum,
                summary.callable_positions,
                self.min_callable_depth,
                minimum,
                summary.maximum_callable_depth
            )?;
            write_ratio(&mut output, summary.callable_depth_sum, length)?;
            for positions in summary.breadth_positions {
                write!(output, "\t")?;
                write_ratio(&mut output, positions, length)?;
            }
            writeln!(output)?;
        }
        output.flush()?;
        Ok(())
    }
}
impl EvidenceAnalyzer for PanelQcAnalyzer {
    fn requirements(&self) -> EvidenceRequirements {
        EvidenceRequirements {
            fields: EvidenceFields::DEPTHS.union(EvidenceFields::QUALITY_SUMS),
            requires_reference: false,
            context_bases: 0,
            retained_bytes: self.additional_memory_bytes(),
        }
    }
    fn on_batch(&mut self, batch: &EvidenceBatch) -> Result<(), EvidenceError> {
        for row in &batch.rows {
            let locus = (batch.contig_id, row.position);
            if self.last_locus.is_some_and(|previous| locus <= previous) {
                return Err(EvidenceError::Analyzer(
                    "panel input loci must be unique and canonically sorted".into(),
                ));
            }
            self.last_locus = Some(locus);
            self.active_targets.retain(|&index| {
                self.summaries[index].target.contig == batch.contig_id
                    && self.summaries[index].target.end > row.position
            });
            while let Some(&index) = self.sorted_targets.get(self.next_target) {
                let target = &self.summaries[index].target;
                if (target.contig, target.start) > locus {
                    break;
                }
                if target.contig == batch.contig_id && target.end > row.position {
                    self.active_targets.push(index);
                }
                self.next_target += 1;
            }
            for &index in &self.active_targets {
                let summary = &mut self.summaries[index];
                checked_add(&mut summary.observed_loci, 1)?;
                checked_add(&mut summary.prefilter_depth_sum, row.prefilter_depth)?;
                checked_add(&mut summary.aligned_depth_sum, row.aligned_depth)?;
                checked_add(&mut summary.callable_depth_sum, row.callable_depth)?;
                checked_add(&mut summary.base_quality_sum, row.base_quality_sum)?;
                checked_add(&mut summary.mapping_quality_sum, row.mapping_quality_sum)?;
                if row.callable_depth >= self.min_callable_depth {
                    checked_add(&mut summary.callable_positions, 1)?;
                }
                summary.minimum_callable_depth =
                    summary.minimum_callable_depth.min(row.callable_depth);
                summary.maximum_callable_depth =
                    summary.maximum_callable_depth.max(row.callable_depth);
                for (i, threshold) in [1, 10, 20, 30].iter().enumerate() {
                    if row.callable_depth >= *threshold {
                        checked_add(&mut summary.breadth_positions[i], 1)?;
                    }
                }
            }
        }
        Ok(())
    }
    fn finish(&mut self) -> Result<(), EvidenceError> {
        for summary in &mut self.summaries {
            if summary.observed_loci != summary.target_length() {
                return Err(EvidenceError::Analyzer(format!(
                    "panel target {} did not receive all requested loci",
                    summary.target.id
                )));
            }
            if summary.minimum_callable_depth == u64::MAX {
                summary.minimum_callable_depth = 0;
            }
        }
        Ok(())
    }
    fn additional_memory_bytes(&self) -> Option<u64> {
        Some(
            self.summaries
                .iter()
                .map(|summary| {
                    std::mem::size_of::<PanelSummary>() as u64 + summary.target.id.len() as u64 + 96
                })
                .sum(),
        )
    }
}
fn write_ratio(out: &mut dyn Write, numerator: u64, denominator: u64) -> Result<(), EvidenceError> {
    let scaled = (numerator as u128 * 1_000_000 + denominator as u128 / 2) / denominator as u128;
    write!(out, "{}.{:06}", scaled / 1_000_000, scaled % 1_000_000)?;
    Ok(())
}

/// Fan one extracted batch to compatible consumers without another BAM pass.
pub struct FusedAnalyzers<'a> {
    analyzers: Vec<&'a mut dyn EvidenceAnalyzer>,
}
impl<'a> FusedAnalyzers<'a> {
    /// Construct the configured value without starting evidence extraction.
    pub fn new(analyzers: Vec<&'a mut dyn EvidenceAnalyzer>) -> Self {
        Self { analyzers }
    }
}
impl std::fmt::Debug for FusedAnalyzers<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("FusedAnalyzers")
            .field("count", &self.analyzers.len())
            .finish()
    }
}
impl EvidenceAnalyzer for FusedAnalyzers<'_> {
    fn requirements(&self) -> EvidenceRequirements {
        let mut requirements = EvidenceRequirements {
            fields: EvidenceFields::from_bits(0).expect("empty capability set"),
            requires_reference: false,
            context_bases: 0,
            retained_bytes: Some(0),
        };
        for analyzer in &self.analyzers {
            let next = analyzer.requirements();
            requirements.fields = requirements.fields.union(next.fields);
            requirements.requires_reference |= next.requires_reference;
            requirements.context_bases = requirements.context_bases.max(next.context_bases);
            requirements.retained_bytes = requirements
                .retained_bytes
                .and_then(|bytes| next.retained_bytes.and_then(|next| bytes.checked_add(next)));
        }
        requirements
    }
    fn on_batch(&mut self, batch: &EvidenceBatch) -> Result<(), EvidenceError> {
        for analyzer in &mut self.analyzers {
            analyzer.on_batch(batch)?;
        }
        Ok(())
    }
    fn finish(&mut self) -> Result<(), EvidenceError> {
        for analyzer in &mut self.analyzers {
            analyzer.finish()?;
        }
        Ok(())
    }
    fn additional_memory_bytes(&self) -> Option<u64> {
        self.analyzers.iter().try_fold(0u64, |total, analyzer| {
            total.checked_add(analyzer.additional_memory_bytes()?)
        })
    }
}
