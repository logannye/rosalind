//! Actionable, read-only preflight diagnostics for a bounded analysis.

use std::path::{Path, PathBuf};

use crate::call::predicted_peak_rss_bytes;
use crate::contract::detected_os_memory_limit_bytes;
use crate::genomics::{AnalysisReference, ReferenceProvider};
use crate::io::bam::{find_bai, inspect_bam_header, StreamingBamSource};
use crate::pileup::ReadSource;
use crate::selection::AnalysisSelection;
use crate::util::rss::peak_rss_bytes;

const DEFAULT_MAX_DEPTH: u32 = 1000;
const DEFAULT_MAX_READ_LEN: u32 = 250;

/// Inputs and optional constraints checked by [`run_doctor`].
#[derive(Debug, Clone)]
pub struct DoctorSpec {
    /// Persisted reference index.
    pub index: PathBuf,
    /// Coordinate-sorted BAM.
    pub alignments: PathBuf,
    /// Planned output path, when known.
    pub output: Option<PathBuf>,
    /// Declared memory budget in MiB, when known.
    pub budget_mb: Option<u64>,
    /// Scan every mapped record to prove order and maximum read length.
    pub deep: bool,
}

/// Machine- and human-readable preflight result with concrete remediation.
#[derive(Debug, Clone)]
pub struct DoctorReport {
    /// Whether all requested checks passed.
    pub ok: bool,
    /// Whether index metadata and its memory-mapped reference could be read.
    pub index_readable: bool,
    /// Whether the BAM header opened and agreed with the index reference dictionary.
    pub alignments_readable: bool,
    /// Sort order declared by the BAM header.
    pub declared_sort_order: Option<String>,
    /// Whether a full scan proved nondecreasing coordinate order.
    pub coordinate_order_proven: Option<bool>,
    /// Largest read observed by a deep scan.
    pub observed_max_read_len: Option<u32>,
    /// Predicted process RSS for the default bounded analyzer contract.
    pub predicted_peak_rss_bytes: u64,
    /// Minimum whole-MiB budget that admits the prediction.
    pub required_budget_mb: u64,
    /// Whether the supplied budget admits the prediction.
    pub budget_feasible: Option<bool>,
    /// Largest default-model depth that fits the supplied budget, when useful.
    pub suggested_max_depth: Option<u32>,
    /// Whether output, partial, and receipt destinations are currently unused.
    pub output_safe: Option<bool>,
    /// Strongest enforcement assurance currently detectable.
    pub available_assurance: String,
    /// Concrete problems found.
    pub issues: Vec<String>,
    /// Suggested commands or configuration changes.
    pub remediation: Vec<String>,
}

impl DoctorReport {
    /// Stable JSON for CI and scheduler preflight.
    pub fn to_json(&self) -> String {
        let strings = |values: &[String]| {
            format!(
                "[{}]",
                values
                    .iter()
                    .map(|value| format!("\"{}\"", escape(value)))
                    .collect::<Vec<_>>()
                    .join(",")
            )
        };
        let optional_bool = |value: Option<bool>| match value {
            Some(value) => value.to_string(),
            None => "null".to_string(),
        };
        let optional_u32 =
            |value: Option<u32>| value.map_or_else(|| "null".to_string(), |v| v.to_string());
        let sort = self
            .declared_sort_order
            .as_ref()
            .map(|value| format!("\"{}\"", escape(value)))
            .unwrap_or_else(|| "null".to_string());
        format!(
            "{{\"schema\":1,\"ok\":{},\"index_readable\":{},\"alignments_readable\":{},\"declared_sort_order\":{sort},\"coordinate_order_proven\":{},\"observed_max_read_len\":{},\"configured_max_read_len\":{DEFAULT_MAX_READ_LEN},\"predicted_peak_rss_bytes\":{},\"required_budget_mb\":{},\"budget_feasible\":{},\"suggested_max_depth\":{},\"output_safe\":{},\"available_assurance\":\"{}\",\"issues\":{},\"remediation\":{}}}",
            self.ok,
            self.index_readable,
            self.alignments_readable,
            optional_bool(self.coordinate_order_proven),
            optional_u32(self.observed_max_read_len),
            self.predicted_peak_rss_bytes,
            self.required_budget_mb,
            optional_bool(self.budget_feasible),
            optional_u32(self.suggested_max_depth),
            optional_bool(self.output_safe),
            escape(&self.available_assurance),
            strings(&self.issues),
            strings(&self.remediation),
        )
    }
}

/// Inspect index/BAM compatibility, output safety, memory feasibility, and optional
/// deep record invariants without creating outputs.
pub fn run_doctor(spec: &DoctorSpec) -> DoctorReport {
    run_doctor_selected(spec, &AnalysisSelection::WholeGenome)
}

/// Inspect preflight requirements for one whole-genome, interval, or shard selection.
pub fn run_doctor_selected(spec: &DoctorSpec, selection: &AnalysisSelection) -> DoctorReport {
    let mut issues = Vec::new();
    let mut remediation = Vec::new();
    let mut report = DoctorReport {
        ok: false,
        index_readable: false,
        alignments_readable: false,
        declared_sort_order: None,
        coordinate_order_proven: None,
        observed_max_read_len: None,
        predicted_peak_rss_bytes: 0,
        required_budget_mb: 0,
        budget_feasible: None,
        suggested_max_depth: None,
        output_safe: None,
        available_assurance: "declared-bound-cooperative".to_string(),
        issues: Vec::new(),
        remediation: Vec::new(),
    };

    let loaded = match AnalysisReference::open(&spec.index) {
        Ok(index) => {
            report.index_readable = true;
            index
        }
        Err(error) => {
            issues.push(format!(
                "analysis reference {} is not readable: {error}",
                spec.index.display()
            ));
            remediation.push(
                "Build an analysis reference with `rosalind reference build --fasta REF --output REF.rref`."
                    .to_string(),
            );
            report.issues = issues;
            report.remediation = remediation;
            return report;
        }
    };
    let contigs = loaded.contigs();

    if !matches!(selection, AnalysisSelection::WholeGenome) && find_bai(&spec.alignments).is_none()
    {
        issues.push(format!(
            "sparse analysis requires a BAM index at {}.bai or {}",
            spec.alignments.display(),
            spec.alignments.with_extension("bai").display()
        ));
        remediation.push("Create a BAI with `samtools index ALIGNMENTS.bam`.".to_string());
    }

    match inspect_bam_header(&spec.alignments, contigs) {
        Ok(header) => {
            report.alignments_readable = true;
            report.declared_sort_order = header.sort_order.clone();
            if header.sort_order.as_deref() != Some("coordinate") && !spec.deep {
                issues.push(format!(
                    "BAM header declares sort order {:?}, not coordinate",
                    header.sort_order.as_deref().unwrap_or("unknown")
                ));
                remediation.push(format!(
                    "Run `rosalind sort --input {} --output sorted.bam` and use the sorted BAM.",
                    spec.alignments.display()
                ));
            }
        }
        Err(error) => {
            issues.push(format!(
                "alignments {} do not match the index: {error}",
                spec.alignments.display()
            ));
            remediation.push("Confirm the BAM was aligned to the exact reference used to build the index, including contig names and lengths.".to_string());
        }
    }

    if spec.deep && report.alignments_readable {
        match deep_scan(&spec.alignments, contigs) {
            Ok(max_len) => {
                report.coordinate_order_proven = Some(true);
                report.observed_max_read_len = Some(max_len);
                if max_len > DEFAULT_MAX_READ_LEN {
                    issues.push(format!(
                        "observed read length {max_len} exceeds the default enforced maximum {DEFAULT_MAX_READ_LEN}"
                    ));
                    remediation.push(format!(
                        "Pass `--max-read-len {max_len}` and re-run planning, or pre-filter longer records."
                    ));
                }
            }
            Err(error) => {
                report.coordinate_order_proven = Some(false);
                issues.push(format!("deep BAM scan failed: {error}"));
                remediation.push(
                    "Coordinate-sort the BAM and repeat `rosalind doctor --deep`.".to_string(),
                );
            }
        }
    }

    let largest = selection.largest_reference_span(contigs);
    let baseline = peak_rss_bytes();
    report.predicted_peak_rss_bytes =
        predicted_peak_rss_bytes(largest, DEFAULT_MAX_DEPTH, DEFAULT_MAX_READ_LEN, baseline);
    report.required_budget_mb = report.predicted_peak_rss_bytes.div_ceil(1 << 20);
    if let Some(budget_mb) = spec.budget_mb {
        let feasible = report.predicted_peak_rss_bytes <= budget_mb.saturating_mul(1 << 20);
        report.budget_feasible = Some(feasible);
        if !feasible {
            report.suggested_max_depth = max_depth_that_fits(largest, baseline, budget_mb);
            issues.push(format!(
                "the default depth contract predicts {} MiB but the budget is {budget_mb} MiB",
                report.required_budget_mb
            ));
            remediation.push(match report.suggested_max_depth {
                Some(depth) if depth > 0 => format!(
                    "Raise `--budget-mb` to at least {} or reduce `--max-depth` to {depth}.",
                    report.required_budget_mb
                ),
                _ => format!(
                    "Raise `--budget-mb` to at least {}; reducing depth alone cannot fit the fixed reference and process cost.",
                    report.required_budget_mb
                ),
            });
        }
    }

    if let Some(output) = &spec.output {
        let destinations = [
            output.clone(),
            PathBuf::from(format!("{}.partial", output.display())),
            PathBuf::from(format!("{}.manifest.json", output.display())),
        ];
        let collisions = destinations
            .iter()
            .filter(|path| path.exists())
            .collect::<Vec<_>>();
        report.output_safe = Some(collisions.is_empty());
        if !collisions.is_empty() {
            issues.push(format!(
                "output reservation collides with {}",
                collisions
                    .iter()
                    .map(|path| path.display().to_string())
                    .collect::<Vec<_>>()
                    .join(", ")
            ));
            remediation.push("Choose a new output path, remove stale files deliberately, or use `--force` for atomic replacement.".to_string());
        }
    }

    if let Some(limit) = detected_os_memory_limit_bytes() {
        report.available_assurance = "cgroup-v2".to_string();
        if let Some(budget) = spec.budget_mb {
            if limit > budget.saturating_mul(1 << 20) {
                remediation.push(format!(
                    "For `--require-os-limit`, lower cgroup v2 memory.max from {} MiB to at most {budget} MiB.",
                    limit / (1 << 20)
                ));
            }
        }
    } else {
        remediation.push("For hard Linux assurance, run inside a cgroup-v2 container or systemd scope with memory.max at or below the declared budget.".to_string());
    }

    report.ok = issues.is_empty();
    report.issues = issues;
    report.remediation = remediation;
    report
}

fn deep_scan(path: &Path, contigs: &crate::core::ContigSet) -> Result<u32, crate::core::CoreError> {
    let mut source = StreamingBamSource::new(path, contigs)?;
    let mut max_len = 0_u32;
    while let Some(read) = source.next_read()? {
        max_len = max_len.max(read.seq.len().try_into().unwrap_or(u32::MAX));
    }
    Ok(max_len)
}

fn max_depth_that_fits(largest: u64, baseline: u64, budget_mb: u64) -> Option<u32> {
    let budget = budget_mb.saturating_mul(1 << 20);
    if predicted_peak_rss_bytes(largest, 1, DEFAULT_MAX_READ_LEN, baseline) > budget {
        return None;
    }
    let mut low = 1_u32;
    let mut high = DEFAULT_MAX_DEPTH;
    while low < high {
        let mid = low + (high - low).div_ceil(2);
        if predicted_peak_rss_bytes(largest, mid, DEFAULT_MAX_READ_LEN, baseline) <= budget {
            low = mid;
        } else {
            high = mid - 1;
        }
    }
    Some(low)
}

fn escape(value: &str) -> String {
    let mut output = String::new();
    for character in value.chars() {
        match character {
            '"' => output.push_str("\\\""),
            '\\' => output.push_str("\\\\"),
            '\n' => output.push_str("\\n"),
            '\r' => output.push_str("\\r"),
            '\t' => output.push_str("\\t"),
            value if (value as u32) < 0x20 => {
                output.push_str(&format!("\\u{:04x}", value as u32));
            }
            value => output.push(value),
        }
    }
    output
}
