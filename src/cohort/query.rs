//! Bounded metadata-only planning for saved-evidence candidate queries.
//!
//! Plans retain one normalized query and per-member/leaf counts, never a dense
//! sample-by-locus matrix. Payload/REF validation remains the strict reader's job.

use super::comparison::{
    check_required_fields, require_named_sample, ComparisonContract, ComparisonMismatch,
    ComparisonMismatchCode,
};
use super::descriptor::CohortLimits;
use super::store::SnapshotHandle;
use super::{CohortError, Result};
use crate::core::ContigSet;
use crate::dataset::{DatasetCoverage, DatasetQuery, DatasetReadPlan};
use crate::evidence::{
    EvidenceAnalyzer, EvidenceBatch, EvidenceError, EvidenceExecution, EvidenceFields,
    EvidenceRequirements, EvidenceSelection, SnvSite,
};
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub(crate) enum MissingPolicy {
    #[default]
    Strict,
    Partial,
}

#[derive(Debug, Clone)]
pub(crate) struct CohortQuery {
    pub selection: EvidenceSelection,
    /// None selects every member; Some([]) selects none. Duplicate/unknown IDs refuse.
    pub member_ids: Option<Vec<String>>,
    pub fields: EvidenceFields,
    pub missing_policy: MissingPolicy,
    pub requirements: EvidenceRequirements,
}
impl CohortQuery {
    pub fn new(selection: EvidenceSelection) -> Self {
        let fields = EvidenceFields::DEPTHS.union(EvidenceFields::ALLELES);
        Self {
            selection,
            member_ids: None,
            fields,
            missing_policy: MissingPolicy::Strict,
            requirements: EvidenceRequirements {
                fields,
                requires_reference: true,
                context_bases: 0,
                retained_bytes: Some(0),
            },
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct CohortQueryLimits {
    pub cohort: CohortLimits,
    /// Bounds the already-owned input before normalization allocates copies.
    pub max_input_sites: usize,
    pub max_query_bytes: u64,
    pub max_plan_bytes: u64,
    /// Retain at most this many explanations; issue_count still records all failures.
    pub max_issues: usize,
}
impl Default for CohortQueryLimits {
    fn default() -> Self {
        Self {
            cohort: CohortLimits::default(),
            max_input_sites: 1_000_000,
            max_query_bytes: 128 << 20,
            max_plan_bytes: 64 << 20,
            max_issues: 1024,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "kebab-case")]
pub(crate) enum CohortQueryIssue {
    Request {
        mismatch: ComparisonMismatch,
    },
    Comparison {
        member_id: String,
        object_id: String,
        mismatch: ComparisonMismatch,
    },
    MissingCoverage {
        member_id: String,
        missing_loci: u64,
        missing_rows: u64,
    },
}

#[derive(Debug, Clone)]
pub(crate) struct LeafQueryPlan {
    pub leaf_index: usize,
    pub object_id: String,
    /// None means the dictionary is incompatible, not that a locus was unmeasured.
    pub covered_loci: Option<u64>,
    pub covered_rows: Option<u64>,
    pub read_plan: Option<DatasetReadPlan>,
}
#[derive(Debug, Clone)]
pub(crate) struct MemberQueryPlan {
    pub member_index: usize,
    pub member_id: String,
    pub covered_loci: Option<u64>,
    pub missing_loci: Option<u64>,
    pub covered_rows: Option<u64>,
    pub missing_rows: Option<u64>,
    pub leaves: Vec<LeafQueryPlan>,
}

#[derive(Debug, Clone)]
pub(crate) struct CohortQueryReservations {
    pub memory_budget_bytes: Option<u64>,
    pub baseline_rss_bytes: u64,
    pub query_bytes: u64,
    pub plan_bytes: u64,
    pub consumer_bytes: u64,
    pub consumer_bound_known: bool,
    pub max_reader_metadata_bytes: u64,
    pub max_source_decoder_bytes: u64,
    pub max_projection_bytes: u64,
    /// Serial maximum, never a sum of per-reader process RSS baselines.
    pub predicted_peak_rss_bytes: u64,
}
impl CohortQueryReservations {
    pub fn retained_bytes(&self) -> u64 {
        self.query_bytes
            .saturating_add(self.plan_bytes)
            .saturating_add(self.consumer_bytes)
    }
}

#[derive(Debug, Clone)]
pub(crate) struct CohortQueryPlan {
    pub snapshot_id: String,
    pub contigs: ContigSet,
    pub sites: Vec<SnvSite>,
    pub fields: EvidenceFields,
    pub missing_policy: MissingPolicy,
    pub requirements: EvidenceRequirements,
    pub members: Vec<MemberQueryPlan>,
    pub requested_loci: u64,
    /// Normalized locus/ALT pairs, before multiplication by selected members.
    pub candidate_rows: u64,
    pub output_rows: u64,
    pub comparison_blake3: Option<String>,
    pub issues: Vec<CohortQueryIssue>,
    pub issue_count: u64,
    pub reservations: CohortQueryReservations,
}
impl CohortQueryPlan {
    pub fn ensure_executable(&self) -> Result<()> {
        if self.issue_count != 0 {
            return Err(CohortError::Incompatible(format!(
                "cohort query has {} blocking issue(s); {} explanation(s) retained",
                self.issue_count,
                self.issues.len()
            )));
        }
        Ok(())
    }

    /// Compute one member's missingness details on demand; never retain this for
    /// all members. Coverage is ownership, not positive read depth. Requested REF
    /// is validated against payload rows only when the reader consumes them.
    pub fn member_coverage(
        &self,
        snapshot: &SnapshotHandle,
        planned_member: usize,
        mut limits: CohortQueryLimits,
    ) -> Result<DatasetCoverage> {
        limits.cohort.validate()?;
        if self.sites.len() > limits.max_input_sites
            || self.reservations.query_bytes > limits.max_query_bytes
        {
            return Err(CohortError::Limit(
                "on-demand coverage exceeds query envelope".into(),
            ));
        }
        if snapshot.id != self.snapshot_id {
            return Err(CohortError::Incompatible(
                "plan belongs to a different snapshot".into(),
            ));
        }
        let member = self
            .members
            .get(planned_member)
            .ok_or_else(|| CohortError::Incompatible("unknown planned member index".into()))?;
        if member.covered_loci.is_none() {
            return Err(CohortError::Incompatible(
                "member reference dictionary is incompatible".into(),
            ));
        }
        limits.cohort.memory_budget_bytes = minimum_budget([
            limits.cohort.memory_budget_bytes,
            limits.cohort.dataset.memory_budget_bytes,
            self.reservations.memory_budget_bytes,
        ]);
        limits.cohort.admit(self.reservations.query_bytes)?;
        snapshot.verify_unchanged()?;
        let mut represented = vec![false; self.sites.len()];
        let source = &snapshot.descriptor.members[member.member_index];
        for leaf in &source.leaves {
            let dataset = snapshot.open_leaf(leaf, limits.cohort)?;
            let intervals = &dataset.descriptor().selection.intervals;
            let mut cursor = 0;
            for (index, site) in self.sites.iter().enumerate() {
                while cursor < intervals.len()
                    && (intervals[cursor].contig < site.contig
                        || (intervals[cursor].contig == site.contig
                            && intervals[cursor].end <= site.position))
                {
                    cursor += 1;
                }
                if intervals.get(cursor).is_some_and(|interval| {
                    interval.contig == site.contig
                        && interval.start <= site.position
                        && site.position < interval.end
                }) {
                    if represented[index] {
                        return Err(CohortError::Corrupt(
                            "overlapping cohort leaf ownership".into(),
                        ));
                    }
                    represented[index] = true;
                }
            }
            dataset.verify_unchanged()?;
        }
        let mut covered = Vec::new();
        let mut missing = Vec::new();
        for (site, present) in self.sites.iter().zip(represented) {
            if present {
                covered.push(site.clone());
            } else {
                missing.push(site.clone());
            }
        }
        snapshot.verify_unchanged()?;
        Ok(DatasetCoverage {
            covered_loci: covered.len() as u64,
            missing_loci: missing.len() as u64,
            covered: EvidenceSelection::Sites(covered),
            missing: EvidenceSelection::Sites(missing),
        })
    }
}

struct PlanConsumer(EvidenceRequirements);
impl EvidenceAnalyzer for PlanConsumer {
    fn requirements(&self) -> EvidenceRequirements {
        self.0.clone()
    }
    fn on_batch(&mut self, _: &EvidenceBatch) -> std::result::Result<(), EvidenceError> {
        Ok(())
    }
}
fn minimum_budget(values: [Option<u64>; 3]) -> Option<u64> {
    values.into_iter().flatten().min()
}
fn count_rows(sites: &[SnvSite]) -> Result<u64> {
    sites.iter().try_fold(0u64, |sum, site| {
        sum.checked_add(site.alternates.len() as u64)
            .ok_or_else(|| CohortError::Limit("candidate row count overflow".into()))
    })
}
fn add_issue(
    plan: &mut CohortQueryPlan,
    issue: CohortQueryIssue,
    limits: CohortQueryLimits,
) -> Result<()> {
    plan.issue_count = plan
        .issue_count
        .checked_add(1)
        .ok_or_else(|| CohortError::Limit("issue count overflow".into()))?;
    if plan.issues.len() < limits.max_issues {
        let bytes = match &issue {
            CohortQueryIssue::Request { mismatch } => {
                mismatch.expected.len() + mismatch.observed.len() + 512
            }
            CohortQueryIssue::Comparison {
                member_id,
                object_id,
                mismatch,
            } => {
                member_id.len()
                    + object_id.len()
                    + mismatch.expected.len()
                    + mismatch.observed.len()
                    + 512
            }
            CohortQueryIssue::MissingCoverage { member_id, .. } => member_id.len() + 256,
        } as u64;
        plan.reservations.plan_bytes = plan
            .reservations
            .plan_bytes
            .checked_add(bytes)
            .ok_or_else(|| CohortError::Limit("plan metadata overflow".into()))?;
        if plan.reservations.plan_bytes > limits.max_plan_bytes {
            return Err(CohortError::Limit(
                "cohort query plan exceeds metadata envelope".into(),
            ));
        }
        limits.cohort.admit(plan.reservations.retained_bytes())?;
        plan.issues.push(issue);
    }
    Ok(())
}

pub(crate) fn plan_query(
    snapshot: &SnapshotHandle,
    query: &CohortQuery,
    execution: &EvidenceExecution,
    mut limits: CohortQueryLimits,
) -> Result<CohortQueryPlan> {
    limits.cohort.validate()?;
    if limits.max_input_sites == 0
        || limits.max_query_bytes == 0
        || limits.max_plan_bytes == 0
        || limits.max_issues == 0
    {
        return Err(CohortError::Limit(
            "query envelopes must be positive".into(),
        ));
    }
    let sites = match &query.selection {
        EvidenceSelection::Sites(sites) => sites,
        _ => {
            return Err(CohortError::Incompatible(
                "cohort candidate queries require an SNV site selection".into(),
            ))
        }
    };
    if sites.len() > limits.max_input_sites {
        return Err(CohortError::Limit(
            "candidate input exceeds site envelope".into(),
        ));
    }
    // Covers normalization maps, temporary per-leaf coverage splits and retained
    // canonical sites. The caller must also bound parsing before constructing input.
    let query_bytes = sites.iter().try_fold(64u64 << 10, |sum, site| {
        sum.checked_add(2048)
            .and_then(|sum| sum.checked_add((site.alternates.len() as u64).saturating_mul(16)))
            .ok_or_else(|| CohortError::Limit("query metadata overflow".into()))
    })?;
    if query_bytes > limits.max_query_bytes {
        return Err(CohortError::Limit(
            "candidate input exceeds query byte envelope".into(),
        ));
    }
    limits.cohort.memory_budget_bytes = minimum_budget([
        limits.cohort.memory_budget_bytes,
        limits.cohort.dataset.memory_budget_bytes,
        execution.memory_budget_bytes,
    ]);
    if query.requirements.context_bases != 0 {
        return Err(CohortError::Incompatible(
            "cohort candidate queries support zero flanking context only".into(),
        ));
    }
    let requirements_fields_ok = query.fields.contains(query.requirements.fields);
    if query.requirements.retained_bytes.is_none() && limits.cohort.memory_budget_bytes.is_some() {
        return Err(CohortError::Incompatible(
            "budgeted cohort consumer requires a declared memory bound".into(),
        ));
    }
    snapshot.verify_unchanged()?;
    let baseline = crate::util::rss::peak_rss_bytes();
    limits.cohort.admit(query_bytes)?;
    let members = &snapshot.descriptor.members;
    let selection_count = query.member_ids.as_ref().map_or(members.len(), Vec::len);
    if selection_count > limits.cohort.max_members {
        return Err(CohortError::Limit(
            "member selection exceeds envelope".into(),
        ));
    }
    let selection_bytes = (selection_count as u64).saturating_mul(512);
    if selection_bytes > limits.max_plan_bytes {
        return Err(CohortError::Limit(
            "member selection exceeds plan metadata envelope".into(),
        ));
    }
    limits
        .cohort
        .admit(query_bytes.saturating_add(selection_bytes))?;
    let indices = match &query.member_ids {
        None => (0..members.len()).collect::<Vec<_>>(),
        Some(ids) => {
            if ids.len() > limits.cohort.max_members {
                return Err(CohortError::Limit(
                    "member selection exceeds envelope".into(),
                ));
            }
            let mut indices = BTreeSet::new();
            for id in ids {
                if id.len() > 256 {
                    return Err(CohortError::Limit("member ID exceeds envelope".into()));
                }
                let index = members
                    .binary_search_by(|member| member.metadata.id.cmp(id))
                    .map_err(|_| {
                        CohortError::Incompatible(format!("unknown cohort member: {id}"))
                    })?;
                if !indices.insert(index) {
                    return Err(CohortError::Incompatible(format!(
                        "duplicate requested member: {id}"
                    )));
                }
            }
            indices.into_iter().collect()
        }
    };
    let leaf_count = indices.iter().try_fold(0usize, |sum, index| {
        sum.checked_add(members[*index].leaves.len())
            .ok_or_else(|| CohortError::Limit("selected leaf count overflow".into()))
    })?;
    if indices.len() > limits.cohort.max_members || leaf_count > limits.cohort.max_leaf_references {
        return Err(CohortError::Limit(
            "selected cohort inventory exceeds envelope".into(),
        ));
    }
    let mut plan_bytes = (indices.len() as u64)
        .saturating_mul(1024)
        .saturating_add((leaf_count as u64).saturating_mul(512))
        .saturating_add(64 << 10);
    if plan_bytes > limits.max_plan_bytes {
        return Err(CohortError::Limit(
            "cohort query plan exceeds metadata envelope".into(),
        ));
    }
    limits
        .cohort
        .admit(query_bytes.saturating_add(plan_bytes))?;
    let dictionary_member = indices
        .first()
        .copied()
        .or_else(|| (!members.is_empty()).then_some(0));
    let contigs = if let Some(index) = dictionary_member {
        let dataset = snapshot.open_leaf(&members[index].leaves[0], limits.cohort)?;
        let bytes = dataset.descriptor().contigs.iter().fold(0u64, |sum, c| {
            sum.saturating_add(c.name.len() as u64).saturating_add(128)
        });
        plan_bytes = plan_bytes.saturating_add(bytes.saturating_mul(16));
        if plan_bytes > limits.max_plan_bytes {
            return Err(CohortError::Limit(
                "comparison dictionary exceeds plan metadata envelope".into(),
            ));
        }
        limits
            .cohort
            .admit(query_bytes.saturating_add(plan_bytes))?;
        dataset.descriptor().contig_set()
    } else {
        if !sites.is_empty() {
            return Err(CohortError::Incompatible(
                "a nonempty candidate query needs a cohort reference dictionary".into(),
            ));
        }
        ContigSet::new()
    };
    let (_, normalized) = query.selection.normalize(&contigs)?;
    let normalized = normalized.into_values().collect::<Vec<_>>();
    let requested_loci = normalized.len() as u64;
    let candidate_rows = count_rows(&normalized)?;
    let output_rows = candidate_rows
        .checked_mul(indices.len() as u64)
        .ok_or_else(|| CohortError::Limit("cohort output row count overflow".into()))?;
    let consumer_bytes = query
        .requirements
        .retained_bytes
        .unwrap_or(0)
        .max(execution.analyzer_bytes);
    let mut plan = CohortQueryPlan {
        snapshot_id: snapshot.id.clone(),
        contigs,
        sites: normalized,
        fields: query.fields,
        missing_policy: query.missing_policy,
        requirements: query.requirements.clone(),
        members: Vec::with_capacity(indices.len()),
        requested_loci,
        candidate_rows,
        output_rows,
        comparison_blake3: None,
        issues: Vec::new(),
        issue_count: 0,
        reservations: CohortQueryReservations {
            memory_budget_bytes: limits.cohort.memory_budget_bytes,
            baseline_rss_bytes: baseline,
            query_bytes,
            plan_bytes,
            consumer_bytes,
            consumer_bound_known: query.requirements.retained_bytes.is_some(),
            max_reader_metadata_bytes: 0,
            max_source_decoder_bytes: 0,
            max_projection_bytes: 0,
            predicted_peak_rss_bytes: baseline
                .saturating_add(query_bytes)
                .saturating_add(plan_bytes)
                .saturating_add(consumer_bytes),
        },
    };
    limits.cohort.admit(plan.reservations.retained_bytes())?;
    if !requirements_fields_ok {
        add_issue(
            &mut plan,
            CohortQueryIssue::Request {
                mismatch: ComparisonMismatch {
                    code: ComparisonMismatchCode::RequiredFields,
                    expected: query.requirements.fields.names().join(","),
                    observed: query.fields.names().join(","),
                },
            },
            limits,
        )?;
    }
    let mut expected: Option<ComparisonContract> = None;
    let mut execution = execution.clone();
    execution.memory_budget_bytes = limits.cohort.memory_budget_bytes;
    execution.analyzer_bytes = 0; // Outer reservation below already contains the caller's reserve.
    for member_index in indices {
        let member = &members[member_index];
        let mut member_plan = MemberQueryPlan {
            member_index,
            member_id: member.metadata.id.clone(),
            covered_loci: Some(0),
            missing_loci: None,
            covered_rows: Some(0),
            missing_rows: None,
            leaves: Vec::with_capacity(member.leaves.len()),
        };
        for (leaf_index, leaf) in member.leaves.iter().enumerate() {
            let dataset = snapshot.open_leaf(leaf, limits.cohort)?;
            let descriptor = dataset.descriptor();
            let mut differences = Vec::new();
            match ComparisonContract::from_descriptor(descriptor) {
                Ok(actual) => {
                    if let Some(reference) = &expected {
                        differences.extend(reference.compare(&actual));
                    } else {
                        plan.comparison_blake3 = Some(actual.digest()?);
                        expected = Some(actual);
                    }
                }
                Err(error) => differences.push(ComparisonMismatch {
                    code: ComparisonMismatchCode::ReferenceContent,
                    expected: "stored effective analysis reference identity".into(),
                    observed: error.to_string(),
                }),
            }
            if let Err(issue) = require_named_sample(descriptor) {
                differences.push(issue);
            }
            // Compare against the normalization dictionary even if the first leaf
            // lacked real reference bases and could not establish a contract.
            let dictionary_matches = descriptor.contigs.len() == plan.contigs.len()
                && descriptor
                    .contigs
                    .iter()
                    .zip(plan.contigs.iter())
                    .all(|(a, b)| {
                        a.id == b.id && a.name.as_str() == b.name.as_ref() && a.length == b.length
                    });
            if !dictionary_matches
                && !differences
                    .iter()
                    .any(|d| d.code == ComparisonMismatchCode::ReferenceDictionary)
            {
                differences.push(ComparisonMismatch {
                    code: ComparisonMismatchCode::ReferenceDictionary,
                    expected: "the query normalization dictionary".into(),
                    observed: "a different ordered dictionary".into(),
                });
            }
            for mismatch in differences {
                add_issue(
                    &mut plan,
                    CohortQueryIssue::Comparison {
                        member_id: member.metadata.id.clone(),
                        object_id: leaf.object_id.clone(),
                        mismatch,
                    },
                    limits,
                )?;
            }
            let mut leaf_plan = LeafQueryPlan {
                leaf_index,
                object_id: leaf.object_id.clone(),
                covered_loci: None,
                covered_rows: None,
                read_plan: None,
            };
            if dictionary_matches {
                let split =
                    dataset.coverage_split(&EvidenceSelection::Sites(plan.sites.clone()))?;
                let covered_rows = match &split.covered {
                    EvidenceSelection::Sites(sites) => count_rows(sites)?,
                    _ => unreachable!(),
                };
                leaf_plan.covered_loci = Some(split.covered_loci);
                leaf_plan.covered_rows = Some(covered_rows);
                if let Some(n) = member_plan.covered_loci {
                    member_plan.covered_loci =
                        Some(n.checked_add(split.covered_loci).ok_or_else(|| {
                            CohortError::Limit("covered locus count overflow".into())
                        })?);
                }
                if let Some(n) = member_plan.covered_rows {
                    member_plan.covered_rows =
                        Some(n.checked_add(covered_rows).ok_or_else(|| {
                            CohortError::Limit("covered row count overflow".into())
                        })?);
                }
                if split.covered_loci != 0 {
                    let fields_ok = match check_required_fields(descriptor, query.fields) {
                        Ok(()) => true,
                        Err(mismatch) => {
                            add_issue(
                                &mut plan,
                                CohortQueryIssue::Comparison {
                                    member_id: member.metadata.id.clone(),
                                    object_id: leaf.object_id.clone(),
                                    mismatch,
                                },
                                limits,
                            )?;
                            false
                        }
                    };
                    if fields_ok && requirements_fields_ok && descriptor.has_reference {
                        let mut requirements = query.requirements.clone();
                        requirements.retained_bytes = Some(plan.reservations.retained_bytes());
                        let read = dataset.plan(
                            &DatasetQuery {
                                selection: split.covered,
                                fields: query.fields,
                            },
                            &PlanConsumer(requirements),
                            &execution,
                        )?;
                        plan.reservations.max_reader_metadata_bytes = plan
                            .reservations
                            .max_reader_metadata_bytes
                            .max(read.metadata_bytes);
                        plan.reservations.max_source_decoder_bytes = plan
                            .reservations
                            .max_source_decoder_bytes
                            .max(read.source_decoder_bytes);
                        plan.reservations.max_projection_bytes = plan
                            .reservations
                            .max_projection_bytes
                            .max(read.projection_bytes);
                        plan.reservations.predicted_peak_rss_bytes = plan
                            .reservations
                            .predicted_peak_rss_bytes
                            .max(read.predicted_peak_rss_bytes);
                        leaf_plan.read_plan = Some(read);
                    }
                }
            } else {
                member_plan.covered_loci = None;
                member_plan.covered_rows = None;
            }
            dataset.verify_unchanged()?;
            member_plan.leaves.push(leaf_plan);
        }
        if let (Some(covered), Some(rows)) = (member_plan.covered_loci, member_plan.covered_rows) {
            let missing = requested_loci
                .checked_sub(covered)
                .ok_or_else(|| CohortError::Corrupt("cohort ownership counts overlap".into()))?;
            let missing_rows = candidate_rows.checked_sub(rows).ok_or_else(|| {
                CohortError::Corrupt("cohort candidate row counts overlap".into())
            })?;
            member_plan.missing_loci = Some(missing);
            member_plan.missing_rows = Some(missing_rows);
            if missing != 0 && query.missing_policy == MissingPolicy::Strict {
                add_issue(
                    &mut plan,
                    CohortQueryIssue::MissingCoverage {
                        member_id: member.metadata.id.clone(),
                        missing_loci: missing,
                        missing_rows,
                    },
                    limits,
                )?;
            }
        }
        plan.members.push(member_plan);
    }
    snapshot.verify_unchanged()?;
    limits.cohort.admit(plan.reservations.retained_bytes())?;
    plan.reservations.predicted_peak_rss_bytes = plan.reservations.predicted_peak_rss_bytes.max(
        baseline
            .saturating_add(plan.reservations.retained_bytes())
            .saturating_add(plan.reservations.max_reader_metadata_bytes)
            .saturating_add(plan.reservations.max_source_decoder_bytes)
            .saturating_add(plan.reservations.max_projection_bytes),
    );
    if limits
        .cohort
        .memory_budget_bytes
        .is_some_and(|budget| plan.reservations.predicted_peak_rss_bytes > budget)
    {
        return Err(CohortError::Limit(
            "cohort query serial working set exceeds admitted budget".into(),
        ));
    }
    Ok(plan)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cohort::store::create_snapshot;
    use crate::cohort::tests::{Fixture, FixtureOptions};
    use crate::evidence::EvidenceProfile;
    use std::path::PathBuf;
    use std::sync::atomic::{AtomicU64, Ordering};

    static ID: AtomicU64 = AtomicU64::new(0);
    struct Root(PathBuf);
    impl Root {
        fn new() -> Self {
            Self(std::env::temp_dir().join(format!(
                "rosalind-cohort-query-{}-{}",
                std::process::id(),
                ID.fetch_add(1, Ordering::Relaxed)
            )))
        }
    }
    impl Drop for Root {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }
    fn site(position: u32, alternates: &[u8]) -> SnvSite {
        SnvSite {
            contig: 0,
            position,
            reference: b'A',
            alternates: alternates.to_vec(),
        }
    }
    fn query(sites: Vec<SnvSite>) -> CohortQuery {
        CohortQuery::new(EvidenceSelection::Sites(sites))
    }
    fn plan(snapshot: &SnapshotHandle, query: &CohortQuery) -> CohortQueryPlan {
        plan_query(
            snapshot,
            query,
            &EvidenceExecution::default(),
            CohortQueryLimits::default(),
        )
        .unwrap()
    }

    #[test]
    fn strict_partial_counts_canonical_members_and_zero_coverage_are_distinct() {
        let a = Fixture::new(FixtureOptions::default());
        let b = Fixture::new(FixtureOptions {
            sample: Some("B".into()),
            fields: EvidenceFields::FULL_V1,
            ..FixtureOptions::default()
        });
        let mut second = b.member("B");
        second
            .manifests
            .push(b.extra_leaf(8, 10, EvidenceFields::FULL_V1));
        let root = Root::new();
        let snapshot = create_snapshot(
            &root.0,
            &[second, a.member("A")],
            None,
            CohortLimits::default(),
        )
        .unwrap();
        drop(a);
        drop(b); // Every query below must work without original source files.
        let mut request = query(vec![
            site(3, b"C"),
            site(1, b"GC"),
            site(1, b"C"),
            site(2, b"C"),
            site(8, b"C"),
        ]);
        request.member_ids = Some(vec!["B".into(), "A".into()]);
        let strict = plan(&snapshot, &request);
        assert_eq!(
            strict
                .members
                .iter()
                .map(|m| m.member_id.as_str())
                .collect::<Vec<_>>(),
            ["A", "B"]
        );
        assert_eq!(
            (
                strict.requested_loci,
                strict.candidate_rows,
                strict.output_rows
            ),
            (4, 5, 10)
        );
        assert_eq!(strict.sites[0].alternates, b"CG");
        assert_eq!(strict.issue_count, 1);
        assert!(strict.ensure_executable().is_err());
        assert!(
            matches!(&strict.issues[0], CohortQueryIssue::MissingCoverage { member_id, missing_loci: 1, missing_rows: 1 } if member_id == "A")
        );
        assert_eq!(
            (
                strict.members[0].covered_loci,
                strict.members[0].covered_rows
            ),
            (Some(3), Some(4))
        );
        request.missing_policy = MissingPolicy::Partial;
        let partial = plan(&snapshot, &request);
        partial.ensure_executable().unwrap();
        assert_eq!(partial.members[0].missing_loci, Some(1));
        assert_eq!(partial.members[1].missing_loci, Some(0));
        let detail = partial
            .member_coverage(&snapshot, 0, CohortQueryLimits::default())
            .unwrap();
        assert_eq!((detail.covered_loci, detail.missing_loci), (3, 1));
        if let EvidenceSelection::Sites(covered) = detail.covered {
            // Fixture position 2 has no reads, but its stored zero is represented.
            assert!(covered.iter().any(|site| site.position == 2));
        } else {
            panic!("expected SNV coverage");
        }
        let readers = partial
            .members
            .iter()
            .flat_map(|m| &m.leaves)
            .filter_map(|leaf| leaf.read_plan.as_ref())
            .collect::<Vec<_>>();
        assert_eq!(
            partial.reservations.max_source_decoder_bytes,
            readers
                .iter()
                .map(|p| p.source_decoder_bytes)
                .max()
                .unwrap()
        );
        assert!(readers
            .iter()
            .any(|p| p.source_fields == EvidenceFields::FULL_V1));
    }

    #[test]
    fn fields_are_required_only_on_consumed_leaves_and_never_partial_nulls() {
        let a = Fixture::new(FixtureOptions::default());
        let mut member = a.member("A");
        member
            .manifests
            .push(a.extra_leaf(8, 10, EvidenceFields::DEPTHS));
        let root = Root::new();
        let snapshot = create_snapshot(&root.0, &[member], None, CohortLimits::default()).unwrap();
        plan(&snapshot, &query(vec![site(1, b"C")]))
            .ensure_executable()
            .unwrap();
        let mut request = query(vec![site(8, b"C")]);
        request.missing_policy = MissingPolicy::Partial;
        let missing_fields = plan(&snapshot, &request);
        assert_eq!(missing_fields.members[0].covered_loci, Some(1));
        assert_eq!(missing_fields.members[0].missing_loci, Some(0));
        assert!(missing_fields.ensure_executable().is_err());
        assert!(missing_fields.issues.iter().any(|issue| matches!(issue,
            CohortQueryIssue::Comparison { mismatch, .. } if mismatch.code == ComparisonMismatchCode::RequiredFields)));
        request.fields = EvidenceFields::DEPTHS;
        let omitted_requirement = plan(&snapshot, &request);
        assert!(
            matches!(&omitted_requirement.issues[0], CohortQueryIssue::Request { mismatch }
            if mismatch.code == ComparisonMismatchCode::RequiredFields)
        );
    }

    #[test]
    fn incompatible_profiles_references_and_unknown_samples_stay_blocking() {
        let a = Fixture::new(FixtureOptions::default());
        let b = Fixture::new(FixtureOptions {
            sample: Some("B".into()),
            reference_base: b'G',
            profile: EvidenceProfile {
                min_base_quality: 36,
                ..EvidenceProfile::default()
            },
            ..FixtureOptions::default()
        });
        let root = Root::new();
        let snapshot = create_snapshot(
            &root.0,
            &[a.member("A"), b.member("B")],
            None,
            CohortLimits::default(),
        )
        .unwrap();
        let mut request = query(vec![site(1, b"C")]);
        request.missing_policy = MissingPolicy::Partial;
        let incompatible = plan(&snapshot, &request);
        assert_eq!(incompatible.issue_count, 2);
        assert!(incompatible.ensure_executable().is_err());
        for code in [
            ComparisonMismatchCode::FilterProfile,
            ComparisonMismatchCode::ReferenceContent,
        ] {
            assert!(incompatible.issues.iter().any(|issue| matches!(issue,
                CohortQueryIssue::Comparison { mismatch, .. } if mismatch.code == code)));
        }
        request.member_ids = Some(vec!["A".into()]);
        plan(&snapshot, &request).ensure_executable().unwrap();
        let unknown = Fixture::new(FixtureOptions {
            sample: None,
            ..FixtureOptions::default()
        });
        let other = Root::new();
        let unnamed = create_snapshot(
            &other.0,
            &[unknown.member("unknown")],
            None,
            CohortLimits::default(),
        )
        .unwrap();
        request.member_ids = None;
        let refused = plan(&unnamed, &request);
        assert!(refused.issues.iter().any(|issue| matches!(issue,
            CohortQueryIssue::Comparison { mismatch, .. } if mismatch.code == ComparisonMismatchCode::SampleScope)));
    }

    #[test]
    fn empty_shapes_member_errors_and_bounded_diagnostics_are_explicit() {
        let a = Fixture::new(FixtureOptions::default());
        let root = Root::new();
        let snapshot =
            create_snapshot(&root.0, &[a.member("A")], None, CohortLimits::default()).unwrap();
        let empty = plan(&snapshot, &query(vec![]));
        empty.ensure_executable().unwrap();
        assert_eq!(
            (empty.output_rows, empty.members[0].covered_loci),
            (0, Some(0))
        );
        let mut request = query(vec![site(1, b"C")]);
        request.member_ids = Some(vec![]);
        let unselected = plan(&snapshot, &request);
        unselected.ensure_executable().unwrap();
        assert_eq!((unselected.requested_loci, unselected.output_rows), (1, 0));
        for ids in [vec!["missing".into()], vec!["A".into(), "A".into()]] {
            request.member_ids = Some(ids);
            assert!(plan_query(
                &snapshot,
                &request,
                &EvidenceExecution::default(),
                CohortQueryLimits::default()
            )
            .is_err());
        }
        let empty_root = Root::new();
        let no_members =
            create_snapshot(&empty_root.0, &[], None, CohortLimits::default()).unwrap();
        plan(&no_members, &query(vec![]))
            .ensure_executable()
            .unwrap();
        assert!(plan_query(
            &no_members,
            &query(vec![site(1, b"C")]),
            &EvidenceExecution::default(),
            CohortQueryLimits::default()
        )
        .is_err());
        let b = Fixture::new(FixtureOptions {
            sample: Some("B".into()),
            reference_base: b'G',
            profile: EvidenceProfile {
                min_base_quality: 36,
                ..EvidenceProfile::default()
            },
            ..FixtureOptions::default()
        });
        let bad_root = Root::new();
        let incompatible = create_snapshot(
            &bad_root.0,
            &[a.member("A"), b.member("B")],
            None,
            CohortLimits::default(),
        )
        .unwrap();
        let bounded = plan_query(
            &incompatible,
            &query(vec![site(1, b"C")]),
            &EvidenceExecution::default(),
            CohortQueryLimits {
                max_issues: 1,
                ..CohortQueryLimits::default()
            },
        )
        .unwrap();
        assert_eq!((bounded.issue_count, bounded.issues.len()), (2, 1));
        assert!(bounded.ensure_executable().is_err());
    }

    #[test]
    fn resource_and_selection_envelopes_refuse_before_reading_payloads() {
        let a = Fixture::new(FixtureOptions::default());
        let root = Root::new();
        let snapshot =
            create_snapshot(&root.0, &[a.member("A")], None, CohortLimits::default()).unwrap();
        let request = query(vec![site(1, b"C"), site(2, b"C")]);
        for limits in [
            CohortQueryLimits {
                max_input_sites: 1,
                ..CohortQueryLimits::default()
            },
            CohortQueryLimits {
                max_query_bytes: 1,
                ..CohortQueryLimits::default()
            },
            CohortQueryLimits {
                max_plan_bytes: 1,
                ..CohortQueryLimits::default()
            },
        ] {
            assert!(matches!(
                plan_query(&snapshot, &request, &EvidenceExecution::default(), limits),
                Err(CohortError::Limit(_))
            ));
        }
        assert!(plan_query(
            &snapshot,
            &request,
            &EvidenceExecution {
                memory_budget_bytes: Some(1),
                ..EvidenceExecution::default()
            },
            CohortQueryLimits::default()
        )
        .is_err());
        let wrong_ref = CohortQuery::new(EvidenceSelection::Sites(vec![SnvSite {
            reference: b'G',
            ..site(1, b"C")
        }]));
        // Planning never decodes payloads or promises REF validation. The strict
        // dataset reader will reject this REF when materialization consumes it.
        plan(&snapshot, &wrong_ref).ensure_executable().unwrap();
    }
}
