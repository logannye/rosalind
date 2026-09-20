//! Saved-only serial windows. Publication and public CLI/Python adapters are
//! separate work; callbacks can observe earlier rows before a later failure.

use super::descriptor::CohortLimits;
use super::pairs::PairPlan;
use super::query::{CohortQueryPlan, MemberQueryPlan, MissingPolicy};
use super::store::SnapshotHandle;
use super::{CohortError, Result};
use crate::dataset::DatasetQuery;
use crate::evidence::{
    EvidenceBatch, EvidenceCallback, EvidenceExecution, EvidenceRowRef, EvidenceSelection, SnvSite,
    CANONICAL_TILE_BASES,
};
use std::collections::BTreeMap;

#[derive(Debug, Clone, Copy)]
pub(crate) struct CohortRow<'a> {
    pub member_index: usize,
    pub member_id: &'a str,
    pub contig: u32,
    pub position: u32,
    pub reference: u8,
    pub alternate: u8,
    /// None means not represented in the saved selection. An observed zero is
    /// Some(row) with exact zero-valued groups, never absent evidence.
    pub evidence: Option<EvidenceRowRef<'a>>,
}

pub(crate) trait CohortConsumer {
    /// Declared maximum retained state, including window reducer/output buffers.
    /// A budgeted run refuses an unspecified bound.
    fn retained_bytes(&self) -> Option<u64>;
    fn begin_window(&mut self, _sites: &[SnvSite]) -> Result<()> {
        Ok(())
    }
    fn on_row(&mut self, row: CohortRow<'_>) -> Result<()>;
    fn end_window(&mut self) -> Result<()> {
        Ok(())
    }
    fn finish(&mut self) -> Result<()> {
        Ok(())
    }
}

/// Explicit ordered comparison. Indices name snapshot members; names or other
/// metadata never choose a direction or imply a pairing.
#[derive(Debug, Clone, Copy)]
pub(crate) struct PairWindow<'a> {
    pub pair_id: &'a str,
    pub left_member_id: &'a str,
    pub right_member_id: &'a str,
    pub left_member_index: usize,
    pub right_member_index: usize,
}

pub(crate) trait PairConsumer: CohortConsumer {
    fn begin_pair_window(&mut self, pair: PairWindow<'_>, sites: &[SnvSite]) -> Result<()>;
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct CohortRunStats {
    pub emitted_rows: u64,
    pub observed_rows: u64,
    pub unmeasured_rows: u64,
    pub execution_window_bases: u32,
    pub largest_window_loci: usize,
    pub predicted_peak_rss_bytes: u64,
    /// Only partition files actually verified by this execution, in path order.
    pub consumed_inputs: BTreeMap<String, String>,
}

/// Canonical long-form order: member, reference dictionary, position, ALT.
pub(crate) fn visit_rows(
    snapshot: &SnapshotHandle,
    plan: &CohortQueryPlan,
    execution: &EvidenceExecution,
    limits: CohortLimits,
    consumer: &mut dyn CohortConsumer,
) -> Result<CohortRunStats> {
    visit(snapshot, plan, execution, limits, consumer, false)
}

/// Candidate-summary order: genomic window then serial members. The consumer
/// retains only current-window reducer state and emits it at end_window.
pub(crate) fn visit_windows(
    snapshot: &SnapshotHandle,
    plan: &CohortQueryPlan,
    execution: &EvidenceExecution,
    limits: CohortLimits,
    consumer: &mut dyn CohortConsumer,
) -> Result<CohortRunStats> {
    visit(snapshot, plan, execution, limits, consumer, true)
}

/// Pair-table order, then canonical candidate order. Only one pair's window is
/// retained; left and right datasets open serially through the ordinary reader.
pub(crate) fn visit_pairs<C: PairConsumer + ?Sized>(
    snapshot: &SnapshotHandle,
    plan: &CohortQueryPlan,
    pairs: &PairPlan,
    execution: &EvidenceExecution,
    limits: CohortLimits,
    consumer: &mut C,
) -> Result<CohortRunStats> {
    plan.ensure_executable()?;
    if plan.snapshot_id != snapshot.id {
        return Err(CohortError::Incompatible(
            "query plan names another snapshot".into(),
        ));
    }
    // Validate all resolved references before the first callback/output.
    for pair in &pairs.pairs {
        if pair.left_member_index == pair.right_member_index
            || pair.left_member_index >= plan.members.len()
            || pair.right_member_index >= plan.members.len()
        {
            return Err(CohortError::Corrupt("invalid resolved pair member".into()));
        }
    }
    let expected = plan
        .candidate_rows
        .checked_mul(pairs.pairs.len() as u64)
        .and_then(|count| count.checked_mul(2))
        .ok_or_else(|| CohortError::Limit("pair side observation count overflow".into()))?;
    snapshot.verify_unchanged()?;
    let retained = consumer
        .retained_bytes()
        .map(|bytes| bytes.saturating_add(pairs.metadata_bytes));
    let admitted = runtime_plan(plan, execution, limits, retained)?;
    let active_limits = CohortLimits {
        memory_budget_bytes: admitted.budget,
        ..limits
    };
    let mut stats = CohortRunStats {
        emitted_rows: 0,
        observed_rows: 0,
        unmeasured_rows: 0,
        execution_window_bases: admitted.width,
        largest_window_loci: 0,
        predicted_peak_rss_bytes: admitted.predicted,
        consumed_inputs: BTreeMap::new(),
    };
    for pair in &pairs.pairs {
        let left = &plan.members[pair.left_member_index];
        let right = &plan.members[pair.right_member_index];
        for sites in Windows::new(&plan.sites, admitted.width) {
            active_limits.admit(0)?;
            consumer.begin_pair_window(
                PairWindow {
                    pair_id: &pair.id,
                    left_member_id: &left.member_id,
                    right_member_id: &right.member_id,
                    left_member_index: left.member_index,
                    right_member_index: right.member_index,
                },
                sites,
            )?;
            for member in [left, right] {
                visit_member_window(
                    snapshot,
                    plan,
                    member,
                    sites,
                    &admitted,
                    active_limits,
                    consumer,
                    &mut stats,
                )?;
            }
            consumer.end_window()?;
        }
    }
    if stats.emitted_rows != expected {
        return Err(CohortError::Corrupt(
            "paired side count differs from requested comparisons".into(),
        ));
    }
    snapshot.verify_unchanged()?;
    active_limits.admit(0)?;
    consumer.finish()?;
    snapshot.verify_unchanged()?;
    active_limits.admit(0)?;
    Ok(stats)
}

fn visit(
    snapshot: &SnapshotHandle,
    plan: &CohortQueryPlan,
    execution: &EvidenceExecution,
    limits: CohortLimits,
    consumer: &mut dyn CohortConsumer,
    window_first: bool,
) -> Result<CohortRunStats> {
    plan.ensure_executable()?;
    if plan.snapshot_id != snapshot.id {
        return Err(CohortError::Incompatible(
            "query plan names another snapshot".into(),
        ));
    }
    snapshot.verify_unchanged()?;
    let admitted = runtime_plan(plan, execution, limits, consumer.retained_bytes())?;
    let mut stats = CohortRunStats {
        emitted_rows: 0,
        observed_rows: 0,
        unmeasured_rows: 0,
        execution_window_bases: admitted.width,
        largest_window_loci: 0,
        predicted_peak_rss_bytes: admitted.predicted,
        consumed_inputs: BTreeMap::new(),
    };
    let active_limits = CohortLimits {
        memory_budget_bytes: admitted.budget,
        ..limits
    };
    if window_first {
        for sites in Windows::new(&plan.sites, admitted.width) {
            active_limits.admit(0)?;
            consumer.begin_window(sites)?;
            for member in &plan.members {
                visit_member_window(
                    snapshot,
                    plan,
                    member,
                    sites,
                    &admitted,
                    active_limits,
                    consumer,
                    &mut stats,
                )?;
            }
            consumer.end_window()?;
        }
    } else {
        for member in &plan.members {
            for sites in Windows::new(&plan.sites, admitted.width) {
                active_limits.admit(0)?;
                consumer.begin_window(sites)?;
                visit_member_window(
                    snapshot,
                    plan,
                    member,
                    sites,
                    &admitted,
                    active_limits,
                    consumer,
                    &mut stats,
                )?;
                consumer.end_window()?;
            }
        }
    }
    if stats.emitted_rows != plan.output_rows {
        return Err(CohortError::Corrupt(
            "cohort row count differs from planned candidate denominator".into(),
        ));
    }
    snapshot.verify_unchanged()?;
    active_limits.admit(0)?;
    consumer.finish()?;
    snapshot.verify_unchanged()?;
    active_limits.admit(0)?;
    Ok(stats)
}

#[derive(Debug)]
struct RuntimePlan {
    width: u32,
    budget: Option<u64>,
    collector_bytes: u64,
    retained_bytes: u64,
    predicted: u64,
}

fn runtime_plan(
    plan: &CohortQueryPlan,
    execution: &EvidenceExecution,
    limits: CohortLimits,
    retained: Option<u64>,
) -> Result<RuntimePlan> {
    let budget = [
        execution.memory_budget_bytes,
        limits.memory_budget_bytes,
        limits.dataset.memory_budget_bytes,
        plan.reservations.memory_budget_bytes,
    ]
    .into_iter()
    .flatten()
    .min();
    let consumer_bytes = match (retained, budget) {
        (Some(bytes), _) => bytes,
        (None, Some(_)) => {
            return Err(CohortError::Incompatible(
                "budgeted cohort consumer requires declared retained bytes".into(),
            ))
        }
        (None, None) => 0,
    };
    if execution.max_microtile_bases == 0 {
        return Err(CohortError::Incompatible(
            "cohort execution window must be positive".into(),
        ));
    }
    let mut width = execution.max_microtile_bases.min(CANONICAL_TILE_BASES);
    // Reader working space includes physical-source decoding even when projection
    // is smaller. Exclude analyzer_bytes: the query planner put the entire outer
    // retention reservation there, and we account for that once below.
    let reader = plan
        .members
        .iter()
        .flat_map(|member| &member.leaves)
        .filter_map(|leaf| leaf.read_plan.as_ref())
        .map(|read| {
            read.metadata_bytes
                .saturating_add(read.source_decoder_bytes)
                .saturating_add(read.projection_bytes)
        })
        .max()
        .unwrap_or(0);
    let baseline = crate::util::rss::peak_rss_bytes();
    // The planned consumer reserve can include finalization buffers beyond this
    // row consumer. Keep that promise, while honoring a larger actual consumer
    // declaration or execution reserve supplied when the plan is run.
    let retained_bytes = plan
        .reservations
        .query_bytes
        .saturating_add(plan.reservations.plan_bytes)
        .saturating_add(plan.reservations.lineage_bytes)
        .saturating_add(
            plan.reservations
                .consumer_bytes
                .max(consumer_bytes)
                .max(execution.analyzer_bytes),
        );
    let fixed = baseline
        .saturating_add(reader)
        .saturating_add(retained_bytes);
    let per_locus = plan.fields.storage_bytes_per_locus().saturating_add(256);
    loop {
        let loci = Windows::new(&plan.sites, width)
            .map(<[SnvSite]>::len)
            .max()
            .unwrap_or(0);
        let collector_bytes = per_locus
            .saturating_mul(loci as u64)
            .saturating_add(64 << 10);
        let predicted = fixed.saturating_add(collector_bytes);
        if budget.is_none_or(|budget| predicted <= budget) {
            return Ok(RuntimePlan {
                width,
                budget,
                collector_bytes,
                retained_bytes,
                predicted,
            });
        }
        if width == 1 {
            return Err(CohortError::Limit(format!(
                "cohort serial window needs {predicted} process bytes; budget {}",
                budget.unwrap()
            )));
        }
        width = (width / 2).max(1);
    }
}

#[allow(clippy::too_many_arguments)]
fn visit_member_window<C: CohortConsumer + ?Sized>(
    snapshot: &SnapshotHandle,
    plan: &CohortQueryPlan,
    member_plan: &MemberQueryPlan,
    sites: &[SnvSite],
    admitted: &RuntimePlan,
    limits: CohortLimits,
    consumer: &mut C,
    stats: &mut CohortRunStats,
) -> Result<()> {
    if sites.is_empty() {
        return Ok(());
    }
    limits.admit(admitted.collector_bytes)?;
    let member = snapshot
        .descriptor
        .members
        .get(member_plan.member_index)
        .filter(|member| member.metadata.id == member_plan.member_id)
        .ok_or_else(|| CohortError::Corrupt("planned member differs from snapshot".into()))?;
    let contig = sites[0].contig;
    let canonical = sites[0].position / CANONICAL_TILE_BASES * CANONICAL_TILE_BASES;
    let name = plan
        .contigs
        .by_id(contig)
        .ok_or_else(|| CohortError::Corrupt("planned contig missing".into()))?
        .name
        .to_string();
    let mut batch = EvidenceBatch::new(contig, name, canonical, plan.fields, Vec::new());
    batch.reserve_exact(sites.len());
    for leaf_plan in &member_plan.leaves {
        if leaf_plan.covered_loci == Some(0) {
            continue;
        }
        let leaf = member
            .leaves
            .get(leaf_plan.leaf_index)
            .filter(|leaf| leaf.object_id == leaf_plan.object_id)
            .ok_or_else(|| CohortError::Corrupt("planned leaf differs from snapshot".into()))?;
        let mut dataset = snapshot.open_leaf(leaf, limits)?;
        let split = dataset.coverage_split(&EvidenceSelection::Sites(sites.to_vec()))?;
        if split.covered_loci == 0 {
            continue;
        }
        let query = DatasetQuery {
            selection: split.covered,
            fields: plan.fields,
        };
        let mut callback = EvidenceCallback::with_fields(
            |source: &EvidenceBatch| {
                if source.contig_id != contig || source.canonical_tile_start != canonical {
                    return Err(crate::evidence::EvidenceError::InvalidInput(
                        "saved batch left its planned window ownership".into(),
                    ));
                }
                for row in source.rows() {
                    if batch.len() >= sites.len() {
                        return Err(crate::evidence::EvidenceError::InvalidInput(
                            "saved leaves emitted overlapping or extra rows".into(),
                        ));
                    }
                    batch.push_row(row)?;
                }
                Ok(())
            },
            admitted
                .collector_bytes
                .saturating_add(admitted.retained_bytes),
            plan.fields,
        );
        let reader_execution = EvidenceExecution {
            memory_budget_bytes: admitted.budget,
            ..EvidenceExecution::default()
        };
        dataset.visit_batches(&query, &mut callback, &reader_execution)?;
        dataset.verify_unchanged()?;
        limits.admit(plan.reservations.lineage_bytes)?;
        for input in dataset.verified_partition_hashes()? {
            if let Some(previous) = stats
                .consumed_inputs
                .insert(input.path, input.blake3.clone())
            {
                if previous != input.blake3 {
                    return Err(CohortError::Corrupt(
                        "consumed partition identity changed between windows".into(),
                    ));
                }
            }
        }
    }
    let mut order: Vec<_> = batch
        .loci()
        .iter()
        .enumerate()
        .map(|(index, locus)| (locus.position, index))
        .collect();
    order.sort_unstable();
    if order.windows(2).any(|pair| pair[0].0 == pair[1].0) {
        return Err(CohortError::Corrupt(
            "multiple leaves supplied one member/locus".into(),
        ));
    }
    stats.largest_window_loci = stats.largest_window_loci.max(sites.len());
    let mut observed_index = 0;
    for site in sites {
        limits.admit(0)?;
        let evidence = if order
            .get(observed_index)
            .is_some_and(|(position, _)| *position == site.position)
        {
            let row = batch
                .row(order[observed_index].1)
                .expect("collected row index");
            observed_index += 1;
            Some(row)
        } else {
            if plan.missing_policy == MissingPolicy::Strict {
                return Err(CohortError::Corrupt(
                    "strict query encountered an unmeasured planned cell".into(),
                ));
            }
            None
        };
        for &alternate in &site.alternates {
            consumer.on_row(CohortRow {
                member_index: member_plan.member_index,
                member_id: &member.metadata.id,
                contig,
                position: site.position,
                reference: site.reference,
                alternate,
                evidence,
            })?;
            stats.emitted_rows = checked_increment(stats.emitted_rows)?;
            if evidence.is_some() {
                stats.observed_rows = checked_increment(stats.observed_rows)?;
            } else {
                stats.unmeasured_rows = checked_increment(stats.unmeasured_rows)?;
            }
        }
    }
    if observed_index != order.len() {
        return Err(CohortError::Corrupt(
            "saved leaf emitted an unselected row".into(),
        ));
    }
    Ok(())
}

fn checked_increment(value: u64) -> Result<u64> {
    value
        .checked_add(1)
        .ok_or_else(|| CohortError::Limit("cohort output count overflow".into()))
}

struct Windows<'a> {
    sites: &'a [SnvSite],
    width: u32,
}
impl<'a> Windows<'a> {
    fn new(sites: &'a [SnvSite], width: u32) -> Self {
        Self { sites, width }
    }
}
impl<'a> Iterator for Windows<'a> {
    type Item = &'a [SnvSite];
    fn next(&mut self) -> Option<Self::Item> {
        let first = self.sites.first()?;
        let canonical = first.position / CANONICAL_TILE_BASES * CANONICAL_TILE_BASES;
        let start = canonical + (first.position - canonical) / self.width * self.width;
        let end = start
            .saturating_add(self.width)
            .min(canonical.saturating_add(CANONICAL_TILE_BASES));
        let count = self
            .sites
            .partition_point(|site| site.contig == first.contig && site.position < end);
        let (window, rest) = self.sites.split_at(count);
        self.sites = rest;
        Some(window)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cohort::query::{plan_query, CohortQuery, CohortQueryLimits};
    use crate::cohort::store::create_snapshot;
    use crate::cohort::summary::{CandidateObservation, CandidateReducer, CandidateSummary};
    use crate::cohort::tests::{Fixture, FixtureOptions};
    use std::fs;
    use std::path::PathBuf;
    use std::sync::atomic::{AtomicU64, Ordering};

    static ID: AtomicU64 = AtomicU64::new(0);
    struct Setup {
        root: PathBuf,
        a: Option<Fixture>,
        b: Option<Fixture>,
        snapshot: SnapshotHandle,
    }
    impl Setup {
        fn new() -> Self {
            let a = Fixture::new(FixtureOptions::default());
            let b = Fixture::new(FixtureOptions {
                sample: Some("B".into()),
                ..FixtureOptions::default()
            });
            let mut second = b.member("B");
            second.manifests = vec![b.extra_leaf(8, 9, crate::evidence::EvidenceFields::ALL)];
            let root = std::env::temp_dir().join(format!(
                "rosalind-cohort-runtime-{}-{}",
                std::process::id(),
                ID.fetch_add(1, Ordering::Relaxed)
            ));
            let snapshot = create_snapshot(
                &root,
                &[a.member("A"), second],
                None,
                CohortLimits::default(),
            )
            .unwrap();
            Self {
                root,
                a: Some(a),
                b: Some(b),
                snapshot,
            }
        }
        fn plan(&self, missing: MissingPolicy) -> CohortQueryPlan {
            let mut query = CohortQuery::new(EvidenceSelection::Sites(vec![
                SnvSite {
                    contig: 0,
                    position: 1,
                    reference: b'A',
                    alternates: vec![b'C', b'T'],
                },
                SnvSite {
                    contig: 0,
                    position: 2,
                    reference: b'A',
                    alternates: vec![b'C'],
                },
                SnvSite {
                    contig: 0,
                    position: 8,
                    reference: b'A',
                    alternates: vec![b'C'],
                },
            ]));
            query.missing_policy = missing;
            plan_query(
                &self.snapshot,
                &query,
                &EvidenceExecution::default(),
                CohortQueryLimits::default(),
            )
            .unwrap()
        }
    }
    impl Drop for Setup {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.root);
        }
    }

    type OutputRow = (String, u32, u8, Option<u64>, Option<u64>);
    #[derive(Default)]
    struct Collect {
        rows: Vec<OutputRow>,
        begins: usize,
        finished: bool,
    }
    impl CohortConsumer for Collect {
        fn retained_bytes(&self) -> Option<u64> {
            Some(64 << 10)
        }
        fn begin_window(&mut self, _: &[SnvSite]) -> Result<()> {
            self.begins += 1;
            Ok(())
        }
        fn on_row(&mut self, row: CohortRow<'_>) -> Result<()> {
            let depth = row
                .evidence
                .and_then(|e| e.depths.map(|d| d.callable_depth));
            let index = b"ACGT"
                .iter()
                .position(|&base| base == row.alternate)
                .unwrap();
            let alt = row
                .evidence
                .and_then(|e| e.alleles.map(|a| a.allele_counts[index]));
            self.rows.push((
                row.member_id.into(),
                row.position,
                row.alternate,
                depth,
                alt,
            ));
            Ok(())
        }
        fn finish(&mut self) -> Result<()> {
            self.finished = true;
            Ok(())
        }
    }

    #[test]
    fn saved_only_partial_rows_distinguish_absent_zero_and_new_alt() {
        let mut setup = Setup::new();
        let plan = setup.plan(MissingPolicy::Partial);
        drop(setup.a.take());
        drop(setup.b.take());
        let mut collect = Collect::default();
        let stats = visit_rows(
            &setup.snapshot,
            &plan,
            &EvidenceExecution::default(),
            CohortLimits::default(),
            &mut collect,
        )
        .unwrap();
        assert_eq!(
            collect.rows,
            vec![
                ("A".into(), 1, b'C', Some(1), Some(1)),
                ("A".into(), 1, b'T', Some(1), Some(0)),
                ("A".into(), 2, b'C', Some(0), Some(0)),
                ("A".into(), 8, b'C', None, None),
                ("B".into(), 1, b'C', None, None),
                ("B".into(), 1, b'T', None, None),
                ("B".into(), 2, b'C', None, None),
                ("B".into(), 8, b'C', Some(1), Some(1)),
            ]
        );
        assert_eq!(
            (
                stats.emitted_rows,
                stats.observed_rows,
                stats.unmeasured_rows
            ),
            (8, 4, 4)
        );
        assert!(collect.finished);
    }

    #[test]
    fn pair_traversal_preserves_explicit_order_missingness_and_saved_only_bytes_across_budgets() {
        use crate::cohort::encoding::CohortOutputFormat;
        use crate::cohort::pair_encoding::PairEncoder;
        use crate::cohort::pairs::{plan_pairs, PairLimits, PairSpec};
        let mut setup = Setup::new();
        let plan = setup.plan(MissingPolicy::Partial);
        let pairs = plan_pairs(
            &[
                PairSpec {
                    id: "first-B-to-A".into(),
                    left: "B".into(),
                    right: "A".into(),
                },
                PairSpec {
                    id: "second-A-to-B".into(),
                    left: "A".into(),
                    right: "B".into(),
                },
            ],
            &plan
                .members
                .iter()
                .map(|member| member.member_id.clone())
                .collect::<Vec<_>>(),
            PairLimits::default(),
        )
        .unwrap();
        // Neither native sources nor original leaf directories remain available.
        drop(setup.a.take());
        drop(setup.b.take());
        let baseline =
            crate::util::rss::peak_rss_bytes().max(plan.reservations.predicted_peak_rss_bytes);
        for format in [CohortOutputFormat::Arrow, CohortOutputFormat::Tsv] {
            let mut expected = None;
            for (width, extra) in [(1, 192u64), (4, 256), (CANONICAL_TILE_BASES, 320)] {
                let execution = EvidenceExecution {
                    max_microtile_bases: width,
                    memory_budget_bytes: Some(baseline + (extra << 20)),
                    ..EvidenceExecution::default()
                };
                let mut encoder =
                    PairEncoder::new(Vec::new(), format, &plan.contigs, 10, 4).unwrap();
                let stats = visit_pairs(
                    &setup.snapshot,
                    &plan,
                    &pairs,
                    &execution,
                    CohortLimits::default(),
                    &mut encoder,
                )
                .unwrap();
                assert_eq!(
                    (
                        stats.emitted_rows,
                        stats.observed_rows,
                        stats.unmeasured_rows
                    ),
                    (16, 8, 8)
                );
                assert!(stats.largest_window_loci <= width as usize);
                // Two leaves, each Arrow payload plus its receipt; shared use in
                // the reversed pair must not duplicate consumed identities.
                assert_eq!(stats.consumed_inputs.len(), 4);
                let output = encoder.into_inner().unwrap();
                if let Some(previous) = &expected {
                    assert_eq!(previous, &output);
                } else {
                    expected = Some(output);
                }
            }
            if format == CohortOutputFormat::Tsv {
                let text = String::from_utf8(expected.unwrap()).unwrap();
                let rows: Vec<_> = text
                    .lines()
                    .map(|line| line.split('\t').collect::<Vec<_>>())
                    .collect();
                let column = |name| rows[0].iter().position(|value| *value == name).unwrap();
                assert_eq!(rows.len(), 9);
                assert_eq!(rows[1][column("pair_id")], "first-B-to-A");
                assert_eq!(rows[5][column("pair_id")], "second-A-to-B");
                assert_eq!(rows[1][column("left_status")], "unmeasured");
                assert_eq!(rows[1][column("right_alt_count")], "1");
                assert_eq!(rows[2][column("right_alt_count")], "0");
                assert_eq!(rows[3][column("right_callable_depth")], "0");
                assert_eq!(rows[3][column("difference_numerator")], ".");
                assert_eq!(rows[4][column("left_alt_count")], "1");
                assert_eq!(rows[4][column("right_status")], "unmeasured");
            }
        }
        let mut encoder =
            PairEncoder::new(Vec::new(), CohortOutputFormat::Tsv, &plan.contigs, 10, 4).unwrap();
        assert!(visit_pairs(
            &setup.snapshot,
            &plan,
            &pairs,
            &EvidenceExecution {
                memory_budget_bytes: Some(1),
                ..EvidenceExecution::default()
            },
            CohortLimits::default(),
            &mut encoder
        )
        .is_err());
        assert!(
            String::from_utf8(encoder.into_inner().unwrap())
                .unwrap()
                .lines()
                .count()
                == 1
        );
    }

    #[test]
    fn strict_coverage_and_changed_inputs_refuse_before_callbacks() {
        let setup = Setup::new();
        let strict = setup.plan(MissingPolicy::Strict);
        let mut collect = Collect::default();
        assert!(visit_rows(
            &setup.snapshot,
            &strict,
            &EvidenceExecution::default(),
            CohortLimits::default(),
            &mut collect
        )
        .is_err());
        assert_eq!(collect.begins, 0);
        let partial = setup.plan(MissingPolicy::Partial);
        let object = &setup.snapshot.descriptor.members[0].leaves[0].object_id;
        let relative = setup.a.as_ref().unwrap().descriptor().partitions[0]
            .arrow
            .path
            .clone();
        fs::write(
            setup.root.join("objects").join(object).join(relative),
            b"changed after planning",
        )
        .unwrap();
        assert!(visit_rows(
            &setup.snapshot,
            &partial,
            &EvidenceExecution::default(),
            CohortLimits::default(),
            &mut collect
        )
        .is_err());
        assert_eq!(collect.begins, 0);
        assert!(!collect.finished);
    }

    #[test]
    fn three_admitted_budgets_and_window_widths_preserve_output_bytes() {
        let setup = Setup::new();
        let plan = setup.plan(MissingPolicy::Partial);
        let mut reference = None;
        let baseline =
            crate::util::rss::peak_rss_bytes().max(plan.reservations.predicted_peak_rss_bytes);
        for (width, extra) in [(1, 128u64), (4, 192), (CANONICAL_TILE_BASES, 256)] {
            let execution = EvidenceExecution {
                max_microtile_bases: width,
                memory_budget_bytes: Some(baseline + (extra << 20)),
                ..EvidenceExecution::default()
            };
            let mut collect = Collect::default();
            let stats = visit_rows(
                &setup.snapshot,
                &plan,
                &execution,
                CohortLimits::default(),
                &mut collect,
            )
            .unwrap();
            assert!(stats.largest_window_loci <= width as usize);
            let bytes = serde_json::to_vec(&collect.rows).unwrap();
            if let Some(expected) = &reference {
                assert_eq!(expected, &bytes);
            } else {
                reference = Some(bytes);
            }
        }
        let execution = EvidenceExecution {
            memory_budget_bytes: Some(1),
            ..EvidenceExecution::default()
        };
        let mut collect = Collect::default();
        assert!(visit_rows(
            &setup.snapshot,
            &plan,
            &execution,
            CohortLimits::default(),
            &mut collect
        )
        .is_err());
        assert_eq!(collect.begins, 0);
    }

    #[test]
    fn serial_admission_counts_outer_retention_once_and_honors_larger_runtime_consumers() {
        let setup = Setup::new();
        let selection = setup.plan(MissingPolicy::Partial).sites;
        let mut query = CohortQuery::new(EvidenceSelection::Sites(selection));
        query.missing_policy = MissingPolicy::Partial;
        // A real reader plan includes this future encoder/finalization promise
        // in analyzer_bytes. It must not become a second reader allocation.
        const RESERVED: u64 = 64 << 20;
        query.requirements.retained_bytes = Some(RESERVED);
        let plan = plan_query(
            &setup.snapshot,
            &query,
            &EvidenceExecution::default(),
            CohortQueryLimits::default(),
        )
        .unwrap();
        let readers: Vec<_> = plan
            .members
            .iter()
            .flat_map(|member| &member.leaves)
            .filter_map(|leaf| leaf.read_plan.as_ref())
            .collect();
        assert!(readers.iter().all(|read| read.analyzer_bytes >= RESERVED));
        let serial_workspace = readers
            .iter()
            .map(|read| read.metadata_bytes + read.source_decoder_bytes + read.projection_bytes)
            .max()
            .unwrap();
        let execution = EvidenceExecution {
            max_microtile_bases: 1,
            memory_budget_bytes: Some(
                crate::util::rss::peak_rss_bytes()
                    + serial_workspace
                    + plan.reservations.retained_bytes()
                    + (8 << 20),
            ),
            ..EvidenceExecution::default()
        };
        let mut collect = Collect::default();
        let stats = visit_rows(
            &setup.snapshot,
            &plan,
            &execution,
            CohortLimits::default(),
            &mut collect,
        )
        .unwrap();
        assert_eq!(stats.emitted_rows, plan.output_rows);
        assert!(collect.finished);
        assert!(runtime_plan(
            &plan,
            &execution,
            CohortLimits::default(),
            Some(RESERVED * 2),
        )
        .is_err());
        let larger_execution = EvidenceExecution {
            analyzer_bytes: RESERVED * 2,
            ..execution
        };
        assert!(runtime_plan(
            &plan,
            &larger_execution,
            CohortLimits::default(),
            collect.retained_bytes(),
        )
        .is_err());
    }

    struct Summaries {
        active: Vec<(u32, u8, CandidateReducer)>,
        output: Vec<(u32, u8, CandidateSummary)>,
    }
    impl CohortConsumer for Summaries {
        fn retained_bytes(&self) -> Option<u64> {
            Some(64 << 10)
        }
        fn begin_window(&mut self, sites: &[SnvSite]) -> Result<()> {
            assert!(self.active.is_empty());
            for site in sites {
                for &alt in &site.alternates {
                    self.active
                        .push((site.position, alt, CandidateReducer::new(2, 1)?));
                }
            }
            Ok(())
        }
        fn on_row(&mut self, row: CohortRow<'_>) -> Result<()> {
            let observation = match row.evidence {
                Some(evidence) => CandidateObservation::from_row(evidence, row.alternate, 1)?,
                None => CandidateObservation::unmeasured(1)?,
            };
            self.active
                .iter_mut()
                .find(|(position, alt, _)| *position == row.position && *alt == row.alternate)
                .unwrap()
                .2
                .push(&observation)
        }
        fn end_window(&mut self) -> Result<()> {
            for (position, alt, reducer) in self.active.drain(..) {
                self.output.push((position, alt, reducer.finish()?));
            }
            Ok(())
        }
    }

    #[test]
    fn window_first_summary_matches_independent_sample_counts() {
        let setup = Setup::new();
        let plan = setup.plan(MissingPolicy::Partial);
        let mut reference = None;
        for width in [1, 2, CANONICAL_TILE_BASES] {
            let execution = EvidenceExecution {
                max_microtile_bases: width,
                ..EvidenceExecution::default()
            };
            let mut consumer = Summaries {
                active: Vec::new(),
                output: Vec::new(),
            };
            visit_windows(
                &setup.snapshot,
                &plan,
                &execution,
                CohortLimits::default(),
                &mut consumer,
            )
            .unwrap();
            let counts: Vec<_> = consumer
                .output
                .iter()
                .map(|(position, alt, summary)| {
                    (
                        *position,
                        *alt,
                        summary.n_requested,
                        summary.n_observed,
                        summary.n_depth_eligible,
                        summary.n_alt_supported,
                        summary.callable_total,
                        summary.alt_total,
                    )
                })
                .collect();
            assert_eq!(
                counts,
                vec![
                    (1, b'C', 2, 1, 1, 1, 1, 1),
                    (1, b'T', 2, 1, 1, 0, 1, 0),
                    (2, b'C', 2, 1, 0, 0, 0, 0),
                    (8, b'C', 2, 1, 1, 1, 1, 1)
                ]
            );
            let bytes = serde_json::to_vec(&consumer.output).unwrap();
            if let Some(expected) = &reference {
                assert_eq!(expected, &bytes);
            } else {
                reference = Some(bytes);
            }
        }
    }

    #[test]
    fn consumer_failure_never_calls_finish() {
        struct Failing {
            finished: bool,
        }
        impl CohortConsumer for Failing {
            fn retained_bytes(&self) -> Option<u64> {
                Some(0)
            }
            fn on_row(&mut self, _: CohortRow<'_>) -> Result<()> {
                Err(CohortError::Incompatible("injected consumer error".into()))
            }
            fn finish(&mut self) -> Result<()> {
                self.finished = true;
                Ok(())
            }
        }
        let setup = Setup::new();
        let plan = setup.plan(MissingPolicy::Partial);
        let mut consumer = Failing { finished: false };
        assert!(visit_rows(
            &setup.snapshot,
            &plan,
            &EvidenceExecution::default(),
            CohortLimits::default(),
            &mut consumer
        )
        .is_err());
        assert!(!consumer.finished);
    }

    #[test]
    fn sparse_queries_and_multiple_leaves_cross_canonical_boundaries() {
        let fixture = Fixture::new(FixtureOptions {
            reference_length: CANONICAL_TILE_BASES + 8,
            ..FixtureOptions::default()
        });
        let mut member = fixture.member("A");
        member.manifests = vec![
            fixture.extra_leaf(
                CANONICAL_TILE_BASES - 2,
                CANONICAL_TILE_BASES,
                crate::evidence::EvidenceFields::ALL,
            ),
            fixture.extra_leaf(
                CANONICAL_TILE_BASES,
                CANONICAL_TILE_BASES + 3,
                crate::evidence::EvidenceFields::DEPTHS
                    .union(crate::evidence::EvidenceFields::ALLELES),
            ),
        ];
        let root = fixture.root.join("cohort");
        let snapshot = create_snapshot(&root, &[member], None, CohortLimits::default()).unwrap();
        let sites: Vec<_> = [
            CANONICAL_TILE_BASES - 2,
            CANONICAL_TILE_BASES - 1,
            CANONICAL_TILE_BASES,
            CANONICAL_TILE_BASES + 2,
        ]
        .into_iter()
        .map(|position| SnvSite {
            contig: 0,
            position,
            reference: b'A',
            alternates: vec![b'C'],
        })
        .collect();
        let plan = plan_query(
            &snapshot,
            &CohortQuery::new(EvidenceSelection::Sites(sites)),
            &EvidenceExecution::default(),
            CohortQueryLimits::default(),
        )
        .unwrap();
        let mut reference = None;
        for width in [1, 3, 7, CANONICAL_TILE_BASES] {
            let mut output = Collect::default();
            let execution = EvidenceExecution {
                max_microtile_bases: width,
                ..EvidenceExecution::default()
            };
            visit_rows(
                &snapshot,
                &plan,
                &execution,
                CohortLimits::default(),
                &mut output,
            )
            .unwrap();
            assert_eq!(output.rows.len(), 4);
            assert!(output
                .rows
                .iter()
                .all(|row| row.3 == Some(0) && row.4 == Some(0)));
            if let Some(expected) = &reference {
                assert_eq!(expected, &output.rows);
            } else {
                reference = Some(output.rows);
            }
        }
    }

    #[test]
    fn empty_selection_finishes_without_rows_or_windows() {
        let setup = Setup::new();
        let plan = plan_query(
            &setup.snapshot,
            &CohortQuery::new(EvidenceSelection::Sites(vec![])),
            &EvidenceExecution::default(),
            CohortQueryLimits::default(),
        )
        .unwrap();
        let mut output = Collect::default();
        let stats = visit_rows(
            &setup.snapshot,
            &plan,
            &EvidenceExecution::default(),
            CohortLimits::default(),
            &mut output,
        )
        .unwrap();
        assert_eq!(stats.emitted_rows, 0);
        assert_eq!(output.begins, 0);
        assert!(output.finished);
    }

    #[test]
    fn cancelled_outer_scope_stops_rows_without_starting_a_nested_runner() {
        const CHILD: &str = "ROSALIND_COHORT_CANCEL_CHILD";
        if std::env::var_os(CHILD).is_none() {
            let status = std::process::Command::new(std::env::current_exe().unwrap())
                .arg("--exact").arg("cohort::runtime::tests::cancelled_outer_scope_stops_rows_without_starting_a_nested_runner")
                .env(CHILD, "1").status().unwrap();
            assert!(status.success());
            return;
        }
        let setup = Setup::new();
        let plan = setup.plan(MissingPolicy::Partial);
        let token = crate::core::cancellation::CancellationToken::new();
        let _scope = crate::core::cancellation::CancellationScope::start(token.clone()).unwrap();
        token.cancel();
        let mut consumer = Collect::default();
        assert!(visit_rows(
            &setup.snapshot,
            &plan,
            &EvidenceExecution::default(),
            CohortLimits::default(),
            &mut consumer
        )
        .is_err());
        assert!(consumer.rows.is_empty());
        assert!(!consumer.finished);
    }
}
