//! Explicit missing-locus extraction. Stored rows are never decoded here: the
//! publication boundary verifies existing objects by their declared byte hashes.

use super::comparison::require_named_sample;
use super::query::{plan_query, CohortQuery, CohortQueryLimits, MissingPolicy};
use super::store::{publish_extension, SnapshotHandle};
use super::{CohortError, Result};
use crate::dataset::{
    publish_evidence_dataset, run_dataset_with_snapshot, DatasetDescriptor, DatasetOptions,
    DescriptorLimits, InputSnapshot, VerifiedInputSession,
};
use crate::evidence::{
    EvidenceBatch, EvidenceCallback, EvidenceEngine, EvidenceExecution, EvidenceRequest,
    EvidenceSampleSelection, EvidenceSelection,
};
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;

/// Every recorded scientific source role must be supplied explicitly. Paths may
/// relocate; original byte identities and resolved sample scope must agree.
#[derive(Debug, Clone)]
pub(crate) struct MemberSourceMapping {
    pub member_id: String,
    pub roles: BTreeMap<String, PathBuf>,
}

#[derive(Debug, Clone)]
pub(crate) struct MemberExtensionStats {
    pub member_id: String,
    pub retained_loci: u64,
    pub computed_loci: u64,
    /// Indexed fetch visits, not unique reads; indexes can visit wider blocks.
    pub native_record_visits: u64,
    /// CRAM's mandatory full-source validation is reported separately.
    pub full_cram_validation_records: u64,
    /// Unique canonical local files hashed by this member's source session.
    pub source_hashed_bytes: u64,
    pub source_verification_ms: u64,
    pub extraction_ms: u64,
}

#[derive(Debug)]
pub(crate) struct ExtensionOutcome {
    pub snapshot: SnapshotHandle,
    pub changed: bool,
    pub members: Vec<MemberExtensionStats>,
    /// Existing unique portable inventories read by byte verification. Excludes
    /// metadata parsing, new-object validation, copying and raw source hashing.
    pub verified_existing_bytes: u64,
    pub publication_ms: u64,
}

/// Prepare every missing-only leaf privately, then atomically publish one child.
/// The caller supplies the outer resource/cancellation lifecycle; this routine
/// never recursively enters the process-wide managed runner.
pub(crate) fn extend_snapshot(
    parent: &SnapshotHandle,
    query: &CohortQuery,
    sources: &[MemberSourceMapping],
    work_parent: &Path,
    execution: &EvidenceExecution,
    limits: CohortQueryLimits,
) -> Result<ExtensionOutcome> {
    extend_snapshot_with_inputs(parent, query, sources, work_parent, execution, limits, None)
}

/// The CLI can additionally retain its candidate-file/source-table guards up to
/// the same final publication boundary as the original scientific sources.
pub(crate) fn extend_snapshot_with_inputs(
    parent: &SnapshotHandle,
    query: &CohortQuery,
    sources: &[MemberSourceMapping],
    work_parent: &Path,
    execution: &EvidenceExecution,
    mut limits: CohortQueryLimits,
    inputs: Option<&InputSnapshot>,
) -> Result<ExtensionOutcome> {
    if let Some(inputs) = inputs {
        inputs.verify()?;
    }
    limits.cohort.memory_budget_bytes = [
        limits.cohort.effective_budget(),
        execution.memory_budget_bytes,
    ]
    .into_iter()
    .flatten()
    .min();
    limits.cohort.validate()?;
    let source_reserve = mapping_reservation(sources, limits)?;
    // Bound the input before cloning to change only the coverage policy.
    bound_query_copy(query, limits)?;
    limits.cohort.admit(source_reserve)?;
    let mut planning_query = query.clone();
    planning_query.missing_policy = MissingPolicy::Partial;
    let plan = plan_query(parent, &planning_query, execution, limits)?;
    drop(planning_query);
    plan.ensure_executable()?;

    let mut mappings = BTreeMap::new();
    for mapping in sources {
        if mappings
            .insert(mapping.member_id.as_str(), mapping)
            .is_some()
        {
            return Err(incompatible("duplicate source mapping member"));
        }
    }
    let affected: BTreeSet<_> = plan
        .members
        .iter()
        .filter(|member| member.missing_loci.is_some_and(|count| count != 0))
        .map(|member| member.member_id.as_str())
        .collect();
    if affected != mappings.keys().copied().collect() {
        return Err(incompatible(
            "source mappings must name exactly the selected members with missing loci; no missing, extra or unselected member mappings",
        ));
    }
    let retained = plan
        .reservations
        .retained_bytes()
        .checked_add(source_reserve)
        // One open original descriptor, source hashing buffer, and staging data.
        .and_then(|n| {
            (limits.cohort.dataset.max_descriptor_bytes as u64)
                .checked_mul(4)
                .and_then(|metadata| n.checked_add(metadata))
        })
        .and_then(|n| n.checked_add(256 << 10))
        .ok_or_else(|| limit("extension memory reservation overflow"))?;
    limits.cohort.admit(retained)?;
    let mut sessions = Vec::with_capacity(affected.len());
    let mut additions = Vec::with_capacity(affected.len());
    let mut stats = Vec::with_capacity(plan.members.len());
    // A no-op does not inspect the staging path or open any raw sources.
    let stage = (!affected.is_empty())
        .then(|| StagingDirectory::create(work_parent, &parent.root))
        .transpose()?;
    for (planned_index, member) in plan.members.iter().enumerate() {
        limits.cohort.admit(retained)?;
        let mut member_stats = MemberExtensionStats {
            member_id: member.member_id.clone(),
            retained_loci: member.covered_loci.unwrap_or(0),
            computed_loci: 0,
            native_record_visits: 0,
            full_cram_validation_records: 0,
            source_hashed_bytes: 0,
            source_verification_ms: 0,
            extraction_ms: 0,
        };
        if member.missing_loci == Some(0) {
            stats.push(member_stats);
            continue;
        }
        let coverage = plan.member_coverage(parent, planned_index, limits)?;
        let original_leaf = parent.descriptor.members[member.member_index]
            .leaves
            .first()
            .ok_or_else(|| incompatible("cannot extend a member without source evidence"))?;
        let original = parent.open_leaf(original_leaf, limits.cohort)?;
        let descriptor = original.descriptor();
        let mapping = mappings[member.member_id.as_str()];
        validate_roles(descriptor, mapping)?;
        let started = Instant::now();
        let session = VerifiedInputSession::open(
            mapping
                .roles
                .iter()
                .map(|(role, path)| (role.clone(), path.clone())),
        )?;
        validate_identities(descriptor, &session)?;
        require_explicit_fasta_indexes(mapping)?;
        member_stats.source_verification_ms = elapsed_ms(started);
        let mut hashed_paths = BTreeSet::new();
        for identity in session.identities() {
            if hashed_paths.insert(&identity.path) {
                member_stats.source_hashed_bytes = member_stats
                    .source_hashed_bytes
                    .checked_add(identity.bytes)
                    .ok_or_else(|| limit("source byte count overflow"))?;
            }
        }
        let mut native_execution = execution.clone();
        native_execution.memory_budget_bytes = limits.cohort.effective_budget();
        native_execution.analyzer_bytes = retained;
        let request = source_request(
            descriptor,
            mapping,
            coverage.missing,
            query,
            native_execution,
        )?;
        let started = Instant::now();
        let mut engine = EvidenceEngine::open(request)?;
        if session.compatibility_key(&engine)? != descriptor.compatibility_blake3 {
            return Err(incompatible(
                "extension engine does not reproduce the original source, reference, filter and named-sample identity",
            ));
        }
        member_stats.full_cram_validation_records = engine
            .plan()
            .cram_envelope
            .as_ref()
            .map_or(0, |envelope| envelope.validated_records);
        session.verify()?;
        original.verify_unchanged()?;
        let namespace = session.dataset_namespace(&engine)?;
        let outcome = run_dataset_with_snapshot(
            &mut engine,
            &namespace,
            &DatasetOptions {
                cache_dir: stage
                    .as_ref()
                    .expect("affected member has staging")
                    .0
                    .join(format!("member-{}", member.member_index)),
                resume: false,
                workers: 1,
            },
            &mut EvidenceCallback::with_fields(|_: &EvidenceBatch| Ok(()), retained, query.fields),
            session.snapshot(),
        )?;
        if outcome.stats.emitted_loci != coverage.missing_loci {
            return Err(CohortError::Corrupt(
                "missing-only extraction emitted an unexpected number of loci".into(),
            ));
        }
        let manifest = publish_evidence_dataset(
            &engine,
            &outcome,
            &session,
            DescriptorLimits {
                max_bytes: limits.cohort.dataset.max_descriptor_bytes,
            },
        )?;
        session.verify()?;
        original.verify_unchanged()?;
        member_stats.computed_loci = outcome.stats.emitted_loci;
        member_stats.native_record_visits = outcome.stats.record_visits;
        member_stats.extraction_ms = elapsed_ms(started);
        additions.push((member.member_index, manifest));
        sessions.push(session);
        stats.push(member_stats);
    }
    let started = Instant::now();
    let publication = publish_extension(parent, &additions, limits.cohort, || {
        if let Some(inputs) = inputs {
            inputs.verify()?;
        }
        for session in &sessions {
            session.verify()?;
        }
        parent.verify_unchanged()
    })?;
    Ok(ExtensionOutcome {
        changed: publication.snapshot.id != parent.id,
        snapshot: publication.snapshot,
        members: stats,
        verified_existing_bytes: publication.verified_existing_bytes,
        publication_ms: elapsed_ms(started),
    })
}

fn source_request(
    descriptor: &DatasetDescriptor,
    mapping: &MemberSourceMapping,
    missing: EvidenceSelection,
    query: &CohortQuery,
    execution: EvidenceExecution,
) -> Result<EvidenceRequest> {
    let sample = require_named_sample(descriptor)
        .map_err(|_| incompatible("extension requires one named analysis sample"))?;
    let mut request = EvidenceRequest::coverage(
        mapping
            .roles
            .get("alignments")
            .ok_or_else(|| incompatible("missing explicit alignment source"))?,
    );
    request.alignment_index = mapping.roles.get("alignment-index").cloned();
    if request.alignment_index.is_none() {
        return Err(incompatible(
            "extension requires an explicit original alignment index",
        ));
    }
    request.reference = mapping.roles.get("reference").cloned();
    request.reference_fai = mapping.roles.get("reference-fai").cloned();
    request.cram_reference = mapping.roles.get("cram-reference").cloned();
    request.cram_reference_fai = mapping.roles.get("cram-reference-fai").cloned();
    request.sample_selection = EvidenceSampleSelection::Named(sample.into());
    request.selection = missing;
    request.fields = query.fields;
    request.profile = descriptor.profile.to_profile();
    request.execution = execution;
    Ok(request)
}

fn scientific_role(role: &str) -> bool {
    !matches!(role, "sites" | "regions")
}

fn validate_roles(descriptor: &DatasetDescriptor, mapping: &MemberSourceMapping) -> Result<()> {
    let expected: BTreeSet<_> = descriptor
        .sources
        .iter()
        .filter(|source| scientific_role(&source.role))
        .map(|source| source.role.as_str())
        .collect();
    if expected != mapping.roles.keys().map(String::as_str).collect() {
        return Err(incompatible("explicit source roles must match every original scientific source role exactly (excluding sites/regions)"));
    }
    Ok(())
}

fn validate_identities(
    descriptor: &DatasetDescriptor,
    session: &VerifiedInputSession,
) -> Result<()> {
    let expected: Vec<_> = descriptor
        .sources
        .iter()
        .filter(|source| scientific_role(&source.role))
        .collect();
    if expected.len() != session.identities().len()
        || expected
            .iter()
            .zip(session.identities())
            .any(|(old, current)| {
                old.role != current.role
                    || old.blake3 != current.blake3
                    || old.bytes != current.bytes
            })
    {
        return Err(incompatible(
            "explicit source content does not match the original role, byte hash and size",
        ));
    }
    Ok(())
}

fn require_explicit_fasta_indexes(mapping: &MemberSourceMapping) -> Result<()> {
    for (reference, index) in [
        ("reference", "reference-fai"),
        ("cram-reference", "cram-reference-fai"),
    ] {
        if let Some(path) = mapping.roles.get(reference) {
            let mut prefix = [0u8; 1];
            let count = fs::File::open(path)?.read(&mut prefix)?;
            if count == 1 && prefix[0] == b'>' && !mapping.roles.contains_key(index) {
                return Err(incompatible("FASTA extension requires its original explicit FAI source role; index discovery is disabled"));
            }
        }
    }
    Ok(())
}

fn bound_query_copy(query: &CohortQuery, limits: CohortQueryLimits) -> Result<()> {
    let EvidenceSelection::Sites(sites) = &query.selection else {
        return Err(incompatible(
            "cohort extension requires candidate SNV sites",
        ));
    };
    if sites.len() > limits.max_input_sites {
        return Err(limit("extension query exceeds input site envelope"));
    }
    let mut bytes = 65536u64;
    for site in sites {
        bytes = bytes
            .saturating_add(2048)
            .saturating_add(site.alternates.len() as u64 * 16);
    }
    if let Some(ids) = &query.member_ids {
        if ids.len() > limits.cohort.max_members {
            return Err(limit("extension member selection exceeds envelope"));
        }
        for id in ids {
            bytes = bytes.saturating_add(id.len() as u64 * 8 + 256);
        }
    }
    if bytes > limits.max_query_bytes {
        return Err(limit("extension query copy exceeds byte envelope"));
    }
    limits.cohort.admit(bytes)
}

fn mapping_reservation(sources: &[MemberSourceMapping], limits: CohortQueryLimits) -> Result<u64> {
    if sources.len() > limits.cohort.max_members {
        return Err(limit("source mappings exceed member envelope"));
    }
    let mut bytes = 65536u64;
    for source in sources {
        if source.member_id.len() > 256 || source.roles.len() > 256 {
            return Err(limit("source mapping identity/role envelope exceeded"));
        }
        bytes = bytes.saturating_add(8192 + source.member_id.len() as u64 * 8);
        for (role, path) in &source.roles {
            // Canonical paths, guards, identities and publication bookkeeping.
            bytes = bytes
                .saturating_add(2048)
                .saturating_add(role.len() as u64 * 8)
                .saturating_add(path.as_os_str().len() as u64 * 16);
        }
    }
    if bytes > limits.max_plan_bytes {
        return Err(limit(
            "source mapping/session metadata exceeds byte envelope",
        ));
    }
    Ok(bytes)
}

static STAGING_ID: AtomicU64 = AtomicU64::new(0);
struct StagingDirectory(PathBuf);
impl StagingDirectory {
    fn create(work_parent: &Path, cohort_root: &Path) -> Result<Self> {
        let work_parent = fs::canonicalize(work_parent)?;
        let cohort_root = fs::canonicalize(cohort_root)?;
        if work_parent.starts_with(&cohort_root) {
            return Err(incompatible(
                "extension staging must be outside the immutable cohort",
            ));
        }
        for _ in 0..128 {
            let path = work_parent.join(format!(
                ".rosalind-extension-{}-{}",
                std::process::id(),
                STAGING_ID.fetch_add(1, Ordering::Relaxed)
            ));
            let mut builder = fs::DirBuilder::new();
            #[cfg(unix)]
            {
                use std::os::unix::fs::DirBuilderExt;
                builder.mode(0o700);
            }
            match builder.create(&path) {
                Ok(()) => return Ok(Self(path)),
                Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
                Err(error) => return Err(error.into()),
            }
        }
        Err(incompatible("could not create private extension staging"))
    }
}
impl Drop for StagingDirectory {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}
fn elapsed_ms(started: Instant) -> u64 {
    started.elapsed().as_millis().try_into().unwrap_or(u64::MAX)
}
fn incompatible(message: &str) -> CohortError {
    CohortError::Incompatible(message.into())
}
fn limit(message: &str) -> CohortError {
    CohortError::Limit(message.into())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cohort::descriptor::CohortLimits;
    use crate::cohort::runtime::{visit_rows, CohortConsumer, CohortRow};
    use crate::cohort::store::{create_snapshot, verify_snapshot};
    use crate::cohort::tests::{Fixture, FixtureOptions};
    use crate::evidence::{EvidenceFields, SnvSite};

    static TEST_ID: AtomicU64 = AtomicU64::new(0);
    struct Work(PathBuf);
    impl Work {
        fn new() -> Self {
            let path = std::env::temp_dir().join(format!(
                "rosalind-extension-test-{}-{}",
                std::process::id(),
                TEST_ID.fetch_add(1, Ordering::Relaxed)
            ));
            fs::create_dir(&path).unwrap();
            Self(path)
        }
    }
    impl Drop for Work {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.0);
        }
    }
    fn query(positions: &[u32]) -> CohortQuery {
        CohortQuery::new(EvidenceSelection::Sites(
            positions
                .iter()
                .map(|position| SnvSite {
                    contig: 0,
                    position: *position,
                    reference: b'A',
                    alternates: vec![b'C'],
                })
                .collect(),
        ))
    }
    fn mapping(fixture: &Fixture, id: &str) -> MemberSourceMapping {
        MemberSourceMapping {
            member_id: id.into(),
            roles: fixture
                .descriptor()
                .sources
                .into_iter()
                .filter(|s| scientific_role(&s.role))
                .map(|s| (s.role, PathBuf::from(s.path)))
                .collect(),
        }
    }
    fn extend(
        parent: &SnapshotHandle,
        query: &CohortQuery,
        mappings: &[MemberSourceMapping],
        work: &Path,
    ) -> Result<ExtensionOutcome> {
        extend_snapshot(
            parent,
            query,
            mappings,
            work,
            &EvidenceExecution::default(),
            CohortQueryLimits::default(),
        )
    }
    #[derive(Default)]
    struct Collect(Vec<(u32, Option<u64>, Option<u64>)>);
    impl CohortConsumer for Collect {
        fn retained_bytes(&self) -> Option<u64> {
            Some(65536)
        }
        fn on_row(&mut self, row: CohortRow<'_>) -> Result<()> {
            self.0.push((
                row.position,
                row.evidence
                    .and_then(|r| r.depths.map(|d| d.callable_depth)),
                row.evidence
                    .and_then(|r| r.alleles.map(|a| a.allele_counts[1])),
            ));
            Ok(())
        }
    }
    fn rows(
        snapshot: &SnapshotHandle,
        query: &CohortQuery,
    ) -> Vec<(u32, Option<u64>, Option<u64>)> {
        let execution = EvidenceExecution::default();
        let plan = plan_query(snapshot, query, &execution, CohortQueryLimits::default()).unwrap();
        let mut collect = Collect::default();
        visit_rows(
            snapshot,
            &plan,
            &execution,
            CohortLimits::default(),
            &mut collect,
        )
        .unwrap();
        collect.0
    }

    #[test]
    fn extends_only_missing_loci_preserves_parent_and_survives_source_removal() {
        let fixture = Fixture::new(FixtureOptions::default());
        let work = Work::new();
        let parent = create_snapshot(
            &work.0.join("cohort"),
            &[fixture.member("A")],
            None,
            CohortLimits::default(),
        )
        .unwrap();
        let original = parent.descriptor.clone();
        let request = query(&[1, 8, 9]);
        let sources = mapping(&fixture, "A");
        let mut fresh = EvidenceEngine::open(
            source_request(
                &fixture.descriptor(),
                &sources,
                request.selection.clone(),
                &request,
                EvidenceExecution::default(),
            )
            .unwrap(),
        )
        .unwrap();
        let mut fresh_rows = Vec::new();
        fresh
            .run(&mut EvidenceCallback::with_fields(
                |batch: &EvidenceBatch| {
                    fresh_rows.extend(batch.rows().map(|row| {
                        (
                            row.position,
                            row.depths.map(|d| d.callable_depth),
                            row.alleles.map(|a| a.allele_counts[1]),
                        )
                    }));
                    Ok(())
                },
                65536,
                request.fields,
            ))
            .unwrap();
        let result = extend(&parent, &request, &[sources], &work.0).unwrap();
        assert!(result.changed);
        assert_eq!(
            result.snapshot.descriptor.parent.as_deref(),
            Some(parent.id.as_str())
        );
        assert_eq!(parent.descriptor, original);
        assert_eq!(result.snapshot.descriptor.members[0].leaves.len(), 2);
        assert_eq!(result.members[0].retained_loci, 1);
        assert_eq!(result.members[0].computed_loci, 2);
        assert_eq!(result.members[0].native_record_visits, 1);
        assert_eq!(result.members[0].full_cram_validation_records, 0);
        assert!(result.members[0].source_hashed_bytes > 0);
        assert!(result.verified_existing_bytes > 0);
        assert!(fs::read_dir(&work.0).unwrap().all(|entry| !entry
            .unwrap()
            .file_name()
            .to_string_lossy()
            .starts_with(".rosalind-extension-")));
        drop(fixture);
        // Independent fixture oracle: C at 1 and 8; position 9 is measured zero.
        assert_eq!(
            rows(&result.snapshot, &request),
            vec![
                (1, Some(1), Some(1)),
                (8, Some(1), Some(1)),
                (9, Some(0), Some(0))
            ]
        );
        assert_eq!(rows(&result.snapshot, &request), fresh_rows);
        let mut partial = request.clone();
        partial.missing_policy = MissingPolicy::Partial;
        assert_eq!(
            rows(&parent, &partial),
            vec![(1, Some(1), Some(1)), (8, None, None), (9, None, None)]
        );
        verify_snapshot(&parent.root, &parent.id, CohortLimits::default()).unwrap();
        verify_snapshot(
            &result.snapshot.root,
            &result.snapshot.id,
            CohortLimits::default(),
        )
        .unwrap();
        let no_op = extend(
            &result.snapshot,
            &request,
            &[],
            &work.0.join("does-not-exist"),
        )
        .unwrap();
        assert!(!no_op.changed);
        assert_eq!(no_op.snapshot.id, result.snapshot.id);
        assert_eq!(no_op.members[0].computed_loci, 0);
        assert_eq!(no_op.members[0].source_hashed_bytes, 0);
    }

    #[test]
    fn explicit_roles_can_relocate_but_missing_extra_and_changed_roles_refuse() {
        let fixture = Fixture::new(FixtureOptions::default());
        let work = Work::new();
        let parent = create_snapshot(
            &work.0.join("cohort"),
            &[fixture.member("A")],
            None,
            CohortLimits::default(),
        )
        .unwrap();
        let mut sources = mapping(&fixture, "A");
        let relocated = work.0.join("raw");
        fs::create_dir(&relocated).unwrap();
        for path in sources.roles.values_mut() {
            let new = relocated.join(path.file_name().unwrap());
            fs::copy(&*path, &new).unwrap();
            *path = new;
        }
        drop(fixture);
        let request = query(&[8]);
        let mut missing = sources.clone();
        missing.roles.remove("alignment-index");
        assert!(extend(&parent, &request, &[missing], &work.0)
            .unwrap_err()
            .to_string()
            .contains("roles"));
        let mut extra = sources.clone();
        extra.roles.insert("sites".into(), work.0.join("unopened"));
        assert!(extend(&parent, &request, &[extra], &work.0)
            .unwrap_err()
            .to_string()
            .contains("roles"));
        let fai = sources.roles["reference-fai"].clone();
        let original = fs::read(&fai).unwrap();
        fs::write(&fai, b"chr1\t32\t7\t32\t33\n").unwrap();
        assert!(extend(&parent, &request, &[sources.clone()], &work.0)
            .unwrap_err()
            .to_string()
            .contains("content"));
        fs::write(&fai, original).unwrap();
        let extended = extend(&parent, &request, &[sources], &work.0).unwrap();
        assert_eq!(
            rows(&extended.snapshot, &request),
            vec![(8, Some(1), Some(1))]
        );
    }

    #[test]
    fn exact_affected_members_and_covered_fields_are_required_before_raw_open() {
        let fixture = Fixture::new(FixtureOptions::default());
        let work = Work::new();
        let parent = create_snapshot(
            &work.0.join("cohort"),
            &[fixture.member("A")],
            None,
            CohortLimits::default(),
        )
        .unwrap();
        let mut sources = mapping(&fixture, "A");
        for path in sources.roles.values_mut() {
            *path = work.0.join("unopened");
        }
        assert!(extend(&parent, &query(&[8]), &[], &work.0)
            .unwrap_err()
            .to_string()
            .contains("exactly"));
        assert!(extend(&parent, &query(&[1]), &[sources.clone()], &work.0)
            .unwrap_err()
            .to_string()
            .contains("exactly"));
        assert!(extend(
            &parent,
            &query(&[8]),
            &[sources.clone(), sources.clone()],
            &work.0
        )
        .unwrap_err()
        .to_string()
        .contains("duplicate"));
        let mut fields = query(&[1, 8]);
        fields.fields = fields.fields.union(EvidenceFields::STRANDS);
        assert!(extend(&parent, &fields, &[sources], &work.0)
            .unwrap_err()
            .to_string()
            .contains("blocking"));
        assert_eq!(
            fs::read_dir(parent.root.join("snapshots")).unwrap().count(),
            1
        );
    }

    #[test]
    fn later_source_failure_does_not_publish_snapshot_or_leave_staging() {
        let a = Fixture::new(FixtureOptions::default());
        let b = Fixture::new(FixtureOptions {
            sample: Some("sample-B".into()),
            ..FixtureOptions::default()
        });
        let work = Work::new();
        let parent = create_snapshot(
            &work.0.join("cohort"),
            &[a.member("A"), b.member("B")],
            None,
            CohortLimits::default(),
        )
        .unwrap();
        let mappings = [mapping(&a, "A"), mapping(&b, "B")];
        fs::write(&mappings[1].roles["alignment-index"], b"corrupted").unwrap();
        assert!(extend(&parent, &query(&[8]), &mappings, &work.0).is_err());
        assert_eq!(
            fs::read_dir(parent.root.join("snapshots")).unwrap().count(),
            1
        );
        assert_eq!(
            fs::read_dir(parent.root.join("objects")).unwrap().count(),
            2
        );
        assert!(fs::read_dir(&work.0).unwrap().all(|entry| !entry
            .unwrap()
            .file_name()
            .to_string_lossy()
            .starts_with(".rosalind-extension-")));
        verify_snapshot(&parent.root, &parent.id, CohortLimits::default()).unwrap();
    }

    #[test]
    fn cram_extension_preserves_named_scope_and_reports_full_validation_cost() {
        let fixture = Fixture::new(FixtureOptions {
            cram: true,
            cram_reference_only: true,
            ..FixtureOptions::default()
        });
        let work = Work::new();
        let parent = create_snapshot(
            &work.0.join("cohort"),
            &[fixture.member("A")],
            None,
            CohortLimits::default(),
        )
        .unwrap();
        let result = extend(&parent, &query(&[8, 9]), &[mapping(&fixture, "A")], &work.0).unwrap();
        assert_eq!(result.members[0].computed_loci, 2);
        assert_eq!(result.members[0].full_cram_validation_records, 3);
        assert_eq!(
            rows(&result.snapshot, &query(&[8, 9])),
            vec![(8, Some(1), Some(1)), (9, Some(0), Some(0))]
        );
    }

    #[test]
    fn bounded_mapping_query_and_staging_admission_precede_output() {
        let fixture = Fixture::new(FixtureOptions::default());
        let work = Work::new();
        let parent = create_snapshot(
            &work.0.join("cohort"),
            &[fixture.member("A")],
            None,
            CohortLimits::default(),
        )
        .unwrap();
        let sources = mapping(&fixture, "A");
        let execution = EvidenceExecution::default();
        let mut limits = CohortQueryLimits {
            max_input_sites: 1,
            ..CohortQueryLimits::default()
        };
        assert!(extend_snapshot(
            &parent,
            &query(&[8, 9]),
            &[sources.clone()],
            &work.0,
            &execution,
            limits
        )
        .unwrap_err()
        .to_string()
        .contains("site envelope"));
        limits = CohortQueryLimits {
            max_plan_bytes: 1,
            ..CohortQueryLimits::default()
        };
        assert!(extend_snapshot(
            &parent,
            &query(&[8]),
            &[sources.clone()],
            &work.0,
            &execution,
            limits
        )
        .is_err());
        assert!(extend(&parent, &query(&[8]), &[sources], &parent.root)
            .unwrap_err()
            .to_string()
            .contains("outside"));
        assert_eq!(
            fs::read_dir(parent.root.join("snapshots")).unwrap().count(),
            1
        );
    }

    #[test]
    fn changed_candidate_or_source_table_guard_refuses_before_publication() {
        let fixture = Fixture::new(FixtureOptions::default());
        let work = Work::new();
        let parent = create_snapshot(
            &work.0.join("cohort"),
            &[fixture.member("A")],
            None,
            CohortLimits::default(),
        )
        .unwrap();
        let input = work.0.join("candidates.vcf");
        fs::write(&input, b"original candidate input").unwrap();
        let inputs = InputSnapshot::capture([input.clone()]).unwrap();
        fs::write(&input, b"changed candidate input").unwrap();
        assert!(extend_snapshot_with_inputs(
            &parent,
            &query(&[8]),
            &[mapping(&fixture, "A")],
            &work.0,
            &EvidenceExecution::default(),
            CohortQueryLimits::default(),
            Some(&inputs)
        )
        .is_err());
        assert_eq!(
            fs::read_dir(parent.root.join("snapshots")).unwrap().count(),
            1
        );
    }
}
