//! One managed lifetime for a complete cohort result. All inputs remain saved
//! evidence; the runner never invokes the single-sample managed runner recursively.

use super::descriptor::{CohortLimits, SNAPSHOT_FILE, SNAPSHOT_RECEIPT};
use super::encoding::{
    memory_bytes, summary_memory_bytes, CohortOutputFormat, ExtractEncoder, SummaryEncoder,
};
use super::query::{plan_query, CohortQuery, CohortQueryLimits, CohortQueryPlan};
use super::runtime::{visit_rows, visit_windows, CohortRunStats};
use super::store::open_snapshot;
use super::summary::required_fields;
use super::{CohortError, Result};
use crate::contract::{EnforcementMode, OutputPolicy};
use crate::core::cancellation::{CancellationScope, CancellationToken, SignalCancellationGuard};
use crate::core::governor::{checkpoint, MemoryGovernor};
use crate::dataset::{canonical_dataset_query, DatasetQuery, InputSnapshot};
use crate::evidence::{
    EvidenceArtifactError, EvidenceError, EvidenceExecution, EvidenceSelection,
    CANONICAL_TILE_BASES,
};
use crate::provenance::{CommandCapture, FileHash, RunManifest};
use crate::util::atomic::{commit_group, ensure_destination, AtomicFile};
use crate::util::rss::peak_rss_bytes;
use crate::variant_io::{parse_snv_record, VariantLimits, VariantReader};
use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::{BufWriter, Read, Seek, Write};
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};

#[derive(Debug, Clone, Copy)]
pub(crate) enum CohortOperation {
    Extract,
    Summarize,
}

#[derive(Debug, Clone)]
pub(crate) struct CohortArtifactSpec {
    pub root: PathBuf,
    pub snapshot_id: String,
    pub query: CohortQuery,
    /// A bounded VCF/VCF.gz/BCF selection parsed inside this managed lifetime.
    pub sites: Option<(PathBuf, VariantLimits)>,
    /// Explicit content-located replay operands; each must match a consumed input.
    pub artifacts: Vec<PathBuf>,
    pub execution: EvidenceExecution,
    pub limits: CohortQueryLimits,
    pub operation: CohortOperation,
    pub format: CohortOutputFormat,
    pub min_callable_depth: u64,
    pub output: PathBuf,
    pub manifest: Option<PathBuf>,
    pub output_policy: OutputPolicy,
    pub enforcement: EnforcementMode,
    pub cancellation: Option<CancellationToken>,
    pub handle_signals: bool,
    pub max_receipt_bytes: usize,
}

#[derive(Debug)]
pub(crate) struct CohortArtifactOutcome {
    pub output: PathBuf,
    pub manifest: PathBuf,
    pub claim_hash: String,
    pub stats: CohortRunStats,
    pub peak_rss_bytes: u64,
}

/// Refusal, cancellation and failed output/receipt publication never create a
/// completed result. Resource failures currently discard staged cohort files;
/// partial cohort artifacts are not a supported product surface.
pub(crate) fn run_cohort_artifact(spec: CohortArtifactSpec) -> Result<CohortArtifactOutcome> {
    match execute(spec, false)? {
        ArtifactExecution::Completed(outcome) => Ok(outcome),
        ArtifactExecution::Planned(_) => unreachable!("execution requested"),
    }
}

pub(crate) fn plan_cohort_artifact(spec: CohortArtifactSpec) -> Result<CohortQueryPlan> {
    match execute(spec, true)? {
        ArtifactExecution::Planned(plan) => Ok(plan),
        ArtifactExecution::Completed(_) => unreachable!("plan requested"),
    }
}

enum ArtifactExecution {
    Planned(CohortQueryPlan),
    Completed(CohortArtifactOutcome),
}

fn execute(mut spec: CohortArtifactSpec, plan_only: bool) -> Result<ArtifactExecution> {
    let started = Instant::now();
    let token = spec.cancellation.clone().unwrap_or_default();
    let cancellation = CancellationScope::start(token.clone())
        .map_err(|error| CohortError::Incompatible(error.to_string()))?;
    let _signals = if spec.handle_signals {
        Some(SignalCancellationGuard::install(&cancellation)?)
    } else {
        None
    };
    checkpoint().map_err(EvidenceError::from)?;
    if spec.min_callable_depth == 0
        || spec.max_receipt_bytes == 0
        || spec.execution.max_microtile_bases == 0
    {
        return Err(CohortError::Incompatible(
            "depth, execution width and receipt envelope must be positive".into(),
        ));
    }
    let declared = [
        spec.execution.memory_budget_bytes,
        spec.limits.cohort.effective_budget(),
    ]
    .into_iter()
    .flatten()
    .min();
    let enforced = spec.enforcement != EnforcementMode::RecordOnly;
    if enforced && declared.is_none() {
        return Err(CohortError::Incompatible(
            "memory enforcement requires a declared budget".into(),
        ));
    }
    let budget = declared.filter(|_| enforced);
    if let Some(budget) = budget {
        if peak_rss_bytes() > budget {
            return Err(EvidenceArtifactError::Refused {
                needed: peak_rss_bytes(),
                budget,
            }
            .into());
        }
    }
    let os_limit = if spec.enforcement == EnforcementMode::RequireOsLimit {
        let limit = crate::contract::detected_os_memory_limit_bytes().ok_or_else(|| {
            EvidenceArtifactError::OsEnforcementUnavailable(
                "require an existing Linux cgroup-v2 memory.max".into(),
            )
        })?;
        if limit > budget.unwrap() {
            return Err(EvidenceArtifactError::OsEnforcementUnavailable(
                "cgroup memory.max exceeds declared budget".into(),
            )
            .into());
        }
        Some(limit)
    } else {
        None
    };
    let _governor = budget
        .map(|bytes| MemoryGovernor::start(bytes, Duration::from_millis(100), peak_rss_bytes))
        .transpose()
        .map_err(|error| CohortError::Incompatible(error.to_string()))?;
    spec.execution.memory_budget_bytes = budget;
    spec.limits.cohort.memory_budget_bytes = budget;
    spec.limits.cohort.dataset.memory_budget_bytes = budget;
    let limits = spec.limits.cohort;
    limits.admit(0)?;
    if spec.artifacts.len() > spec.limits.cohort.max_leaf_references {
        return Err(CohortError::Limit(
            "replay operands exceed metadata envelope".into(),
        ));
    }
    limits.admit((spec.artifacts.len() as u64).saturating_mul(16_384))?;
    spec.artifacts = spec
        .artifacts
        .iter()
        .map(fs::canonicalize)
        .collect::<std::io::Result<_>>()?;
    let replay_guard = InputSnapshot::capture(spec.artifacts.iter().cloned())?;
    let snapshot = open_snapshot(&spec.root, &spec.snapshot_id, limits)?;
    spec.query.requirements.fields = required_fields();
    spec.query.requirements.requires_reference = true;
    spec.query.requirements.context_bases = 0;
    // The first metadata plan establishes dictionary and coverage; no output has
    // been created. A second plan admits the complete encoder/finalization bound.
    spec.query.requirements.retained_bytes = Some(0);
    let mut dictionary_query = spec.query.clone();
    if spec.sites.is_some() {
        dictionary_query.selection = EvidenceSelection::Sites(Vec::new());
    }
    let mut plan = plan_query(&snapshot, &dictionary_query, &spec.execution, spec.limits)?;
    let candidate_input = if let Some((path, variant_limits)) = &spec.sites {
        let path = fs::canonicalize(path)?;
        let guard = InputSnapshot::capture([path.clone()])?;
        let hash = hash_file(&path, limits)?;
        let parser_bytes = (variant_limits.max_header_bytes as u64)
            .saturating_add(variant_limits.max_record_bytes as u64)
            .saturating_mul(4)
            .saturating_add(128 << 10);
        limits.admit(parser_bytes)?;
        let mut reader = VariantReader::open(&path, *variant_limits)?;
        let mut sites = Vec::new();
        while let Some(record) = reader.read()? {
            if sites.len() >= spec.limits.max_input_sites {
                return Err(CohortError::Limit(
                    "candidate file exceeds record envelope".into(),
                ));
            }
            let bytes = (sites.len() as u64 + 1).saturating_mul(4096);
            if bytes > spec.limits.max_query_bytes {
                return Err(CohortError::Limit(
                    "candidate file exceeds query byte envelope".into(),
                ));
            }
            limits.admit(parser_bytes.saturating_add(bytes))?;
            sites.push(parse_snv_record(record, &plan.contigs)?);
        }
        drop(reader);
        guard.verify()?;
        spec.query.selection = EvidenceSelection::Sites(sites);
        plan = plan_query(&snapshot, &spec.query, &spec.execution, spec.limits)?;
        Some((path, hash, guard))
    } else {
        None
    };
    if !plan_only {
        plan.ensure_executable()?;
    }
    let max_text = plan
        .contigs
        .iter()
        .map(|contig| contig.name.len())
        .max()
        .unwrap_or(0)
        .max(256);
    if max_text > 4096 {
        return Err(CohortError::Limit(
            "cohort report text exceeds 4096-byte envelope".into(),
        ));
    }
    let max_window_candidates = (spec.execution.max_microtile_bases.min(CANONICAL_TILE_BASES)
        as usize)
        .checked_mul(3)
        .ok_or_else(|| CohortError::Limit("window candidate count overflow".into()))?;
    let encoder_bytes = match spec.operation {
        CohortOperation::Extract => memory_bytes(plan.fields, max_text),
        CohortOperation::Summarize => summary_memory_bytes(max_window_candidates, max_text),
    };
    let leaf_count = plan.members.iter().try_fold(0u64, |sum, member| {
        sum.checked_add(member.leaves.len() as u64)
            .ok_or_else(|| CohortError::Limit("lineage count overflow".into()))
    })?;
    // Reserve before duplicating canonical query/member identities or formatting
    // lineage. Receipt JSON escaping and finalization use the same admitted bound.
    let receipt_reserve = plan
        .reservations
        .query_bytes
        .saturating_mul(4)
        .saturating_add(leaf_count.saturating_mul(16_384))
        .saturating_add((plan.members.len() as u64).saturating_mul(4096))
        .saturating_add(2 << 20);
    let reserve = encoder_bytes
        .checked_add(receipt_reserve)
        .ok_or_else(|| CohortError::Limit("encoder/receipt reservation overflow".into()))?;
    limits.admit(reserve)?;
    spec.query.requirements.retained_bytes = Some(reserve);
    plan = plan_query(&snapshot, &spec.query, &spec.execution, spec.limits)?;
    if plan_only {
        if let Some((_, _, guard)) = &candidate_input {
            guard.verify()?;
        }
        return Ok(ArtifactExecution::Planned(plan));
    }
    plan.ensure_executable()?;
    let receipt_path = spec
        .manifest
        .clone()
        .unwrap_or_else(|| suffix(&spec.output, ".manifest.json"));
    let output = absolute_destination(&spec.output)?;
    let manifest = absolute_destination(&receipt_path)?;
    if output == manifest
        || [&output, &manifest]
            .iter()
            .any(|path| path.starts_with(&snapshot.root))
        || candidate_input
            .as_ref()
            .is_some_and(|(path, _, _)| path == &output || path == &manifest)
        || spec
            .artifacts
            .iter()
            .any(|path| path == &output || path == &manifest)
    {
        return Err(CohortError::Incompatible(
            "cohort outputs must be distinct and outside the immutable store".into(),
        ));
    }
    let replace = spec.output_policy == OutputPolicy::ReplaceAtomic;
    ensure_destination(&output, replace)?;
    ensure_destination(&manifest, replace)?;
    let normalized = canonical_dataset_query(
        &DatasetQuery {
            selection: EvidenceSelection::Sites(plan.sites.clone()),
            fields: plan.fields,
        },
        &plan.contigs,
    )?;
    if candidate_input.is_none() && normalized.len() > crate::evidence::MAX_ARTIFACT_QUERY_BYTES {
        return Err(CohortError::Limit(
            "inline replay query exceeds 32 KiB; supply a candidate file instead".into(),
        ));
    }
    let selected = serde_json::to_string(
        &plan
            .members
            .iter()
            .map(|member| member.member_id.as_str())
            .collect::<Vec<_>>(),
    )
    .map_err(|error| CohortError::Corrupt(error.to_string()))?;
    let operation = match spec.operation {
        CohortOperation::Extract => "extract",
        CohortOperation::Summarize => "summarize",
    };
    let format = match spec.format {
        CohortOutputFormat::Arrow => "arrow-ipc",
        CohortOutputFormat::Tsv => "tsv",
    };
    let mut receipt = RunManifest::new(format!("cohort {operation}"));
    for (key, value) in [
        ("run_status", "completed".into()),
        ("cohort.result_semantics", "cohort-candidate-v1".into()),
        ("cohort.snapshot_blake3", snapshot.id.clone()),
        (
            "cohort.comparison_blake3",
            plan.comparison_blake3.clone().unwrap_or_default(),
        ),
        ("cohort.normalized_query", normalized.clone()),
        (
            "cohort.query_blake3",
            blake3::hash(normalized.as_bytes()).to_hex().to_string(),
        ),
        ("cohort.members", selected),
        (
            "cohort.missing_policy",
            match plan.missing_policy {
                super::query::MissingPolicy::Strict => "strict",
                super::query::MissingPolicy::Partial => "partial",
            }
            .into(),
        ),
        (
            "cohort.min_callable_depth",
            spec.min_callable_depth.to_string(),
        ),
        ("cohort.counting_unit", "read".into()),
        ("cohort.operation", operation.into()),
        ("cohort.format", format.into()),
        (
            "contract.assurance",
            match spec.enforcement {
                EnforcementMode::RecordOnly => "observed-only",
                EnforcementMode::Cooperative => "declared-bound-cooperative",
                EnforcementMode::RequireOsLimit => "cgroup-v2",
            }
            .into(),
        ),
        ("contract_verdict", "within".into()),
        ("artifact.output.0.role", "cohort-candidate-result".into()),
        (
            "resource.peak_sampling_phase",
            "post-encoding-sync-and-source-validation-before-atomic-commit".into(),
        ),
    ] {
        receipt.params.insert(key.into(), value);
    }
    if let Some(bytes) = declared {
        receipt
            .params
            .insert("memory_budget_bytes".into(), bytes.to_string());
    }
    let mut lineage = BTreeMap::new();
    if let Some((path, hash, guard)) = &candidate_input {
        guard.verify()?;
        lineage.insert(path.clone(), hash.clone());
    }
    let snapshot_path = snapshot
        .root
        .join("snapshots")
        .join(&snapshot.id)
        .join(SNAPSHOT_FILE);
    lineage.insert(snapshot_path.clone(), snapshot.id.clone());
    let snapshot_receipt = snapshot
        .root
        .join("snapshots")
        .join(&snapshot.id)
        .join(SNAPSHOT_RECEIPT);
    lineage.insert(
        snapshot_receipt.clone(),
        hash_file(&snapshot_receipt, limits)?,
    );
    for member in &plan.members {
        for leaf_plan in &member.leaves {
            if leaf_plan.covered_loci.unwrap_or(0) == 0 {
                continue;
            }
            let leaf =
                &snapshot.descriptor.members[member.member_index].leaves[leaf_plan.leaf_index];
            let dataset = snapshot.open_leaf(leaf, limits)?;
            for identity in dataset.source_hashes() {
                lineage.insert(PathBuf::from(&identity.path), identity.blake3.clone());
            }
            dataset.verify_unchanged()?;
        }
    }
    let lineage_guard = InputSnapshot::capture(lineage.keys().cloned())?;
    receipt.inputs = lineage
        .into_iter()
        .map(|(path, blake3)| FileHash {
            path: path.display().to_string(),
            blake3,
        })
        .collect();
    let estimate = receipt.to_canonical_json().len();
    if estimate.saturating_add(8192) > spec.max_receipt_bytes {
        return Err(CohortError::Limit(
            "cohort receipt exceeds envelope; reduce the query or member selection".into(),
        ));
    }
    snapshot.verify_unchanged()?;
    limits.admit(0)?;
    let mut pending = AtomicFile::create(&output)?;
    let stats = {
        let mut writer = BufWriter::new(pending.file_mut());
        let stats = match spec.operation {
            CohortOperation::Extract => {
                let mut encoder = ExtractEncoder::new(
                    &mut writer,
                    spec.format,
                    &plan.contigs,
                    plan.fields,
                    spec.min_callable_depth,
                )?;
                visit_rows(&snapshot, &plan, &spec.execution, limits, &mut encoder)?
            }
            CohortOperation::Summarize => {
                let mut encoder = SummaryEncoder::new(
                    &mut writer,
                    spec.format,
                    &plan.contigs,
                    plan.members.len() as u64,
                    spec.min_callable_depth,
                    max_window_candidates,
                )?;
                visit_windows(&snapshot, &plan, &spec.execution, limits, &mut encoder)?
            }
        };
        writer.flush()?;
        stats
    };
    pending.file_mut().sync_all()?;
    for (path, blake3) in &stats.consumed_inputs {
        receipt.inputs.push(FileHash {
            path: path.clone(),
            blake3: blake3.clone(),
        });
    }
    for path in &spec.artifacts {
        let hash = hash_file(path, limits)?;
        if !receipt.inputs.iter().any(|input| input.blake3 == hash) {
            return Err(CohortError::Incompatible(
                "replay operand does not match a verified cohort dependency".into(),
            ));
        }
    }
    replay_guard.verify()?;
    receipt.outputs.push(FileHash {
        path: output.display().to_string(),
        blake3: hash_file(pending.temporary_path(), limits)?,
    });
    let mut capture = CommandCapture::from_argv_prefix(["cohort", operation]);
    for input in &receipt.inputs {
        let flag = if Path::new(&input.path) == snapshot_path {
            "--cohort-snapshot-manifest"
        } else if candidate_input
            .as_ref()
            .is_some_and(|(path, _, _)| Path::new(&input.path) == path)
        {
            "--sites"
        } else {
            "--cohort-artifact"
        };
        capture.input_hashed(flag, &input.path, &input.blake3);
    }
    if candidate_input.is_none() {
        capture.opt("--cohort-query", &normalized);
    }
    capture
        .opt("--snapshot", &snapshot.id)
        .opt("--cohort-members", &receipt.params["cohort.members"])
        .opt("--fields", plan.fields.names().join(","))
        .opt("--missing", &receipt.params["cohort.missing_policy"])
        .opt("--min-callable-depth", spec.min_callable_depth)
        .opt("--format", format)
        .opt("--max-microtile-bases", spec.execution.max_microtile_bases)
        .opt("--max-receipt-bytes", spec.max_receipt_bytes)
        .opt(
            "--max-snapshot-bytes",
            spec.limits.cohort.max_snapshot_bytes,
        )
        .opt(
            "--max-dataset-metadata-bytes",
            spec.limits
                .cohort
                .dataset
                .max_manifest_bytes
                .max(spec.limits.cohort.dataset.max_descriptor_bytes),
        )
        .opt("--max-candidate-sites", spec.limits.max_input_sites)
        .flag_if(enforced, "--enforce")
        .flag_if(
            spec.enforcement == EnforcementMode::RequireOsLimit,
            "--require-os-limit",
        );
    if let Some((_, variant_limits)) = &spec.sites {
        capture
            .opt(
                "--max-variant-header-bytes",
                variant_limits.max_header_bytes,
            )
            .opt(
                "--max-variant-record-bytes",
                variant_limits.max_record_bytes,
            );
    }
    if let Some(bytes) = declared {
        capture.opt("--memory-budget-bytes", bytes);
    }
    for output in &receipt.outputs {
        capture.output_hashed("--output", &output.path, &output.blake3);
    }
    capture.record_into(&mut receipt);
    receipt.measurements.extend([
        ("peak_rss_bytes".into(), "00000000000000000000".into()),
        (
            "elapsed_ms".into(),
            started.elapsed().as_millis().to_string(),
        ),
        (
            "execution.predicted_peak_rss_bytes".into(),
            stats.predicted_peak_rss_bytes.to_string(),
        ),
        (
            "execution.window_bases".into(),
            stats.execution_window_bases.to_string(),
        ),
        (
            "execution.original_alignment_records_decoded".into(),
            "0".into(),
        ),
        (
            "outcome.sample_candidate_rows".into(),
            stats.emitted_rows.to_string(),
        ),
        (
            "outcome.observed_rows".into(),
            stats.observed_rows.to_string(),
        ),
        (
            "outcome.unmeasured_rows".into(),
            stats.unmeasured_rows.to_string(),
        ),
    ]);
    if let Some(limit) = os_limit {
        receipt
            .measurements
            .insert("resource.os_limit_bytes".into(), limit.to_string());
    }
    receipt.finalize();
    let claim_hash = receipt.params["manifest_blake3"].clone();
    let mut bytes = receipt.to_canonical_json();
    if bytes.len() > spec.max_receipt_bytes {
        return Err(CohortError::Limit(
            "cohort receipt exceeds final envelope".into(),
        ));
    }
    let mut pending_receipt = AtomicFile::create(&manifest)?;
    pending_receipt.file_mut().write_all(bytes.as_bytes())?;
    pending_receipt.file_mut().sync_all()?;
    test_before_publication()?;
    snapshot.verify_unchanged()?;
    lineage_guard.verify()?;
    if let Some((_, _, guard)) = &candidate_input {
        guard.verify()?;
    }
    limits.admit(0)?;
    let peak = peak_rss_bytes();
    crate::evidence::reseal_observed_measurements(&mut bytes, peak, declared)?;
    pending_receipt.file_mut().rewind()?;
    pending_receipt.file_mut().set_len(0)?;
    pending_receipt.file_mut().write_all(bytes.as_bytes())?;
    pending_receipt.file_mut().sync_all()?;
    snapshot.verify_unchanged()?;
    lineage_guard.verify()?;
    if let Some((_, _, guard)) = &candidate_input {
        guard.verify()?;
    }
    limits.admit(0)?;
    token.check().map_err(EvidenceError::from)?;
    replay_guard.verify()?;
    commit_group(
        vec![
            (pending, output.clone()),
            (pending_receipt, manifest.clone()),
        ],
        replace,
    )?;
    Ok(ArtifactExecution::Completed(CohortArtifactOutcome {
        output,
        manifest,
        claim_hash,
        stats,
        peak_rss_bytes: peak,
    }))
}

fn absolute_destination(path: &Path) -> Result<PathBuf> {
    let filename = path
        .file_name()
        .filter(|name| !name.is_empty())
        .ok_or_else(|| CohortError::Incompatible("output needs a file name".into()))?;
    let parent = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    Ok(fs::canonicalize(parent)?.join(filename))
}
fn suffix(path: &Path, ending: &str) -> PathBuf {
    let mut name = path.as_os_str().to_owned();
    name.push(ending);
    name.into()
}
fn hash_file(path: &Path, limits: CohortLimits) -> Result<String> {
    let mut file = File::open(path)?;
    let mut hash = blake3::Hasher::new();
    let mut buffer = [0u8; 64 << 10];
    loop {
        limits.admit(0)?;
        let count = file.read(&mut buffer)?;
        if count == 0 {
            break;
        }
        hash.update(&buffer[..count]);
    }
    Ok(hash.finalize().to_hex().to_string())
}

#[cfg(not(test))]
fn test_before_publication() -> Result<()> {
    Ok(())
}
#[cfg(test)]
type BeforePublicationHook = Box<dyn FnOnce() -> Result<()>>;
#[cfg(test)]
thread_local! { static BEFORE_PUBLICATION: std::cell::RefCell<Option<BeforePublicationHook>> = const { std::cell::RefCell::new(None) }; }
#[cfg(test)]
fn test_before_publication() -> Result<()> {
    BEFORE_PUBLICATION
        .with(|hook| hook.borrow_mut().take())
        .map_or(Ok(()), |hook| hook())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cohort::query::MissingPolicy;
    use crate::cohort::store::create_snapshot;
    use crate::cohort::tests::{Fixture, FixtureOptions};
    use crate::evidence::SnvSite;

    fn specification(root: &Path, snapshot: &str, output: &Path) -> CohortArtifactSpec {
        let mut query = CohortQuery::new(EvidenceSelection::Sites(
            [1, 2, 8]
                .into_iter()
                .map(|position| SnvSite {
                    contig: 0,
                    position,
                    reference: b'A',
                    alternates: vec![b'C'],
                })
                .collect(),
        ));
        query.missing_policy = MissingPolicy::Partial;
        CohortArtifactSpec {
            root: root.to_owned(),
            snapshot_id: snapshot.into(),
            query,
            sites: None,
            artifacts: Vec::new(),
            execution: EvidenceExecution {
                memory_budget_bytes: Some(512 << 20),
                ..EvidenceExecution::default()
            },
            limits: CohortQueryLimits::default(),
            operation: CohortOperation::Extract,
            format: CohortOutputFormat::Tsv,
            min_callable_depth: 10,
            output: output.to_owned(),
            manifest: None,
            output_policy: OutputPolicy::CreateNewAtomic,
            enforcement: EnforcementMode::Cooperative,
            cancellation: None,
            handle_signals: false,
            max_receipt_bytes: 32 << 20,
        }
    }

    // Scope/governor tests run in a dedicated process so cancellation cannot
    // interfere with unrelated parallel evidence tests in this Rust test binary.
    #[test]
    fn managed_cohort_outputs_and_failure_boundaries() {
        const CHILD: &str = "ROSALIND_COHORT_ARTIFACT_TEST_CHILD";
        if std::env::var_os(CHILD).is_none() {
            let result = std::process::Command::new(std::env::current_exe().unwrap())
                .args([
                    "--exact",
                    "cohort::artifact::tests::managed_cohort_outputs_and_failure_boundaries",
                    "--test-threads=1",
                    "--nocapture",
                ])
                .env(CHILD, "1")
                .output()
                .unwrap();
            assert!(
                result.status.success(),
                "{}\n{}",
                String::from_utf8_lossy(&result.stdout),
                String::from_utf8_lossy(&result.stderr)
            );
            return;
        }
        struct Temporary(PathBuf);
        impl Temporary {
            fn path(&self) -> &Path {
                &self.0
            }
        }
        impl Drop for Temporary {
            fn drop(&mut self) {
                let _ = fs::remove_dir_all(&self.0);
            }
        }
        let temporary = Temporary(std::env::temp_dir().join(format!(
            "rosalind-cohort-artifact-test-{}",
            std::process::id()
        )));
        fs::create_dir(temporary.path()).unwrap();
        let store = temporary.path().join("cohort");
        let a = Fixture::new(FixtureOptions::default());
        let b = Fixture::new(FixtureOptions {
            sample: Some("second".into()),
            ..FixtureOptions::default()
        });
        let snapshot = create_snapshot(
            &store,
            &[a.member("a"), b.member("b")],
            None,
            CohortLimits::default(),
        )
        .unwrap();
        drop(a);
        drop(b);
        for operation in [CohortOperation::Extract, CohortOperation::Summarize] {
            for format in [CohortOutputFormat::Tsv, CohortOutputFormat::Arrow] {
                let mut expected = None;
                for (index, (budget, width)) in [(384u64, 1u32), (512, 3), (768, 16_384)]
                    .into_iter()
                    .enumerate()
                {
                    let output = temporary
                        .path()
                        .join(format!("{operation:?}-{format:?}-{index}"));
                    let mut spec = specification(&store, &snapshot.id, &output);
                    spec.operation = operation;
                    spec.format = format;
                    spec.execution.memory_budget_bytes = Some(budget << 20);
                    spec.execution.max_microtile_bases = width;
                    let outcome = run_cohort_artifact(spec).unwrap();
                    assert_eq!(outcome.stats.emitted_rows, 6);
                    assert_eq!(outcome.stats.observed_rows, 4);
                    assert_eq!(outcome.stats.unmeasured_rows, 2);
                    assert!(outcome.peak_rss_bytes <= budget << 20);
                    assert!(outcome.stats.predicted_peak_rss_bytes >= outcome.peak_rss_bytes);
                    let bytes = fs::read(&outcome.output).unwrap();
                    if let Some(expected) = &expected {
                        assert_eq!(&bytes, expected);
                    } else {
                        expected = Some(bytes);
                    }
                    let manifest = RunManifest::from_canonical_json(
                        &fs::read_to_string(&outcome.manifest).unwrap(),
                    )
                    .unwrap();
                    assert_eq!(manifest.self_hash_ok(), Some(true));
                    assert_eq!(manifest.measurement_hash_ok(), Some(true));
                    assert_eq!(manifest.params["manifest_blake3"], outcome.claim_hash);
                    assert_eq!(manifest.params["cohort.snapshot_blake3"], snapshot.id);
                    assert_eq!(manifest.params["replay_schema"], "3");
                    assert!(manifest.params["command_argv"].contains("--cohort-snapshot-manifest"));
                    assert!(manifest.params["command_argv"].contains("--cohort-query"));
                    assert_eq!(outcome.stats.consumed_inputs.len(), 4);
                    for (path, hash) in &outcome.stats.consumed_inputs {
                        assert!(manifest
                            .inputs
                            .iter()
                            .any(|input| &input.path == path && &input.blake3 == hash));
                    }
                    assert_eq!(
                        manifest.measurements["execution.original_alignment_records_decoded"],
                        "0"
                    );
                    assert_eq!(
                        manifest.outputs[0].blake3,
                        blake3::hash(&fs::read(&output).unwrap())
                            .to_hex()
                            .to_string()
                    );
                }
            }
        }
        let output = temporary.path().join("refused");
        let mut strict = specification(&store, &snapshot.id, &output);
        strict.query.missing_policy = MissingPolicy::Strict;
        assert!(run_cohort_artifact(strict).is_err());
        assert!(!output.exists());
        assert!(!suffix(&output, ".manifest.json").exists());
        let mut tiny = specification(&store, &snapshot.id, &output);
        tiny.execution.memory_budget_bytes = Some(1);
        assert!(run_cohort_artifact(tiny).is_err());
        assert!(!output.exists());
        let cancelled = CancellationToken::new();
        cancelled.cancel();
        let mut pre_cancelled = specification(&store, &snapshot.id, &output);
        pre_cancelled.cancellation = Some(cancelled);
        assert!(run_cohort_artifact(pre_cancelled).is_err());
        assert!(!output.exists());
        let final_token = CancellationToken::new();
        let mut final_cancel = specification(&store, &snapshot.id, &output);
        final_cancel.cancellation = Some(final_token.clone());
        BEFORE_PUBLICATION.with(|hook| {
            *hook.borrow_mut() = Some(Box::new(move || {
                final_token.cancel();
                Ok(())
            }))
        });
        assert!(run_cohort_artifact(final_cancel).is_err());
        assert!(!output.exists());
        assert!(!suffix(&output, ".manifest.json").exists());
        let mut protected = specification(&store, &snapshot.id, &store.join("inside"));
        assert!(run_cohort_artifact(protected.clone()).is_err());
        protected.output = temporary.path().join("outside");
        protected.manifest = Some(protected.output.clone());
        assert!(run_cohort_artifact(protected).is_err());
        assert!(!temporary.path().join("outside").exists());

        // File selection is parsed inside the managed lifetime; a strict plan
        // explains missing loci without creating an output. The file itself and
        // every consumed partition remain receipt inputs for relocated replay.
        let candidates = temporary.path().join("candidates.vcf");
        let vcf = "##fileformat=VCFv4.3\n##contig=<ID=chr1,length=32>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t2\t.\tA\tC\t.\tPASS\t.\nchr1\t3\t.\tA\tC\t.\tPASS\t.\nchr1\t9\t.\tA\tC\t.\tPASS\t.\n";
        fs::write(&candidates, vcf).unwrap();
        let mut file_spec = specification(&store, &snapshot.id, &output);
        file_spec.sites = Some((candidates.clone(), VariantLimits::default()));
        file_spec.query.missing_policy = MissingPolicy::Strict;
        let plan = plan_cohort_artifact(file_spec.clone()).unwrap();
        assert_eq!(plan.issue_count, 2);
        assert!(!output.exists());
        file_spec.query.missing_policy = MissingPolicy::Partial;
        let outcome = run_cohort_artifact(file_spec.clone()).unwrap();
        let receipt =
            RunManifest::from_canonical_json(&fs::read_to_string(&outcome.manifest).unwrap())
                .unwrap();
        assert!(receipt.params["command_argv"].contains("--sites"));
        assert!(!receipt.params["command_argv"].contains("--cohort-query"));
        assert!(receipt.inputs.iter().any(
            |input| input.path == fs::canonicalize(&candidates).unwrap().display().to_string()
        ));
        fs::remove_file(&output).unwrap();
        fs::remove_file(&outcome.manifest).unwrap();

        let mut overwrite = file_spec.clone();
        overwrite.output = candidates.clone();
        overwrite.output_policy = OutputPolicy::ReplaceAtomic;
        assert!(run_cohort_artifact(overwrite).is_err());
        assert_eq!(fs::read_to_string(&candidates).unwrap(), vcf);
        let copy = temporary.path().join("dependency-copy");
        fs::copy(&candidates, &copy).unwrap();
        let mut replay_overwrite = file_spec.clone();
        replay_overwrite.artifacts = vec![copy.clone()];
        replay_overwrite.manifest = Some(copy.clone());
        replay_overwrite.output_policy = OutputPolicy::ReplaceAtomic;
        assert!(run_cohort_artifact(replay_overwrite).is_err());
        assert_eq!(fs::read_to_string(&copy).unwrap(), vcf);
        let candidates_for_hook = candidates.clone();
        BEFORE_PUBLICATION.with(|hook| {
            *hook.borrow_mut() = Some(Box::new(move || {
                fs::write(candidates_for_hook, b"changed candidate input")?;
                Ok(())
            }))
        });
        assert!(run_cohort_artifact(file_spec).is_err());
        assert!(!output.exists());
        assert!(!suffix(&output, ".manifest.json").exists());

        // Failure at the second publication rolls back the first output. Existing
        // sentinels survive; this uses the same group helper as single-sample artifacts.
        let bad_receipt = suffix(&output, ".manifest.json");
        let bad_receipt_for_hook = bad_receipt.clone();
        BEFORE_PUBLICATION.with(|hook| {
            *hook.borrow_mut() = Some(Box::new(move || {
                fs::create_dir(&bad_receipt_for_hook)?;
                Ok(())
            }))
        });
        assert!(run_cohort_artifact(specification(&store, &snapshot.id, &output)).is_err());
        assert!(!output.exists());
        assert!(bad_receipt.is_dir());
        fs::remove_dir(&bad_receipt).unwrap();
        // A late dataset mutation cannot result in a successful artifact/receipt.
        let object = &snapshot.descriptor.members[0].leaves[0];
        let payload = snapshot
            .open_leaf(object, CohortLimits::default())
            .unwrap()
            .descriptor()
            .partitions[0]
            .arrow
            .path
            .clone();
        let payload = store.join("objects").join(&object.object_id).join(payload);
        BEFORE_PUBLICATION.with(|hook| {
            *hook.borrow_mut() = Some(Box::new(move || {
                fs::write(payload, b"changed")?;
                Ok(())
            }))
        });
        assert!(run_cohort_artifact(specification(&store, &snapshot.id, &output)).is_err());
        assert!(!output.exists());
        assert!(!bad_receipt.exists());
    }
}
