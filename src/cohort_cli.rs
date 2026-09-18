//! Private implementation bridge for the experimental cohort command family.
//! This hidden module is callable by the binary, not a stable Rust cohort SDK.
#![allow(missing_docs)]

use crate::cohort::artifact::{
    plan_cohort_artifact, run_cohort_artifact, CohortArtifactSpec, CohortOperation,
};
use crate::cohort::descriptor::{CohortLimits, MemberMetadata};
use crate::cohort::encoding::CohortOutputFormat;
use crate::cohort::pair_table::PairInput;
use crate::cohort::pairs::PairLimits;
use crate::cohort::query::{
    plan_query, CohortQuery, CohortQueryLimits, CohortQueryPlan, MissingPolicy,
};
use crate::cohort::store::{create_snapshot, open_snapshot, verify_snapshot, ImportMember};
use crate::cohort::summary::{required_fields, DEFAULT_MIN_CALLABLE_DEPTH};
use crate::cohort::CohortError;
use crate::contract::{EnforcementMode, OutputPolicy};
use crate::core::cancellation::{CancellationScope, CancellationToken, SignalCancellationGuard};
use crate::core::governor::{checkpoint, MemoryGovernor};
use crate::core::CoreError;
use crate::dataset::{InputSnapshot, VerifiedEvidenceDataset};
use crate::evidence::{
    EvidenceArtifactError, EvidenceError, EvidenceExecution, EvidenceFields, EvidenceSelection,
};
use crate::util::rss::peak_rss_bytes;
use crate::variant_io::{parse_snv_record, VariantLimits, VariantReader};
use anyhow::{anyhow, bail, Context, Result};
use clap::{Args, Subcommand, ValueEnum};
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet};
use std::fs::{self, File};
use std::io::Read;
use std::path::{Path, PathBuf};
use std::time::Duration;

#[derive(Subcommand, Debug)]
pub enum CohortAction {
    /// Bounded internal transport for validated cohort receipt recipes.
    #[command(hide = true)]
    Replay(ReplayArgs),
    /// Import verified portable datasets and publish one immutable snapshot (preview).
    Create(CreateArgs),
    /// Inspect snapshot and dataset metadata without decoding payload rows.
    Inspect(SnapshotArgs),
    /// Verify the current snapshot's complete stored evidence, without original sources.
    Verify(SnapshotArgs),
    /// Extract exact sample-by-candidate evidence from a saved snapshot.
    Extract(QueryArgs),
    /// Summarize observed, depth-eligible and ALT-supported samples for supplied SNVs.
    Summarize(QueryArgs),
    /// Compare explicitly ordered sample pairs; report exact right-minus-left read fractions.
    ComparePairs(ComparePairsArgs),
    /// Explicitly extend missing loci using local source mappings (preview).
    Extend(ExtendArgs),
}

#[derive(Args, Debug)]
pub struct ReplayArgs {
    #[arg(long)]
    request: PathBuf,
    #[arg(long)]
    memory_budget_bytes: Option<u64>,
    #[arg(long)]
    enforce: bool,
    #[arg(long, requires = "enforce")]
    require_os_limit: bool,
}

#[derive(clap::Parser, Debug)]
struct CohortInvocation {
    #[command(subcommand)]
    action: CohortAction,
}

#[derive(Args, Debug, Clone)]
pub struct ResourceArgs {
    /// Declared whole-process memory budget in MiB. Add --enforce to honor it.
    #[arg(long)]
    memory_budget_mb: Option<u64>,
    /// Exact-byte replay operand.
    #[arg(long, hide = true, conflicts_with = "memory_budget_mb")]
    memory_budget_bytes: Option<u64>,
    /// Cooperative admission and monitoring; requires --memory-budget-mb.
    #[arg(long)]
    enforce: bool,
    /// Require an existing Linux cgroup-v2 memory.max at or below the budget.
    #[arg(long, requires = "enforce")]
    require_os_limit: bool,
    /// Maximum serialized snapshot metadata bytes.
    #[arg(long, default_value_t = 8_388_608)]
    max_snapshot_bytes: usize,
    /// Maximum portable dataset descriptor or manifest bytes.
    #[arg(long, default_value_t = 33_554_432)]
    max_dataset_metadata_bytes: usize,
}
impl ResourceArgs {
    fn declared_budget(&self) -> Result<Option<u64>> {
        if let Some(bytes) = self.memory_budget_bytes {
            if bytes == 0 {
                bail!("memory budget must be positive");
            }
            return Ok(Some(bytes));
        }
        self.memory_budget_mb
            .map(|mb| {
                mb.checked_mul(1 << 20)
                    .filter(|bytes| *bytes > 0)
                    .ok_or_else(|| anyhow!("memory budget must be positive and fit uint64"))
            })
            .transpose()
    }
    fn enforcement(&self) -> EnforcementMode {
        if self.require_os_limit {
            EnforcementMode::RequireOsLimit
        } else if self.enforce {
            EnforcementMode::Cooperative
        } else {
            EnforcementMode::RecordOnly
        }
    }
    fn limits(&self, active_budget: Option<u64>) -> Result<CohortLimits> {
        if self.max_dataset_metadata_bytes == 0 || self.max_snapshot_bytes == 0 {
            bail!("metadata envelopes must be positive");
        }
        let mut limits = CohortLimits {
            memory_budget_bytes: active_budget,
            max_snapshot_bytes: self.max_snapshot_bytes,
            ..CohortLimits::default()
        };
        limits.dataset.memory_budget_bytes = active_budget;
        limits.dataset.max_manifest_bytes = self.max_dataset_metadata_bytes;
        limits.dataset.max_descriptor_bytes = self.max_dataset_metadata_bytes;
        Ok(limits)
    }
}

#[derive(Args, Debug)]
pub struct SnapshotArgs {
    /// Local cohort directory; copying this directory preserves saved-only usability.
    #[arg(
        long = "cohort",
        visible_alias = "cohort-root",
        required_unless_present = "cohort_snapshot_manifest",
        conflicts_with = "cohort_snapshot_manifest"
    )]
    cohort: Option<PathBuf>,
    /// File-addressed immutable snapshot operand used for relocated replay.
    #[arg(long, hide = true)]
    cohort_snapshot_manifest: Option<PathBuf>,
    /// Immutable 64-character snapshot content identity.
    #[arg(long)]
    snapshot: String,
    #[command(flatten)]
    resources: ResourceArgs,
}

#[derive(Args, Debug)]
pub struct CreateArgs {
    /// Destination cohort directory; imports and snapshots are immutable and deduplicated.
    #[arg(long = "cohort", visible_alias = "cohort-root")]
    cohort: PathBuf,
    /// TSV columns: id, manifest; optional group, subject, timepoint. Paths are table-relative.
    #[arg(long)]
    members: PathBuf,
    /// Existing immutable parent snapshot identity in this cohort.
    #[arg(long)]
    parent: Option<String>,
    /// Verify source metadata and show the import inventory without copying or publishing.
    #[arg(long)]
    plan: bool,
    /// Maximum members-table bytes; also bounded by member and leaf cardinality.
    #[arg(long, default_value_t = 8_388_608)]
    max_table_bytes: usize,
    #[command(flatten)]
    resources: ResourceArgs,
}

#[derive(ValueEnum, Debug, Clone, Copy)]
pub enum MissingMode {
    Strict,
    Partial,
}
#[derive(ValueEnum, Debug, Clone, Copy)]
pub enum CohortFormat {
    #[value(alias = "arrow")]
    ArrowIpc,
    Tsv,
}

#[derive(Args, Debug)]
pub struct CandidateArgs {
    #[command(flatten)]
    input: SnapshotArgs,
    /// Local candidate-SNV VCF, VCF.gz or BCF, with A/C/G/T REF and ALT alleles.
    #[arg(
        long,
        required_unless_present = "cohort_query",
        conflicts_with = "cohort_query"
    )]
    sites: Option<PathBuf>,
    /// Canonical inline dataset query used for replay.
    #[arg(long, hide = true)]
    cohort_query: Option<String>,
    /// Select a member ID; repeat. Omitted means all members in canonical ID order.
    #[arg(long = "member")]
    members: Vec<String>,
    /// JSON selected member IDs, including an explicit empty selection.
    #[arg(long, hide = true, conflicts_with = "members")]
    cohort_members: Option<String>,
    /// Content-matched verified replay input operands.
    #[arg(long = "cohort-artifact", hide = true)]
    artifacts: Vec<PathBuf>,
    /// Strict refuses absent loci; partial emits unmeasured cells with null metrics.
    #[arg(long, value_enum, default_value_t = MissingMode::Strict)]
    missing: MissingMode,
    /// Exact evidence groups, comma separated; depths and alleles are required.
    #[arg(long, value_delimiter = ',', default_value = "depths,alleles")]
    fields: Vec<String>,
    /// Technical depth screen; not a confidence estimate or genotype threshold.
    #[arg(long, default_value_t = DEFAULT_MIN_CALLABLE_DEPTH)]
    min_callable_depth: u64,
    /// Print requested samples, coverage, incompatibilities and resource reservations as JSON.
    #[arg(long)]
    plan: bool,
    /// Execution genomic window width; output batching and scientific results are unchanged.
    #[arg(long, alias = "max-microtile-bases", default_value_t = 16_384)]
    tile_bases: u32,
    /// Maximum input candidate records before normalization.
    #[arg(long, default_value_t = 1_000_000)]
    max_candidate_sites: usize,
    /// Cooperative maximum variant header/formatting envelope.
    #[arg(long, default_value_t = 8_388_608)]
    max_variant_header_bytes: usize,
    /// Cooperative maximum decoded variant record/formatting envelope.
    #[arg(long, default_value_t = 1_048_576)]
    max_variant_record_bytes: usize,
    /// Maximum encoded result receipt bytes.
    #[arg(long, default_value_t = 33_554_432)]
    max_receipt_bytes: usize,
}

#[derive(Args, Debug)]
pub struct QueryArgs {
    #[command(flatten)]
    candidate: CandidateArgs,
    /// Result path; required unless --plan. Receipt defaults to <output>.manifest.json.
    #[arg(short, long, required_unless_present = "plan")]
    output: Option<PathBuf>,
    /// Explicit result receipt path.
    #[arg(long, requires = "output")]
    manifest: Option<PathBuf>,
    /// Atomically replace existing output and receipt after successful computation.
    #[arg(long)]
    force: bool,
    #[arg(long, value_enum, default_value_t = CohortFormat::Tsv)]
    format: CohortFormat,
}

#[derive(Args, Debug)]
pub struct ComparePairsArgs {
    #[command(flatten)]
    query: QueryArgs,
    /// Explicit ordered TSV with exactly id, left, right columns. Never inferred from metadata.
    #[arg(long)]
    pairs: PathBuf,
    /// Maximum serialized pair-table bytes.
    #[arg(long, default_value_t = 8_388_608)]
    max_pair_table_bytes: usize,
    /// Maximum explicit pair records; also limited by the metadata envelope.
    #[arg(long, default_value_t = 65_536)]
    max_pairs: usize,
}

#[derive(Args, Debug)]
pub struct ExtendArgs {
    #[command(flatten)]
    candidate: CandidateArgs,
    /// Explicit source mapping table; no discovery of original BAM/CRAM inputs.
    #[arg(long)]
    sources: PathBuf,
    /// Existing scratch directory outside the cohort; temporary extraction data are cleaned up.
    #[arg(long, default_value = ".")]
    work_dir: PathBuf,
    /// Maximum source mapping table bytes.
    #[arg(long, default_value_t = 8_388_608)]
    max_table_bytes: usize,
}

/// Binary bridge with the established CLI exit categories.
pub fn run(action: CohortAction) -> Result<()> {
    if let Err(error) = run_inner(action) {
        let code = error
            .chain()
            .find_map(|cause| cause.downcast_ref::<CohortError>().map(cohort_exit_code))
            .unwrap_or(2);
        eprintln!("{error:#}");
        std::process::exit(code);
    }
    Ok(())
}
fn cohort_exit_code(error: &CohortError) -> i32 {
    match error {
        CohortError::Artifact(error) => error.exit_code(),
        CohortError::Dataset(error) => error.exit_code(),
        CohortError::Limit(_) => 3,
        CohortError::Corrupt(_) => 5,
        CohortError::Evidence(EvidenceError::Core(CoreError::Cancelled)) => 130,
        CohortError::Evidence(EvidenceError::Refused { .. }) => 3,
        CohortError::Evidence(
            EvidenceError::Core(CoreError::BudgetExceeded { .. }) | EvidenceError::RecordLimit(_),
        ) => 4,
        CohortError::Evidence(EvidenceError::CounterOverflow | EvidenceError::Analyzer(_)) => 1,
        _ => 2,
    }
}
fn run_inner(action: CohortAction) -> Result<()> {
    match action {
        CohortAction::Replay(args) => run_inner(read_replay_action(args)?),
        CohortAction::Create(args) => create(args),
        CohortAction::Inspect(args) => inspect(args, false),
        CohortAction::Verify(args) => inspect(args, true),
        CohortAction::Extract(args) => query(args, CohortOperation::Extract, None),
        CohortAction::Summarize(args) => query(args, CohortOperation::Summarize, None),
        CohortAction::ComparePairs(args) => query(
            args.query,
            CohortOperation::ComparePairs,
            Some(PairInput {
                path: args.pairs,
                max_table_bytes: args.max_pair_table_bytes,
                limits: PairLimits {
                    max_pairs: args.max_pairs,
                    ..PairLimits::default()
                },
            }),
        ),
        CohortAction::Extend(args) => extend(args),
    }
}

fn read_replay_action(args: ReplayArgs) -> Result<CohortAction> {
    const MAX_BYTES: usize = 32 << 20;
    const MAX_TOKENS: usize = 1_000_000;
    if args.memory_budget_bytes == Some(0) || (args.enforce && args.memory_budget_bytes.is_none()) {
        bail!("replay transport enforcement requires a positive byte budget");
    }
    let limits = CohortLimits {
        memory_budget_bytes: args.memory_budget_bytes.filter(|_| args.enforce),
        ..CohortLimits::default()
    };
    let path = fs::canonicalize(&args.request)?;
    let guard = InputSnapshot::capture([path.clone()])?;
    let length = usize::try_from(fs::metadata(&path)?.len())
        .context("replay request size exceeds address space")?;
    if length > MAX_BYTES {
        bail!("cohort replay request exceeds 32 MiB");
    }
    // The request file carries only already-resolved argument values. Start no
    // second managed scope: reserve expansion before parsing, then let the actual
    // extract/summarize runner own cancellation, inputs and publication.
    limits.admit((length as u64).saturating_mul(24).saturating_add(2 << 20))?;
    let mut bytes = Vec::with_capacity(length);
    File::open(&path)?
        .take((length as u64).saturating_add(1))
        .read_to_end(&mut bytes)?;
    guard.verify()?;
    if bytes.len() != length {
        bail!("cohort replay request changed while reading");
    }
    struct TokensVisitor;
    impl<'de> serde::de::Visitor<'de> for TokensVisitor {
        type Value = Vec<String>;
        fn expecting(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            formatter.write_str("a bounded array of cohort extract/summarize arguments")
        }
        fn visit_seq<A: serde::de::SeqAccess<'de>>(
            self,
            mut sequence: A,
        ) -> std::result::Result<Self::Value, A::Error> {
            let mut tokens = Vec::new();
            while let Some(token) = sequence.next_element::<String>()? {
                if tokens.len() >= MAX_TOKENS {
                    return Err(serde::de::Error::custom("too many replay tokens"));
                }
                if token.contains('\0') {
                    return Err(serde::de::Error::custom(
                        "replay arguments cannot contain NUL",
                    ));
                }
                tokens.push(token);
            }
            Ok(tokens)
        }
    }
    use clap::Parser;
    use serde::Deserializer;
    let mut reader = serde_json::Deserializer::from_slice(&bytes);
    let tokens = reader
        .deserialize_seq(TokensVisitor)
        .context("invalid cohort replay request")?;
    reader.end()?;
    if tokens.first().map(String::as_str) != Some("cohort")
        || !matches!(
            tokens.get(1).map(String::as_str),
            Some("extract" | "summarize" | "compare-pairs")
        )
    {
        bail!("cohort replay request must select cohort extract, summarize or compare-pairs");
    }
    let invocation = CohortInvocation::try_parse_from(tokens)?;
    let target = match &invocation.action {
        CohortAction::Extract(query) | CohortAction::Summarize(query) => {
            &query.candidate.input.resources
        }
        CohortAction::ComparePairs(args) => &args.query.candidate.input.resources,
        _ => bail!("nested or unrelated replay commands are forbidden"),
    };
    if target.declared_budget()? != args.memory_budget_bytes
        || target.enforce != args.enforce
        || target.require_os_limit != args.require_os_limit
    {
        bail!("replay transport budget and enforcement must exactly match the target command");
    }
    guard.verify()?;
    limits.admit(0)?;
    Ok(invocation.action)
}

fn with_scope(
    resources: &ResourceArgs,
    operation: impl FnOnce(CohortLimits) -> Result<()>,
) -> Result<()> {
    let scope = CancellationScope::start(CancellationToken::default())?;
    let _signals = SignalCancellationGuard::install(&scope)?;
    let declared = resources.declared_budget()?;
    let active = declared.filter(|_| resources.enforce);
    if resources.enforce && active.is_none() {
        bail!("--enforce requires --memory-budget-mb");
    }
    if resources.require_os_limit {
        let os_limit = crate::contract::detected_os_memory_limit_bytes().ok_or_else(|| {
            CohortError::Artifact(EvidenceArtifactError::OsEnforcementUnavailable(
                "require an existing Linux cgroup-v2 memory.max".into(),
            ))
        })?;
        if os_limit > active.ok_or_else(|| anyhow!("--require-os-limit requires --enforce"))? {
            return Err(
                CohortError::Artifact(EvidenceArtifactError::OsEnforcementUnavailable(
                    "cgroup memory.max exceeds declared budget".into(),
                ))
                .into(),
            );
        }
    }
    let _governor = active
        .map(|bytes| MemoryGovernor::start(bytes, Duration::from_millis(100), peak_rss_bytes))
        .transpose()?;
    checkpoint()
        .map_err(EvidenceError::from)
        .map_err(CohortError::from)?;
    let limits = resources.limits(active)?;
    limits.admit(0)?;
    operation(limits)?;
    limits.admit(0)?;
    Ok(())
}

fn create(args: CreateArgs) -> Result<()> {
    with_scope(&args.resources, |limits| {
        let members = read_members(&args.members, args.max_table_bytes, limits)?;
        if args.plan {
            if let Some(parent) = &args.parent {
                open_snapshot(&args.cohort, parent, limits)?.verify_unchanged()?;
            }
            let mut inventories = Vec::new();
            for member in &members {
                let mut objects = Vec::new();
                for manifest in &member.manifests {
                    let dataset = VerifiedEvidenceDataset::open(manifest, limits.dataset_limits())?;
                    objects.push(json!({
                        "manifest": manifest,
                        "object_id": dataset.source_hashes()[0].blake3,
                        "fields": dataset.fields().names(),
                        "stored_loci": dataset.descriptor().selection.intervals.iter().map(|interval| (interval.end - interval.start) as u64).sum::<u64>(),
                        "partition_count": dataset.descriptor().partitions.len(),
                    }));
                    dataset.verify_unchanged()?;
                }
                inventories.push(json!({"metadata":member.metadata,"datasets":objects}));
            }
            print_json(
                &json!({"status":"planned","operation":"create","parent":args.parent,"members":inventories,"payloads_verified":false,"published":false}),
            )?;
        } else {
            let snapshot = create_snapshot(&args.cohort, &members, args.parent.as_deref(), limits)?;
            print_json(
                &json!({"status":"created","cohort":snapshot.root,"snapshot_id":snapshot.id,"members":snapshot.descriptor.members.len(),"parent":snapshot.descriptor.parent}),
            )?;
        }
        Ok(())
    })
}
fn inspect(args: SnapshotArgs, verify: bool) -> Result<()> {
    let root = snapshot_root(&args)?;
    with_scope(&args.resources, |limits| {
        let snapshot = if verify {
            verify_snapshot(&root, &args.snapshot, limits)?
        } else {
            open_snapshot(&root, &args.snapshot, limits)?
        };
        snapshot.verify_unchanged()?;
        print_json(&json!({
            "status": if verify {"verified"} else {"inspected"},
            "snapshot_id":snapshot.id,"descriptor":snapshot.descriptor,
            "current_payloads_verified":verify,"ancestor_payloads_verified":false,
            "original_sources_opened":false,
        }))
    })
}
fn query(args: QueryArgs, operation: CohortOperation, pairs: Option<PairInput>) -> Result<()> {
    let planning = args.candidate.plan;
    let mut spec = artifact_spec(args, operation)?;
    spec.pairs = pairs;
    if planning {
        let plan = plan_cohort_artifact(spec)?;
        let mut report = plan_json(&plan.query);
        if let Some(pairs) = &plan.pairs {
            let output_rows = plan
                .candidate_rows
                .checked_mul(pairs.pairs.len() as u64)
                .ok_or_else(|| anyhow!("paired row count overflow"))?;
            report["source_member_candidate_rows"] = json!(plan.output_rows);
            report["output_rows"] = json!(output_rows);
            report["paired_candidate_rows"] = json!(output_rows);
            report["pair_direction"] = json!("right-minus-left");
            report["pairs"] = json!(pairs
                .pairs
                .iter()
                .map(|pair| json!({
                    "id":pair.id,"left":plan.members[pair.left_member_index].member_id,
                    "right":plan.members[pair.right_member_index].member_id
                }))
                .collect::<Vec<_>>());
            report["resources"]["pair_metadata_bytes"] = json!(pairs.metadata_bytes);
        }
        print_json(&report)?;
    } else {
        let outcome = run_cohort_artifact(spec)?;
        let mut report = json!({"status":"completed","output":outcome.output,"manifest":outcome.manifest,"claim_blake3":outcome.claim_hash,
            "sample_candidate_rows":outcome.stats.emitted_rows,"observed_rows":outcome.stats.observed_rows,
            "unmeasured_rows":outcome.stats.unmeasured_rows,"original_alignment_records_decoded":0});
        if let Some(rows) = outcome.paired_candidate_rows {
            report["paired_candidate_rows"] = json!(rows);
            report["pair_direction"] = json!("right-minus-left");
        }
        print_json(&report)?;
    }
    Ok(())
}
fn snapshot_root(args: &SnapshotArgs) -> Result<PathBuf> {
    if let Some(root) = &args.cohort {
        return Ok(root.clone());
    }
    let file = fs::canonicalize(
        args.cohort_snapshot_manifest
            .as_ref()
            .ok_or_else(|| anyhow!("a cohort root or snapshot manifest is required"))?,
    )?;
    let directory = file
        .parent()
        .ok_or_else(|| anyhow!("invalid snapshot layout"))?;
    let snapshots = directory
        .parent()
        .ok_or_else(|| anyhow!("invalid snapshot layout"))?;
    if file.file_name().and_then(|name| name.to_str()) != Some("snapshot.json")
        || directory.file_name().and_then(|name| name.to_str()) != Some(args.snapshot.as_str())
        || snapshots.file_name().and_then(|name| name.to_str()) != Some("snapshots")
    {
        bail!("snapshot manifest must be <cohort>/snapshots/<snapshot ID>/snapshot.json");
    }
    Ok(snapshots
        .parent()
        .ok_or_else(|| anyhow!("invalid snapshot root"))?
        .to_owned())
}
fn selected_members(args: &CandidateArgs) -> Result<Option<Vec<String>>> {
    let active_budget = args
        .input
        .resources
        .declared_budget()?
        .filter(|_| args.input.resources.enforce);
    let limits = args.input.resources.limits(active_budget)?;
    if let Some(text) = &args.cohort_members {
        if text.len() > 32 << 20 {
            bail!("selected member JSON exceeds envelope");
        }
        // Raw CLI values already exist; admit all expansion before deserializing.
        limits.admit(
            (text.len() as u64)
                .saturating_mul(24)
                .saturating_add(128 << 10),
        )?;
        struct MembersVisitor(usize);
        impl<'de> serde::de::Visitor<'de> for MembersVisitor {
            type Value = Vec<String>;
            fn expecting(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                formatter.write_str("a bounded JSON array of valid member IDs")
            }
            fn visit_seq<A: serde::de::SeqAccess<'de>>(
                self,
                mut sequence: A,
            ) -> std::result::Result<Self::Value, A::Error> {
                let mut members = Vec::new();
                while let Some(id) = sequence.next_element::<String>()? {
                    if members.len() >= self.0 {
                        return Err(serde::de::Error::custom("too many selected members"));
                    }
                    validate_member_id(&id).map_err(serde::de::Error::custom)?;
                    members.push(id);
                }
                Ok(members)
            }
        }
        use serde::Deserializer;
        let mut reader = serde_json::Deserializer::from_str(text);
        let members = reader
            .deserialize_seq(MembersVisitor(limits.max_members))
            .context("invalid selected-member JSON")?;
        reader.end()?;
        return Ok(Some(members));
    }
    if args.members.len() > limits.max_members {
        bail!("too many selected members");
    }
    for member in &args.members {
        validate_member_id(member)?;
    }
    limits.admit(
        (args.members.len() as u64)
            .saturating_mul(512)
            .saturating_add(128 << 10),
    )?;
    Ok((!args.members.is_empty()).then(|| args.members.clone()))
}
fn validate_member_id(id: &str) -> Result<()> {
    if id.is_empty() || id.len() > 256 || id.chars().any(char::is_control) {
        bail!("member IDs must be nonempty, at most 256 UTF-8 bytes, without control characters");
    }
    Ok(())
}

fn artifact_spec(args: QueryArgs, operation: CohortOperation) -> Result<CohortArtifactSpec> {
    let QueryArgs {
        candidate: args,
        output,
        manifest,
        force,
        format,
    } = args;
    let root = snapshot_root(&args.input)?;
    let budget = args.input.resources.declared_budget()?;
    let enforcement = args.input.resources.enforcement();
    let fields = parse_fields(&args.fields)?;
    let mut query = CohortQuery::new(EvidenceSelection::Sites(Vec::new()));
    query.fields = fields;
    query.requirements.fields = required_fields();
    query.member_ids = selected_members(&args)?;
    if let Some(inline) = &args.cohort_query {
        let parsed = crate::evidence::parse_artifact_query(inline)?;
        if parsed.fields != fields {
            bail!("inline query fields differ from --fields");
        }
        query.selection = parsed.selection;
    }
    query.missing_policy = match args.missing {
        MissingMode::Strict => MissingPolicy::Strict,
        MissingMode::Partial => MissingPolicy::Partial,
    };
    Ok(CohortArtifactSpec {
        root,
        snapshot_id: args.input.snapshot,
        query,
        sites: args.sites.map(|path| {
            (
                path,
                VariantLimits {
                    max_header_bytes: args.max_variant_header_bytes,
                    max_record_bytes: args.max_variant_record_bytes,
                },
            )
        }),
        pairs: None,
        artifacts: args.artifacts,
        execution: EvidenceExecution {
            memory_budget_bytes: budget,
            max_microtile_bases: args.tile_bases,
            ..EvidenceExecution::default()
        },
        limits: CohortQueryLimits {
            cohort: args.input.resources.limits(budget)?,
            max_input_sites: args.max_candidate_sites,
            ..CohortQueryLimits::default()
        },
        operation,
        format: match format {
            CohortFormat::ArrowIpc => CohortOutputFormat::Arrow,
            CohortFormat::Tsv => CohortOutputFormat::Tsv,
        },
        min_callable_depth: args.min_callable_depth,
        output: output.unwrap_or_default(),
        manifest,
        output_policy: if force {
            OutputPolicy::ReplaceAtomic
        } else {
            OutputPolicy::CreateNewAtomic
        },
        enforcement,
        cancellation: None,
        handle_signals: true,
        max_receipt_bytes: args.max_receipt_bytes,
    })
}
fn parse_fields(names: &[String]) -> Result<EvidenceFields> {
    let mut fields = EvidenceFields::from_bits(0)?;
    for name in names {
        fields = fields.union(match name.as_str() {
            "all" if names.len() == 1 => EvidenceFields::ALL,
            "all-supported" if names.len() == 1 => EvidenceFields::ALL_SUPPORTED,
            "depths" => EvidenceFields::DEPTHS,
            "alleles" => EvidenceFields::ALLELES,
            "strands" => EvidenceFields::STRANDS,
            "quality-sums" => EvidenceFields::QUALITY_SUMS,
            "quality-histograms" => EvidenceFields::QUALITY_HISTOGRAMS,
            "read-position" => EvidenceFields::READ_POSITION,
            "allele-quality" => EvidenceFields::ALLELE_QUALITY,
            _ => bail!("invalid evidence field group {name:?}"),
        });
    }
    if !fields.contains(required_fields()) {
        bail!("cohort candidates require depths and alleles");
    }
    Ok(fields)
}
fn plan_json(plan: &CohortQueryPlan) -> Value {
    json!({
        "status":if plan.issue_count == 0 {"ready"} else {"blocked"},
        "snapshot_id":plan.snapshot_id,"comparison_blake3":plan.comparison_blake3,
        "requested_loci":plan.requested_loci,"candidate_rows":plan.candidate_rows,"output_rows":plan.output_rows,
        "fields":plan.fields.names(),"missing_policy":plan.missing_policy,
        "issues":plan.issues,"issue_count":plan.issue_count,
        "members":plan.members.iter().map(|member| json!({
            "id":member.member_id,"covered_loci":member.covered_loci,"missing_loci":member.missing_loci,
            "covered_rows":member.covered_rows,"missing_rows":member.missing_rows,
            "datasets":member.leaves.iter().map(|leaf| json!({"object_id":leaf.object_id,"covered_loci":leaf.covered_loci,"covered_rows":leaf.covered_rows})).collect::<Vec<_>>()
        })).collect::<Vec<_>>(),
        "resources":{
            "memory_budget_bytes":plan.reservations.memory_budget_bytes,"baseline_rss_bytes":plan.reservations.baseline_rss_bytes,
            "query_bytes":plan.reservations.query_bytes,"plan_bytes":plan.reservations.plan_bytes,
            "consumer_bytes":plan.reservations.consumer_bytes,"lineage_bytes":plan.reservations.lineage_bytes,"consumer_bound_known":plan.reservations.consumer_bound_known,
            "max_reader_metadata_bytes":plan.reservations.max_reader_metadata_bytes,
            "max_source_decoder_bytes":plan.reservations.max_source_decoder_bytes,
            "max_projection_bytes":plan.reservations.max_projection_bytes,
            "predicted_peak_rss_bytes":plan.reservations.predicted_peak_rss_bytes,
        },
        "original_alignment_records_decoded":0,"published":false,
    })
}
fn print_json(value: &Value) -> Result<()> {
    use std::io::Write;
    let stdout = std::io::stdout();
    let mut out = stdout.lock();
    serde_json::to_writer_pretty(&mut out, value)?;
    writeln!(out)?;
    Ok(())
}

/// The table is read once into immutable request values. Its bytes are guarded
/// during parsing; snapshot identity commits the resulting metadata and objects.
fn read_members(path: &Path, max_bytes: usize, limits: CohortLimits) -> Result<Vec<ImportMember>> {
    let (path, text) = read_table(path, max_bytes, limits)?;
    let mut lines = text.lines();
    let header: Vec<_> = lines
        .next()
        .ok_or_else(|| anyhow!("members table is empty"))?
        .split('\t')
        .collect();
    let allowed = ["id", "manifest", "group", "subject", "timepoint"];
    let mut seen = BTreeSet::new();
    for name in &header {
        if !allowed.contains(name) || !seen.insert(*name) {
            bail!("unknown or duplicate members-table column {name:?}");
        }
    }
    let column = |name: &str| header.iter().position(|value| *value == name);
    let id_col = column("id").ok_or_else(|| anyhow!("members table requires id column"))?;
    let manifest_col =
        column("manifest").ok_or_else(|| anyhow!("members table requires manifest column"))?;
    let parent = path.parent().unwrap();
    let mut members = BTreeMap::<String, ImportMember>::new();
    let mut leaf_count = 0usize;
    let mut supplied = BTreeSet::new();
    for (line_index, line) in lines.enumerate() {
        limits.admit(0)?;
        if line.is_empty() {
            bail!("blank members-table row {}", line_index + 2);
        }
        let values: Vec<_> = line.split('\t').collect();
        if values.len() != header.len() {
            bail!(
                "members-table row {} has the wrong column count",
                line_index + 2
            );
        }
        let optional = |name| {
            column(name)
                .map(|index| values[index])
                .filter(|value| !value.is_empty())
                .map(str::to_owned)
        };
        let metadata = MemberMetadata {
            id: values[id_col].into(),
            group: optional("group"),
            subject: optional("subject"),
            timepoint: optional("timepoint"),
        };
        metadata.validate()?;
        let manifest = values[manifest_col];
        if manifest.is_empty() || manifest.len() > 4096 || manifest.chars().any(char::is_control) {
            bail!("members-table manifest path is empty, too long or contains a control character");
        }
        leaf_count = leaf_count
            .checked_add(1)
            .ok_or_else(|| anyhow!("table leaf count overflow"))?;
        if leaf_count > limits.max_leaf_references {
            bail!("members table exceeds leaf-reference envelope");
        }
        let manifest = parent.join(manifest);
        if !supplied.insert((metadata.id.clone(), manifest.clone())) {
            bail!("duplicate member/manifest row for {:?}", metadata.id);
        }
        match members.entry(metadata.id.clone()) {
            std::collections::btree_map::Entry::Vacant(entry) => {
                entry.insert(ImportMember {
                    metadata,
                    manifests: vec![manifest],
                });
            }
            std::collections::btree_map::Entry::Occupied(mut entry) => {
                if entry.get().metadata != metadata {
                    bail!("repeated member {:?} has conflicting metadata", metadata.id);
                }
                entry.get_mut().manifests.push(manifest);
            }
        }
        if members.len() > limits.max_members {
            bail!("members table exceeds member envelope");
        }
    }
    if members.is_empty() {
        bail!("members table must contain at least one member");
    }
    Ok(members.into_values().collect())
}
fn read_table(path: &Path, max_bytes: usize, limits: CohortLimits) -> Result<(PathBuf, String)> {
    if max_bytes == 0 {
        bail!("table byte envelope must be positive");
    }
    let path =
        fs::canonicalize(path).with_context(|| format!("cannot open table {}", path.display()))?;
    let guard = InputSnapshot::capture([path.clone()])?;
    let metadata = fs::metadata(&path)?;
    let length =
        usize::try_from(metadata.len()).map_err(|_| anyhow!("table size exceeds address space"))?;
    if length > max_bytes {
        bail!("table exceeds --max-table-bytes");
    }
    limits.admit((length as u64).saturating_mul(24).saturating_add(128 << 10))?;
    let mut bytes = Vec::with_capacity(length);
    File::open(&path)?
        .take((length as u64).saturating_add(1))
        .read_to_end(&mut bytes)?;
    guard.verify()?;
    if bytes.len() != length {
        bail!("table changed while reading");
    }
    let text = String::from_utf8(bytes).context("table must be UTF-8")?;
    if text.lines().any(|line| line.len() > 8192) {
        bail!("table row exceeds 8192-byte envelope");
    }
    Ok((path, text))
}

fn extend(args: ExtendArgs) -> Result<()> {
    use crate::cohort::extend::{extend_snapshot_with_inputs, MemberSourceMapping};
    let candidate = &args.candidate;
    let root = snapshot_root(&candidate.input)?;
    if !candidate.artifacts.is_empty() {
        bail!("cohort extension has no replay recipe; artifact operands are not accepted");
    }
    with_scope(&candidate.input.resources, |limits| {
        let snapshot = open_snapshot(&root, &candidate.input.snapshot, limits)?;
        let execution = EvidenceExecution {
            memory_budget_bytes: limits.effective_budget(),
            max_microtile_bases: candidate.tile_bases,
            ..EvidenceExecution::default()
        };
        let query_limits = CohortQueryLimits {
            cohort: limits,
            max_input_sites: candidate.max_candidate_sites,
            ..CohortQueryLimits::default()
        };
        let fields = parse_fields(&candidate.fields)?;
        let mut query = CohortQuery::new(EvidenceSelection::Sites(Vec::new()));
        query.fields = fields;
        query.member_ids = selected_members(candidate)?;
        query.missing_policy = MissingPolicy::Partial;
        let mut guarded_paths = vec![fs::canonicalize(&args.sources)?];
        if let Some(path) = &candidate.sites {
            guarded_paths.push(fs::canonicalize(path)?);
        }
        let guard = InputSnapshot::capture(guarded_paths)?;
        let dictionary = plan_query(&snapshot, &query, &execution, query_limits)?;
        if let Some(text) = &candidate.cohort_query {
            let parsed = crate::evidence::parse_artifact_query(text)?;
            if parsed.fields != fields {
                bail!("inline query fields differ from --fields");
            }
            query.selection = parsed.selection;
        } else {
            let path = candidate
                .sites
                .as_ref()
                .ok_or_else(|| anyhow!("extension requires candidate SNVs"))?;
            let parser_bytes = (candidate.max_variant_header_bytes as u64)
                .saturating_add(candidate.max_variant_record_bytes as u64)
                .saturating_mul(4)
                .saturating_add(128 << 10);
            limits.admit(parser_bytes)?;
            let mut reader = VariantReader::open(
                path,
                VariantLimits {
                    max_header_bytes: candidate.max_variant_header_bytes,
                    max_record_bytes: candidate.max_variant_record_bytes,
                },
            )?;
            let mut sites = Vec::new();
            while let Some(record) = reader.read()? {
                if sites.len() >= query_limits.max_input_sites {
                    bail!("candidate file exceeds record envelope");
                }
                let bytes = (sites.len() as u64 + 1).saturating_mul(4096);
                if bytes > query_limits.max_query_bytes {
                    bail!("candidate file exceeds query byte envelope");
                }
                limits.admit(parser_bytes.saturating_add(bytes))?;
                sites.push(parse_snv_record(record, &dictionary.contigs)?);
            }
            query.selection = EvidenceSelection::Sites(sites);
        }
        drop(dictionary);
        let mappings: Vec<MemberSourceMapping> =
            read_sources(&args.sources, args.max_table_bytes, limits)?;
        guard.verify()?;
        if candidate.plan {
            let plan = plan_query(&snapshot, &query, &execution, query_limits)?;
            let affected: BTreeSet<_> = plan
                .members
                .iter()
                .filter(|member| member.missing_loci.is_some_and(|count| count > 0))
                .map(|member| member.member_id.as_str())
                .collect();
            let supplied: BTreeSet<_> = mappings
                .iter()
                .map(|mapping| mapping.member_id.as_str())
                .collect();
            let matches = affected == supplied;
            let mut report = plan_json(&plan);
            report["operation"] = json!("extend");
            report["affected_members"] = json!(affected);
            report["source_mapping_members"] = json!(supplied);
            report["source_mapping_members_match"] = json!(matches);
            report["source_identities_verified"] = json!(false);
            report["raw_sources_opened"] = json!(false);
            if !matches {
                report["status"] = json!("blocked");
            }
            guard.verify()?;
            snapshot.verify_unchanged()?;
            print_json(&report)?;
        } else {
            let outcome = extend_snapshot_with_inputs(
                &snapshot,
                &query,
                &mappings,
                &args.work_dir,
                &execution,
                query_limits,
                Some(&guard),
            )?;
            print_json(&json!({
                "status":"completed","snapshot_id":outcome.snapshot.id,"parent":outcome.snapshot.descriptor.parent,
                "changed":outcome.changed,"verified_existing_bytes":outcome.verified_existing_bytes,"publication_ms":outcome.publication_ms,
                "members":outcome.members.iter().map(|member| json!({
                    "id":member.member_id,"retained_loci":member.retained_loci,"computed_loci":member.computed_loci,
                    "native_record_visits":member.native_record_visits,"full_cram_validation_records":member.full_cram_validation_records,
                    "source_hashed_bytes":member.source_hashed_bytes,"source_verification_ms":member.source_verification_ms,
                    "extraction_ms":member.extraction_ms,
                })).collect::<Vec<_>>(),
            }))?;
        }
        Ok(())
    })
}
fn read_sources(
    path: &Path,
    max_bytes: usize,
    limits: CohortLimits,
) -> Result<Vec<crate::cohort::extend::MemberSourceMapping>> {
    use crate::cohort::extend::MemberSourceMapping;
    let (path, text) = read_table(path, max_bytes, limits)?;
    let mut lines = text.lines();
    if lines.next() != Some("id\trole\tpath") {
        bail!("source table header must be id, role, path (tab-separated, in that order)");
    }
    let mut mappings = BTreeMap::<String, MemberSourceMapping>::new();
    let mut count = 0usize;
    for line in lines {
        let values: Vec<_> = line.split('\t').collect();
        if values.len() != 3
            || values
                .iter()
                .any(|value| value.is_empty() || value.chars().any(char::is_control))
        {
            bail!("source table rows require id, role and path");
        }
        MemberMetadata {
            id: values[0].into(),
            group: None,
            subject: None,
            timepoint: None,
        }
        .validate()?;
        if values[1].len() > 256 || values[2].len() > 4096 {
            bail!("source role or path exceeds envelope");
        }
        count += 1;
        if count > limits.max_leaf_references || mappings.len() > limits.max_members {
            bail!("source mapping inventory exceeds envelope");
        }
        if !mappings.contains_key(values[0]) && mappings.len() >= limits.max_members {
            bail!("source mapping exceeds member envelope");
        }
        let entry = mappings
            .entry(values[0].into())
            .or_insert_with(|| MemberSourceMapping {
                member_id: values[0].into(),
                roles: BTreeMap::new(),
            });
        if entry
            .roles
            .insert(values[1].into(), path.parent().unwrap().join(values[2]))
            .is_some()
        {
            bail!("duplicate source role for member {:?}", values[0]);
        }
    }
    Ok(mappings.into_values().collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::Parser;
    #[derive(Parser)]
    struct TestCli {
        #[command(subcommand)]
        action: CohortAction,
    }
    fn candidate(extra: &[&str]) -> CandidateArgs {
        let mut args = vec![
            "cohort",
            "extract",
            "--cohort",
            "unused",
            "--snapshot",
            "unused",
            "--sites",
            "unused",
            "--plan",
        ];
        args.extend_from_slice(extra);
        match TestCli::try_parse_from(args).unwrap().action {
            CohortAction::Extract(query) => query.candidate,
            _ => unreachable!(),
        }
    }
    #[test]
    fn member_selection_is_admitted_before_bounded_json_expansion() {
        let empty = candidate(&["--cohort-members", "[]"]);
        assert_eq!(selected_members(&empty).unwrap(), Some(vec![]));
        let over_count = format!(
            "[{}]",
            std::iter::repeat_n("\"a\"", 65_537)
                .collect::<Vec<_>>()
                .join(",")
        );
        assert!(selected_members(&candidate(&["--cohort-members", &over_count])).is_err());
        for malformed in ["[\"\"]", "[\"line\\nbreak\"]", "[1]", "[] garbage"] {
            assert!(selected_members(&candidate(&["--cohort-members", malformed])).is_err());
        }
        let huge_id = serde_json::to_string(&vec!["a".repeat(257)]).unwrap();
        assert!(selected_members(&candidate(&["--cohort-members", &huge_id])).is_err());
        let under_budget = candidate(&[
            "--cohort-members",
            "[\"A\"]",
            "--memory-budget-bytes",
            "1",
            "--enforce",
        ]);
        let error = selected_members(&under_budget).unwrap_err();
        assert!(matches!(
            error.downcast_ref::<CohortError>(),
            Some(CohortError::Limit(_))
        ));
    }
    #[test]
    fn clap_families_and_projection_requirements_are_consistent() {
        use clap::CommandFactory;
        TestCli::command().debug_assert();
        assert!(parse_fields(&["depths".into()]).is_err());
        assert_eq!(
            parse_fields(&["depths".into(), "alleles".into()]).unwrap(),
            required_fields()
        );
        assert!(parse_fields(&["made-up".into()]).is_err());
        assert!(
            TestCli::try_parse_from([
                "cohort",
                "extract",
                "--cohort",
                "unused",
                "--snapshot",
                "unused",
                "--sites",
                "unused"
            ])
            .is_err(),
            "materialization requires output"
        );
        assert!(
            TestCli::try_parse_from([
                "cohort",
                "extend",
                "--cohort",
                "unused",
                "--snapshot",
                "unused",
                "--sites",
                "unused",
                "--sources",
                "mapping.tsv"
            ])
            .is_ok(),
            "extension publishes a snapshot, not a result file"
        );
    }
}
