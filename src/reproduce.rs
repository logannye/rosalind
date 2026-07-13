//! `reproduce` — re-derive a recorded result and compare it byte-for-byte. Standalone:
//! filesystem + blake3 + subprocess only (no htslib). The verdict is over OUTPUT bytes;
//! code / inputs / resource are diagnostic context. v1 supports the deterministic text
//! outputs Rosalind emits (VCF, TSV); BAM/bgzf is reported INCONCLUSIVE, never a false
//! DIVERGED.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use anyhow::{anyhow, Context, Result};

use crate::provenance::{blake3_file, CommandCapture, FileHash, RunManifest};

/// The output-byte verdict. Exit codes are stable: a stranger's CI can branch on them.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Verdict {
    /// Every recorded output re-derived byte-identically.
    Reproduced,
    /// At least one output differs from what the receipt recorded.
    Diverged,
    /// Could not run the comparison (pre-schema-5, input not located, unsupported output).
    Inconclusive,
}

impl Verdict {
    /// The process exit code for this verdict (0 / 6 / 7).
    pub fn exit_code(self) -> i32 {
        match self {
            Verdict::Reproduced => 0,
            Verdict::Diverged => 6,
            Verdict::Inconclusive => 7,
        }
    }
    /// The human-readable label.
    pub fn label(self) -> &'static str {
        match self {
            Verdict::Reproduced => "REPRODUCED",
            Verdict::Diverged => "DIVERGED",
            Verdict::Inconclusive => "INCONCLUSIVE",
        }
    }
}

/// The result of comparing produced outputs against the recorded ones.
#[derive(Debug)]
pub struct Outcome {
    /// REPRODUCED iff every recorded output matched; DIVERGED otherwise.
    pub verdict: Verdict,
    /// One human-readable line per mismatch (empty when REPRODUCED).
    pub diffs: Vec<String>,
}

/// Pure: compare produced `(role, hash)` pairs against the recorded outputs. REPRODUCED
/// iff every recorded output has a positional produced match; DIVERGED otherwise, naming
/// each mismatch.
pub fn classify_outputs(recorded: &[FileHash], produced: &[(String, String)]) -> Outcome {
    let mut diffs = Vec::new();
    for (i, rec) in recorded.iter().enumerate() {
        match produced.get(i) {
            Some((_, got)) if *got == rec.blake3 => {}
            Some((role, got)) => {
                diffs.push(format!("output {role}: recorded {} got {got}", rec.blake3))
            }
            None => diffs.push(format!(
                "output {i}: recorded {} but not produced",
                rec.blake3
            )),
        }
    }
    let verdict = if diffs.is_empty() {
        Verdict::Reproduced
    } else {
        Verdict::Diverged
    };
    Outcome { verdict, diffs }
}

/// Per-output comparison detail (carried for the reproduction certificate).
#[derive(Debug, Clone)]
pub struct OutputCmp {
    /// Output role label (e.g. `output[0]`).
    pub role: String,
    /// The hash the receipt recorded.
    pub recorded: String,
    /// The hash the re-derivation produced.
    pub observed: String,
    /// Whether they match.
    pub matched: bool,
}

/// Machine-local resource observation of the re-run (explicitly NOT cross-machine).
#[derive(Debug, Clone, Default)]
pub struct ResourceHere {
    /// Realized peak RSS of the re-run on this machine, if its receipt was readable.
    pub peak_rss_bytes: Option<u64>,
    /// The budget the original declared (MiB), if any.
    pub declared_budget_mb: Option<u64>,
    /// Whether the re-run's peak fit the declared budget here.
    pub fit: Option<bool>,
}

/// The full result of a `reproduce` run: the printable report, the exit code, and the
/// structured data the reproduction certificate is minted from.
#[derive(Debug)]
pub struct ReproReport {
    /// `REPRODUCED` | `DIVERGED` | `INCONCLUSIVE` | `TAMPERED`.
    pub verdict_label: String,
    /// Process exit code: 0 | 6 | 7 | 5 (tamper).
    pub exit_code: i32,
    /// Aligned, human-readable report lines.
    pub lines: Vec<String>,
    /// The parent receipt's `content_hash()` — the cross-machine id the certificate chains to.
    pub parent_claim: String,
    /// The parent receipt's subcommand.
    pub parent_subcommand: String,
    /// Per-output comparison (populated for REPRODUCED / DIVERGED).
    pub outputs: Vec<OutputCmp>,
    /// Machine-local resource observation of the re-run.
    pub resource_here: ResourceHere,
    /// Build identity recorded by the original receipt.
    pub original_code: BTreeMap<String, String>,
    /// The build-identity of the binary that performed this reproduction.
    pub reproducer_code: BTreeMap<String, String>,
    /// Executable explicitly selected for the replay (or the current executable).
    pub execution_binary: PathBuf,
    /// True only for REPRODUCED / DIVERGED (an actual byte comparison happened).
    pub compared: bool,
}

impl ReproReport {
    /// Stable machine-readable summary for `rosalind reproduce --json`.
    pub fn to_json(&self) -> String {
        let outputs = self
            .outputs
            .iter()
            .map(|o| {
                format!(
                    "{{\"role\":\"{}\",\"recorded_blake3\":\"{}\",\"observed_blake3\":\"{}\",\"matched\":{}}}",
                    json_escape(&o.role),
                    json_escape(&o.recorded),
                    json_escape(&o.observed),
                    o.matched
                )
            })
            .collect::<Vec<_>>()
            .join(",");
        let encode_identity = |identity: &BTreeMap<String, String>| {
            identity
                .iter()
                .map(|(key, value)| format!("\"{}\":\"{}\"", json_escape(key), json_escape(value)))
                .collect::<Vec<_>>()
                .join(",")
        };
        let number = |value: Option<u64>| {
            value.map_or_else(|| "null".to_string(), |value| value.to_string())
        };
        let boolean = |value: Option<bool>| match value {
            Some(true) => "true",
            Some(false) => "false",
            None => "null",
        };
        format!(
            "{{\"schema\":1,\"verdict\":\"{}\",\"exit_code\":{},\"compared\":{},\"parent_claim\":\"{}\",\"parent_subcommand\":\"{}\",\"execution_binary\":\"{}\",\"code_identity_matches\":{},\"original_code\":{{{}}},\"reproducer_code\":{{{}}},\"resource_here\":{{\"peak_rss_bytes\":{},\"declared_budget_mb\":{},\"fit\":{}}},\"outputs\":[{}]}}",
            json_escape(&self.verdict_label),
            self.exit_code,
            self.compared,
            json_escape(&self.parent_claim),
            json_escape(&self.parent_subcommand),
            json_escape(&self.execution_binary.display().to_string()),
            self.original_code == self.reproducer_code,
            encode_identity(&self.original_code),
            encode_identity(&self.reproducer_code),
            number(self.resource_here.peak_rss_bytes),
            number(self.resource_here.declared_budget_mb),
            boolean(self.resource_here.fit),
            outputs
        )
    }
}

/// Whether replay uses the trusted upstream command surface or an explicitly supplied analyzer.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProducerKind {
    /// A receipt produced by the upstream Rosalind binary (legacy receipts included).
    Rosalind,
    /// A third-party analyzer executable explicitly selected by the caller.
    ExternalAnalyzer,
}

/// One input located by its recorded content hash.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LocatedArtifact {
    /// Stable positional role such as `input[0]`.
    pub role: String,
    /// Recorded BLAKE3 digest.
    pub blake3: String,
    /// Local path whose bytes match the digest.
    pub path: PathBuf,
}

/// One output path reserved inside the reproduction working directory.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PlannedOutput {
    /// Stable positional role such as `output[0]`.
    pub role: String,
    /// Digest expected after replay.
    pub recorded_blake3: String,
    /// Fresh local destination used only for the re-run.
    pub path: PathBuf,
}

/// A validated receipt-driven execution plan. Construct only through
/// [`plan_reproduction`]; every content marker and executable policy has already
/// been checked when this value is returned.
#[derive(Debug, Clone)]
pub struct ReproductionPlan {
    /// Executable selected by policy, never by an untrusted recorded path.
    pub binary: PathBuf,
    /// Fully substituted argv passed directly to the executable (never a shell).
    pub argv: Vec<String>,
    /// Fresh isolated working directory.
    pub work_dir: PathBuf,
    /// Upstream or explicitly selected external producer.
    pub producer_kind: ProducerKind,
    /// Content-located inputs.
    pub inputs: Vec<LocatedArtifact>,
    /// Temporary output destinations.
    pub outputs: Vec<PlannedOutput>,
    manifest: RunManifest,
    explicit_binary: bool,
}

impl ReproductionPlan {
    /// Stable JSON suitable for `rosalind reproduce --dry-run --json`.
    pub fn to_json(&self) -> String {
        let argv = self
            .argv
            .iter()
            .map(|value| format!("\"{}\"", json_escape(value)))
            .collect::<Vec<_>>()
            .join(",");
        format!(
            "{{\"schema\":1,\"validated\":true,\"binary\":\"{}\",\"producer_kind\":\"{}\",\"work_dir\":\"{}\",\"argv\":[{}],\"inputs\":{},\"outputs\":{}}}",
            json_escape(&self.binary.display().to_string()),
            match self.producer_kind {
                ProducerKind::Rosalind => "rosalind",
                ProducerKind::ExternalAnalyzer => "external-analyzer",
            },
            json_escape(&self.work_dir.display().to_string()),
            argv,
            self.inputs.len(),
            self.outputs.len(),
        )
    }
}

/// Receipt recipes rejected before process execution.
#[derive(Debug)]
pub enum ReplaySafetyError {
    /// Receipt integrity is absent or invalid.
    Integrity(String),
    /// Recipe shape or producer policy is unsafe or unsupported.
    UnsafeRecipe(String),
    /// A required input cannot be content-located.
    InputNotLocated(String),
    /// Filesystem or parsing failure while constructing the plan.
    Io(String),
}

impl std::fmt::Display for ReplaySafetyError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Integrity(message) => write!(f, "receipt integrity failed: {message}"),
            Self::UnsafeRecipe(message) => write!(f, "unsafe replay recipe: {message}"),
            Self::InputNotLocated(message) => write!(f, "input not located: {message}"),
            Self::Io(message) => write!(f, "cannot plan reproduction: {message}"),
        }
    }
}

impl std::error::Error for ReplaySafetyError {}

/// Deterministic first-party output types `reproduce` can byte-compare. BAM/bgzf
/// remains out of scope because its C compression implementation is not captured
/// by `deps_lock_blake3`.
fn output_is_byte_comparable(path: &str) -> bool {
    let p = path.to_ascii_lowercase();
    p.ends_with(".vcf")
        || p.ends_with(".tsv")
        || p.ends_with(".txt")
        || p.ends_with(".gvcf")
        || p.ends_with(".arrow")
        || p.ends_with(".rref")
}

/// Build a content-hash → path index of every file directly under `inputs_dir`.
fn index_inputs(inputs_dir: &Path) -> Result<BTreeMap<String, PathBuf>> {
    let mut idx = BTreeMap::new();
    let rd = std::fs::read_dir(inputs_dir)
        .with_context(|| format!("failed to read --inputs dir {}", inputs_dir.display()))?;
    for entry in rd {
        let p = entry?.path();
        if p.is_file() {
            if let Ok(h) = blake3_file(&p) {
                let absolute = std::fs::canonicalize(&p).unwrap_or(p);
                idx.entry(h).or_insert(absolute);
            }
        }
    }
    Ok(idx)
}

/// A fresh, unique temp directory for the re-run's outputs.
fn make_temp_dir() -> Result<PathBuf> {
    static COUNTER: AtomicU64 = AtomicU64::new(0);
    let n = COUNTER.fetch_add(1, Ordering::Relaxed);
    let d = std::env::temp_dir().join(format!("rosalind-reproduce-{}-{n}", std::process::id()));
    std::fs::create_dir_all(&d)
        .with_context(|| format!("failed to create temp dir {}", d.display()))?;
    Ok(d)
}

/// The build-identity of the *current* binary (the one performing the reproduction),
/// read by finalizing a throwaway manifest (which stamps build-identity from `build.rs`).
fn current_build_identity() -> BTreeMap<String, String> {
    let mut m = RunManifest::new("reproduce-probe");
    m.tool_version = env!("CARGO_PKG_VERSION").to_string();
    m.finalize();
    let mut out = BTreeMap::new();
    for k in [
        "code_git_sha",
        "code_dirty",
        "rustc_version",
        "target_triple",
        "deps_lock_blake3",
    ] {
        if let Some(v) = m.params.get(k) {
            out.insert(k.to_string(), v.clone());
        }
    }
    out
}

/// Read `peak_rss_bytes` from the receipt the re-run wrote next to `output_temp`.
fn read_rerun_manifest(output_temp: &Path) -> Option<RunManifest> {
    let mut name = output_temp.as_os_str().to_os_string();
    name.push(".manifest.json");
    let text = std::fs::read_to_string(PathBuf::from(name)).ok()?;
    RunManifest::from_canonical_json(&text).ok()
}

/// Validate a receipt and construct a shell-free execution plan without running it.
pub fn plan_reproduction(
    manifest_path: &Path,
    inputs_dir: &Path,
    binary: Option<&Path>,
) -> std::result::Result<ReproductionPlan, ReplaySafetyError> {
    let text = std::fs::read_to_string(manifest_path)
        .map_err(|error| ReplaySafetyError::Io(error.to_string()))?;
    let manifest = RunManifest::from_canonical_json(&text)
        .map_err(|error| ReplaySafetyError::Integrity(error.to_string()))?;
    build_reproduction_plan(manifest, inputs_dir, binary)
}

fn build_reproduction_plan(
    manifest: RunManifest,
    inputs_dir: &Path,
    binary: Option<&Path>,
) -> std::result::Result<ReproductionPlan, ReplaySafetyError> {
    if manifest.self_hash_ok() != Some(true) {
        return Err(ReplaySafetyError::Integrity(
            "claim self-hash is absent or does not match".to_string(),
        ));
    }
    if manifest.measurement_hash_ok() == Some(false) {
        return Err(ReplaySafetyError::Integrity(
            "measurement self-hash does not match".to_string(),
        ));
    }
    if manifest.claims_measurements() && manifest.measurement_hash_ok() != Some(true) {
        return Err(ReplaySafetyError::Integrity(
            "claimed measurement block is missing or unverifiable".to_string(),
        ));
    }
    if manifest.params.get("run_status").map(String::as_str) == Some("breached") {
        return Err(ReplaySafetyError::UnsafeRecipe(
            "breach receipts preserve partial evidence and are not executable success recipes"
                .to_string(),
        ));
    }
    if manifest.outputs.is_empty() {
        return Err(ReplaySafetyError::UnsafeRecipe(
            "receipt records no output artifact".to_string(),
        ));
    }
    if let Some(output) = manifest
        .outputs
        .iter()
        .find(|output| !output_is_byte_comparable(&output.path))
    {
        return Err(ReplaySafetyError::UnsafeRecipe(format!(
            "output {} is not byte-comparable",
            output.path
        )));
    }

    let tokens = crate::provenance::command_template_tokens(&manifest)
        .map_err(ReplaySafetyError::UnsafeRecipe)?;
    validate_recipe_tokens(&manifest, &tokens, binary.is_some())?;

    let producer_kind = producer_kind(&manifest);
    let execution_binary = match (producer_kind, binary) {
        (ProducerKind::Rosalind, Some(path)) | (ProducerKind::ExternalAnalyzer, Some(path)) => {
            std::fs::canonicalize(path).map_err(|error| {
                ReplaySafetyError::Io(format!("cannot resolve binary {}: {error}", path.display()))
            })?
        }
        (ProducerKind::Rosalind, None) => {
            std::env::current_exe().map_err(|error| ReplaySafetyError::Io(error.to_string()))?
        }
        (ProducerKind::ExternalAnalyzer, None) => {
            return Err(ReplaySafetyError::UnsafeRecipe(
                "external analyzer receipts require --binary PATH".to_string(),
            ))
        }
    };

    let index =
        index_inputs(inputs_dir).map_err(|error| ReplaySafetyError::Io(error.to_string()))?;
    let mut inputs = Vec::with_capacity(manifest.inputs.len());
    for (position, input) in manifest.inputs.iter().enumerate() {
        let path = index
            .get(&input.blake3)
            .cloned()
            .ok_or_else(|| ReplaySafetyError::InputNotLocated(format!("@in:{}", input.blake3)))?;
        inputs.push(LocatedArtifact {
            role: format!("input[{position}]"),
            blake3: input.blake3.clone(),
            path,
        });
    }

    let work_dir = make_temp_dir().map_err(|error| ReplaySafetyError::Io(error.to_string()))?;
    let mut output_by_hash = BTreeMap::new();
    let mut outputs = Vec::with_capacity(manifest.outputs.len());
    for (position, output) in manifest.outputs.iter().enumerate() {
        let extension = Path::new(&output.path)
            .extension()
            .and_then(|extension| extension.to_str())
            .unwrap_or("out");
        let path = work_dir.join(format!("out_{position}.{extension}"));
        output_by_hash.insert(output.blake3.clone(), path.clone());
        outputs.push(PlannedOutput {
            role: format!("output[{position}]"),
            recorded_blake3: output.blake3.clone(),
            path,
        });
    }
    let input_by_hash = inputs
        .iter()
        .map(|input| (input.blake3.clone(), input.path.clone()))
        .collect::<BTreeMap<_, _>>();
    let argv = match CommandCapture::argv_from_manifest(
        &manifest,
        &|hash| {
            input_by_hash
                .get(hash)
                .map(|path| path.display().to_string())
        },
        &|hash| {
            output_by_hash
                .get(hash)
                .map(|path| path.display().to_string())
                .unwrap_or_else(|| work_dir.join("out_unknown").display().to_string())
        },
    ) {
        Ok(argv) => argv,
        Err(error) => {
            let _ = std::fs::remove_dir_all(&work_dir);
            return Err(ReplaySafetyError::UnsafeRecipe(error));
        }
    };

    Ok(ReproductionPlan {
        binary: execution_binary,
        argv,
        work_dir,
        producer_kind,
        inputs,
        outputs,
        manifest,
        explicit_binary: binary.is_some(),
    })
}

fn producer_kind(manifest: &RunManifest) -> ProducerKind {
    match manifest.params.get("replay.kind").map(String::as_str) {
        Some("external-analyzer") => ProducerKind::ExternalAnalyzer,
        Some("rosalind") => ProducerKind::Rosalind,
        _ if manifest.params.get("producer.name").map(String::as_str) == Some("rosalind") => {
            ProducerKind::Rosalind
        }
        _ if manifest.params.contains_key("producer.name") => ProducerKind::ExternalAnalyzer,
        _ => ProducerKind::Rosalind,
    }
}

fn validate_recipe_tokens(
    manifest: &RunManifest,
    tokens: &[String],
    explicit_binary: bool,
) -> std::result::Result<(), ReplaySafetyError> {
    let prefix = manifest
        .subcommand
        .split_whitespace()
        .map(str::to_string)
        .collect::<Vec<_>>();
    if tokens.len() < prefix.len() || tokens[..prefix.len()] != prefix {
        return Err(ReplaySafetyError::UnsafeRecipe(
            "command prefix does not match receipt subcommand".to_string(),
        ));
    }
    if tokens.iter().any(|token| token == "--manifest") {
        return Err(ReplaySafetyError::UnsafeRecipe(
            "recorded --manifest paths are never replayed".to_string(),
        ));
    }
    if tokens.iter().any(|token| {
        token.starts_with("--manifest=")
            || token.starts_with("--output=")
            || token.starts_with("-o=")
    }) {
        return Err(ReplaySafetyError::UnsafeRecipe(
            "embedded manifest or output paths are never replayed".to_string(),
        ));
    }
    if let Some(command) = manifest.params.get("command") {
        if command != &tokens.join(" ") {
            return Err(ReplaySafetyError::UnsafeRecipe(
                "command and command_argv disagree".to_string(),
            ));
        }
    }

    let kind = producer_kind(manifest);
    if kind == ProducerKind::ExternalAnalyzer {
        let replay_schema = manifest
            .params
            .get("replay_schema")
            .and_then(|value| value.parse::<u32>().ok());
        if replay_schema != Some(3) {
            return Err(ReplaySafetyError::UnsafeRecipe(
                "external analyzers require replay_schema=3".to_string(),
            ));
        }
        if !explicit_binary {
            return Err(ReplaySafetyError::UnsafeRecipe(
                "external analyzer receipts require --binary PATH".to_string(),
            ));
        }
    } else if !built_in_prefix_allowed(&prefix) {
        return Err(ReplaySafetyError::UnsafeRecipe(format!(
            "Rosalind subcommand {:?} is not in the replay allowlist",
            manifest.subcommand
        )));
    }

    let markers = |prefix: &str| {
        let mut values = tokens
            .iter()
            .enumerate()
            .filter_map(|(index, token)| {
                token
                    .strip_prefix(prefix)
                    .map(|hash| (index, hash.to_string()))
            })
            .collect::<Vec<_>>();
        values.sort_by(|a, b| a.1.cmp(&b.1));
        values
    };
    let expected = |files: &[FileHash]| {
        let mut hashes = files
            .iter()
            .map(|file| file.blake3.clone())
            .collect::<Vec<_>>();
        hashes.sort();
        hashes
    };
    let inputs = markers("@in:");
    let outputs = markers("@out:");
    if inputs
        .iter()
        .map(|(_, hash)| hash.clone())
        .collect::<Vec<_>>()
        != expected(&manifest.inputs)
    {
        return Err(ReplaySafetyError::UnsafeRecipe(
            "input markers do not match receipt inputs exactly".to_string(),
        ));
    }
    if outputs
        .iter()
        .map(|(_, hash)| hash.clone())
        .collect::<Vec<_>>()
        != expected(&manifest.outputs)
    {
        return Err(ReplaySafetyError::UnsafeRecipe(
            "output markers do not match receipt outputs exactly".to_string(),
        ));
    }
    for (index, _) in inputs.iter().chain(outputs.iter()) {
        if *index == 0 || !tokens[index - 1].starts_with('-') {
            return Err(ReplaySafetyError::UnsafeRecipe(
                "every content marker must be the value of an option flag".to_string(),
            ));
        }
    }
    for (index, token) in tokens.iter().enumerate() {
        if (token == "-o" || token == "--output")
            && !tokens
                .get(index + 1)
                .is_some_and(|value| value.starts_with("@out:"))
        {
            return Err(ReplaySafetyError::UnsafeRecipe(
                "output flags must target a recorded @out marker".to_string(),
            ));
        }
    }
    let mut index = prefix.len();
    while index < tokens.len() {
        let flag = &tokens[index];
        if !flag.starts_with('-') {
            return Err(ReplaySafetyError::UnsafeRecipe(format!(
                "unrecorded positional operand {flag:?}"
            )));
        }
        let next = tokens.get(index + 1);
        if next.is_some_and(|value| value.starts_with("@in:") || value.starts_with("@out:")) {
            index += 2;
            continue;
        }
        let key = flag.trim_start_matches('-').replace('-', "_");
        let recorded = manifest.params.get(&key).ok_or_else(|| {
            ReplaySafetyError::UnsafeRecipe(format!(
                "command option {flag} has no matching claim parameter"
            ))
        })?;
        if recorded == "true" {
            index += 1;
            continue;
        }
        let value = next.ok_or_else(|| {
            ReplaySafetyError::UnsafeRecipe(format!("command option {flag} has no value"))
        })?;
        if value != recorded {
            return Err(ReplaySafetyError::UnsafeRecipe(format!(
                "command option {flag} disagrees with claim parameter {key}"
            )));
        }
        index += 2;
    }
    Ok(())
}

fn built_in_prefix_allowed(prefix: &[String]) -> bool {
    matches!(
        prefix.first().map(String::as_str),
        Some("variants" | "features" | "somatic")
    ) || (prefix.first().map(String::as_str) == Some("analyze")
        && matches!(
            prefix.get(1).map(String::as_str),
            Some("features" | "coverage")
        ))
        || (prefix.first().map(String::as_str) == Some("reference")
            && matches!(prefix.get(1).map(String::as_str), Some("build" | "convert")))
}

/// Re-derive the result recorded in `manifest_path` from inputs content-located under
/// `inputs_dir`, and compare byte-for-byte. See the module docs for scope.
pub fn reproduce(manifest_path: &Path, inputs_dir: &Path) -> Result<ReproReport> {
    reproduce_with_binary(manifest_path, inputs_dir, None)
}

/// Re-derive a result using an explicitly supplied executable. Receipt paths are
/// never executed automatically: callers must opt in to a third-party binary.
pub fn reproduce_with_binary(
    manifest_path: &Path,
    inputs_dir: &Path,
    binary: Option<&Path>,
) -> Result<ReproReport> {
    let text = std::fs::read_to_string(manifest_path)
        .with_context(|| format!("failed to read receipt {}", manifest_path.display()))?;
    let manifest = RunManifest::from_canonical_json(&text)
        .map_err(|e| anyhow!("malformed receipt {}: {e}", manifest_path.display()))?;

    let parent_claim = manifest.content_hash();
    let parent_subcommand = manifest.subcommand.clone();
    let original_code = build_identity_from_manifest(&manifest);
    let reproducer_code = current_build_identity();
    let execution_binary = match binary {
        Some(path) => path.to_path_buf(),
        None => std::env::current_exe().context("locating the rosalind binary to re-run")?,
    };

    let report =
        |verdict_label: &str, exit_code: i32, lines: Vec<String>, compared: bool| ReproReport {
            verdict_label: verdict_label.to_string(),
            exit_code,
            lines,
            parent_claim: parent_claim.clone(),
            parent_subcommand: parent_subcommand.clone(),
            outputs: Vec::new(),
            resource_here: ResourceHere::default(),
            original_code: original_code.clone(),
            reproducer_code: reproducer_code.clone(),
            execution_binary: execution_binary.clone(),
            compared,
        };

    // 1. Integrity gate (tamper-evidence). Budget-fit is verify's job, NOT reproduce's —
    //    a receipt that honestly records an over-budget run is still reproducible.
    if manifest.self_hash_ok() == Some(false) || manifest.measurement_hash_ok() == Some(false) {
        return Ok(report(
            "TAMPERED",
            5,
            vec![
                "  the receipt's self-hash does not match — it was modified after it was written"
                    .to_string(),
                "  VERDICT     : TAMPERED".to_string(),
            ],
            false,
        ));
    }

    // 2. Need a recorded command (schema >= 5) to replay.
    if crate::provenance::command_template_tokens(&manifest).is_err() {
        return Ok(report(
            "INCONCLUSIVE",
            7,
            vec![
                "  pre-schema-5 receipt: no recorded command to replay".to_string(),
                "  VERDICT     : INCONCLUSIVE".to_string(),
            ],
            false,
        ));
    }

    // 3. Supported, deterministic outputs (at least one to compare).
    if manifest.outputs.is_empty() {
        return Ok(report(
            "INCONCLUSIVE",
            7,
            vec![
                "  receipt records no outputs to compare (e.g. a stdout run)".to_string(),
                "  VERDICT     : INCONCLUSIVE".to_string(),
            ],
            false,
        ));
    }
    if let Some(o) = manifest
        .outputs
        .iter()
        .find(|o| !output_is_byte_comparable(&o.path))
    {
        return Ok(report(
            "INCONCLUSIVE",
            7,
            vec![
                format!(
                    "  output {} is not byte-comparable (BAM/bgzf rests on a C zlib outside the contract)",
                    o.path
                ),
                "  VERDICT     : INCONCLUSIVE".to_string(),
            ],
            false,
        ));
    }

    let plan = match build_reproduction_plan(manifest, inputs_dir, binary) {
        Ok(plan) => plan,
        Err(error) => {
            return Ok(report(
                "INCONCLUSIVE",
                7,
                vec![
                    format!("  cannot reproduce safely: {error}"),
                    "  VERDICT     : INCONCLUSIVE".to_string(),
                ],
                false,
            ))
        }
    };
    execute_reproduction(plan)
}

/// Execute a previously validated plan and compare every produced output byte-for-byte.
pub fn execute_reproduction(plan: ReproductionPlan) -> Result<ReproReport> {
    let ReproductionPlan {
        binary: execution_binary,
        argv,
        work_dir: work,
        producer_kind: _,
        inputs,
        outputs,
        manifest,
        explicit_binary,
    } = plan;
    let parent_claim = manifest.content_hash();
    let parent_subcommand = manifest.subcommand.clone();
    let original_code = build_identity_from_manifest(&manifest);
    let reproducer_code = current_build_identity();
    let report =
        |verdict_label: &str, exit_code: i32, lines: Vec<String>, compared: bool| ReproReport {
            verdict_label: verdict_label.to_string(),
            exit_code,
            lines,
            parent_claim: parent_claim.clone(),
            parent_subcommand: parent_subcommand.clone(),
            outputs: Vec::new(),
            resource_here: ResourceHere::default(),
            original_code: original_code.clone(),
            reproducer_code: reproducer_code.clone(),
            execution_binary: execution_binary.clone(),
            compared,
        };

    let mut command = std::process::Command::new(&execution_binary);
    command
        .args(&argv)
        .current_dir(&work)
        .env_clear()
        .env("HOME", &work)
        .env("TMPDIR", &work)
        .env("LC_ALL", "C")
        .env("TZ", "UTC");
    if let Some(path) = std::env::var_os("PATH") {
        command.env("PATH", path);
    }
    let child = command
        .output()
        .context("re-running the validated recorded command")?;
    if !child.status.success() {
        let _ = std::fs::remove_dir_all(&work);
        return Ok(report(
            "INCONCLUSIVE",
            7,
            vec![
                format!(
                    "  the re-run did not complete cleanly (exit {:?})",
                    child.status.code()
                ),
                format!(
                    "  stderr: {}",
                    String::from_utf8_lossy(&child.stderr).trim()
                ),
                "  VERDICT     : INCONCLUSIVE".to_string(),
            ],
            false,
        ));
    }

    // 7. Hash produced outputs and classify.
    let mut produced: Vec<(String, String)> = Vec::new();
    let mut cmps: Vec<OutputCmp> = Vec::new();
    for (i, (recorded, planned)) in manifest.outputs.iter().zip(outputs.iter()).enumerate() {
        let observed = blake3_file(&planned.path).with_context(|| {
            format!(
                "the re-run did not produce the expected output {}",
                planned.path.display()
            )
        })?;
        produced.push((format!("output[{i}]"), observed.clone()));
        cmps.push(OutputCmp {
            role: format!("output[{i}]"),
            recorded: recorded.blake3.clone(),
            observed: observed.clone(),
            matched: observed == recorded.blake3,
        });
    }
    let outcome = classify_outputs(&manifest.outputs, &produced);

    // 8. Machine-local resource + code context (diagnostic; never fatal).
    let declared_budget_mb = manifest
        .get_recorded("memory_budget_mb")
        .and_then(|v| v.parse().ok());
    let rerun_manifest = manifest
        .outputs
        .first()
        .and_then(|_| outputs.first().map(|output| &output.path))
        .and_then(|p| read_rerun_manifest(p));
    let peak_here = rerun_manifest
        .as_ref()
        .and_then(|m| m.get_recorded("peak_rss_bytes"))
        .and_then(|v| v.parse().ok());
    let fit = match (peak_here, declared_budget_mb) {
        (Some(peak), Some(mb)) => Some(crate::core::MemoryBudget::from_mb(mb).admits(peak)),
        _ => None,
    };
    let resource_here = ResourceHere {
        peak_rss_bytes: peak_here,
        declared_budget_mb,
        fit,
    };

    // 9. Build the report.
    let recorded_sha = manifest.params.get("code_git_sha").map(String::as_str);
    let here_sha = rerun_manifest
        .as_ref()
        .and_then(|m| m.params.get("code_git_sha"))
        .map(String::as_str);
    let code_line = match (recorded_sha, here_sha) {
        (Some(r), Some(h)) if r == h => format!("  code        : git {} — matches", short(h)),
        (Some(r), Some(h)) => format!(
            "  code        : git {} — DIFFERS from recorded {} (byte-match still counts)",
            short(h),
            short(r)
        ),
        _ => "  code        : (build-identity unavailable)".to_string(),
    };

    let mut lines = vec![
        format!("  claim       : {} (re-derived)", short(&parent_claim)),
        code_line,
        format!(
            "  inputs      : {} file(s) located by content hash",
            inputs.len()
        ),
    ];
    for c in &cmps {
        if c.matched {
            lines.push(format!(
                "  output      : {}  OK byte-identical (blake3 {})",
                c.role,
                short(&c.observed)
            ));
        } else {
            lines.push(format!(
                "  output      : {}  DIVERGED (recorded {} got {})",
                c.role,
                short(&c.recorded),
                short(&c.observed)
            ));
        }
    }
    if let (Some(peak), Some(mb)) = (
        resource_here.peak_rss_bytes,
        resource_here.declared_budget_mb,
    ) {
        lines.push(format!(
            "  resource    : peak {} MiB vs declared {} MiB (here: this machine)",
            peak / (1 << 20),
            mb
        ));
    }
    lines.push(format!("  VERDICT     : {}", outcome.verdict.label()));

    // Clean up the temp work dir (best-effort).
    std::fs::remove_dir_all(&work).ok();

    Ok(ReproReport {
        verdict_label: outcome.verdict.label().to_string(),
        exit_code: outcome.verdict.exit_code(),
        lines,
        parent_claim,
        parent_subcommand,
        outputs: cmps,
        resource_here,
        original_code,
        reproducer_code: rerun_manifest
            .as_ref()
            .map(build_identity_from_manifest)
            .unwrap_or_else(|| {
                if explicit_binary {
                    BTreeMap::new()
                } else {
                    reproducer_code
                }
            }),
        execution_binary,
        compared: true,
    })
}

fn build_identity_from_manifest(manifest: &RunManifest) -> BTreeMap<String, String> {
    [
        "code_git_sha",
        "code_dirty",
        "rustc_version",
        "target_triple",
        "deps_lock_blake3",
    ]
    .into_iter()
    .filter_map(|key| {
        manifest
            .params
            .get(key)
            .map(|value| (key.to_string(), value.clone()))
    })
    .collect()
}

fn json_escape(value: &str) -> String {
    let mut out = String::new();
    for ch in value.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out
}

/// First 10 hex chars of a digest for compact display.
fn short(hex: &str) -> &str {
    &hex[..hex.len().min(10)]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::provenance::{CommandCapture, FileHash};

    fn fh(hash: &str) -> FileHash {
        FileHash {
            path: "x".into(),
            blake3: hash.into(),
        }
    }

    #[test]
    fn classify_reproduced_when_all_outputs_match() {
        let recorded = vec![fh("h1")];
        let produced = vec![("o0".to_string(), "h1".to_string())];
        let v = classify_outputs(&recorded, &produced);
        assert!(matches!(v.verdict, Verdict::Reproduced));
        assert!(v.diffs.is_empty());
    }

    #[test]
    fn classify_diverged_names_the_first_mismatch() {
        let recorded = vec![fh("h1")];
        let produced = vec![("o0".to_string(), "DIFFERENT".to_string())];
        let v = classify_outputs(&recorded, &produced);
        assert!(matches!(v.verdict, Verdict::Diverged));
        assert!(v
            .diffs
            .iter()
            .any(|d| d.contains("h1") && d.contains("DIFFERENT")));
    }

    #[test]
    fn exit_codes_are_stable() {
        assert_eq!(Verdict::Reproduced.exit_code(), 0);
        assert_eq!(Verdict::Diverged.exit_code(), 6);
        assert_eq!(Verdict::Inconclusive.exit_code(), 7);
    }

    fn replayable_manifest() -> (RunManifest, Vec<String>) {
        let mut manifest = RunManifest::new("variants");
        manifest
            .params
            .insert("producer.name".to_string(), "rosalind".to_string());
        manifest
            .params
            .insert("replay.kind".to_string(), "rosalind".to_string());
        let mut capture = CommandCapture::new("variants");
        capture.input_hashed("--index", "path with spaces/参考.idx", "index-hash");
        capture.input_hashed("--alignments", "reads.bam", "bam-hash");
        capture.opt("--quality-threshold", "value with spaces 和");
        capture.output_hashed("-o", "calls.vcf", "output-hash");
        capture.record_into(&mut manifest);
        manifest.finalize();
        let tokens = crate::provenance::command_template_tokens(&manifest).unwrap();
        (manifest, tokens)
    }

    #[test]
    fn strict_recipe_accepts_tokenized_spaces_and_unicode() {
        let (manifest, tokens) = replayable_manifest();
        validate_recipe_tokens(&manifest, &tokens, false).unwrap();
    }

    #[test]
    fn strict_recipe_rejects_prefix_duplicate_and_unrecorded_output_operands() {
        let (manifest, tokens) = replayable_manifest();

        let mut wrong_prefix = tokens.clone();
        wrong_prefix[0] = "sort".to_string();
        assert!(validate_recipe_tokens(&manifest, &wrong_prefix, false).is_err());

        let mut duplicate = tokens.clone();
        duplicate.extend(["--index".to_string(), "@in:index-hash".to_string()]);
        assert!(validate_recipe_tokens(&manifest, &duplicate, false).is_err());

        let mut unrecorded_output = tokens;
        unrecorded_output.extend(["--output".to_string(), "/tmp/escape.vcf".to_string()]);
        assert!(validate_recipe_tokens(&manifest, &unrecorded_output, false).is_err());
    }

    #[test]
    fn strict_recipe_rejects_manifest_paths_unsupported_builtins_and_legacy_externals() {
        let (manifest, mut tokens) = replayable_manifest();
        tokens.push("--manifest=/tmp/escape.json".to_string());
        assert!(validate_recipe_tokens(&manifest, &tokens, false).is_err());

        let (mut unsupported, unsupported_tokens) = replayable_manifest();
        unsupported.subcommand = "index".to_string();
        let mut unsupported_tokens = unsupported_tokens;
        unsupported_tokens[0] = "index".to_string();
        unsupported
            .params
            .insert("command".to_string(), unsupported_tokens.join(" "));
        assert!(validate_recipe_tokens(&unsupported, &unsupported_tokens, false).is_err());

        let (mut external, external_tokens) = replayable_manifest();
        external
            .params
            .insert("replay.kind".to_string(), "external-analyzer".to_string());
        external
            .params
            .insert("replay_schema".to_string(), "2".to_string());
        assert!(validate_recipe_tokens(&external, &external_tokens, true).is_err());
    }
}
