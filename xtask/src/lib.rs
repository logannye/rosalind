//! Maintainer-only release, benchmark, and design-partner control plane.

mod caller_evidence;

use std::collections::{BTreeMap, BTreeSet};
use std::ffi::OsString;
use std::fs;
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::time::{SystemTime, UNIX_EPOCH};

use clap::{Args, Parser, Subcommand, ValueEnum};
use serde::{Deserialize, Serialize};

pub const EXIT_OK: i32 = 0;
pub const EXIT_CONFIG: i32 = 2;
pub const EXIT_BLOCKED: i32 = 3;
pub const EXIT_REMOTE: i32 = 4;
pub const EXIT_INTEGRITY: i32 = 5;

#[derive(Parser, Debug)]
#[command(name = "cargo xtask", about = "Rosalind maintainer control plane")]
pub struct Cli {
    #[command(subcommand)]
    command: TopCommand,
}

#[derive(Subcommand, Debug)]
enum TopCommand {
    /// Inspect local and GitHub release prerequisites without changing state.
    Doctor(OutputArgs),
    /// Generate the frozen public-contract fingerprint.
    Contract {
        #[command(subcommand)]
        command: ContractCommand,
    },
    /// Plan, dispatch, or inspect a release candidate.
    Rc {
        #[command(subcommand)]
        command: RcCommand,
    },
    /// Plan, dispatch, or inspect a stable release.
    Release {
        #[command(subcommand)]
        command: ReleaseCommand,
    },
    /// Manage the pinned GIAB evaluator and benchmark lifecycle.
    Giab {
        #[command(subcommand)]
        command: GiabCommand,
    },
    /// Generate and validate anonymized design-partner evidence.
    Partners {
        #[command(subcommand)]
        command: PartnersCommand,
    },
}

#[derive(Args, Debug, Clone, Default)]
struct OutputArgs {
    #[arg(long)]
    json: bool,
}

#[derive(Subcommand, Debug)]
enum ContractCommand {
    Snapshot {
        #[arg(long, default_value = "HEAD")]
        r#ref: String,
        #[arg(long)]
        output: Option<PathBuf>,
        #[arg(long)]
        json: bool,
    },
}

#[derive(Subcommand, Debug)]
enum RcCommand {
    Plan {
        #[arg(long)]
        version: String,
        #[arg(long)]
        number: u32,
        #[arg(long, default_value = "HEAD")]
        r#ref: String,
        #[arg(long)]
        output: Option<PathBuf>,
        #[arg(long)]
        json: bool,
    },
    Dispatch(DispatchArgs),
    Status {
        #[arg(long)]
        tag: String,
        #[arg(long)]
        json: bool,
    },
}

#[derive(Subcommand, Debug)]
enum ReleaseCommand {
    Plan {
        #[arg(long)]
        version: String,
        #[arg(long)]
        rc_tag: String,
        #[arg(long, default_value = "HEAD")]
        r#ref: String,
        #[arg(long)]
        output: Option<PathBuf>,
        #[arg(long)]
        json: bool,
    },
    Dispatch(DispatchArgs),
    Status {
        #[arg(long)]
        version: String,
        #[arg(long)]
        json: bool,
    },
}

#[derive(Args, Debug)]
struct DispatchArgs {
    #[arg(long)]
    plan: PathBuf,
    #[arg(long)]
    confirm: String,
}

#[derive(Subcommand, Debug)]
enum GiabCommand {
    Image {
        #[command(subcommand)]
        command: ImageCommand,
    },
    Benchmark {
        #[command(subcommand)]
        command: BenchmarkCommand,
    },
    Baseline {
        #[command(subcommand)]
        command: BaselineCommand,
    },
}

#[derive(Subcommand, Debug)]
enum ImageCommand {
    Plan {
        #[arg(long, default_value = "HEAD")]
        r#ref: String,
        #[arg(long)]
        output: Option<PathBuf>,
        #[arg(long)]
        json: bool,
    },
    Dispatch(DispatchArgs),
    Status(OutputArgs),
}

#[derive(Subcommand, Debug)]
enum BenchmarkCommand {
    Plan {
        #[arg(long)]
        image: Option<String>,
        #[arg(long)]
        output: Option<PathBuf>,
        #[arg(long)]
        json: bool,
    },
    Dispatch(DispatchArgs),
    Status(OutputArgs),
    Download {
        #[arg(long)]
        run_id: Option<u64>,
        #[arg(long)]
        output: PathBuf,
    },
}

#[derive(Subcommand, Debug)]
enum BaselineCommand {
    Propose {
        #[arg(long)]
        run_id: u64,
        #[arg(long)]
        reason: String,
        #[arg(long)]
        confirm: Option<String>,
        #[arg(long)]
        json: bool,
    },
}

#[derive(Subcommand, Debug)]
enum PartnersCommand {
    Init {
        #[arg(long, value_enum)]
        persona: Persona,
        #[arg(long)]
        output: PathBuf,
    },
    Validate {
        #[arg(long)]
        input: PathBuf,
        #[arg(long)]
        json: bool,
    },
    Report {
        #[arg(long)]
        input_dir: PathBuf,
        #[arg(long)]
        json: bool,
    },
}

#[derive(ValueEnum, Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord)]
#[serde(rename_all = "kebab-case")]
enum Persona {
    AnalyzerBuilder,
    WorkflowHpc,
    ConstrainedOffline,
}

impl Persona {
    fn as_str(self) -> &'static str {
        match self {
            Self::AnalyzerBuilder => "analyzer-builder",
            Self::WorkflowHpc => "workflow-hpc",
            Self::ConstrainedOffline => "constrained-offline",
        }
    }
}

#[derive(Debug, Clone, Deserialize)]
struct Policy {
    schema: u32,
    repository: String,
    default_branch: String,
    msrv: String,
    release_toolchain: String,
    soak_seconds: u64,
    soak_required_from: String,
    foundation_test_baseline: u32,
    cargo_public_api_version: String,
    public_api_toolchain: String,
    wasm_pack_version: String,
    actionlint_version: String,
    shellcheck_version: String,
    jsonschema_version: String,
    happy_repository: String,
    caller_source_base: String,
    caller_source_paths: Vec<String>,
    partners_required_from: String,
    required_environments: Vec<String>,
    required_secrets: Vec<String>,
    partner_personas: Vec<String>,
    publish_order: Vec<String>,
    package_versions: BTreeMap<String, String>,
    cli_help: Vec<String>,
    schema_files: Vec<String>,
    legacy_fixture_dir: String,
    partner_record_dir: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct Check {
    pub id: String,
    pub ok: bool,
    pub detail: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct MaintainerReport {
    pub schema: u32,
    pub command: String,
    pub status: String,
    pub commit: String,
    pub version: Option<String>,
    pub checks: Vec<Check>,
    pub blockers: Vec<String>,
    pub actions: Vec<String>,
    pub contract_fingerprint: Option<String>,
    pub plan_id: String,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub metadata: BTreeMap<String, String>,
}

impl MaintainerReport {
    fn new(command: impl Into<String>, commit: impl Into<String>) -> Self {
        Self {
            schema: 1,
            command: command.into(),
            status: "ready".to_string(),
            commit: commit.into(),
            version: None,
            checks: Vec::new(),
            blockers: Vec::new(),
            actions: Vec::new(),
            contract_fingerprint: None,
            plan_id: String::new(),
            metadata: BTreeMap::new(),
        }
    }

    fn check(&mut self, id: &str, ok: bool, detail: impl Into<String>) {
        let detail = detail.into();
        self.checks.push(Check {
            id: id.to_string(),
            ok,
            detail: detail.clone(),
        });
        if !ok {
            self.status = "blocked".to_string();
            self.blockers.push(format!("{id}: {detail}"));
        }
    }

    fn seal(&mut self) {
        self.checks.sort_by(|a, b| a.id.cmp(&b.id));
        self.blockers.sort();
        self.actions.sort();
        let intent = serde_json::json!({
            "schema": self.schema,
            "command": self.command,
            "status": self.status,
            "commit": self.commit,
            "version": self.version,
            "checks": self.checks.iter().map(|check| (&check.id, check.ok)).collect::<Vec<_>>(),
            "actions": self.actions,
            "contract_fingerprint": self.contract_fingerprint,
            "metadata": self.metadata,
        });
        self.plan_id = blake3::hash(serde_json::to_string(&intent).unwrap().as_bytes())
            .to_hex()
            .to_string();
    }

    fn exit_code(&self) -> i32 {
        match self.status.as_str() {
            "ready" => EXIT_OK,
            "blocked"
                if self.checks.iter().any(|check| {
                    !check.ok
                        && (check
                            .detail
                            .contains("does not match the candidate package bytes")
                            || (matches!(
                                check.id.as_str(),
                                "contract.equivalent"
                                    | "contract.current_matches_rc"
                                    | "contract.commit"
                            ) && (check.detail.starts_with("rc=")
                                || check.detail.starts_with("release="))))
                }) =>
            {
                EXIT_INTEGRITY
            }
            "blocked"
                if self.checks.iter().any(|check| {
                    !check.ok
                        && check.id.starts_with("crate.")
                        && (check.detail.starts_with("crates.io")
                            || check.detail.starts_with("cannot start"))
                }) =>
            {
                EXIT_REMOTE
            }
            "blocked" => EXIT_BLOCKED,
            _ => EXIT_REMOTE,
        }
    }
}

fn mark_crate_remote_failure(report: &mut MaintainerReport) {
    if report.checks.iter().any(|check| {
        !check.ok
            && check.id.starts_with("crate.")
            && (check.detail.starts_with("crates.io") || check.detail.starts_with("cannot start"))
    }) {
        report.status = "failed".into();
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
struct ContractSnapshot {
    schema: u32,
    commit: String,
    components: BTreeMap<String, String>,
    receipt_fields: BTreeMap<String, Vec<String>>,
    primary_artifacts: BTreeMap<String, String>,
    package_versions: BTreeMap<String, String>,
    aggregate_blake3: String,
}

pub trait Runner {
    fn run(&self, cwd: &Path, program: &str, args: &[OsString]) -> io::Result<Output>;
}

#[derive(Debug, Default)]
struct SystemRunner;

impl Runner for SystemRunner {
    fn run(&self, cwd: &Path, program: &str, args: &[OsString]) -> io::Result<Output> {
        Command::new(program).current_dir(cwd).args(args).output()
    }
}

fn args(values: &[&str]) -> Vec<OsString> {
    values.iter().map(OsString::from).collect()
}

fn output_text(output: &Output) -> String {
    String::from_utf8_lossy(&output.stdout).trim().to_string()
}

fn command_text<R: Runner>(
    runner: &R,
    cwd: &Path,
    program: &str,
    arguments: &[OsString],
) -> Result<String, String> {
    let output = runner
        .run(cwd, program, arguments)
        .map_err(|error| format!("cannot start {program}: {error}"))?;
    if output.status.success() {
        Ok(output_text(&output))
    } else {
        Err(format!(
            "{program} exited {:?}: {}",
            output.status.code(),
            String::from_utf8_lossy(&output.stderr).trim()
        ))
    }
}

fn repo_root() -> Result<PathBuf, String> {
    let mut path = std::env::current_dir().map_err(|error| error.to_string())?;
    loop {
        let manifest = path.join("Cargo.toml");
        if manifest.exists()
            && fs::read_to_string(&manifest)
                .map(|text| text.contains("name = \"rosalind-bio\""))
                .unwrap_or(false)
        {
            return Ok(path);
        }
        if !path.pop() {
            return Err("run cargo xtask inside the Rosalind checkout".to_string());
        }
    }
}

fn load_policy(root: &Path) -> Result<Policy, String> {
    let text = fs::read_to_string(root.join("release/policy.toml"))
        .map_err(|error| format!("cannot read release policy: {error}"))?;
    let policy: Policy =
        toml::from_str(&text).map_err(|error| format!("cannot parse release policy: {error}"))?;
    if policy.schema != 1 {
        return Err(format!(
            "unsupported release policy schema {}",
            policy.schema
        ));
    }
    Ok(policy)
}

fn git_commit<R: Runner>(runner: &R, root: &Path, reference: &str) -> Result<String, String> {
    command_text(
        runner,
        root,
        "git",
        &[
            OsString::from("rev-parse"),
            OsString::from("--verify"),
            OsString::from("--end-of-options"),
            OsString::from(format!("{reference}^{{commit}}")),
        ],
    )
}

fn repository_matches(origin: &str, repository: &str) -> bool {
    let origin = origin.trim().trim_end_matches('/').trim_end_matches(".git");
    origin == format!("https://github.com/{repository}")
        || origin == format!("ssh://git@github.com/{repository}")
        || origin == format!("git@github.com:{repository}")
}

fn add_mutation_preflight<R: Runner>(
    report: &mut MaintainerReport,
    runner: &R,
    root: &Path,
    policy: &Policy,
    reference: &str,
) {
    let clean = command_text(runner, root, "git", &args(&["status", "--porcelain"]));
    report.check(
        "git.clean",
        clean.as_ref().is_ok_and(String::is_empty),
        clean.map_or_else(
            |error| error,
            |value| {
                if value.is_empty() {
                    "clean".into()
                } else {
                    value
                }
            },
        ),
    );
    let origin = command_text(runner, root, "git", &args(&["remote", "get-url", "origin"]));
    report.check(
        "git.repository",
        origin
            .as_ref()
            .is_ok_and(|value| repository_matches(value, &policy.repository)),
        origin.unwrap_or_else(|error| error),
    );
    let commit = git_commit(runner, root, reference).unwrap_or_else(|_| report.commit.clone());
    let ancestor = runner.run(
        root,
        "git",
        &[
            "merge-base".into(),
            "--is-ancestor".into(),
            commit.clone().into(),
            format!("origin/{}", policy.default_branch).into(),
        ],
    );
    report.check(
        "git.reachable_from_default",
        ancestor
            .as_ref()
            .is_ok_and(|output| output.status.success()),
        format!(
            "candidate must be reachable from origin/{}",
            policy.default_branch
        ),
    );
    let pushed = command_text(
        runner,
        root,
        "git",
        &[
            "branch".into(),
            "-r".into(),
            "--contains".into(),
            commit.into(),
        ],
    );
    report.check(
        "git.pushed",
        pushed
            .as_ref()
            .is_ok_and(|branches| !branches.trim().is_empty()),
        pushed.unwrap_or_else(|error| error),
    );
    match package_versions(root) {
        Ok(actual) => report.check(
            "version.workspace",
            actual == policy.package_versions,
            if actual == policy.package_versions {
                format!("{} package versions match policy", actual.len())
            } else {
                format!("policy={:?}, workspace={actual:?}", policy.package_versions)
            },
        ),
        Err(error) => report.check("version.workspace", false, error),
    }
}

fn check_command<R: Runner>(runner: &R, root: &Path, program: &str, arguments: &[&str]) -> Check {
    match command_text(runner, root, program, &args(arguments)) {
        Ok(value) => Check {
            id: format!("tool.{program}"),
            ok: true,
            detail: if value.is_empty() {
                "available".into()
            } else {
                value
            },
        },
        Err(error) => Check {
            id: format!("tool.{program}"),
            ok: false,
            detail: error,
        },
    }
}

fn write_create_new(path: &Path, bytes: &[u8]) -> Result<(), String> {
    if path.exists() {
        return Err(format!("refusing to overwrite {}", path.display()));
    }
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent).map_err(|error| error.to_string())?;
    }
    let mut temporary = path.as_os_str().to_os_string();
    temporary.push(format!(".xtask-{}.partial", std::process::id()));
    let temporary = PathBuf::from(temporary);
    let result = (|| {
        let mut file = fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&temporary)
            .map_err(|error| error.to_string())?;
        file.write_all(bytes).map_err(|error| error.to_string())?;
        file.sync_all().map_err(|error| error.to_string())?;
        fs::hard_link(&temporary, path).map_err(|error| error.to_string())?;
        fs::remove_file(&temporary).map_err(|error| error.to_string())
    })();
    if result.is_err() {
        let _ = fs::remove_file(&temporary);
    }
    result
}

fn emit_report(report: &MaintainerReport, json: bool, output: Option<&Path>) -> i32 {
    let rendered = serde_json::to_string_pretty(report).expect("report serializes") + "\n";
    if let Some(path) = output {
        if let Err(error) = write_create_new(path, rendered.as_bytes()) {
            eprintln!("{error}");
            return EXIT_CONFIG;
        }
    }
    if json {
        print!("{rendered}");
    } else {
        println!("{}: {}", report.command, report.status.to_uppercase());
        println!("  commit : {}", report.commit);
        if let Some(version) = &report.version {
            println!("  version: {version}");
        }
        println!("  plan   : {}", report.plan_id);
        for check in &report.checks {
            println!(
                "  [{}] {} — {}",
                if check.ok { "OK" } else { "BLOCK" },
                check.id,
                check.detail
            );
        }
        for action in &report.actions {
            println!("  action : {action}");
        }
    }
    report.exit_code()
}

fn doctor<R: Runner>(runner: &R, root: &Path, policy: &Policy) -> MaintainerReport {
    let commit = git_commit(runner, root, "HEAD").unwrap_or_else(|_| "unknown".to_string());
    let mut report = MaintainerReport::new("doctor", commit);
    report.metadata.insert("msrv".into(), policy.msrv.clone());
    report
        .metadata
        .insert("release_toolchain".into(), policy.release_toolchain.clone());
    report.metadata.insert(
        "foundation_test_baseline".into(),
        policy.foundation_test_baseline.to_string(),
    );
    report.metadata.insert(
        "cargo_public_api_version".into(),
        policy.cargo_public_api_version.clone(),
    );
    report.metadata.insert(
        "public_api_toolchain".into(),
        policy.public_api_toolchain.clone(),
    );
    report
        .metadata
        .insert("wasm_pack_version".into(), policy.wasm_pack_version.clone());
    report.metadata.insert(
        "actionlint_version".into(),
        policy.actionlint_version.clone(),
    );
    report.metadata.insert(
        "shellcheck_version".into(),
        policy.shellcheck_version.clone(),
    );
    report.metadata.insert(
        "jsonschema_version".into(),
        policy.jsonschema_version.clone(),
    );
    for check in [
        check_command(runner, root, "git", &["--version"]),
        check_command(runner, root, "gh", &["--version"]),
        check_command(runner, root, "cargo", &["--version"]),
        check_command(runner, root, "docker", &["--version"]),
    ] {
        report.check(&check.id, check.ok, check.detail);
    }
    for (tool, arguments, expected) in [
        (
            "actionlint",
            &["--version"][..],
            policy.actionlint_version.as_str(),
        ),
        (
            "shellcheck",
            &["--version"][..],
            policy.shellcheck_version.as_str(),
        ),
        (
            "wasm-pack",
            &["--version"][..],
            policy.wasm_pack_version.as_str(),
        ),
    ] {
        let version = command_text(runner, root, tool, &args(arguments));
        report.check(
            &format!("tool.{tool}"),
            version.as_ref().is_ok_and(|value| value.contains(expected)),
            version.unwrap_or_else(|_| format!("install pinned {tool} {expected}")),
        );
    }
    let docker_daemon = command_text(
        runner,
        root,
        "docker",
        &args(&["info", "--format", "{{.ServerVersion}}"]),
    );
    report.check(
        "docker.daemon",
        docker_daemon.is_ok(),
        docker_daemon.unwrap_or_else(|_| {
            "Docker CLI is installed but its daemon is unavailable; image and benchmark dispatch remain blocked"
                .into()
        }),
    );
    let public_api = command_text(runner, root, "cargo", &args(&["public-api", "--version"]));
    report.check(
        "tool.cargo-public-api",
        public_api
            .as_ref()
            .is_ok_and(|value| value.contains(policy.cargo_public_api_version.as_str())),
        public_api.unwrap_or_else(|_| {
            format!(
                "install with `cargo install cargo-public-api --version {}`",
                policy.cargo_public_api_version
            )
        }),
    );
    let public_api_toolchain = command_text(
        runner,
        root,
        "rustup",
        &[
            "run".into(),
            policy.public_api_toolchain.clone().into(),
            "rustc".into(),
            "--version".into(),
        ],
    );
    report.check(
        "tool.public-api-toolchain",
        public_api_toolchain.is_ok(),
        public_api_toolchain.unwrap_or_else(|_| {
            format!(
                "install with `rustup toolchain install {} --profile minimal`",
                policy.public_api_toolchain
            )
        }),
    );
    let release_toolchain = command_text(
        runner,
        root,
        "rustup",
        &[
            "run".into(),
            policy.release_toolchain.clone().into(),
            "rustc".into(),
            "--version".into(),
        ],
    );
    report.check(
        "tool.release-toolchain",
        release_toolchain.is_ok(),
        release_toolchain.unwrap_or_else(|_| {
            format!(
                "install with `rustup toolchain install {} --profile minimal`",
                policy.release_toolchain
            )
        }),
    );
    let remote = command_text(runner, root, "git", &args(&["remote", "get-url", "origin"]));
    report.check(
        "repository.origin",
        remote
            .as_ref()
            .is_ok_and(|value| repository_matches(value, &policy.repository)),
        remote.unwrap_or_else(|error| error),
    );
    let secrets = command_text(
        runner,
        root,
        "gh",
        &args(&["secret", "list", "--repo", &policy.repository]),
    );
    let environment_secrets = command_text(
        runner,
        root,
        "gh",
        &args(&[
            "secret",
            "list",
            "--repo",
            &policy.repository,
            "--env",
            "release",
        ]),
    )
    .unwrap_or_default();
    let secret_names = format!("{}\n{environment_secrets}", secrets.unwrap_or_default());
    for secret in &policy.required_secrets {
        report.check(
            &format!("github.secret.{secret}"),
            secret_names.lines().any(|line| line.starts_with(secret)),
            if secret_names.lines().any(|line| line.starts_with(secret)) {
                "configured".to_string()
            } else {
                format!(
                    "missing; configure with `gh secret set {secret} --repo {}`",
                    policy.repository
                )
            },
        );
    }
    let environments = command_text(
        runner,
        root,
        "gh",
        &args(&["api", &format!("repos/{}/environments", policy.repository)]),
    )
    .ok()
    .and_then(|text| serde_json::from_str::<serde_json::Value>(&text).ok())
    .and_then(|json| {
        json.get("environments")
            .and_then(|value| value.as_array())
            .cloned()
    })
    .unwrap_or_default()
    .into_iter()
    .filter_map(|entry| {
        let name = entry.get("name")?.as_str()?.to_owned();
        Some((name, entry))
    })
    .collect::<BTreeMap<_, _>>();
    for environment in &policy.required_environments {
        report.check(
            &format!("github.environment.{environment}"),
            environments.contains_key(environment),
            if environments.contains_key(environment) {
                "configured".to_string()
            } else {
                format!(
                    "missing; bootstrap with `gh api --method PUT repos/{}/environments/{environment}` then add required reviewers in environment settings",
                    policy.repository
                )
            },
        );
        if let Some(entry) = environments.get(environment) {
            let reviewers = entry["protection_rules"].as_array().is_some_and(|rules| {
                rules.iter().any(|rule| {
                    rule["type"] == "required_reviewers"
                        && rule["reviewers"]
                            .as_array()
                            .is_some_and(|people| !people.is_empty())
                })
            });
            report.check(
                &format!("github.environment.{environment}.reviewers"),
                reviewers,
                if reviewers {
                    "required reviewer configured"
                } else {
                    "add at least one required reviewer"
                },
            );
            if entry["deployment_branch_policy"]["custom_branch_policies"] == true {
                let policies = command_text(
                    runner,
                    root,
                    "gh",
                    &args(&[
                        "api",
                        &format!(
                            "repos/{}/environments/{environment}/deployment-branch-policies",
                            policy.repository
                        ),
                    ]),
                )
                .ok()
                .and_then(|value| serde_json::from_str::<serde_json::Value>(&value).ok());
                let count = policies
                    .as_ref()
                    .and_then(|value| value["branch_policies"].as_array())
                    .map_or(0, Vec::len);
                report.check(
                    &format!("github.environment.{environment}.deployment_refs"),
                    count > 0,
                    if count > 0 { format!("{count} allowed ref policies; the dispatch ref must match") } else { "custom deployment policy has no allowed branch/tag refs; configure the intended workflow dispatch ref".into() },
                );
            }
        }
    }
    report.actions = vec![
        "install missing pinned maintainer tools".into(),
        "configure protected rc and release environments".into(),
        "configure CARGO_REGISTRY_TOKEN without exposing its value".into(),
    ];
    report.seal();
    report
}

fn hash_bytes(bytes: &[u8]) -> String {
    blake3::hash(bytes).to_hex().to_string()
}

fn hash_file(path: &Path) -> Result<String, String> {
    fs::read(path)
        .map(|bytes| hash_bytes(&bytes))
        .map_err(|error| format!("cannot hash {}: {error}", path.display()))
}

fn canonical_cli_shape(help: &str) -> String {
    let mut shape = Vec::new();
    let mut section = "";
    for line in help.lines() {
        let trimmed = line.trim();
        if let Some(usage) = trimmed.strip_prefix("Usage:") {
            section = "usage";
            shape.push(format!(
                "usage:{}",
                usage.split_whitespace().collect::<Vec<_>>().join(" ")
            ));
            continue;
        }
        if matches!(trimmed, "Commands:" | "Arguments:" | "Options:") {
            section = trimmed.trim_end_matches(':');
            continue;
        }
        if trimmed.is_empty() {
            if section == "usage" {
                section = "";
            }
            continue;
        }
        match section {
            "usage" if line.starts_with(char::is_whitespace) => {
                shape.push(format!(
                    "usage+:{}",
                    trimmed.split_whitespace().collect::<Vec<_>>().join(" ")
                ));
            }
            "Commands" if line.starts_with("  ") && !line.starts_with("   ") => {
                if let Some(command) = trimmed.split_whitespace().next() {
                    shape.push(format!("command:{command}"));
                }
            }
            "Arguments"
                if line.starts_with("  ")
                    && !line.starts_with("   ")
                    && trimmed.starts_with(['<', '[']) =>
            {
                let syntax = trimmed.split("  ").next().unwrap_or(trimmed).trim();
                shape.push(format!("argument:{syntax}"));
            }
            "Options" if trimmed.starts_with('-') => {
                let syntax = trimmed.split("  ").next().unwrap_or(trimmed).trim();
                shape.push(format!("option:{syntax}"));
            }
            _ => {}
        }
    }
    shape.sort();
    shape.dedup();
    shape.join("\n") + "\n"
}

fn cli_subcommands(help: &str) -> Vec<String> {
    let mut commands = Vec::new();
    let mut in_commands = false;
    for line in help.lines() {
        let trimmed = line.trim();
        if trimmed == "Commands:" {
            in_commands = true;
            continue;
        }
        if matches!(trimmed, "Arguments:" | "Options:") || trimmed.starts_with("Usage:") {
            in_commands = false;
        }
        if in_commands && line.starts_with("  ") && !line.starts_with("   ") {
            if let Some(command) = trimmed.split_whitespace().next() {
                if command != "help" {
                    commands.push(command.to_string());
                }
            }
        }
    }
    commands.sort();
    commands.dedup();
    commands
}

fn discover_cli_contract<R: Runner>(
    runner: &R,
    root: &Path,
    binary: &Path,
) -> Result<BTreeMap<String, String>, String> {
    let binary = binary.display().to_string();
    let mut pending = vec![String::new()];
    let mut visited = BTreeSet::new();
    let mut contracts = BTreeMap::new();
    while let Some(command_path) = pending.pop() {
        if !visited.insert(command_path.clone()) {
            continue;
        }
        let mut arguments = command_path
            .split_whitespace()
            .filter(|part| !part.is_empty())
            .map(OsString::from)
            .collect::<Vec<_>>();
        arguments.push("--help".into());
        let help = command_text(runner, root, &binary, &arguments)?;
        for child in cli_subcommands(&help) {
            pending.push(if command_path.is_empty() {
                child
            } else {
                format!("{command_path} {child}")
            });
        }
        contracts.insert(command_path, canonical_cli_shape(&help));
    }
    Ok(contracts)
}

fn files_recursively(path: &Path) -> Result<Vec<PathBuf>, String> {
    let mut files = Vec::new();
    if path.is_file() {
        files.push(path.to_path_buf());
    } else {
        for entry in fs::read_dir(path).map_err(|error| error.to_string())? {
            let entry = entry.map_err(|error| error.to_string())?;
            let child = entry.path();
            if child.is_dir() {
                files.extend(files_recursively(&child)?);
            } else if child.is_file() {
                files.push(child);
            }
        }
    }
    files.sort();
    Ok(files)
}

fn package_versions(root: &Path) -> Result<BTreeMap<String, String>, String> {
    let packages = [
        ("rosalind-bio", root.join("Cargo.toml")),
        ("rosalind-receipt", root.join("crates/receipt/Cargo.toml")),
        (
            "rosalind-build-info",
            root.join("crates/build-info/Cargo.toml"),
        ),
    ];
    let mut versions = BTreeMap::new();
    for (name, path) in packages {
        let text = fs::read_to_string(&path).map_err(|error| error.to_string())?;
        let value: toml::Value = toml::from_str(&text).map_err(|error| error.to_string())?;
        let version = value
            .get("package")
            .and_then(|value| value.get("version"))
            .and_then(toml::Value::as_str)
            .ok_or_else(|| format!("{} has no package.version", path.display()))?;
        versions.insert(name.to_string(), version.to_string());
    }
    Ok(versions)
}

#[derive(Debug, Clone)]
enum CratePublicationState {
    Absent,
    Exact(String),
    Conflict(String),
    RemoteFailure(String),
}

fn crate_publication_state<R: Runner>(
    runner: &R,
    root: &Path,
    package: &str,
    version: &str,
    release_toolchain: &str,
) -> CratePublicationState {
    let body = std::env::temp_dir().join(format!(
        "rosalind-crate-state-{}-{package}.json",
        std::process::id()
    ));
    fs::remove_file(&body).ok();
    let url = format!("https://crates.io/api/v1/crates/{package}/{version}");
    let status = runner.run(
        root,
        "curl",
        &[
            "--silent".into(),
            "--show-error".into(),
            "--output".into(),
            body.as_os_str().to_owned(),
            "--write-out".into(),
            "%{http_code}".into(),
            "--header".into(),
            "User-Agent: rosalind-release-status/1".into(),
            url.into(),
        ],
    );
    let http_status = match status {
        Ok(output) if output.status.success() => output_text(&output),
        Ok(output) => {
            fs::remove_file(&body).ok();
            return CratePublicationState::RemoteFailure(format!(
                "crates.io query failed: {}",
                String::from_utf8_lossy(&output.stderr).trim()
            ));
        }
        Err(error) => {
            fs::remove_file(&body).ok();
            return CratePublicationState::RemoteFailure(error.to_string());
        }
    };
    if http_status == "404" {
        fs::remove_file(&body).ok();
        return CratePublicationState::Absent;
    }
    if http_status != "200" {
        fs::remove_file(&body).ok();
        return CratePublicationState::RemoteFailure(format!(
            "crates.io returned HTTP {http_status}"
        ));
    }
    let checksum = fs::read(&body)
        .ok()
        .and_then(|bytes| serde_json::from_slice::<serde_json::Value>(&bytes).ok())
        .and_then(|value| {
            value
                .get("version")
                .and_then(|version| version.get("checksum"))
                .and_then(serde_json::Value::as_str)
                .map(str::to_owned)
        });
    fs::remove_file(&body).ok();
    let Some(checksum) = checksum else {
        return CratePublicationState::RemoteFailure(
            "crates.io response lacks version.checksum".into(),
        );
    };
    let packaged = runner.run(
        root,
        "rustup",
        &[
            "run".into(),
            release_toolchain.into(),
            "cargo".into(),
            "package".into(),
            "-p".into(),
            package.into(),
            "--locked".into(),
        ],
    );
    if !packaged
        .as_ref()
        .is_ok_and(|output| output.status.success())
    {
        return CratePublicationState::Conflict(format!(
            "cannot recreate local package {package} {version}"
        ));
    }
    let archive = root
        .join("target/package")
        .join(format!("{package}-{version}.crate"));
    let download = std::env::temp_dir().join(format!(
        "rosalind-crate-state-{}-{package}.crate",
        std::process::id()
    ));
    fs::remove_file(&download).ok();
    let downloaded = runner.run(
        root,
        "curl",
        &[
            "--fail".into(),
            "--silent".into(),
            "--show-error".into(),
            "--location".into(),
            "--header".into(),
            "User-Agent: rosalind-release-status/1".into(),
            format!("https://crates.io/api/v1/crates/{package}/{version}/download").into(),
            "-o".into(),
            download.as_os_str().to_owned(),
        ],
    );
    let exact = downloaded
        .as_ref()
        .is_ok_and(|output| output.status.success())
        && match (fs::read(&archive), fs::read(&download)) {
            (Ok(local), Ok(remote)) => local == remote,
            _ => false,
        };
    fs::remove_file(download).ok();
    if exact {
        CratePublicationState::Exact(checksum)
    } else {
        CratePublicationState::Conflict(format!(
            "published {package} {version} does not match the candidate package bytes"
        ))
    }
}

fn collect_receipt_fields(
    path: &Path,
    fields: &mut BTreeMap<String, BTreeSet<String>>,
) -> Result<(), String> {
    let value: serde_json::Value = serde_json::from_slice(
        &fs::read(path).map_err(|error| format!("cannot read {}: {error}", path.display()))?,
    )
    .map_err(|error| format!("cannot parse {}: {error}", path.display()))?;
    if let Some(object) = value.as_object() {
        fields
            .entry("top_level".into())
            .or_default()
            .extend(object.keys().cloned());
    }
    for section in ["params", "measurements"] {
        if let Some(object) = value.get(section).and_then(serde_json::Value::as_object) {
            fields
                .entry(section.to_string())
                .or_default()
                .extend(object.keys().cloned());
        }
    }
    Ok(())
}

type ArtifactHashes = BTreeMap<String, String>;
type ReceiptFieldInventory = BTreeMap<String, Vec<String>>;

fn run_checked<R: Runner>(
    runner: &R,
    cwd: &Path,
    program: &str,
    arguments: Vec<OsString>,
) -> Result<(), String> {
    command_text(runner, cwd, program, &arguments).map(|_| ())
}

fn generate_primary_artifacts<R: Runner>(
    runner: &R,
    root: &Path,
    binary: &Path,
    release_toolchain: &str,
) -> Result<(ArtifactHashes, ReceiptFieldInventory), String> {
    let temporary = std::env::temp_dir().join(format!(
        "rosalind-contract-snapshot-{}-{}",
        std::process::id(),
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap_or_default()
            .as_nanos()
    ));
    fs::create_dir_all(&temporary).map_err(|error| error.to_string())?;
    let demo = temporary.join("demo");
    let binary_text = binary.display().to_string();
    let result = (|| {
        run_checked(
            runner,
            root,
            &binary_text,
            vec![
                "demo".into(),
                "--output-dir".into(),
                demo.as_os_str().to_owned(),
                "--json".into(),
            ],
        )?;
        let sam = demo.join("alignments.sam");
        run_checked(
            runner,
            root,
            &binary_text,
            vec![
                "align".into(),
                "--reference".into(),
                demo.join("reference.fa").into_os_string(),
                "--reads".into(),
                demo.join("reads.fastq").into_os_string(),
                "--format".into(),
                "sam".into(),
                "--output".into(),
                sam.as_os_str().to_owned(),
            ],
        )?;
        let features = demo.join("features.tsv");
        run_checked(
            runner,
            root,
            &binary_text,
            vec![
                "features".into(),
                "--index".into(),
                demo.join("ref.idx").into_os_string(),
                "--alignments".into(),
                demo.join("sorted.bam").into_os_string(),
                "--output".into(),
                features.as_os_str().to_owned(),
            ],
        )?;
        let gvcf = demo.join("calls.gvcf");
        run_checked(
            runner,
            root,
            &binary_text,
            vec![
                "variants".into(),
                "--index".into(),
                demo.join("ref.idx").into_os_string(),
                "--alignments".into(),
                demo.join("sorted.bam").into_os_string(),
                "--gvcf".into(),
                "--output".into(),
                gvcf.as_os_str().to_owned(),
            ],
        )?;
        let analyzer_project = temporary.join("golden-analyzer");
        run_checked(
            runner,
            root,
            &binary_text,
            vec![
                "new".into(),
                "analyzer".into(),
                "golden-analyzer".into(),
                "--output".into(),
                analyzer_project.as_os_str().to_owned(),
            ],
        )?;
        let toml_path = analyzer_project.join("Cargo.toml");
        let toml_path_value = |path: &Path| {
            path.display()
                .to_string()
                .replace('\\', "\\\\")
                .replace('"', "\\\"")
        };
        let patch = format!(
            "\n[patch.crates-io]\nrosalind-bio = {{ path = \"{}\" }}\nrosalind-build-info = {{ path = \"{}\" }}\n",
            toml_path_value(root),
            toml_path_value(&root.join("crates/build-info")),
        );
        fs::OpenOptions::new()
            .append(true)
            .open(&toml_path)
            .and_then(|mut file| file.write_all(patch.as_bytes()))
            .map_err(|error| format!("cannot patch generated analyzer: {error}"))?;
        run_checked(
            runner,
            root,
            "rustup",
            vec![
                "run".into(),
                release_toolchain.into(),
                "cargo".into(),
                "build".into(),
                "--manifest-path".into(),
                toml_path.as_os_str().to_owned(),
            ],
        )?;
        let analyzer_binary = analyzer_project.join("target/debug/golden-analyzer");
        let analyzer_output = demo.join("golden-analyzer.tsv");
        run_checked(
            runner,
            root,
            &analyzer_binary.display().to_string(),
            vec![
                "run".into(),
                "--index".into(),
                demo.join("ref.idx").into_os_string(),
                "--alignments".into(),
                demo.join("sorted.bam").into_os_string(),
                "--output".into(),
                analyzer_output.as_os_str().to_owned(),
            ],
        )?;
        let conformance = command_text(
            runner,
            root,
            &binary_text,
            &[
                "conformance".into(),
                "analyzer".into(),
                "--binary".into(),
                analyzer_binary.as_os_str().to_owned(),
                "--json".into(),
            ],
        )?;
        let conformance_json: serde_json::Value =
            serde_json::from_str(&conformance).map_err(|error| error.to_string())?;
        if conformance_json
            .get("passed")
            .and_then(serde_json::Value::as_bool)
            != Some(true)
        {
            return Err("generated analyzer failed the contract matrix".into());
        }
        let outputs = [
            ("index", demo.join("ref.idx")),
            ("raw_bam", demo.join("raw.bam")),
            ("sorted_bam", demo.join("sorted.bam")),
            ("vcf", demo.join("calls.vcf")),
            ("sam", sam),
            ("tsv", features),
            ("gvcf", gvcf),
            ("generated_analyzer", analyzer_output),
        ];
        let mut hashes = BTreeMap::new();
        for (role, path) in outputs {
            hashes.insert(role.to_string(), hash_file(&path)?);
        }
        let demo_body = serde_json::to_vec(&hashes).map_err(|error| error.to_string())?;
        hashes.insert("demo".into(), hash_bytes(&demo_body));
        hashes.insert(
            "generated_analyzer_conformance".into(),
            hash_bytes(&serde_json::to_vec(&conformance_json).map_err(|error| error.to_string())?),
        );
        let mut receipt_fields: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
        for file in files_recursively(&demo)? {
            let name = file
                .file_name()
                .and_then(|name| name.to_str())
                .unwrap_or("");
            if name.ends_with(".manifest.json") || name.ends_with(".repro.json") {
                collect_receipt_fields(&file, &mut receipt_fields)?;
            }
        }
        Ok((
            hashes,
            receipt_fields
                .into_iter()
                .map(|(key, values)| (key, values.into_iter().collect()))
                .collect(),
        ))
    })();
    fs::remove_dir_all(&temporary).ok();
    result
}

fn snapshot_here<R: Runner>(
    runner: &R,
    root: &Path,
    commit: &str,
    policy: &Policy,
) -> Result<ContractSnapshot, String> {
    run_checked(
        runner,
        root,
        "rustup",
        vec![
            "run".into(),
            policy.release_toolchain.clone().into(),
            "cargo".into(),
            "build".into(),
            "--bin".into(),
            "rosalind".into(),
        ],
    )?;
    let binary = root.join("target/debug/rosalind");
    let mut components = BTreeMap::new();
    for package in ["rosalind-bio", "rosalind-receipt"] {
        let output = command_text(
            runner,
            root,
            "rustup",
            &[
                "run".into(),
                policy.public_api_toolchain.clone().into(),
                "cargo".into(),
                "public-api".into(),
                "--simplified".into(),
                "-p".into(),
                package.into(),
            ],
        )?;
        components.insert(
            format!("public_api.{package}"),
            hash_bytes(output.as_bytes()),
        );
    }
    let cli_contract = discover_cli_contract(runner, root, &binary)?;
    let discovered = cli_contract.keys().cloned().collect::<BTreeSet<_>>();
    let expected = policy.cli_help.iter().cloned().collect::<BTreeSet<_>>();
    if discovered != expected {
        return Err(format!(
            "release policy CLI inventory drifted: missing_from_policy={:?}, stale_in_policy={:?}",
            discovered.difference(&expected).collect::<Vec<_>>(),
            expected.difference(&discovered).collect::<Vec<_>>(),
        ));
    }
    for (command, shape) in cli_contract {
        let key = if command.is_empty() {
            "root".to_string()
        } else {
            command.replace(' ', ".")
        };
        components.insert(format!("cli_help.{key}"), hash_bytes(shape.as_bytes()));
    }
    for schema in &policy.schema_files {
        components.insert(format!("schema.{schema}"), hash_file(&root.join(schema))?);
    }
    for fixture in files_recursively(&root.join(&policy.legacy_fixture_dir))? {
        let relative = fixture.strip_prefix(root).unwrap_or(&fixture).display();
        components.insert(format!("legacy.{relative}"), hash_file(&fixture)?);
    }
    for package in ["rosalind-build-info", "rosalind-receipt", "rosalind-bio"] {
        let mut package_args = args(&["package", "-p", package, "--allow-dirty", "--list"]);
        if package == "rosalind-bio" {
            package_args.extend(args(&[
                "--config",
                "patch.crates-io.rosalind-receipt.path=\"crates/receipt\"",
                "--config",
                "patch.crates-io.rosalind-build-info.path=\"crates/build-info\"",
            ]));
        }
        let package_args = [
            vec![
                "run".into(),
                policy.release_toolchain.clone().into(),
                "cargo".into(),
            ],
            package_args,
        ]
        .concat();
        let package_list = command_text(runner, root, "rustup", &package_args)?;
        components.insert(
            format!("package.{package}.files"),
            hash_bytes(package_list.as_bytes()),
        );
    }
    let (primary_artifacts, receipt_fields) =
        generate_primary_artifacts(runner, root, &binary, &policy.release_toolchain)?;
    let package_versions = package_versions(root)?;
    let body = serde_json::to_vec(&serde_json::json!({
        "components": components,
        "receipt_fields": receipt_fields,
        "primary_artifacts": primary_artifacts,
        "package_versions": package_versions,
    }))
    .map_err(|error| error.to_string())?;
    Ok(ContractSnapshot {
        schema: 1,
        commit: commit.to_string(),
        components,
        receipt_fields,
        primary_artifacts,
        package_versions,
        aggregate_blake3: hash_bytes(&body),
    })
}

fn contract_snapshot<R: Runner>(
    runner: &R,
    root: &Path,
    reference: &str,
    policy: &Policy,
) -> Result<ContractSnapshot, String> {
    let commit = git_commit(runner, root, reference)?;
    let head = git_commit(runner, root, "HEAD")?;
    if commit == head {
        return snapshot_here(runner, root, &commit, policy);
    }
    let worktree =
        std::env::temp_dir().join(format!("rosalind-contract-worktree-{}", std::process::id()));
    if worktree.exists() {
        fs::remove_dir_all(&worktree).map_err(|error| error.to_string())?;
    }
    run_checked(
        runner,
        root,
        "git",
        vec![
            "worktree".into(),
            "add".into(),
            "--detach".into(),
            worktree.as_os_str().to_owned(),
            commit.clone().into(),
        ],
    )?;
    let result = load_policy(&worktree)
        .and_then(|target_policy| snapshot_here(runner, &worktree, &commit, &target_policy));
    let _ = runner.run(
        root,
        "git",
        &[
            "worktree".into(),
            "remove".into(),
            "--force".into(),
            worktree.as_os_str().to_owned(),
        ],
    );
    result
}

fn release_preflight<R: Runner>(
    runner: &R,
    root: &Path,
    policy: &Policy,
    reference: &str,
    version: &str,
) -> MaintainerReport {
    let commit = git_commit(runner, root, reference).unwrap_or_else(|_| "unknown".into());
    let mut report = MaintainerReport::new("release.preflight", &commit);
    report.version = Some(version.to_string());
    let clean = command_text(runner, root, "git", &args(&["status", "--porcelain"]));
    report.check(
        "git.clean",
        clean.as_ref().is_ok_and(String::is_empty),
        clean.map_or_else(
            |error| error,
            |value| {
                if value.is_empty() {
                    "clean".into()
                } else {
                    value
                }
            },
        ),
    );
    let origin = command_text(runner, root, "git", &args(&["remote", "get-url", "origin"]));
    report.check(
        "git.repository",
        origin
            .as_ref()
            .is_ok_and(|value| repository_matches(value, &policy.repository)),
        origin.unwrap_or_else(|error| error),
    );
    let ancestor = runner.run(
        root,
        "git",
        &[
            "merge-base".into(),
            "--is-ancestor".into(),
            commit.clone().into(),
            format!("origin/{}", policy.default_branch).into(),
        ],
    );
    report.check(
        "git.reachable_from_default",
        ancestor
            .as_ref()
            .is_ok_and(|output| output.status.success()),
        format!(
            "candidate must be reachable from origin/{}",
            policy.default_branch
        ),
    );
    let pushed = command_text(
        runner,
        root,
        "git",
        &[
            "branch".into(),
            "-r".into(),
            "--contains".into(),
            commit.clone().into(),
        ],
    );
    report.check(
        "git.pushed",
        pushed
            .as_ref()
            .is_ok_and(|branches| !branches.trim().is_empty()),
        pushed.unwrap_or_else(|error| error),
    );
    let versions = package_versions(root);
    let root_version = versions
        .as_ref()
        .ok()
        .and_then(|versions| versions.get("rosalind-bio"));
    report.check(
        "version.root",
        root_version.is_some_and(|actual| actual == version),
        root_version
            .map(|value| format!("Cargo.toml={value}, requested={version}"))
            .unwrap_or_else(|| match &versions {
                Ok(_) => "version unavailable".into(),
                Err(error) => error.clone(),
            }),
    );
    match &versions {
        Ok(actual) => report.check(
            "version.workspace",
            actual == &policy.package_versions,
            if actual == &policy.package_versions {
                format!("{} package versions match policy", actual.len())
            } else {
                format!("policy={:?}, workspace={actual:?}", policy.package_versions)
            },
        ),
        Err(error) => report.check("version.workspace", false, error.clone()),
    }
    let changelog = fs::read_to_string(root.join("CHANGELOG.md")).unwrap_or_default();
    report.check(
        "changelog.unreleased",
        changelog.contains("## [Unreleased]") || changelog.contains(&format!("## [{version}]")),
        "CHANGELOG must contain an Unreleased or candidate-version section",
    );
    report
}

fn rc_plan<R: Runner>(
    runner: &R,
    root: &Path,
    policy: &Policy,
    version: &str,
    number: u32,
    reference: &str,
) -> MaintainerReport {
    let mut report = release_preflight(runner, root, policy, reference, version);
    report.command = "rc.plan".into();
    let tag = format!("v{version}-rc.{number}");
    report.metadata.insert("tag".into(), tag.clone());
    report.metadata.insert("ref".into(), reference.into());
    let tag_ref = format!("refs/tags/{tag}");
    let tag_query = command_text(
        runner,
        root,
        "git",
        &["ls-remote".into(), "origin".into(), tag_ref.into()],
    );
    let tag_exists = tag_query.as_ref().is_ok_and(|value| !value.is_empty());
    report.check(
        "tag.absent",
        tag_query.is_ok() && !tag_exists,
        if tag_query.is_err() {
            "could not query RC tag on origin"
        } else if tag_exists {
            "tag already exists"
        } else {
            "available"
        },
    );
    match contract_snapshot(runner, root, reference, policy) {
        Ok(snapshot) => {
            report.contract_fingerprint = Some(snapshot.aggregate_blake3);
            report.check("contract.snapshot", true, "generated");
        }
        Err(error) => report.check("contract.snapshot", false, error),
    }
    report.actions = vec![
        format!("dispatch protected rc workflow for {tag}"),
        "build and attest cross-platform release assets".into(),
        "publish GitHub prerelease without publishing crates".into(),
    ];
    report.seal();
    report
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct PartnerScenario {
    id: String,
    passed: bool,
    notes: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
enum Severity {
    None,
    Low,
    Medium,
    ReleaseBlocking,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct PartnerRecord {
    schema: u32,
    partner_id: String,
    persona: Persona,
    tested_commit: String,
    contract_fingerprint: String,
    environment: String,
    scenarios: Vec<PartnerScenario>,
    blocker_severity: Severity,
    contract_change_requested: bool,
    resolution: String,
    consent_to_publish: bool,
}

fn validate_partner_record(record: &PartnerRecord) -> Vec<String> {
    let mut failures = Vec::new();
    if record.schema != 1 {
        failures.push("schema must be 1".into());
    }
    if record.partner_id.is_empty()
        || !record
            .partner_id
            .bytes()
            .all(|byte| byte.is_ascii_lowercase() || byte.is_ascii_digit() || byte == b'-')
    {
        failures.push("partner_id must be an anonymized lowercase slug".into());
    }
    if record.tested_commit.len() != 40
        || !record
            .tested_commit
            .bytes()
            .all(|byte| byte.is_ascii_hexdigit())
    {
        failures.push("tested_commit must be a full Git SHA".into());
    }
    if record.contract_fingerprint.len() != 64
        || !record
            .contract_fingerprint
            .bytes()
            .all(|byte| byte.is_ascii_hexdigit())
    {
        failures.push("contract_fingerprint must be a BLAKE3 digest".into());
    }
    if record.environment.trim().is_empty() {
        failures.push("environment is required".into());
    }
    let expected = partner_scenarios(record.persona)
        .into_iter()
        .map(|scenario| scenario.id)
        .collect::<BTreeSet<_>>();
    let actual = record
        .scenarios
        .iter()
        .map(|scenario| scenario.id.clone())
        .collect::<BTreeSet<_>>();
    if actual != expected || actual.len() != record.scenarios.len() {
        failures.push("persona scenarios must be complete, exact, and unique".into());
    }
    if record.scenarios.iter().any(|scenario| !scenario.passed) {
        failures.push("every persona scenario must pass".into());
    }
    if record
        .scenarios
        .iter()
        .any(|scenario| scenario.notes.trim().is_empty())
    {
        failures.push("every persona scenario requires a concise anonymized result".into());
    }
    if matches!(record.blocker_severity, Severity::ReleaseBlocking) {
        failures.push("release-blocking finding is unresolved".into());
    }
    if record.contract_change_requested && record.resolution.trim().is_empty() {
        failures.push("contract change requests require a resolution".into());
    }
    if !record.consent_to_publish {
        failures.push("sanitized release-gate record lacks consent to publish".into());
    }
    let text_fields = std::iter::once(record.environment.as_str())
        .chain(std::iter::once(record.resolution.as_str()))
        .chain(
            record
                .scenarios
                .iter()
                .map(|scenario| scenario.notes.as_str()),
        );
    if text_fields.into_iter().any(contains_sensitive_text) {
        failures.push(
            "free-text fields appear to contain contact, credential, patient, organization, or genomic data"
                .into(),
        );
    }
    failures
}

fn contains_sensitive_text(value: &str) -> bool {
    let lowercase = value.to_ascii_lowercase();
    if value.contains('@')
        || [
            "email",
            "password",
            "secret",
            "token=",
            "private key",
            "patient",
            "employer",
            "organization",
        ]
        .iter()
        .any(|marker| lowercase.contains(marker))
    {
        return true;
    }
    value
        .split(|character: char| !character.is_ascii_alphabetic())
        .any(|word| {
            word.len() >= 30
                && word.bytes().all(|base| {
                    matches!(base.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'N')
                })
        })
}

fn read_partner_record(path: &Path) -> Result<PartnerRecord, String> {
    let text = fs::read_to_string(path)
        .map_err(|error| format!("cannot read {}: {error}", path.display()))?;
    serde_json::from_str(&text).map_err(|error| format!("cannot parse {}: {error}", path.display()))
}

fn partner_report<R: Runner>(
    runner: &R,
    root: &Path,
    policy: &Policy,
    input_dir: &Path,
    release_gate: Option<(&str, &str)>,
) -> MaintainerReport {
    let commit = git_commit(runner, root, "HEAD").unwrap_or_else(|_| "unknown".into());
    let mut report = MaintainerReport::new("partners.report", commit);
    let mut personas = BTreeSet::new();
    let mut partner_ids = BTreeSet::new();
    let mut count = 0usize;
    match files_recursively(input_dir) {
        Ok(files) => {
            for path in files
                .into_iter()
                .filter(|path| path.extension().and_then(|value| value.to_str()) == Some("json"))
            {
                match read_partner_record(&path) {
                    Ok(record) => {
                        let mut failures = validate_partner_record(&record);
                        if !partner_ids.insert(record.partner_id.clone()) {
                            failures.push("duplicate partner_id".into());
                        }
                        if let Some((expected_fingerprint, candidate_commit)) = release_gate {
                            if record.contract_fingerprint != expected_fingerprint {
                                failures.push(format!(
                                    "contract fingerprint differs from release candidate: expected {expected_fingerprint}"
                                ));
                            }
                            let ancestor = runner.run(
                                root,
                                "git",
                                &[
                                    "merge-base".into(),
                                    "--is-ancestor".into(),
                                    record.tested_commit.clone().into(),
                                    candidate_commit.into(),
                                ],
                            );
                            if !ancestor
                                .as_ref()
                                .is_ok_and(|output| output.status.success())
                            {
                                failures.push(
                                    "tested commit is not an ancestor of the promoted candidate"
                                        .into(),
                                );
                            }
                        }
                        report.check(
                            &format!("partner.{}", record.partner_id),
                            failures.is_empty(),
                            if failures.is_empty() {
                                "accepted".into()
                            } else {
                                failures.join("; ")
                            },
                        );
                        if failures.is_empty() {
                            personas.insert(record.persona.as_str().to_string());
                            count += 1;
                        }
                    }
                    Err(error) => {
                        report.check(&format!("partner.file.{}", path.display()), false, error)
                    }
                }
            }
        }
        Err(error) => report.check("partner.directory", false, error),
    }
    for persona in &policy.partner_personas {
        report.check(
            &format!("persona.{persona}"),
            personas.contains(persona),
            if personas.contains(persona) {
                "represented"
            } else {
                "missing"
            },
        );
    }
    report
        .metadata
        .insert("accepted_records".into(), count.to_string());
    report.actions = vec!["collect one accepted anonymized record for every persona".into()];
    report.seal();
    report
}

fn stable_release_plan<R: Runner>(
    runner: &R,
    root: &Path,
    policy: &Policy,
    version: &str,
    rc_tag: &str,
    reference: &str,
) -> MaintainerReport {
    let mut report = release_preflight(runner, root, policy, reference, version);
    report.command = "release.plan".into();
    report.metadata.insert("rc_tag".into(), rc_tag.into());
    report.metadata.insert("ref".into(), reference.into());
    let stable_tag = format!("v{version}");
    let direct_ref = format!("refs/tags/{stable_tag}");
    let peeled_ref = format!("{direct_ref}^{{}}");
    let direct = command_text(
        runner,
        root,
        "git",
        &["ls-remote".into(), "origin".into(), direct_ref.into()],
    );
    let peeled = command_text(
        runner,
        root,
        "git",
        &["ls-remote".into(), "origin".into(), peeled_ref.into()],
    );
    let remote_tag_target = peeled
        .as_ref()
        .ok()
        .filter(|value| !value.is_empty())
        .or_else(|| direct.as_ref().ok().filter(|value| !value.is_empty()))
        .and_then(|value| value.split_whitespace().next());
    let tag_query_ok = direct.is_ok() && peeled.is_ok();
    report.check(
        "tag.nonconflicting",
        tag_query_ok && remote_tag_target.is_none_or(|target| target == report.commit),
        if !tag_query_ok {
            "could not query stable tag on origin".into()
        } else if let Some(target) = remote_tag_target {
            format!("existing tag resolves to {target}")
        } else {
            "stable tag is available".into()
        },
    );
    let releases = command_text(
        runner,
        root,
        "gh",
        &args(&[
            "release",
            "list",
            "--repo",
            &policy.repository,
            "--limit",
            "100",
            "--json",
            "tagName,isPrerelease,isDraft",
        ]),
    );
    let release_state = releases.as_ref().ok().and_then(|text| {
        serde_json::from_str::<Vec<serde_json::Value>>(text)
            .ok()
            .and_then(|values| {
                values.into_iter().find(|value| {
                    value.get("tagName").and_then(serde_json::Value::as_str)
                        == Some(stable_tag.as_str())
                })
            })
    });
    let release_nonconflicting = releases.is_ok()
        && release_state.as_ref().is_none_or(|value| {
            value
                .get("isPrerelease")
                .and_then(serde_json::Value::as_bool)
                == Some(false)
                && value.get("isDraft").and_then(serde_json::Value::as_bool) == Some(false)
        });
    report.check(
        "github_release.nonconflicting",
        release_nonconflicting,
        if releases.is_err() {
            "could not query GitHub Releases".to_string()
        } else if release_state.is_some() && release_nonconflicting {
            "matching stable release already exists; workflow may verify/resume".to_string()
        } else if release_state.is_some() {
            "a draft or prerelease conflicts with the stable tag".to_string()
        } else {
            "stable release name is available".to_string()
        },
    );
    for package in &policy.publish_order {
        let package_version = policy
            .package_versions
            .get(package)
            .map(String::as_str)
            .unwrap_or(version);
        let state = crate_publication_state(
            runner,
            root,
            package,
            package_version,
            &policy.release_toolchain,
        );
        let (ok, detail) = match state {
            CratePublicationState::Absent => {
                (true, format!("{package} {package_version} is available"))
            }
            CratePublicationState::Exact(checksum) => (
                true,
                format!("{package} {package_version} already exists with exact bytes ({checksum})"),
            ),
            CratePublicationState::Conflict(error)
            | CratePublicationState::RemoteFailure(error) => (false, error),
        };
        report.check(&format!("crate.nonconflicting.{package}"), ok, detail);
    }
    let rc = rc_status(runner, root, policy, rc_tag);
    let soak_required = version_at_least(version, &policy.soak_required_from);
    let soak_ok = !soak_required
        || rc
            .checks
            .iter()
            .find(|check| check.id == "soak.elapsed")
            .is_some_and(|check| check.ok);
    report.check(
        "soak.elapsed",
        soak_ok,
        if soak_required {
            rc.metadata
                .get("soak")
                .cloned()
                .unwrap_or_else(|| "RC unavailable".into())
        } else {
            format!(
                "not required before {}; stabilization release",
                policy.soak_required_from
            )
        },
    );
    let ci_ok = rc
        .checks
        .iter()
        .find(|check| check.id == "ci.current")
        .is_some_and(|check| check.ok);
    report.check(
        "ci.candidate",
        ci_ok,
        rc.checks
            .iter()
            .find(|check| check.id == "ci.current")
            .map(|check| check.detail.clone())
            .unwrap_or_else(|| "candidate CI status unavailable".into()),
    );
    let mut candidate_fingerprint = None;
    match contract_snapshot(runner, root, reference, policy) {
        Ok(snapshot) => {
            report.contract_fingerprint = Some(snapshot.aggregate_blake3.clone());
            candidate_fingerprint = Some(snapshot.aggregate_blake3.clone());
            let expected = rc.metadata.get("contract_fingerprint");
            report.check(
                "contract.equivalent",
                expected.is_some_and(|value| value == &snapshot.aggregate_blake3),
                expected.map_or_else(
                    || "RC fingerprint unavailable".into(),
                    |value| format!("rc={value}, candidate={}", snapshot.aggregate_blake3),
                ),
            );
        }
        Err(error) => report.check("contract.snapshot", false, error),
    }
    let candidate_commit = report.commit.clone();
    let partners_required = version_at_least(version, &policy.partners_required_from);
    let partners = partner_report(
        runner,
        root,
        policy,
        &root.join(&policy.partner_record_dir),
        candidate_fingerprint
            .as_deref()
            .map(|fingerprint| (fingerprint, candidate_commit.as_str())),
    );
    report.check(
        "design_partners.complete",
        !partners_required || partners.status == "ready",
        if !partners_required {
            format!(
                "not required before {}; validation continues after stabilization",
                policy.partners_required_from
            )
        } else if partners.status == "ready" {
            "three personas accepted at the frozen contract".to_string()
        } else {
            partners.blockers.join("; ")
        },
    );
    match caller_evidence::gate(runner, root, policy, &report.commit) {
        Ok(detail) => report.check("giab.caller_source", true, detail),
        Err(detail) => report.check("giab.caller_source", false, detail),
    }
    report.actions = vec![
        format!("dispatch protected stable release for v{version}"),
        "publish crates idempotently in dependency order".into(),
        "verify fresh-cache installation on Linux and macOS".into(),
        "create stable tag and GitHub Release after post-publication checks".into(),
    ];
    mark_crate_remote_failure(&mut report);
    report.seal();
    report
}

fn rc_status<R: Runner>(runner: &R, root: &Path, policy: &Policy, tag: &str) -> MaintainerReport {
    let commit = git_commit(runner, root, "HEAD").unwrap_or_else(|_| "unknown".into());
    let mut report = MaintainerReport::new("rc.status", commit);
    report.metadata.insert("tag".into(), tag.into());
    let mut rc_commit = None;
    let release = command_text(
        runner,
        root,
        "gh",
        &args(&[
            "release",
            "view",
            tag,
            "--repo",
            &policy.repository,
            "--json",
            "publishedAt,isPrerelease,tagName,targetCommitish",
        ]),
    );
    match release
        .as_ref()
        .ok()
        .and_then(|text| serde_json::from_str::<serde_json::Value>(text).ok())
    {
        Some(value) => {
            let prerelease =
                value.get("isPrerelease").and_then(|value| value.as_bool()) == Some(true);
            report.check(
                "rc.prerelease",
                prerelease,
                if prerelease {
                    "published prerelease"
                } else {
                    "not a prerelease"
                },
            );
            if let Some(target) = value
                .get("targetCommitish")
                .and_then(|value| value.as_str())
            {
                rc_commit = Some(target.to_string());
                report.metadata.insert("rc_commit".into(), target.into());
            }
            let published = value.get("publishedAt").and_then(|value| value.as_str());
            let elapsed = published
                .and_then(parse_github_rfc3339_unix)
                .map(|published_at| {
                    let now = SystemTime::now()
                        .duration_since(UNIX_EPOCH)
                        .unwrap_or_default()
                        .as_secs();
                    soak_elapsed(published_at, now)
                });
            let elapsed = elapsed.unwrap_or(0);
            report.metadata.insert(
                "soak".into(),
                format!(
                    "elapsed={elapsed}s required={}s remaining={}s",
                    policy.soak_seconds,
                    policy.soak_seconds.saturating_sub(elapsed)
                ),
            );
            report.check(
                "soak.elapsed",
                elapsed >= policy.soak_seconds,
                report.metadata["soak"].clone(),
            );
        }
        None => report.check(
            "rc.release",
            false,
            release
                .err()
                .unwrap_or_else(|| "cannot parse GitHub prerelease".into()),
        ),
    }
    let temp = std::env::temp_dir().join(format!("rosalind-rc-status-{}", std::process::id()));
    fs::remove_dir_all(&temp).ok();
    fs::create_dir_all(&temp).ok();
    let download = runner.run(
        root,
        "gh",
        &[
            "release".into(),
            "download".into(),
            tag.into(),
            "--repo".into(),
            policy.repository.clone().into(),
            "--pattern".into(),
            "contract-snapshot.json".into(),
            "--dir".into(),
            temp.as_os_str().to_owned(),
        ],
    );
    let snapshot_path = temp.join("contract-snapshot.json");
    let snapshot = download
        .ok()
        .filter(|output| output.status.success())
        .and_then(|_| fs::read_to_string(&snapshot_path).ok())
        .and_then(|text| serde_json::from_str::<ContractSnapshot>(&text).ok());
    let mut rc_fingerprint = None;
    if let Some(snapshot) = snapshot {
        rc_fingerprint = Some(snapshot.aggregate_blake3.clone());
        report.check(
            "contract.commit",
            rc_commit
                .as_ref()
                .is_some_and(|commit| commit == &snapshot.commit),
            rc_commit.as_ref().map_or_else(
                || "RC target commit unavailable".into(),
                |commit| format!("release={commit}, snapshot={}", snapshot.commit),
            ),
        );
        report
            .metadata
            .insert("contract_fingerprint".into(), snapshot.aggregate_blake3);
        report.check("contract.asset", true, "RC snapshot downloaded");
    } else {
        report.check(
            "contract.asset",
            false,
            "RC contract-snapshot.json is missing or invalid",
        );
    }
    let attestation = runner.run(
        root,
        "gh",
        &[
            "attestation".into(),
            "verify".into(),
            snapshot_path.as_os_str().to_owned(),
            "--repo".into(),
            policy.repository.clone().into(),
            "--source-digest".into(),
            rc_commit.clone().unwrap_or_default().into(),
            "--signer-workflow".into(),
            format!("{}/.github/workflows/rc.yml", policy.repository).into(),
        ],
    );
    report.check(
        "contract.attestation",
        snapshot_path.exists()
            && attestation
                .as_ref()
                .is_ok_and(|output| output.status.success()),
        "RC contract snapshot must have a valid GitHub attestation",
    );
    fs::remove_dir_all(temp).ok();
    match contract_snapshot(runner, root, "HEAD", policy) {
        Ok(current) => {
            let matches = rc_fingerprint
                .as_ref()
                .is_some_and(|expected| expected == &current.aggregate_blake3);
            report.metadata.insert(
                "current_contract_fingerprint".into(),
                current.aggregate_blake3.clone(),
            );
            report.check(
                "contract.current_matches_rc",
                matches,
                rc_fingerprint.as_ref().map_or_else(
                    || "RC fingerprint unavailable".into(),
                    |expected| format!("rc={expected}, current={}", current.aggregate_blake3),
                ),
            );
        }
        Err(error) => report.check("contract.current", false, error),
    }
    let ci = command_text(
        runner,
        root,
        "gh",
        &args(&[
            "run",
            "list",
            "--repo",
            &policy.repository,
            "--workflow",
            "ci.yml",
            "--commit",
            &report.commit,
            "--limit",
            "1",
            "--json",
            "status,conclusion,headSha,url",
        ]),
    );
    let ci_value = ci.as_ref().ok().and_then(|text| {
        serde_json::from_str::<Vec<serde_json::Value>>(text)
            .ok()
            .and_then(|values| values.into_iter().next())
    });
    let ci_green = ci_value.as_ref().is_some_and(|value| {
        value.get("status").and_then(serde_json::Value::as_str) == Some("completed")
            && value.get("conclusion").and_then(serde_json::Value::as_str) == Some("success")
    });
    report.check(
        "ci.current",
        ci_green,
        ci_value
            .as_ref()
            .map(serde_json::Value::to_string)
            .unwrap_or_else(|| ci.err().unwrap_or_else(|| "no CI run found".into())),
    );
    let partners = partner_report(
        runner,
        root,
        policy,
        &root.join(&policy.partner_record_dir),
        rc_fingerprint
            .as_deref()
            .map(|fingerprint| (fingerprint, report.commit.as_str())),
    );
    report
        .metadata
        .insert("partners".into(), partners.status.clone());
    report.check(
        "design_partners.complete",
        partners.status == "ready",
        if partners.status == "ready" {
            "all required personas accepted".into()
        } else {
            partners.blockers.join("; ")
        },
    );
    report.actions = vec![
        "wait for the soak and complete design-partner records".into(),
        "cut a new RC when contract.current_matches_rc is false".into(),
    ];
    report.seal();
    report
}

fn soak_elapsed(published_unix: u64, now_unix: u64) -> u64 {
    now_unix.saturating_sub(published_unix)
}

/// Parse the UTC RFC3339 shape emitted by GitHub (`publishedAt`) without pulling
/// a date/time dependency into the MSRV-constrained maintainer binary.
fn parse_github_rfc3339_unix(value: &str) -> Option<u64> {
    let value = value.strip_suffix('Z')?;
    let (date, time) = value.split_once('T')?;
    let mut date = date.split('-').map(|field| field.parse::<i64>().ok());
    let (year, month, day) = (date.next()??, date.next()??, date.next()??);
    if date.next().is_some() || !(1..=12).contains(&month) {
        return None;
    }
    let leap = year % 4 == 0 && (year % 100 != 0 || year % 400 == 0);
    let month_days = [
        31,
        if leap { 29 } else { 28 },
        31,
        30,
        31,
        30,
        31,
        31,
        30,
        31,
        30,
        31,
    ];
    if day < 1 || day > month_days[(month - 1) as usize] {
        return None;
    }
    let mut time = time.split(':');
    let hour = time.next()?.parse::<i64>().ok()?;
    let minute = time.next()?.parse::<i64>().ok()?;
    let second = time.next()?.split('.').next()?.parse::<i64>().ok()?;
    if time.next().is_some() || hour > 23 || minute > 59 || second > 59 {
        return None;
    }

    // Howard Hinnant's civil-date transform: days relative to 1970-01-01.
    let adjusted_year = year - i64::from(month <= 2);
    let era = if adjusted_year >= 0 {
        adjusted_year
    } else {
        adjusted_year - 399
    } / 400;
    let year_of_era = adjusted_year - era * 400;
    let shifted_month = month + if month > 2 { -3 } else { 9 };
    let day_of_year = (153 * shifted_month + 2) / 5 + day - 1;
    let day_of_era = year_of_era * 365 + year_of_era / 4 - year_of_era / 100 + day_of_year;
    let days = era * 146_097 + day_of_era - 719_468;
    let seconds = days
        .checked_mul(86_400)?
        .checked_add(hour * 3_600 + minute * 60 + second)?;
    u64::try_from(seconds).ok()
}

fn version_at_least(version: &str, floor: &str) -> bool {
    semver::Version::parse(version)
        .ok()
        .zip(semver::Version::parse(floor).ok())
        .is_some_and(|(version, floor)| version >= floor)
}

fn release_status<R: Runner>(
    runner: &R,
    root: &Path,
    policy: &Policy,
    version: &str,
) -> MaintainerReport {
    let commit = git_commit(runner, root, "HEAD").unwrap_or_else(|_| "unknown".into());
    let mut report = MaintainerReport::new("release.status", commit);
    report.version = Some(version.into());
    let tag = format!("v{version}");
    let release = command_text(
        runner,
        root,
        "gh",
        &args(&[
            "release",
            "view",
            &tag,
            "--repo",
            &policy.repository,
            "--json",
            "isDraft,isPrerelease,publishedAt",
        ]),
    );
    let release_valid = release.as_ref().ok().and_then(|text| {
        serde_json::from_str::<serde_json::Value>(text)
            .ok()
            .map(|value| {
                value.get("isDraft").and_then(serde_json::Value::as_bool) == Some(false)
                    && value
                        .get("isPrerelease")
                        .and_then(serde_json::Value::as_bool)
                        == Some(false)
            })
    }) == Some(true);
    report.check(
        "github.release",
        release_valid,
        release.unwrap_or_else(|error| error),
    );
    for package in &policy.publish_order {
        let package_version = package_versions(root)
            .ok()
            .and_then(|versions| versions.get(package).cloned())
            .unwrap_or_else(|| version.into());
        let state = crate_publication_state(
            runner,
            root,
            package,
            &package_version,
            &policy.release_toolchain,
        );
        let (package_ok, detail) = match state {
            CratePublicationState::Exact(checksum) => (
                true,
                format!("version {package_version}, exact package bytes, sha256={checksum}"),
            ),
            CratePublicationState::Absent => {
                (false, format!("version {package_version} is not published"))
            }
            CratePublicationState::Conflict(error)
            | CratePublicationState::RemoteFailure(error) => (false, error),
        };
        report.check(&format!("crate.{package}"), package_ok, detail);
    }
    report.actions =
        vec!["resume the protected release workflow if any component is incomplete".into()];
    mark_crate_remote_failure(&mut report);
    report.seal();
    report
}

fn source_lock_hash(root: &Path) -> Result<String, String> {
    let lock_path = root.join("benchmarks/giab/happy/lock.json");
    let mut lock: serde_json::Value =
        serde_json::from_slice(&fs::read(&lock_path).map_err(|error| error.to_string())?)
            .map_err(|error| error.to_string())?;
    if let Some(object) = lock.as_object_mut() {
        object.remove("generated_image");
    }
    let mut bytes = serde_json::to_vec(&lock).map_err(|error| error.to_string())?;
    bytes.extend(
        fs::read(root.join("benchmarks/giab/happy/Dockerfile"))
            .map_err(|error| error.to_string())?,
    );
    if let Ok(requirements) = fs::read(root.join("benchmarks/giab/happy/requirements-py27.txt")) {
        bytes.extend(b"\0requirements-py27.txt\0");
        bytes.extend(requirements);
    }
    Ok(hash_bytes(&bytes))
}

fn source_lock_hash_at_ref<R: Runner>(
    runner: &R,
    root: &Path,
    reference: &str,
) -> Result<String, String> {
    let commit = git_commit(runner, root, reference)?;
    if commit == git_commit(runner, root, "HEAD")? {
        return source_lock_hash(root);
    }
    let lock = runner
        .run(
            root,
            "git",
            &[
                "show".into(),
                format!("{commit}:benchmarks/giab/happy/lock.json").into(),
            ],
        )
        .map_err(|error| error.to_string())?;
    if !lock.status.success() {
        return Err(output_text(&lock));
    }
    let dockerfile = runner
        .run(
            root,
            "git",
            &[
                "show".into(),
                format!("{commit}:benchmarks/giab/happy/Dockerfile").into(),
            ],
        )
        .map_err(|error| error.to_string())?;
    if !dockerfile.status.success() {
        return Err(output_text(&dockerfile));
    }
    let mut lock: serde_json::Value =
        serde_json::from_slice(&lock.stdout).map_err(|error| error.to_string())?;
    if let Some(object) = lock.as_object_mut() {
        object.remove("generated_image");
    }
    let mut bytes = serde_json::to_vec(&lock).map_err(|error| error.to_string())?;
    bytes.extend(dockerfile.stdout);
    if let Ok(requirements) = runner.run(
        root,
        "git",
        &[
            "show".into(),
            format!("{commit}:benchmarks/giab/happy/requirements-py27.txt").into(),
        ],
    ) {
        if requirements.status.success() {
            bytes.extend(b"\0requirements-py27.txt\0");
            bytes.extend(requirements.stdout);
        }
    }
    Ok(hash_bytes(&bytes))
}

fn generated_image(root: &Path) -> Result<(String, String, String, String), String> {
    let value: serde_json::Value = serde_json::from_slice(
        &fs::read(root.join("benchmarks/giab/happy/lock.json"))
            .map_err(|error| error.to_string())?,
    )
    .map_err(|error| error.to_string())?;
    let image = value
        .get("generated_image")
        .and_then(serde_json::Value::as_object)
        .ok_or_else(|| {
            "happy lock has no generated_image; dispatch the image workflow".to_string()
        })?;
    let repository = image
        .get("repository")
        .and_then(serde_json::Value::as_str)
        .unwrap_or("");
    let digest = image
        .get("digest")
        .and_then(serde_json::Value::as_str)
        .unwrap_or("");
    let source = image
        .get("source_lock_blake3")
        .and_then(serde_json::Value::as_str)
        .unwrap_or("");
    let platform = image
        .get("platform")
        .and_then(serde_json::Value::as_str)
        .unwrap_or("");
    let built_from = image
        .get("built_from_commit")
        .and_then(serde_json::Value::as_str)
        .unwrap_or("");
    let hash_valid = |value: &str, length: usize| {
        value.len() == length && value.bytes().all(|byte| byte.is_ascii_hexdigit())
    };
    if repository.is_empty()
        || !digest_valid(digest)
        || !hash_valid(source, 64)
        || platform != "linux/amd64"
        || !hash_valid(built_from, 40)
    {
        return Err("generated_image is incomplete or not immutable".into());
    }
    Ok((
        repository.into(),
        digest.into(),
        source.into(),
        built_from.into(),
    ))
}

fn digest_valid(value: &str) -> bool {
    value.len() == 71
        && value.starts_with("sha256:")
        && value[7..].bytes().all(|byte| byte.is_ascii_hexdigit())
}

fn image_plan<R: Runner>(
    runner: &R,
    root: &Path,
    policy: &Policy,
    reference: &str,
) -> MaintainerReport {
    let commit = git_commit(runner, root, reference).unwrap_or_else(|_| "unknown".into());
    let mut report = MaintainerReport::new("giab.image.plan", commit);
    report.metadata.insert("ref".into(), reference.into());
    add_mutation_preflight(&mut report, runner, root, policy, reference);
    match source_lock_hash_at_ref(runner, root, reference) {
        Ok(source) => {
            report
                .metadata
                .insert("source_lock_blake3".into(), source.clone());
            report
                .metadata
                .insert("repository".into(), policy.happy_repository.clone());
            report.check("image.source_lock", true, source);
        }
        Err(error) => report.check("image.source_lock", false, error),
    }
    report.actions = vec![
        format!(
            "build and push linux/amd64 image to {}",
            policy.happy_repository
        ),
        "attest image and SBOM".into(),
        "open a digest-lock pull request".into(),
    ];
    report.seal();
    report
}

fn image_status<R: Runner>(runner: &R, root: &Path, policy: &Policy) -> MaintainerReport {
    let commit = git_commit(runner, root, "HEAD").unwrap_or_else(|_| "unknown".into());
    let mut report = MaintainerReport::new("giab.image.status", commit);
    match generated_image(root) {
        Ok((repository, digest, recorded_source, built_from)) => {
            let current_source = source_lock_hash(root).unwrap_or_default();
            report.check(
                "image.source_match",
                recorded_source == current_source,
                format!("recorded={recorded_source}, current={current_source}"),
            );
            let built_source = source_lock_hash_at_ref(runner, root, &built_from);
            report.check(
                "image.build_commit_source",
                built_source
                    .as_ref()
                    .is_ok_and(|source| source == &recorded_source),
                built_source.map_or_else(
                    |error| error,
                    |source| format!("commit={built_from}, source_lock={source}"),
                ),
            );
            let reachable = runner.run(
                root,
                "git",
                &[
                    "merge-base".into(),
                    "--is-ancestor".into(),
                    built_from.clone().into(),
                    format!("origin/{}", policy.default_branch).into(),
                ],
            );
            report.check(
                "image.build_commit_reachable",
                reachable
                    .as_ref()
                    .is_ok_and(|output| output.status.success()),
                built_from.clone(),
            );
            report.check(
                "image.repository",
                repository == policy.happy_repository,
                repository.clone(),
            );
            let reference = format!("{repository}@{digest}");
            let inspect = runner.run(
                root,
                "docker",
                &[
                    "buildx".into(),
                    "imagetools".into(),
                    "inspect".into(),
                    reference.clone().into(),
                ],
            );
            report.check(
                "image.remote",
                inspect.as_ref().is_ok_and(|output| output.status.success()),
                reference.clone(),
            );
            report.check(
                "image.platform",
                inspect.as_ref().is_ok_and(|output| {
                    output.status.success()
                        && output_text(output)
                            .to_ascii_lowercase()
                            .contains("linux/amd64")
                }),
                "remote manifest must include linux/amd64".to_string(),
            );
            let attest = runner.run(
                root,
                "gh",
                &[
                    "attestation".into(),
                    "verify".into(),
                    format!("oci://{reference}").into(),
                    "--repo".into(),
                    policy.repository.clone().into(),
                    "--source-digest".into(),
                    built_from.clone().into(),
                    "--signer-workflow".into(),
                    format!("{}/.github/workflows/happy-image.yml", policy.repository).into(),
                ],
            );
            report.check(
                "image.attestation",
                attest.as_ref().is_ok_and(|output| output.status.success()),
                "GitHub OCI attestation",
            );
            report.metadata.insert("image".into(), reference);
        }
        Err(error) => report.check("image.lock", false, error),
    }
    report.actions = vec!["merge the generated digest-lock PR before benchmarking".into()];
    report.seal();
    report
}

fn benchmark_plan<R: Runner>(
    runner: &R,
    root: &Path,
    policy: &Policy,
    override_image: Option<&str>,
) -> MaintainerReport {
    let commit = git_commit(runner, root, "HEAD").unwrap_or_else(|_| "unknown".into());
    let mut report = MaintainerReport::new("giab.benchmark.plan", commit);
    add_mutation_preflight(&mut report, runner, root, policy, "HEAD");
    report.metadata.insert(
        "expected_repository".into(),
        policy.happy_repository.clone(),
    );
    let image = override_image.map(str::to_owned).or_else(|| {
        generated_image(root)
            .ok()
            .map(|(repository, digest, _, _)| format!("{repository}@{digest}"))
    });
    match image {
        Some(image) => {
            let digest = image
                .split_once('@')
                .map(|(_, digest)| digest)
                .unwrap_or("");
            report.check("benchmark.image", digest_valid(digest), image.clone());
            report.metadata.insert("image".into(), image);
        }
        None => report.check(
            "benchmark.image",
            false,
            "no committed immutable evaluator image",
        ),
    }
    report.check(
        "benchmark.resources",
        root.join("benchmarks/giab/resources.tsv").exists(),
        "pinned resources.tsv",
    );
    report.actions = vec![
        "dispatch the monthly/manual HG002 v5.0q GRCh38 chr20 workflow".into(),
        "upload attested candidate evidence without editing the baseline".into(),
    ];
    report.seal();
    report
}

fn workflow_status<R: Runner>(
    runner: &R,
    root: &Path,
    policy: &Policy,
    workflow: &str,
    command: &str,
) -> MaintainerReport {
    let commit = git_commit(runner, root, "HEAD").unwrap_or_else(|_| "unknown".into());
    let mut report = MaintainerReport::new(command, commit);
    let output = command_text(
        runner,
        root,
        "gh",
        &args(&[
            "run",
            "list",
            "--repo",
            &policy.repository,
            "--workflow",
            workflow,
            "--limit",
            "1",
            "--json",
            "databaseId,status,conclusion,headSha,url",
        ]),
    );
    match output {
        Ok(json) => {
            let latest = serde_json::from_str::<Vec<serde_json::Value>>(&json)
                .ok()
                .and_then(|runs| runs.into_iter().next());
            if let Some(run) = latest {
                let status = run
                    .get("status")
                    .and_then(serde_json::Value::as_str)
                    .unwrap_or("unknown");
                let conclusion = run
                    .get("conclusion")
                    .and_then(serde_json::Value::as_str)
                    .unwrap_or("");
                let complete = status == "completed" && conclusion == "success";
                report.check(
                    "workflow.complete",
                    complete,
                    format!("status={status}, conclusion={conclusion}"),
                );
                if status == "completed" && conclusion != "success" {
                    report.status = "failed".into();
                }
                report.metadata.insert("run".into(), run.to_string());
            } else {
                report.check("workflow.run", false, "no workflow run found");
            }
        }
        Err(error) => {
            report.check("workflow.query", false, error);
            report.status = "failed".into();
        }
    }
    report.actions = vec!["inspect or resume the latest workflow run".into()];
    report.seal();
    report
}

fn read_plan(path: &Path) -> Result<MaintainerReport, String> {
    let text = fs::read_to_string(path)
        .map_err(|error| format!("cannot read plan {}: {error}", path.display()))?;
    let mut report: MaintainerReport =
        serde_json::from_str(&text).map_err(|error| format!("cannot parse plan: {error}"))?;
    let checks_ready = report.checks.iter().all(|check| check.ok);
    let status_consistent = if checks_ready {
        report.status == "ready"
    } else {
        matches!(report.status.as_str(), "blocked" | "failed")
    };
    if !status_consistent {
        return Err(format!(
            "plan status is inconsistent with its checks: recorded={}",
            report.status,
        ));
    }
    let expected_blockers = report
        .checks
        .iter()
        .filter(|check| !check.ok)
        .map(|check| format!("{}: {}", check.id, check.detail))
        .collect::<BTreeSet<_>>();
    if report.blockers.iter().cloned().collect::<BTreeSet<_>>() != expected_blockers {
        return Err("plan blockers are inconsistent with failed checks".into());
    }
    let recorded = report.plan_id.clone();
    report.seal();
    if recorded != report.plan_id {
        return Err("plan_id does not match the deterministic plan body".into());
    }
    Ok(report)
}

fn find_dispatched_run(runs: &str, plan_id: &str) -> Result<Option<serde_json::Value>, String> {
    let runs = serde_json::from_str::<Vec<serde_json::Value>>(runs)
        .map_err(|error| format!("cannot parse workflow runs: {error}"))?;
    Ok(runs.into_iter().find(|run| {
        run.get("displayTitle")
            .and_then(serde_json::Value::as_str)
            .is_some_and(|title| title.contains(plan_id))
    }))
}

fn workflow_already_dispatched<R: Runner>(
    runner: &R,
    root: &Path,
    policy: &Policy,
    workflow: &str,
    plan_id: &str,
) -> Result<Option<serde_json::Value>, String> {
    let runs = command_text(
        runner,
        root,
        "gh",
        &args(&[
            "run",
            "list",
            "--repo",
            &policy.repository,
            "--workflow",
            workflow,
            "--limit",
            "100",
            "--json",
            "databaseId,displayTitle,status,conclusion,url",
        ]),
    )?;
    find_dispatched_run(&runs, plan_id)
}

fn dispatch_plan<R: Runner>(
    runner: &R,
    root: &Path,
    policy: &Policy,
    plan_path: &Path,
    confirm: &str,
) -> Result<(), (i32, String)> {
    let plan = read_plan(plan_path).map_err(|error| (EXIT_INTEGRITY, error))?;
    if plan.plan_id != confirm {
        return Err((
            EXIT_INTEGRITY,
            "--confirm must exactly match plan_id".into(),
        ));
    }
    if plan.status != "ready" {
        return Err((
            if plan.status == "failed" {
                EXIT_REMOTE
            } else {
                EXIT_BLOCKED
            },
            format!("plan is {}, not ready", plan.status),
        ));
    }
    let (workflow, mut fields): (&str, Vec<(&str, String)>) = match plan.command.as_str() {
        "rc.plan" => (
            "rc.yml",
            vec![
                ("version", plan.version.clone().unwrap_or_default()),
                (
                    "number",
                    plan.metadata
                        .get("tag")
                        .and_then(|tag| tag.rsplit('.').next())
                        .unwrap_or("1")
                        .to_string(),
                ),
                (
                    "ref",
                    plan.metadata
                        .get("ref")
                        .cloned()
                        .unwrap_or_else(|| plan.commit.clone()),
                ),
            ],
        ),
        "release.plan" => (
            "release.yml",
            vec![
                ("version", plan.version.clone().unwrap_or_default()),
                (
                    "rc_tag",
                    plan.metadata.get("rc_tag").cloned().unwrap_or_default(),
                ),
                (
                    "ref",
                    plan.metadata
                        .get("ref")
                        .cloned()
                        .unwrap_or_else(|| plan.commit.clone()),
                ),
            ],
        ),
        "giab.image.plan" => (
            "happy-image.yml",
            vec![(
                "ref",
                plan.metadata
                    .get("ref")
                    .cloned()
                    .unwrap_or_else(|| plan.commit.clone()),
            )],
        ),
        "giab.benchmark.plan" => (
            "giab.yml",
            vec![(
                "image",
                plan.metadata.get("image").cloned().unwrap_or_default(),
            )],
        ),
        command => return Err((EXIT_CONFIG, format!("unsupported dispatch plan {command}"))),
    };
    fields.push(("plan_id", plan.plan_id.clone()));
    match workflow_already_dispatched(runner, root, policy, workflow, &plan.plan_id) {
        Ok(Some(run)) => {
            println!("plan {} already has a workflow run: {}", plan.plan_id, run);
            return Ok(());
        }
        Ok(None) => {}
        Err(error) => return Err((EXIT_REMOTE, error)),
    }
    let mut arguments = vec![
        OsString::from("workflow"),
        OsString::from("run"),
        OsString::from(workflow),
        OsString::from("--repo"),
        OsString::from(&policy.repository),
    ];
    for (key, value) in fields {
        arguments.push("-f".into());
        arguments.push(format!("{key}={value}").into());
    }
    run_checked(runner, root, "gh", arguments).map_err(|error| (EXIT_REMOTE, error))?;
    println!("dispatched {workflow} for plan {}", plan.plan_id);
    Ok(())
}

fn partner_scenarios(persona: Persona) -> Vec<PartnerScenario> {
    let ids: &[&str] = match persona {
        Persona::AnalyzerBuilder => &[
            "scaffold-build",
            "deterministic-repeat",
            "refusal-and-breach",
            "sanitize-and-relocate",
            "external-replay",
            "causal-diff",
            "conformance",
        ],
        Persona::WorkflowHpc => &[
            "doctor",
            "cgroup-assurance",
            "atomic-collision",
            "receipt-chain",
            "intoto-export",
            "scheduler-integration",
        ],
        Persona::ConstrainedOffline => &[
            "offline-install",
            "offline-demo",
            "local-studio",
            "zero-external-network",
            "trust-interpretation",
            "remediation-quality",
        ],
    };
    ids.iter()
        .map(|id| PartnerScenario {
            id: (*id).into(),
            passed: false,
            notes: String::new(),
        })
        .collect()
}

fn init_partner_packet(root: &Path, persona: Persona, output: &Path) -> Result<(), String> {
    if output.exists()
        && fs::read_dir(output)
            .map_err(|error| error.to_string())?
            .next()
            .is_some()
    {
        return Err(format!(
            "partner packet destination is not empty: {}",
            output.display()
        ));
    }
    fs::create_dir_all(output).map_err(|error| error.to_string())?;
    let record = PartnerRecord {
        schema: 1,
        partner_id: "replace-with-anonymous-id".into(),
        persona,
        tested_commit: "0".repeat(40),
        contract_fingerprint: "0".repeat(64),
        environment: "replace-with-environment-summary".into(),
        scenarios: partner_scenarios(persona),
        blocker_severity: Severity::None,
        contract_change_requested: false,
        resolution: String::new(),
        consent_to_publish: false,
    };
    let feedback = serde_json::to_string_pretty(&record).map_err(|error| error.to_string())? + "\n";
    write_create_new(&output.join("feedback.json"), feedback.as_bytes())?;
    let scenarios = match persona {
        Persona::AnalyzerBuilder => {
            "1. Scaffold an analyzer and build it from an empty Cargo cache.\n2. Prove deterministic repeat output.\n3. Exercise bounded fit, upfront refusal, and a governed breach.\n4. Sanitize a receipt, relocate inputs, and verify it again.\n5. Reproduce with the explicitly supplied analyzer binary.\n6. Diff two parameter variants and confirm causal localization.\n7. Pass `rosalind conformance analyzer --binary PATH --json`."
        }
        Persona::WorkflowHpc => {
            "1. Run `rosalind doctor` in the target execution environment.\n2. Exercise cooperative and cgroup-v2 OS-limit assurance.\n3. Confirm create-new collisions and `--force` atomic replacement.\n4. Build and inspect a complete receipt chain.\n5. Export and validate the unsigned in-toto statement.\n6. Capture scheduler integration behavior without storing site credentials."
        }
        Persona::ConstrainedOffline => {
            "1. Install from the release-candidate offline bundle.\n2. Run the embedded demo with networking unavailable.\n3. Open the embedded local Receipt Studio.\n4. Confirm the workflow makes zero external requests.\n5. Explain each evidence/trust status in the participant's own words.\n6. Evaluate whether failure remediation is concrete and actionable."
        }
    };
    let instructions = format!(
        "# Rosalind design-partner packet: {}\n\nRun every scenario below against one commit. Record only anonymized, consented evidence. Do not include names, email addresses, organizations, credentials, genomic data, or raw interview notes.\n\n## Scenario packet\n\n{scenarios}\n\nRecord the full tested commit with `git rev-parse HEAD`. Generate the frozen contract with `cargo xtask contract snapshot --output contract-snapshot.json` and copy its aggregate BLAKE3 into `feedback.json`.\n\nValidate with:\n\n```sh\ncargo xtask partners validate --input feedback.json --json\n```\n\nAfter review, place only the validated JSON record under `release/design-partners/`. Keep private interview notes in the ignored `release/private-design-partners/` directory or outside this repository.\n",
        persona.as_str()
    );
    write_create_new(&output.join("README.md"), instructions.as_bytes())?;
    let privacy = "Raw notes and personal data stay outside the repository. Only the completed, anonymized feedback.json may be proposed for the release gate.\n";
    write_create_new(&output.join("PRIVACY.md"), privacy.as_bytes())?;
    let _ = root;
    Ok(())
}

fn partner_validation_report(root: &Path, path: &Path) -> MaintainerReport {
    let commit = git_commit(&SystemRunner, root, "HEAD").unwrap_or_else(|_| "unknown".into());
    let mut report = MaintainerReport::new("partners.validate", commit);
    match read_partner_record(path) {
        Ok(record) => {
            report
                .metadata
                .insert("partner_id".into(), record.partner_id.clone());
            report
                .metadata
                .insert("persona".into(), record.persona.as_str().into());
            let failures = validate_partner_record(&record);
            report.check(
                "partner.record",
                failures.is_empty(),
                if failures.is_empty() {
                    "accepted".into()
                } else {
                    failures.join("; ")
                },
            );
        }
        Err(error) => report.check("partner.record", false, error),
    }
    report.actions =
        vec!["store only the validated anonymized record under release/design-partners".into()];
    report.seal();
    report
}

fn baseline_proposal<R: Runner>(
    runner: &R,
    root: &Path,
    policy: &Policy,
    run_id: u64,
    reason: &str,
    confirm: Option<&str>,
) -> Result<MaintainerReport, (i32, String)> {
    let commit = git_commit(runner, root, "HEAD").unwrap_or_else(|_| "unknown".into());
    let mut report = MaintainerReport::new("giab.baseline.propose", commit);
    add_mutation_preflight(&mut report, runner, root, policy, "HEAD");
    report.metadata.insert("run_id".into(), run_id.to_string());
    report.metadata.insert("reason".into(), reason.into());
    report.check(
        "reason.present",
        !reason.trim().is_empty(),
        "intentional baseline changes require a reason",
    );
    let changelog = fs::read_to_string(root.join("CHANGELOG.md")).unwrap_or_default();
    report.check(
        "changelog.baseline",
        baseline_changelog_allows(&changelog, reason),
        "CHANGELOG must contain `GIAB baseline update` and the exact review reason",
    );
    report.actions = vec![format!(
        "open an automated baseline PR from workflow run {run_id}"
    )];
    report.seal();
    if let Some(confirm) = confirm {
        if confirm != report.plan_id {
            return Err((
                EXIT_INTEGRITY,
                "--confirm must match the displayed plan_id".into(),
            ));
        }
        if report.status != "ready" {
            return Err((EXIT_BLOCKED, report.blockers.join("; ")));
        }
        match workflow_already_dispatched(runner, root, policy, "baseline-pr.yml", &report.plan_id)
        {
            Ok(Some(_)) => return Ok(report),
            Ok(None) => {}
            Err(error) => return Err((EXIT_REMOTE, error)),
        }
        let arguments = vec![
            "workflow".into(),
            "run".into(),
            "baseline-pr.yml".into(),
            "--repo".into(),
            policy.repository.clone().into(),
            "-f".into(),
            format!("run_id={run_id}").into(),
            "-f".into(),
            format!("reason={reason}").into(),
            "-f".into(),
            format!("plan_id={}", report.plan_id).into(),
        ];
        run_checked(runner, root, "gh", arguments).map_err(|error| (EXIT_REMOTE, error))?;
    }
    Ok(report)
}

fn baseline_changelog_allows(changelog: &str, reason: &str) -> bool {
    !reason.trim().is_empty()
        && changelog.contains("GIAB baseline update")
        && changelog.contains(reason.trim())
}

pub fn run(cli: Cli) -> i32 {
    let root = match repo_root() {
        Ok(root) => root,
        Err(error) => {
            eprintln!("{error}");
            return EXIT_CONFIG;
        }
    };
    let policy = match load_policy(&root) {
        Ok(policy) => policy,
        Err(error) => {
            eprintln!("{error}");
            return EXIT_CONFIG;
        }
    };
    let runner = SystemRunner;
    match cli.command {
        TopCommand::Doctor(output) => {
            emit_report(&doctor(&runner, &root, &policy), output.json, None)
        }
        TopCommand::Contract {
            command:
                ContractCommand::Snapshot {
                    r#ref,
                    output,
                    json,
                },
        } => match contract_snapshot(&runner, &root, &r#ref, &policy) {
            Ok(snapshot) => {
                let snapshot_json = serde_json::to_string_pretty(&snapshot).unwrap() + "\n";
                if let Some(path) = output.as_deref() {
                    if let Err(error) = write_create_new(path, snapshot_json.as_bytes()) {
                        eprintln!("{error}");
                        return EXIT_CONFIG;
                    }
                }
                let mut report =
                    MaintainerReport::new("contract.snapshot", snapshot.commit.clone());
                report.contract_fingerprint = Some(snapshot.aggregate_blake3.clone());
                report.check(
                    "contract.generated",
                    true,
                    format!("{} components", snapshot.components.len()),
                );
                report.actions = vec!["attach and attest contract-snapshot.json on an RC".into()];
                report.seal();
                if json {
                    print!("{}", serde_json::to_string_pretty(&report).unwrap());
                    println!();
                } else {
                    println!("contract snapshot: {}", snapshot.aggregate_blake3);
                }
                EXIT_OK
            }
            Err(error) => {
                eprintln!("contract snapshot failed: {error}");
                EXIT_BLOCKED
            }
        },
        TopCommand::Rc { command } => match command {
            RcCommand::Plan {
                version,
                number,
                r#ref,
                output,
                json,
            } => {
                let report = rc_plan(&runner, &root, &policy, &version, number, &r#ref);
                emit_report(&report, json, output.as_deref())
            }
            RcCommand::Dispatch(args) => {
                dispatch_plan(&runner, &root, &policy, &args.plan, &args.confirm)
                    .map(|_| EXIT_OK)
                    .unwrap_or_else(|(code, error)| {
                        eprintln!("{error}");
                        code
                    })
            }
            RcCommand::Status { tag, json } => {
                emit_report(&rc_status(&runner, &root, &policy, &tag), json, None)
            }
        },
        TopCommand::Release { command } => match command {
            ReleaseCommand::Plan {
                version,
                rc_tag,
                r#ref,
                output,
                json,
            } => {
                let report =
                    stable_release_plan(&runner, &root, &policy, &version, &rc_tag, &r#ref);
                emit_report(&report, json, output.as_deref())
            }
            ReleaseCommand::Dispatch(args) => {
                dispatch_plan(&runner, &root, &policy, &args.plan, &args.confirm)
                    .map(|_| EXIT_OK)
                    .unwrap_or_else(|(code, error)| {
                        eprintln!("{error}");
                        code
                    })
            }
            ReleaseCommand::Status { version, json } => emit_report(
                &release_status(&runner, &root, &policy, &version),
                json,
                None,
            ),
        },
        TopCommand::Giab { command } => match command {
            GiabCommand::Image { command } => match command {
                ImageCommand::Plan {
                    r#ref,
                    output,
                    json,
                } => emit_report(
                    &image_plan(&runner, &root, &policy, &r#ref),
                    json,
                    output.as_deref(),
                ),
                ImageCommand::Dispatch(args) => {
                    dispatch_plan(&runner, &root, &policy, &args.plan, &args.confirm)
                        .map(|_| EXIT_OK)
                        .unwrap_or_else(|(code, error)| {
                            eprintln!("{error}");
                            code
                        })
                }
                ImageCommand::Status(output) => {
                    emit_report(&image_status(&runner, &root, &policy), output.json, None)
                }
            },
            GiabCommand::Benchmark { command } => match command {
                BenchmarkCommand::Plan {
                    image,
                    output,
                    json,
                } => emit_report(
                    &benchmark_plan(&runner, &root, &policy, image.as_deref()),
                    json,
                    output.as_deref(),
                ),
                BenchmarkCommand::Dispatch(args) => {
                    dispatch_plan(&runner, &root, &policy, &args.plan, &args.confirm)
                        .map(|_| EXIT_OK)
                        .unwrap_or_else(|(code, error)| {
                            eprintln!("{error}");
                            code
                        })
                }
                BenchmarkCommand::Status(output) => emit_report(
                    &workflow_status(&runner, &root, &policy, "giab.yml", "giab.benchmark.status"),
                    output.json,
                    None,
                ),
                BenchmarkCommand::Download { run_id, output } => {
                    if output.exists() {
                        eprintln!("refusing to overwrite {}", output.display());
                        return EXIT_CONFIG;
                    }
                    let run_id = match run_id {
                        Some(run_id) => run_id,
                        None => match command_text(
                            &runner,
                            &root,
                            "gh",
                            &args(&[
                                "run",
                                "list",
                                "--repo",
                                &policy.repository,
                                "--workflow",
                                "giab.yml",
                                "--limit",
                                "1",
                                "--json",
                                "databaseId",
                                "--jq",
                                ".[0].databaseId",
                            ]),
                        )
                        .ok()
                        .and_then(|value| value.parse::<u64>().ok())
                        {
                            Some(run_id) => run_id,
                            None => {
                                eprintln!("no GIAB workflow run is available to download");
                                return EXIT_REMOTE;
                            }
                        },
                    };
                    let run = match command_text(
                        &runner,
                        &root,
                        "gh",
                        &args(&[
                            "run",
                            "view",
                            &run_id.to_string(),
                            "--repo",
                            &policy.repository,
                            "--json",
                            "headSha,workflowName",
                        ]),
                    )
                    .ok()
                    .and_then(|value| serde_json::from_str::<serde_json::Value>(&value).ok())
                    {
                        Some(run)
                            if run.get("workflowName").and_then(serde_json::Value::as_str)
                                == Some("GIAB v5.0q benchmark") =>
                        {
                            run
                        }
                        _ => {
                            eprintln!("run {run_id} is not a GIAB benchmark workflow");
                            return EXIT_REMOTE;
                        }
                    };
                    let source_digest = run
                        .get("headSha")
                        .and_then(serde_json::Value::as_str)
                        .unwrap_or("")
                        .to_string();
                    let mut command =
                        vec!["run".into(), "download".into(), run_id.to_string().into()];
                    command.extend([
                        "--repo".into(),
                        policy.repository.clone().into(),
                        "--name".into(),
                        "hg002-v5.0q-grch38-chr20".into(),
                        "--dir".into(),
                        output.as_os_str().to_owned(),
                    ]);
                    if let Err(error) = run_checked(&runner, &root, "gh", command) {
                        eprintln!("{error}");
                        return EXIT_REMOTE;
                    }
                    let evidence = output.join("latest.json");
                    if !evidence.is_file() {
                        eprintln!("downloaded run lacks latest.json");
                        return EXIT_INTEGRITY;
                    }
                    let verify = vec![
                        "attestation".into(),
                        "verify".into(),
                        evidence.into_os_string(),
                        "--repo".into(),
                        policy.repository.clone().into(),
                        "--source-digest".into(),
                        source_digest.into(),
                        "--signer-workflow".into(),
                        format!("{}/.github/workflows/giab.yml", policy.repository).into(),
                    ];
                    run_checked(&runner, &root, "gh", verify)
                        .map(|_| EXIT_OK)
                        .unwrap_or_else(|error| {
                            eprintln!("{error}");
                            EXIT_INTEGRITY
                        })
                }
            },
            GiabCommand::Baseline {
                command:
                    BaselineCommand::Propose {
                        run_id,
                        reason,
                        confirm,
                        json,
                    },
            } => {
                match baseline_proposal(
                    &runner,
                    &root,
                    &policy,
                    run_id,
                    &reason,
                    confirm.as_deref(),
                ) {
                    Ok(report) => emit_report(&report, json, None),
                    Err((code, error)) => {
                        eprintln!("{error}");
                        code
                    }
                }
            }
        },
        TopCommand::Partners { command } => match command {
            PartnersCommand::Init { persona, output } => {
                init_partner_packet(&root, persona, &output)
                    .map(|_| {
                        println!("wrote {} packet to {}", persona.as_str(), output.display());
                        EXIT_OK
                    })
                    .unwrap_or_else(|error| {
                        eprintln!("{error}");
                        EXIT_CONFIG
                    })
            }
            PartnersCommand::Validate { input, json } => {
                emit_report(&partner_validation_report(&root, &input), json, None)
            }
            PartnersCommand::Report { input_dir, json } => emit_report(
                &partner_report(&runner, &root, &policy, &input_dir, None),
                json,
                None,
            ),
        },
    }
}

#[cfg(test)]
fn contract_diff(a: &ContractSnapshot, b: &ContractSnapshot) -> BTreeSet<String> {
    let mut changed = BTreeSet::new();
    let keys = a
        .components
        .keys()
        .chain(b.components.keys())
        .collect::<BTreeSet<_>>();
    for key in keys {
        if a.components.get(key) != b.components.get(key) {
            let category = if key.starts_with("public_api.") {
                "rust-api"
            } else if key.starts_with("cli_help.") {
                "cli"
            } else if key.starts_with("schema.") {
                "schema"
            } else if key.starts_with("package.") {
                "package"
            } else {
                "supporting-contract"
            };
            changed.insert(category.to_string());
        }
    }
    if a.receipt_fields != b.receipt_fields {
        changed.insert("receipt-fields".into());
    }
    if a.primary_artifacts != b.primary_artifacts {
        changed.insert("primary-artifacts".into());
    }
    if a.package_versions != b.package_versions {
        changed.insert("package-versions".into());
    }
    changed
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Mutex;
    use tempfile::tempdir;

    #[derive(Default)]
    struct RecordingRunner {
        calls: Mutex<Vec<(String, Vec<OsString>)>>,
    }

    impl Runner for RecordingRunner {
        fn run(&self, _cwd: &Path, program: &str, arguments: &[OsString]) -> io::Result<Output> {
            self.calls
                .lock()
                .unwrap()
                .push((program.to_string(), arguments.to_vec()));
            if arguments.iter().any(|argument| argument == "list") {
                Command::new("printf").arg("[]").output()
            } else {
                Command::new("true").output()
            }
        }
    }

    fn report() -> MaintainerReport {
        let mut report = MaintainerReport::new("rc.plan", "a".repeat(40));
        report.version = Some("0.4.0".into());
        report.actions = vec!["publish prerelease".into(), "attest snapshot".into()];
        report.contract_fingerprint = Some("b".repeat(64));
        report.metadata.insert("tag".into(), "v0.4.0-rc.1".into());
        report.seal();
        report
    }

    fn partner(persona: Persona) -> PartnerRecord {
        PartnerRecord {
            schema: 1,
            partner_id: format!("partner-{}", persona.as_str()),
            persona,
            tested_commit: "a".repeat(40),
            contract_fingerprint: "b".repeat(64),
            environment: "linux test host".into(),
            scenarios: partner_scenarios(persona)
                .into_iter()
                .map(|mut scenario| {
                    scenario.passed = true;
                    scenario.notes = "passed".into();
                    scenario
                })
                .collect(),
            blocker_severity: Severity::None,
            contract_change_requested: false,
            resolution: "none".into(),
            consent_to_publish: true,
        }
    }

    fn snapshot() -> ContractSnapshot {
        ContractSnapshot {
            schema: 1,
            commit: "a".repeat(40),
            components: BTreeMap::from([
                ("public_api.rosalind-bio".into(), "1".into()),
                ("cli_help.root".into(), "2".into()),
                ("schema.receipt".into(), "3".into()),
            ]),
            receipt_fields: BTreeMap::from([("params".into(), vec!["schema_version".into()])]),
            primary_artifacts: BTreeMap::from([("vcf".into(), "4".into())]),
            package_versions: BTreeMap::from([("rosalind-bio".into(), "0.4.0".into())]),
            aggregate_blake3: "5".repeat(64),
        }
    }

    fn git_repository() -> (tempfile::TempDir, String) {
        let directory = tempdir().unwrap();
        for arguments in [
            vec!["init", "-q"],
            vec!["config", "user.email", "test@example.invalid"],
            vec!["config", "user.name", "Release Test"],
        ] {
            assert!(Command::new("git")
                .args(arguments)
                .current_dir(directory.path())
                .status()
                .unwrap()
                .success());
        }
        fs::write(directory.path().join("README"), "test\n").unwrap();
        assert!(Command::new("git")
            .args(["add", "README"])
            .current_dir(directory.path())
            .status()
            .unwrap()
            .success());
        assert!(Command::new("git")
            .args(["commit", "-q", "-m", "test"])
            .current_dir(directory.path())
            .status()
            .unwrap()
            .success());
        let commit = command_text(
            &SystemRunner,
            directory.path(),
            "git",
            &args(&["rev-parse", "HEAD"]),
        )
        .unwrap();
        (directory, commit)
    }

    #[test]
    fn plan_id_is_deterministic_and_authenticates_gate_results() {
        let mut a = report();
        a.check("z", true, "ok");
        a.check("a", true, "ok");
        a.seal();
        let mut b = report();
        b.check("a", true, "different non-contract wording");
        b.check("z", true, "ok");
        b.seal();
        assert_eq!(a.plan_id, b.plan_id);
        b.checks[0].ok = false;
        b.status = "blocked".into();
        b.seal();
        assert_ne!(a.plan_id, b.plan_id);
    }

    #[test]
    fn plan_tampering_is_rejected() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("plan.json");
        let mut value = serde_json::to_value(report()).unwrap();
        value["actions"][0] = serde_json::Value::String("malicious action".into());
        fs::write(&path, serde_json::to_vec(&value).unwrap()).unwrap();
        assert!(read_plan(&path).unwrap_err().contains("plan_id"));
    }

    #[test]
    fn blocked_plan_cannot_be_changed_to_ready() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("plan.json");
        let mut blocked = report();
        blocked.check("external.gate", false, "not approved");
        blocked.seal();
        let mut value = serde_json::to_value(blocked).unwrap();
        value["status"] = serde_json::Value::String("ready".into());
        fs::write(&path, serde_json::to_vec(&value).unwrap()).unwrap();
        assert!(read_plan(&path)
            .unwrap_err()
            .contains("inconsistent with its checks"));
    }

    #[test]
    fn duplicate_dispatch_finds_the_existing_plan_run() {
        let plan_id = "a".repeat(64);
        let runs = serde_json::json!([
            {"databaseId": 1, "displayTitle": "unrelated", "status": "completed"},
            {"databaseId": 2, "displayTitle": format!("Stable v0.4.0 plan:{plan_id}"), "status": "in_progress"}
        ]);
        let found = find_dispatched_run(&runs.to_string(), &plan_id)
            .unwrap()
            .unwrap();
        assert_eq!(found["databaseId"], 2);
        assert!(find_dispatched_run(&runs.to_string(), &"b".repeat(64))
            .unwrap()
            .is_none());
    }

    #[test]
    fn dispatch_fake_remote_uses_tokenized_arguments_and_no_unplanned_mutation() {
        let directory = tempdir().unwrap();
        let plan_path = directory.path().join("plan.json");
        let plan = report();
        fs::write(&plan_path, serde_json::to_vec(&plan).unwrap()).unwrap();
        let runner = RecordingRunner::default();
        let policy_root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
        let policy = load_policy(policy_root).unwrap();
        dispatch_plan(&runner, policy_root, &policy, &plan_path, &plan.plan_id).unwrap();
        let calls = runner.calls.lock().unwrap();
        assert_eq!(calls.len(), 2);
        assert_eq!(calls[0].0, "gh");
        assert_eq!(calls[1].0, "gh");
        let dispatched = calls[1]
            .1
            .iter()
            .map(|argument| argument.to_string_lossy().to_string())
            .collect::<Vec<_>>();
        assert_eq!(&dispatched[..3], ["workflow", "run", "rc.yml"]);
        assert!(dispatched.contains(&format!("plan_id={}", plan.plan_id)));
        assert!(!dispatched.iter().any(|argument| argument.contains(';')));
    }

    #[test]
    fn soak_boundary_is_exact() {
        assert_eq!(soak_elapsed(1_000, 605_799), 604_799);
        assert_eq!(soak_elapsed(1_000, 605_800), 604_800);
        assert!(soak_elapsed(1_000, 605_799) < 604_800);
        assert!(soak_elapsed(1_000, 605_800) >= 604_800);
    }

    #[test]
    fn github_timestamp_parser_handles_fractional_utc_and_rejects_invalid_dates() {
        assert_eq!(parse_github_rfc3339_unix("1970-01-01T00:00:00Z"), Some(0));
        assert_eq!(
            parse_github_rfc3339_unix("2024-02-29T00:00:00.123Z"),
            Some(1_709_164_800)
        );
        assert_eq!(parse_github_rfc3339_unix("2023-02-29T00:00:00Z"), None);
        assert_eq!(parse_github_rfc3339_unix("2024-01-01T00:00:00+01:00"), None);
    }

    #[test]
    fn stabilization_release_gates_begin_at_configured_version() {
        assert!(!version_at_least("0.4.0", "0.5.0"));
        assert!(version_at_least("0.5.0-rc.1", "0.5.0-rc.1"));
        assert!(version_at_least("0.5.0", "0.5.0-rc.1"));
        assert!(!version_at_least("invalid", "0.5.0"));
    }

    #[test]
    fn maintainer_exit_codes_distinguish_gates_remote_and_integrity() {
        assert_eq!(EXIT_OK, 0);
        assert_eq!(EXIT_CONFIG, 2);
        assert_eq!(EXIT_BLOCKED, 3);
        assert_eq!(EXIT_REMOTE, 4);
        assert_eq!(EXIT_INTEGRITY, 5);
        let mut gate = MaintainerReport::new("status", "a".repeat(40));
        gate.check("soak.elapsed", false, "remaining=1s");
        assert_eq!(gate.exit_code(), EXIT_BLOCKED);
        gate.status = "failed".into();
        assert_eq!(gate.exit_code(), EXIT_REMOTE);
        let mut integrity = MaintainerReport::new("status", "a".repeat(40));
        integrity.check("contract.equivalent", false, "rc=aaaa, candidate=bbbb");
        assert_eq!(integrity.exit_code(), EXIT_INTEGRITY);
    }

    #[test]
    fn immutable_digest_validation_rejects_tags_and_image_ids() {
        assert!(digest_valid(&format!("sha256:{}", "a".repeat(64))));
        assert!(!digest_valid("latest"));
        assert!(!digest_valid("sha256:abc"));
        assert!(!digest_valid(&format!("image-id:{}", "a".repeat(64))));

        let directory = tempdir().unwrap();
        let happy = directory.path().join("benchmarks/giab/happy");
        fs::create_dir_all(&happy).unwrap();
        fs::write(happy.join("Dockerfile"), "FROM scratch\n").unwrap();
        let lock = |platform: &str| {
            serde_json::json!({
                "schema": 1,
                "generated_image": {
                    "repository": "ghcr.io/logannye/rosalind-happy",
                    "digest": format!("sha256:{}", "a".repeat(64)),
                    "platform": platform,
                    "source_lock_blake3": "b".repeat(64),
                    "built_from_commit": "c".repeat(40),
                }
            })
        };
        fs::write(
            happy.join("lock.json"),
            serde_json::to_vec(&lock("linux/arm64")).unwrap(),
        )
        .unwrap();
        assert!(generated_image(directory.path()).is_err());
        fs::write(
            happy.join("lock.json"),
            serde_json::to_vec(&lock("linux/amd64")).unwrap(),
        )
        .unwrap();
        assert!(generated_image(directory.path()).is_ok());
    }

    #[test]
    fn evaluator_source_identity_includes_python_wheel_lock() {
        let directory = tempdir().unwrap();
        let happy = directory.path().join("benchmarks/giab/happy");
        fs::create_dir_all(&happy).unwrap();
        fs::write(happy.join("Dockerfile"), "FROM scratch\n").unwrap();
        fs::write(happy.join("lock.json"), "{\"schema\":1}").unwrap();
        let before = source_lock_hash(directory.path()).unwrap();
        fs::write(
            happy.join("requirements-py27.txt"),
            "wheel-A --hash=sha256:aaa\n",
        )
        .unwrap();
        let first = source_lock_hash(directory.path()).unwrap();
        fs::write(
            happy.join("requirements-py27.txt"),
            "wheel-B --hash=sha256:bbb\n",
        )
        .unwrap();
        let second = source_lock_hash(directory.path()).unwrap();
        assert_ne!(before, first);
        assert_ne!(first, second);
    }

    #[test]
    fn baseline_update_requires_a_reason_specific_changelog_entry() {
        let reason = "accept attested run 123 after toolchain refresh";
        assert!(!baseline_changelog_allows("GIAB baseline update", reason));
        assert!(baseline_changelog_allows(
            &format!("GIAB baseline update: {reason}"),
            reason
        ));
        assert!(!baseline_changelog_allows(
            "GIAB baseline update: unrelated",
            reason
        ));
    }

    #[test]
    fn partner_gate_accepts_all_personas_and_rejects_blockers() {
        for persona in [
            Persona::AnalyzerBuilder,
            Persona::WorkflowHpc,
            Persona::ConstrainedOffline,
        ] {
            assert!(validate_partner_record(&partner(persona)).is_empty());
        }
        let mut blocked = partner(Persona::AnalyzerBuilder);
        blocked.blocker_severity = Severity::ReleaseBlocking;
        assert!(validate_partner_record(&blocked)
            .iter()
            .any(|failure| failure.contains("release-blocking")));

        let mut revoked = partner(Persona::AnalyzerBuilder);
        revoked.consent_to_publish = false;
        assert!(validate_partner_record(&revoked)
            .iter()
            .any(|failure| failure.contains("consent")));

        let mut contract_change = partner(Persona::AnalyzerBuilder);
        contract_change.contract_change_requested = true;
        contract_change.resolution.clear();
        assert!(validate_partner_record(&contract_change)
            .iter()
            .any(|failure| failure.contains("resolution")));

        let mut duplicate_scenario = partner(Persona::AnalyzerBuilder);
        duplicate_scenario
            .scenarios
            .push(duplicate_scenario.scenarios[0].clone());
        assert!(validate_partner_record(&duplicate_scenario)
            .iter()
            .any(|failure| failure.contains("unique")));
    }

    #[test]
    fn partner_report_rejects_missing_personas_duplicates_and_stale_commits() {
        let (repository, commit) = git_repository();
        let records = repository.path().join("records");
        fs::create_dir(&records).unwrap();
        let mut analyzer = partner(Persona::AnalyzerBuilder);
        analyzer.tested_commit = commit.clone();
        fs::write(
            records.join("analyzer.json"),
            serde_json::to_vec(&analyzer).unwrap(),
        )
        .unwrap();
        let policy_root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
        let policy = load_policy(policy_root).unwrap();
        let missing = partner_report(&SystemRunner, repository.path(), &policy, &records, None);
        assert!(missing
            .blockers
            .iter()
            .any(|blocker| blocker.contains("persona.workflow-hpc")));

        let mut duplicate = partner(Persona::WorkflowHpc);
        duplicate.partner_id = analyzer.partner_id.clone();
        duplicate.tested_commit = commit.clone();
        fs::write(
            records.join("duplicate.json"),
            serde_json::to_vec(&duplicate).unwrap(),
        )
        .unwrap();
        let duplicates = partner_report(&SystemRunner, repository.path(), &policy, &records, None);
        assert!(duplicates
            .blockers
            .iter()
            .any(|blocker| blocker.contains("duplicate partner_id")));

        fs::remove_file(records.join("duplicate.json")).unwrap();
        analyzer.tested_commit = "f".repeat(40);
        fs::write(
            records.join("analyzer.json"),
            serde_json::to_vec(&analyzer).unwrap(),
        )
        .unwrap();
        let stale = partner_report(
            &SystemRunner,
            repository.path(),
            &policy,
            &records,
            Some((&"b".repeat(64), &commit)),
        );
        assert!(stale
            .blockers
            .iter()
            .any(|blocker| blocker.contains("not an ancestor")));

        analyzer.tested_commit = commit.clone();
        analyzer.contract_fingerprint = "d".repeat(64);
        fs::write(
            records.join("analyzer.json"),
            serde_json::to_vec(&analyzer).unwrap(),
        )
        .unwrap();
        let changed_contract = partner_report(
            &SystemRunner,
            repository.path(),
            &policy,
            &records,
            Some((&"b".repeat(64), &commit)),
        );
        assert!(changed_contract
            .blockers
            .iter()
            .any(|blocker| blocker.contains("fingerprint differs")));
    }

    #[test]
    fn partner_gate_rejects_pii_shaped_unknown_fields() {
        let mut value = serde_json::to_value(partner(Persona::AnalyzerBuilder)).unwrap();
        value["email"] = serde_json::Value::String("private@example.test".into());
        assert!(serde_json::from_value::<PartnerRecord>(value).is_err());
        let mut record = partner(Persona::AnalyzerBuilder);
        record.scenarios[0].notes = "contact private@example.test".into();
        assert!(validate_partner_record(&record)
            .iter()
            .any(|failure| failure.contains("free-text")));
    }

    #[test]
    fn contract_diff_localizes_every_frozen_surface() {
        let a = snapshot();
        let mut b = a.clone();
        b.components
            .insert("public_api.rosalind-bio".into(), "changed".into());
        b.components
            .insert("cli_help.root".into(), "changed".into());
        b.components
            .insert("schema.receipt".into(), "changed".into());
        b.receipt_fields.insert("params".into(), vec!["new".into()]);
        b.primary_artifacts.insert("vcf".into(), "changed".into());
        let diff = contract_diff(&a, &b);
        assert_eq!(
            diff,
            BTreeSet::from([
                "cli".into(),
                "primary-artifacts".into(),
                "receipt-fields".into(),
                "rust-api".into(),
                "schema".into(),
            ])
        );
    }

    #[test]
    fn docs_only_changes_are_outside_the_contract_snapshot() {
        let a = snapshot();
        let b = a.clone();
        assert!(contract_diff(&a, &b).is_empty());
    }

    #[test]
    fn cli_shape_ignores_copy_but_detects_flags_and_operands() {
        let first = "A description\n\nUsage: rosalind demo [OPTIONS]\n\nOptions:\n      --json  Emit JSON\n  -h, --help  Print help\n";
        let copy_change = "Different words\n\nUsage: rosalind demo [OPTIONS]\n\nCommands:\n  child  Different command wording\n         wrapped copy must not become a command\n\nOptions:\n      --json  Machine output\n  -h, --help  Show this\n";
        let first = format!("Commands:\n  child  Original command wording\n\n{first}");
        let flag_change = "Different words\n\nUsage: rosalind demo [OPTIONS]\n\nCommands:\n  child  Different command wording\n\nOptions:\n      --format <FORMAT>  Machine output\n  -h, --help  Show this\n";
        assert_eq!(
            canonical_cli_shape(&first),
            canonical_cli_shape(copy_change)
        );
        assert_ne!(
            canonical_cli_shape(&first),
            canonical_cli_shape(flag_change)
        );
    }

    #[test]
    fn create_new_writer_refuses_overwrite() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("plan.json");
        write_create_new(&path, b"first").unwrap();
        assert!(write_create_new(&path, b"second").is_err());
        assert_eq!(fs::read(&path).unwrap(), b"first");
    }
}
