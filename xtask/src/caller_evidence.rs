//! Change-based scientific release gate. Feature-only work does not inherit a
//! caller benchmark requirement; shared caller dependencies remain covered.

use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

use quote::ToTokens;
use syn::visit::Visit;

use crate::{args, command_text, Policy, Runner};

#[derive(Default)]
struct Dependencies(BTreeSet<String>);

impl<'ast> Visit<'ast> for Dependencies {
    fn visit_path(&mut self, path: &'ast syn::Path) {
        if let Some(segment) = path.segments.first() {
            self.0.insert(segment.ident.to_string());
        }
        syn::visit::visit_path(self, path);
    }
}

#[derive(Default)]
struct VariantDispatch(Vec<syn::Arm>);

impl<'ast> Visit<'ast> for VariantDispatch {
    fn visit_arm(&mut self, arm: &'ast syn::Arm) {
        if let syn::Pat::Struct(pattern) = &arm.pat {
            if pattern.path.to_token_stream().to_string() == "Commands :: Variants" {
                self.0.push(arm.clone());
            }
        }
        syn::visit::visit_arm(self, arm);
    }
}

fn item_name(item: &syn::Item) -> Option<String> {
    match item {
        syn::Item::Fn(item) => Some(item.sig.ident.to_string()),
        syn::Item::Struct(item) => Some(item.ident.to_string()),
        syn::Item::Enum(item) => Some(item.ident.to_string()),
        syn::Item::Type(item) => Some(item.ident.to_string()),
        syn::Item::Const(item) => Some(item.ident.to_string()),
        syn::Item::Static(item) => Some(item.ident.to_string()),
        _ => None,
    }
}

/// Parse Rust, including nested braces/comments/raw strings. Hash the Variants
/// defaults and dispatch plus transitively referenced local helpers and types.
/// Adding an unrelated CLI command therefore does not change the caller identity.
fn caller_cli(source: &str) -> Result<String, String> {
    let parsed = syn::parse_file(source).map_err(|error| format!("caller CLI parse: {error}"))?;
    let items = parsed
        .items
        .iter()
        .filter_map(|item| item_name(item).map(|name| (name, item)))
        .collect::<BTreeMap<_, _>>();
    let Some(syn::Item::Enum(commands)) = items.get("Commands").copied() else {
        return Err("caller CLI lacks Commands enum".into());
    };
    let variant = commands
        .variants
        .iter()
        .find(|variant| variant.ident == "Variants")
        .ok_or("caller CLI lacks Variants arguments")?;
    let mut dispatch = VariantDispatch::default();
    dispatch.visit_file(&parsed);
    if dispatch.0.len() != 1 {
        return Err("caller CLI must have exactly one Variants dispatch arm".into());
    }
    let mut dependencies = Dependencies::default();
    dependencies.visit_variant(variant);
    dependencies.visit_arm(&dispatch.0[0]);
    let mut output = vec![
        variant.to_token_stream().to_string(),
        dispatch.0[0].to_token_stream().to_string(),
    ];
    let mut visited = BTreeSet::from(["Commands".to_string(), "main".to_string()]);
    loop {
        let pending = dependencies
            .0
            .difference(&visited)
            .cloned()
            .collect::<Vec<_>>();
        if pending.is_empty() {
            break;
        }
        for name in pending {
            visited.insert(name.clone());
            if let Some(item) = items.get(&name) {
                output.push(item.to_token_stream().to_string());
                dependencies.visit_item(item);
            }
        }
    }
    output.sort();
    Ok(output.join("\n"))
}

fn normalized_lock(bytes: &str) -> Result<String, String> {
    let mut lock: toml::Value =
        toml::from_str(bytes).map_err(|error| format!("Cargo.lock: {error}"))?;
    let packages = lock
        .get_mut("package")
        .and_then(toml::Value::as_array_mut)
        .ok_or("Cargo.lock lacks packages")?;
    packages.retain(|package| package.get("name").and_then(toml::Value::as_str) != Some("xtask"));
    for package in packages {
        if package.get("name").and_then(toml::Value::as_str) == Some("rosalind-bio") {
            package.as_table_mut().unwrap().remove("version");
        }
    }
    Ok(lock.to_string())
}

pub(crate) fn fingerprint<R: Runner>(
    runner: &R,
    root: &Path,
    reference: &str,
    paths: &[String],
) -> Result<String, String> {
    // Resolve before constructing git object names; arbitrary option/ref strings
    // never become arguments to git show.
    let commit = crate::git_commit(runner, root, reference)?;
    if commit.len() != 40 || !commit.bytes().all(|byte| byte.is_ascii_hexdigit()) {
        return Err("caller evidence reference did not resolve to a full commit".into());
    }
    let mut arguments = args(&["ls-tree", "-r", "--name-only", &commit, "--"]);
    arguments.extend(paths.iter().map(Into::into));
    let mut files = command_text(runner, root, "git", &arguments)?
        .lines()
        .map(str::to_owned)
        .collect::<BTreeSet<_>>();
    files.insert("src/main.rs".into());
    files.insert("Cargo.lock".into());
    let mut hash = blake3::Hasher::new();
    hash.update(b"rosalind-caller-source-v1\0");
    for path in files {
        let source = command_text(
            runner,
            root,
            "git",
            &args(&["show", &format!("{commit}:{path}")]),
        )?;
        let content = match path.as_str() {
            "src/main.rs" => caller_cli(&source)?,
            "Cargo.lock" => normalized_lock(&source)?,
            _ => source,
        };
        hash.update(path.as_bytes());
        hash.update(&[0]);
        hash.update(content.as_bytes());
        hash.update(&[0]);
    }
    Ok(hash.finalize().to_hex().to_string())
}

pub(crate) fn gate<R: Runner>(
    runner: &R,
    root: &Path,
    policy: &Policy,
    candidate: &str,
) -> Result<String, String> {
    let candidate_hash = fingerprint(runner, root, candidate, &policy.caller_source_paths)?;
    let baseline_path = root.join("benchmarks/giab/baseline.json");
    let baseline: serde_json::Value = serde_json::from_slice(
        &std::fs::read(&baseline_path)
            .map_err(|error| format!("GIAB baseline unavailable: {error}"))?,
    )
    .map_err(|error| format!("invalid GIAB baseline: {error}"))?;
    if baseline.get("status").and_then(serde_json::Value::as_str) == Some("not-yet-established") {
        let bootstrap = fingerprint(
            runner,
            root,
            &policy.caller_source_base,
            &policy.caller_source_paths,
        )?;
        if candidate_hash == bootstrap {
            return Ok("experimental caller unchanged from declared bootstrap; no scientific score claimed".into());
        }
        return Err("caller or shared caller dependencies changed: establish and review GIAB evidence for this candidate".into());
    }
    for field in [
        "calls_filter_all",
        "calls_filter_pass",
        "external_happy_vcfeval",
        "data_manifest",
    ] {
        if baseline
            .get(field)
            .and_then(serde_json::Value::as_object)
            .is_none_or(|value| value.is_empty())
        {
            return Err(format!("accepted GIAB baseline lacks {field}"));
        }
    }
    let evaluator = baseline["external_happy_vcfeval"]
        .get("container")
        .and_then(serde_json::Value::as_str)
        .unwrap_or_default();
    if !evaluator
        .split_once('@')
        .is_some_and(|(_, digest)| crate::digest_valid(digest))
    {
        return Err("accepted GIAB baseline lacks an immutable evaluator image".into());
    }
    let identity = baseline
        .get("producer_identity")
        .ok_or("GIAB baseline lacks producer identity")?;
    if identity
        .get("code_dirty")
        .and_then(serde_json::Value::as_str)
        != Some("false")
    {
        return Err("GIAB baseline must come from a clean committed build".into());
    }
    let tested = identity
        .get("code_git_sha")
        .and_then(serde_json::Value::as_str)
        .ok_or("GIAB baseline lacks tested source commit")?;
    let tested_hash = fingerprint(runner, root, tested, &policy.caller_source_paths)?;
    if tested_hash != candidate_hash {
        return Err("caller or shared caller dependencies changed since accepted GIAB evidence; review a new baseline".into());
    }
    Ok(format!(
        "caller source matches reviewed GIAB build {tested}"
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    const CLI: &str = r##"
        enum Commands { Variants { #[arg(default_value_t = 30)] quality: u8 }, Features {} }
        fn main() { match c { Commands::Variants { quality } => call(quality), Commands::Features {} => feature() } }
        fn call(q: u8) { helper(q); let example = r#"brace }"#; }
        fn helper(q: u8) { consume(q); }
        fn feature() { do_features(); }
    "##;

    #[test]
    fn feature_only_cli_changes_do_not_require_caller_evidence() {
        let original = caller_cli(CLI).unwrap();
        assert_eq!(
            original,
            caller_cli(&CLI.replace("do_features()", "a_new_feature()")).unwrap()
        );
        assert_eq!(
            original,
            caller_cli(&CLI.replace("Features {}", "Features { extra: bool }")).unwrap()
        );
    }

    #[test]
    fn defaults_dispatch_and_transitive_helpers_change_caller_identity() {
        let original = caller_cli(CLI).unwrap();
        for edited in [
            CLI.replace("= 30", "= 20"),
            CLI.replace("call(quality)", "call(quality + 1)"),
            CLI.replace("consume(q)", "consume(q + 1)"),
        ] {
            assert_ne!(original, caller_cli(&edited).unwrap());
        }
    }

    #[test]
    fn invalid_or_missing_caller_source_fails_closed() {
        assert!(caller_cli("invalid Rust").is_err());
        assert!(caller_cli("enum Commands { Features {} }").is_err());
    }

    #[test]
    fn release_version_and_xtask_do_not_change_dependency_identity() {
        let lock = "version = 3\n[[package]]\nname = 'rosalind-bio'\nversion = '0.4.0'\n[[package]]\nname = 'xtask'\nversion = '0.1.0'\n[[package]]\nname = 'rust-htslib'\nversion = '0.44.1'\n";
        assert_eq!(
            normalized_lock(lock).unwrap(),
            normalized_lock(&lock.replace("0.4.0", "0.4.1").replace("0.1.0", "0.2.0")).unwrap()
        );
        assert_ne!(
            normalized_lock(lock).unwrap(),
            normalized_lock(&lock.replace("0.44.1", "0.44.2")).unwrap()
        );
    }

    #[test]
    fn shared_caller_changes_require_fresh_reviewed_evidence_but_feature_changes_do_not() {
        use crate::SystemRunner;
        let directory = tempfile::tempdir().unwrap();
        let root = directory.path();
        let git = |arguments: &[&str]| {
            command_text(&SystemRunner, root, "git", &args(arguments)).unwrap()
        };
        git(&["init", "-q"]);
        git(&["config", "user.email", "test@example.invalid"]);
        git(&["config", "user.name", "Caller evidence test"]);
        std::fs::create_dir_all(root.join("src/pileup")).unwrap();
        std::fs::create_dir_all(root.join("benchmarks/giab")).unwrap();
        std::fs::write(root.join("src/main.rs"), CLI).unwrap();
        std::fs::write(root.join("src/pileup/engine.rs"), "fn pileup() {}\n").unwrap();
        std::fs::write(
            root.join("Cargo.lock"),
            "version = 3\n[[package]]\nname = 'rosalind-bio'\nversion = '0.4.0'\n",
        )
        .unwrap();
        let baseline = root.join("benchmarks/giab/baseline.json");
        std::fs::write(&baseline, r#"{"status":"not-yet-established"}"#).unwrap();
        git(&["add", "."]);
        git(&["commit", "-qm", "bootstrap"]);
        let base = git(&["rev-parse", "HEAD"]);
        let mut policy =
            crate::load_policy(Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap()).unwrap();
        policy.caller_source_base = base.clone();
        policy.caller_source_paths = vec!["src/pileup/".into()];
        assert!(gate(&SystemRunner, root, &policy, &base).is_ok());

        std::fs::write(
            root.join("src/main.rs"),
            CLI.replace("do_features()", "improved_features()"),
        )
        .unwrap();
        git(&["add", "."]);
        git(&["commit", "-qm", "feature only"]);
        assert!(gate(&SystemRunner, root, &policy, "HEAD").is_ok());

        std::fs::write(
            root.join("src/pileup/engine.rs"),
            "fn pileup() { changed(); }\n",
        )
        .unwrap();
        git(&["add", "."]);
        git(&["commit", "-qm", "shared pileup change"]);
        let changed = git(&["rev-parse", "HEAD"]);
        assert!(gate(&SystemRunner, root, &policy, &changed)
            .unwrap_err()
            .contains("shared caller"));
        let accepted = |source: &str| {
            serde_json::json!({
                "calls_filter_all": {"precision": 0.1}, "calls_filter_pass": {"precision": 0.1},
                "external_happy_vcfeval": {"container": format!("image@sha256:{}", "a".repeat(64))},
                "data_manifest": {"schema": 1},
                "producer_identity": {"code_git_sha": source, "code_dirty": "false"}
            })
        };
        std::fs::write(&baseline, accepted(&base).to_string()).unwrap();
        assert!(gate(&SystemRunner, root, &policy, &changed).is_err());
        std::fs::write(&baseline, accepted(&changed).to_string()).unwrap();
        assert!(gate(&SystemRunner, root, &policy, &changed).is_ok());
        std::fs::write(&baseline, "{}").unwrap();
        assert!(gate(&SystemRunner, root, &policy, &changed).is_err());
        std::fs::remove_file(&baseline).unwrap();
        assert!(gate(&SystemRunner, root, &policy, &changed).is_err());
    }
}
