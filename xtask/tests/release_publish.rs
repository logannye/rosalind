#![cfg(unix)]

use std::fs;
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use tempfile::TempDir;

const PACKAGES: [(&str, &str); 3] = [
    ("rosalind-build-info", "0.1.0"),
    ("rosalind-receipt", "0.3.0"),
    ("rosalind-bio", "0.4.0"),
];

struct FakeRegistry {
    temp: TempDir,
    repo: PathBuf,
    registry: PathBuf,
    bin: PathBuf,
}

impl FakeRegistry {
    fn new() -> Self {
        let temp = tempfile::tempdir().unwrap();
        let repo = temp.path().join("repo");
        let registry = temp.path().join("registry");
        let bin = temp.path().join("bin");
        fs::create_dir_all(repo.join("scripts")).unwrap();
        fs::create_dir_all(&registry).unwrap();
        fs::create_dir_all(&bin).unwrap();
        fs::copy(
            Path::new(env!("CARGO_MANIFEST_DIR")).join("../scripts/release-publish.sh"),
            repo.join("scripts/release-publish.sh"),
        )
        .unwrap();
        write_executable(&bin.join("cargo"), FAKE_CARGO);
        write_executable(&bin.join("curl"), FAKE_CURL);
        Self {
            temp,
            repo,
            registry,
            bin,
        }
    }

    fn archive(&self, package: &str, version: &str, contents: &str) {
        fs::write(
            self.registry.join(format!("{package}-{version}.crate")),
            contents,
        )
        .unwrap();
    }

    fn exact_archive(&self, package: &str, version: &str) {
        self.archive(package, version, &format!("{package}-{version}\n"));
    }

    fn run(&self, race_package: Option<&str>) -> Output {
        let inherited_path = std::env::var_os("PATH").unwrap_or_default();
        let path = std::env::join_paths(
            std::iter::once(self.bin.clone()).chain(std::env::split_paths(&inherited_path)),
        )
        .unwrap();
        let mut command = Command::new("bash");
        command
            .arg(self.repo.join("scripts/release-publish.sh"))
            .current_dir(&self.repo)
            .env("PATH", path)
            .env("CARGO_REGISTRY_TOKEN", "fake-token-never-printed")
            .env("CARGO_REGISTRY_API", "https://registry.invalid/api/v1")
            .env("FAKE_REGISTRY", &self.registry)
            .env("FAKE_PUBLISH_LOG", self.temp.path().join("publish.log"));
        if let Some(package) = race_package {
            command.env("FAKE_PUBLISH_RACE_PACKAGE", package);
        }
        command.output().unwrap()
    }

    fn state(&self) -> serde_json::Value {
        serde_json::from_slice(&fs::read(self.repo.join("publish-state.json")).unwrap()).unwrap()
    }

    fn publish_log(&self) -> String {
        fs::read_to_string(self.temp.path().join("publish.log")).unwrap_or_default()
    }
}

fn write_executable(path: &Path, body: &str) {
    fs::write(path, body).unwrap();
    let mut permissions = fs::metadata(path).unwrap().permissions();
    permissions.set_mode(0o755);
    fs::set_permissions(path, permissions).unwrap();
}

const FAKE_CARGO: &str = r#"#!/usr/bin/env bash
set -euo pipefail
case "${1:-}" in
  metadata)
    printf '%s\n' '{"packages":[{"name":"rosalind-build-info","version":"0.1.0"},{"name":"rosalind-receipt","version":"0.3.0"},{"name":"rosalind-bio","version":"0.4.0"}]}'
    ;;
  package|publish)
    action="$1"; shift
    package=""
    while [ "$#" -gt 0 ]; do
      if [ "$1" = -p ]; then package="$2"; shift 2; else shift; fi
    done
    case "$package" in
      rosalind-build-info) version=0.1.0 ;;
      rosalind-receipt) version=0.3.0 ;;
      rosalind-bio) version=0.4.0 ;;
      *) exit 2 ;;
    esac
    mkdir -p target/package "$FAKE_REGISTRY"
    archive="target/package/$package-$version.crate"
    if [ "$action" = package ]; then
      printf '%s\n' "$package-$version" > "$archive"
    else
      printf '%s\n' "$package" >> "$FAKE_PUBLISH_LOG"
      cp "$archive" "$FAKE_REGISTRY/$package-$version.crate"
      if [ "${FAKE_PUBLISH_RACE_PACKAGE:-}" = "$package" ]; then exit 1; fi
    fi
    ;;
  *) exit 2 ;;
esac
"#;

const FAKE_CURL: &str = r#"#!/usr/bin/env bash
set -euo pipefail
url=""; output=""
while [ "$#" -gt 0 ]; do
  case "$1" in
    -o) output="$2"; shift 2 ;;
    -H) shift 2 ;;
    --*) shift ;;
    *) url="$1"; shift ;;
  esac
done
trimmed="${url%/download}"
version="${trimmed##*/}"
parent="${trimmed%/*}"
package="${parent##*/}"
archive="$FAKE_REGISTRY/$package-$version.crate"
test -f "$archive" || exit 22
if [ -n "$output" ]; then
  cp "$archive" "$output"
else
  checksum=$(sha256sum "$archive" | awk '{print $1}')
  printf '{"version":{"checksum":"%s"}}\n' "$checksum"
fi
"#;

#[test]
fn absent_versions_publish_in_dependency_order() {
    let fake = FakeRegistry::new();
    let output = fake.run(None);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(!String::from_utf8_lossy(&output.stdout).contains("fake-token"));
    assert!(!String::from_utf8_lossy(&output.stderr).contains("fake-token"));
    assert_eq!(
        fake.publish_log().lines().collect::<Vec<_>>(),
        PACKAGES.iter().map(|(name, _)| *name).collect::<Vec<_>>()
    );
    let state = fake.state();
    for (package, _) in PACKAGES {
        assert_eq!(state["packages"][package]["disposition"], "published");
    }
}

#[test]
fn exact_existing_versions_resume_without_upload() {
    let fake = FakeRegistry::new();
    for (package, version) in PACKAGES {
        fake.exact_archive(package, version);
    }
    let output = fake.run(None);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(fake.publish_log().is_empty());
    let state = fake.state();
    for (package, _) in PACKAGES {
        assert_eq!(
            state["packages"][package]["disposition"],
            "verified-existing"
        );
    }
}

#[test]
fn checksum_conflict_fails_permanently() {
    let fake = FakeRegistry::new();
    fake.archive("rosalind-build-info", "0.1.0", "different bytes\n");
    let output = fake.run(None);
    assert_eq!(output.status.code(), Some(5));
    assert!(String::from_utf8_lossy(&output.stderr).contains("conflict"));
    assert!(fake.publish_log().is_empty());
}

#[test]
fn partial_publication_resumes_remaining_packages() {
    let fake = FakeRegistry::new();
    fake.exact_archive("rosalind-build-info", "0.1.0");
    let output = fake.run(None);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        fake.publish_log().lines().collect::<Vec<_>>(),
        vec!["rosalind-receipt", "rosalind-bio"]
    );
    let state = fake.state();
    assert_eq!(
        state["packages"]["rosalind-build-info"]["disposition"],
        "verified-existing"
    );
}

#[test]
fn concurrent_upload_failure_is_accepted_only_after_exact_requery() {
    let fake = FakeRegistry::new();
    let output = fake.run(Some("rosalind-receipt"));
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    for (package, version) in PACKAGES {
        assert!(fake
            .registry
            .join(format!("{package}-{version}.crate"))
            .exists());
    }
}
