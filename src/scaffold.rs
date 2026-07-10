//! Project scaffolding for downstream, statically linked analyzer binaries.

use std::path::{Path, PathBuf};

/// Files written by [`create_analyzer_project`].
#[derive(Debug, Clone)]
pub struct ScaffoldReport {
    /// Project root.
    pub root: PathBuf,
    /// Created paths, relative to `root`.
    pub files: Vec<PathBuf>,
}

/// Create a standalone analyzer crate. The destination must be absent or an
/// empty directory; existing files are never overwritten.
pub fn create_analyzer_project(name: &str, output: &Path) -> std::io::Result<ScaffoldReport> {
    validate_name(name)?;
    if output.exists() && std::fs::read_dir(output)?.next().is_some() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::AlreadyExists,
            format!("destination is not empty: {}", output.display()),
        ));
    }
    std::fs::create_dir_all(output.join("src"))?;
    std::fs::create_dir_all(output.join("tests"))?;
    std::fs::create_dir_all(output.join("scripts"))?;
    std::fs::create_dir_all(output.join(".github/workflows"))?;

    let crate_ident = name.replace('-', "_");
    let version = env!("CARGO_PKG_VERSION");
    let replacements = |template: &str| {
        template
            .replace("__PACKAGE_NAME__", name)
            .replace("__CRATE_IDENT__", &crate_ident)
            .replace("__ROSALIND_VERSION__", version)
    };
    let files = [
        ("Cargo.toml", replacements(CARGO_TEMPLATE)),
        ("build.rs", BUILD_TEMPLATE.to_string()),
        ("src/main.rs", replacements(MAIN_TEMPLATE)),
        ("tests/contract.rs", TEST_TEMPLATE.to_string()),
        ("scripts/contract-check.sh", replacements(CHECK_TEMPLATE)),
        (".github/workflows/ci.yml", CI_TEMPLATE.to_string()),
        ("README.md", replacements(README_TEMPLATE)),
    ];
    let mut created = Vec::new();
    for (relative, body) in files {
        let path = output.join(relative);
        crate::util::atomic::write_atomic(&path, body.as_bytes(), false)?;
        created.push(PathBuf::from(relative));
    }
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        let path = output.join("scripts/contract-check.sh");
        let mut permissions = std::fs::metadata(&path)?.permissions();
        permissions.set_mode(0o755);
        std::fs::set_permissions(path, permissions)?;
    }
    Ok(ScaffoldReport {
        root: output.to_path_buf(),
        files: created,
    })
}

fn validate_name(name: &str) -> std::io::Result<()> {
    let valid = !name.is_empty()
        && name
            .bytes()
            .all(|b| b.is_ascii_lowercase() || b.is_ascii_digit() || b == b'-')
        && name.as_bytes()[0].is_ascii_lowercase()
        && !name.ends_with('-')
        && !name.contains("--");
    if valid {
        Ok(())
    } else {
        Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "analyzer name must be lowercase kebab-case and start with a letter",
        ))
    }
}

const CARGO_TEMPLATE: &str = include_str!("../assets/scaffold/Cargo.toml");
const BUILD_TEMPLATE: &str = include_str!("../assets/scaffold/build.rs");
const MAIN_TEMPLATE: &str = include_str!("../assets/scaffold/src/main.rs");
const TEST_TEMPLATE: &str = include_str!("../assets/scaffold/tests/contract.rs");
const CHECK_TEMPLATE: &str = include_str!("../assets/scaffold/scripts/contract-check.sh");
const CI_TEMPLATE: &str = include_str!("../assets/scaffold/.github/workflows/ci.yml");
const README_TEMPLATE: &str = include_str!("../assets/scaffold/README.md");

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicU64, Ordering};

    fn temporary(name: &str) -> PathBuf {
        static COUNTER: AtomicU64 = AtomicU64::new(0);
        std::env::temp_dir().join(format!(
            "rosalind-scaffold-{name}-{}-{}",
            std::process::id(),
            COUNTER.fetch_add(1, Ordering::Relaxed)
        ))
    }

    #[test]
    fn creates_the_complete_analyzer_project_and_refuses_overwrite() {
        let destination = temporary("complete");
        let report = create_analyzer_project("depth-track", &destination).unwrap();
        for required in [
            "Cargo.toml",
            "build.rs",
            "src/main.rs",
            "tests/contract.rs",
            "scripts/contract-check.sh",
            ".github/workflows/ci.yml",
            "README.md",
        ] {
            assert!(report.files.contains(&PathBuf::from(required)));
            assert!(destination.join(required).is_file());
        }
        let cargo = std::fs::read_to_string(destination.join("Cargo.toml")).unwrap();
        assert!(cargo.contains(&format!("version = \"={}\"", env!("CARGO_PKG_VERSION"))));
        assert!(cargo.contains("contract-testkit"));
        let error = create_analyzer_project("depth-track", &destination).unwrap_err();
        assert_eq!(error.kind(), std::io::ErrorKind::AlreadyExists);
        std::fs::remove_dir_all(destination).ok();
    }

    #[test]
    fn rejects_names_that_are_not_cargo_kebab_case() {
        let destination = temporary("invalid");
        let error = create_analyzer_project("Not Valid", &destination).unwrap_err();
        assert_eq!(error.kind(), std::io::ErrorKind::InvalidInput);
        assert!(!destination.exists());
    }
}
