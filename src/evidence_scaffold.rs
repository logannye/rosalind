//! Standalone batch evidence analyzers using Rosalind's artifact lifecycle.

use crate::scaffold::ScaffoldReport;
use std::path::{Path, PathBuf};

/// Create a standalone evidence analyzer. Existing files are never overwritten.
/// The legacy column scaffold remains [`crate::scaffold::create_analyzer_project`].
pub fn create_evidence_analyzer_project(
    name: &str,
    output: &Path,
) -> std::io::Result<ScaffoldReport> {
    let valid = !name.is_empty()
        && name.as_bytes()[0].is_ascii_lowercase()
        && name
            .bytes()
            .all(|b| b.is_ascii_lowercase() || b.is_ascii_digit() || b == b'-')
        && !name.ends_with('-')
        && !name.contains("--");
    if !valid {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "analyzer name must be lowercase kebab-case and start with a letter",
        ));
    }
    if output.exists() && std::fs::read_dir(output)?.next().is_some() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::AlreadyExists,
            format!("destination is not empty: {}", output.display()),
        ));
    }
    let files = [
        (
            "Cargo.toml",
            include_str!("../assets/evidence-scaffold/Cargo.toml.template"),
        ),
        ("build.rs", include_str!("../assets/scaffold/build.rs")),
        (
            "src/main.rs",
            include_str!("../assets/evidence-scaffold/src/main.rs"),
        ),
        (
            "tests/contract.rs",
            include_str!("../assets/evidence-scaffold/tests/contract.rs"),
        ),
        (
            "scripts/contract-check.sh",
            include_str!("../assets/evidence-scaffold/scripts/contract-check.sh"),
        ),
        (
            ".github/workflows/ci.yml",
            include_str!("../assets/evidence-scaffold/.github/workflows/ci.yml"),
        ),
        (
            "README.md",
            include_str!("../assets/evidence-scaffold/README.md"),
        ),
    ];
    for (relative, template) in files {
        let path = output.join(relative);
        std::fs::create_dir_all(path.parent().expect("scaffold file has parent"))?;
        let body = template
            .replace("__PACKAGE_NAME__", name)
            .replace("__ROSALIND_VERSION__", env!("CARGO_PKG_VERSION"));
        crate::util::atomic::write_atomic(&path, body.as_bytes(), false)?;
    }
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        std::fs::set_permissions(
            output.join("scripts/contract-check.sh"),
            std::fs::Permissions::from_mode(0o755),
        )?;
    }
    Ok(ScaffoldReport {
        root: output.to_path_buf(),
        files: files.iter().map(|(p, _)| PathBuf::from(p)).collect(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn generates_separate_artifact_project_without_overwriting() {
        let root =
            std::env::temp_dir().join(format!("rosalind-evidence-scaffold-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&root);
        let result = create_evidence_analyzer_project("candidate-qc", &root).unwrap();
        assert_eq!(result.files.len(), 7);
        for path in result.files {
            assert!(root.join(path).is_file());
        }
        let main = std::fs::read_to_string(root.join("src/main.rs")).unwrap();
        assert!(main.contains("run_evidence_artifact"));
        assert!(!main.contains("EvidenceEngine::open"));
        assert!(!main.contains("__PACKAGE_NAME__"));
        let manifest = std::fs::read_to_string(root.join("Cargo.toml")).unwrap();
        assert!(manifest.contains(&format!("version = \"={}\"", env!("CARGO_PKG_VERSION"))));
        assert_eq!(
            create_evidence_analyzer_project("candidate-qc", &root)
                .unwrap_err()
                .kind(),
            std::io::ErrorKind::AlreadyExists
        );
        std::fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn rejects_invalid_names_before_creating_paths() {
        let root = std::env::temp_dir().join(format!(
            "rosalind-evidence-scaffold-invalid-{}",
            std::process::id()
        ));
        assert_eq!(
            create_evidence_analyzer_project("../bad", &root)
                .unwrap_err()
                .kind(),
            std::io::ErrorKind::InvalidInput
        );
        assert!(!root.exists());
    }
}
