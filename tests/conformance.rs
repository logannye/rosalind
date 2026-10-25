use std::path::PathBuf;
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

fn rosalind_bin() -> &'static str {
    env!("CARGO_BIN_EXE_rosalind")
}

#[test]
fn generated_analyzer_passes_the_embedded_conformance_harness() {
    static COUNTER: AtomicU64 = AtomicU64::new(0);
    let root = std::env::temp_dir().join(format!(
        "rosalind-generated-conformance-{}-{}",
        std::process::id(),
        COUNTER.fetch_add(1, Ordering::Relaxed)
    ));
    let project = root.join("depth-check");
    rosalind::scaffold::create_analyzer_project("depth-check", &project).unwrap();

    // Published scaffolds pin crates.io exactly. This repository acceptance test
    // substitutes only local paths so it can run before publication and offline.
    let manifest_path = project.join("Cargo.toml");
    let manifest = std::fs::read_to_string(&manifest_path).unwrap();
    let workspace = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let local = manifest
        .replace(
            &format!(
                "rosalind-bio = {{ version = \"={}\", features = [\"contract-testkit\"] }}",
                env!("CARGO_PKG_VERSION")
            ),
            &format!(
                "rosalind-bio = {{ path = {:?}, features = [\"contract-testkit\"] }}",
                workspace
            ),
        )
        .replace(
            "rosalind-build-info = \"=0.1.0\"",
            &format!(
                "rosalind-build-info = {{ path = {:?} }}",
                workspace.join("crates/build-info")
            ),
        );
    std::fs::write(&manifest_path, local).unwrap();
    // Reuse the acceptance job's target directory. Nested Cargo otherwise
    // honors CARGO_TARGET_DIR while the binary lookup below assumes a private
    // scaffold target, making configured hosted builds fail after compiling.
    let target = std::env::var_os("CARGO_TARGET_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|| workspace.join("target"));
    let target = if target.is_absolute() {
        target
    } else {
        std::env::current_dir().unwrap().join(target)
    };
    let build = Command::new("cargo")
        .args(["build", "--offline"])
        .arg("--target-dir")
        .arg(&target)
        .current_dir(&project)
        .output()
        .unwrap();
    assert!(
        build.status.success(),
        "generated analyzer did not build:\n{}",
        String::from_utf8_lossy(&build.stderr)
    );

    let analyzer = target
        .join("debug")
        .join(format!("depth-check{}", std::env::consts::EXE_SUFFIX));
    let conformance = Command::new(rosalind_bin())
        .args(["conformance", "analyzer", "--binary"])
        .arg(&analyzer)
        .arg("--json")
        .output()
        .unwrap();
    let stdout = String::from_utf8_lossy(&conformance.stdout);
    assert!(
        conformance.status.success(),
        "conformance failed:\nstdout={stdout}\nstderr={}",
        String::from_utf8_lossy(&conformance.stderr)
    );
    assert!(stdout.contains("\"passed\":true"), "{stdout}");
    assert!(
        stdout.contains("\"external_reproduction\":true"),
        "{stdout}"
    );
    assert!(
        stdout.contains("\"forced_breach_partial\":true"),
        "{stdout}"
    );
    std::fs::remove_dir_all(root).ok();
}
