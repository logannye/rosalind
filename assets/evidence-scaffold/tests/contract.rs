use std::process::Command;

#[test]
fn binary_supports_native_and_persisted_inputs() {
    let result = Command::new(env!("CARGO_BIN_EXE___PACKAGE_NAME__"))
        .args(["run", "--help"])
        .output()
        .unwrap();
    assert!(result.status.success());
    let help = String::from_utf8(result.stdout).unwrap();
    for flag in [
        "--alignments",
        "--dataset",
        "--sites",
        "--output",
        "--enforce",
    ] {
        assert!(help.contains(flag), "missing {flag}");
    }
}

#[test]
fn conflicting_sources_are_rejected_before_creating_output() {
    let output = std::env::temp_dir().join(format!(
        "__PACKAGE_NAME__-invalid-{}.tsv",
        std::process::id()
    ));
    let result = Command::new(env!("CARGO_BIN_EXE___PACKAGE_NAME__"))
        .args([
            "run",
            "--alignments",
            "missing.bam",
            "--dataset",
            "missing.json",
            "--output",
        ])
        .arg(&output)
        .output()
        .unwrap();
    assert_eq!(result.status.code(), Some(2));
    assert!(!output.exists());
}
