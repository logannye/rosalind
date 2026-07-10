use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

use rosalind::provenance::{CommandCapture, ReproOutput, ReproReceipt, RunManifest, TrustState};

fn unique_dir() -> PathBuf {
    static COUNTER: AtomicU64 = AtomicU64::new(0);
    let id = COUNTER.fetch_add(1, Ordering::Relaxed);
    let path = std::env::temp_dir().join(format!(
        "rosalind-receipt-tools-{}-{id}",
        std::process::id()
    ));
    std::fs::create_dir_all(&path).unwrap();
    path
}

fn receipt(dir: &std::path::Path) -> (PathBuf, PathBuf, PathBuf, RunManifest) {
    let input = dir.join("private patient input.bam");
    let output = dir.join("private result.vcf");
    let manifest_path = dir.join("run.manifest.json");
    std::fs::write(&input, b"input bytes\n").unwrap();
    std::fs::write(&output, b"output bytes\n").unwrap();
    let mut manifest = RunManifest::new("variants");
    manifest.params.insert(
        "artifact.input.0.role".to_string(),
        "alignments".to_string(),
    );
    manifest
        .params
        .insert("artifact.output.0.role".to_string(), "calls".to_string());
    manifest
        .params
        .insert("memory_budget_mb".to_string(), "128".to_string());
    manifest.record_measurement("peak_rss_bytes", (16_u64 << 20).to_string());
    let mut command = CommandCapture::new("variants");
    command.input("--alignments", &input).unwrap();
    command.output("--output", &output).unwrap();
    command.record_into(&mut manifest);
    manifest.finalize();
    std::fs::write(&manifest_path, manifest.to_canonical_json()).unwrap();
    (manifest_path, input, output, manifest)
}

#[test]
fn inspect_matches_artifacts_and_a_linked_certificate_by_hash() {
    let dir = unique_dir();
    let (path, input, output, manifest) = receipt(&dir);
    let certificate = ReproReceipt::build(
        &manifest.content_hash(),
        "variants",
        "REPRODUCED",
        1,
        &[ReproOutput {
            role: "output[0]".to_string(),
            recorded_blake3: manifest.outputs[0].blake3.clone(),
            observed_blake3: manifest.outputs[0].blake3.clone(),
            matched: true,
        }],
        None,
        Some(128),
    );
    let certificate_path = dir.join("run.repro.json");
    std::fs::write(&certificate_path, certificate.to_canonical_json()).unwrap();
    let report =
        rosalind::inspect_receipt(&path, &[input, output], Some(certificate_path.as_path()))
            .unwrap();
    assert!(report.missing_artifacts.is_empty());
    assert_eq!(
        report.trust.artifact_completeness.state,
        TrustState::Satisfied
    );
    assert_eq!(
        report.trust.reproduction_evidence.state,
        TrustState::Satisfied
    );
    std::fs::remove_dir_all(dir).ok();
}

#[test]
fn sanitize_removes_paths_without_changing_claim() {
    let dir = unique_dir();
    let (path, _, _, manifest) = receipt(&dir);
    let sanitized = dir.join("sanitized.json");
    let claim = rosalind::sanitize_receipt(&path, &sanitized).unwrap();
    let parsed =
        RunManifest::from_canonical_json(&std::fs::read_to_string(&sanitized).unwrap()).unwrap();
    assert_eq!(claim, manifest.content_hash());
    assert_eq!(parsed.content_hash(), manifest.content_hash());
    assert_eq!(parsed.inputs[0].path, "role://input/alignments");
    assert_eq!(parsed.outputs[0].path, "role://output/calls");
    assert!(rosalind::sanitize_receipt(&path, &sanitized).is_err());
    std::fs::remove_dir_all(dir).ok();
}

#[test]
fn export_intoto_is_unsigned_and_carries_native_evidence() {
    let dir = unique_dir();
    let (path, _, _, manifest) = receipt(&dir);
    let statement = dir.join("run.intoto.json");
    rosalind::export_intoto(&path, &statement).unwrap();
    let text = std::fs::read_to_string(&statement).unwrap();
    for expected in [
        "https://in-toto.io/Statement/v1",
        "https://logannye.github.io/rosalind/schema/run-attestation-v1",
        &manifest.content_hash(),
        "\"nativeReceiptSchema\":5",
        "\"argv\":[\"variants\"",
        "\"role\":\"alignments\"",
        "\"role\":\"calls\"",
    ] {
        assert!(text.contains(expected), "missing {expected}: {text}");
    }
    assert!(!text.contains("signature"));
    std::fs::remove_dir_all(dir).ok();
}
