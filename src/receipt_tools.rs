//! Builder-facing receipt inspection, portable sanitization, and standards export.

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

use anyhow::{anyhow, bail, Context, Result};

use crate::provenance::{
    blake3_file, command_template_tokens, ArtifactEvidence, CertificateEvidence, ReproReceipt,
    RunManifest, TrustReport,
};
use crate::util::atomic::write_atomic;

/// Shared inspection result used by `receipt inspect` and local tooling.
#[derive(Debug, Clone)]
pub struct ReceiptInspection {
    /// Independent trust dimensions.
    pub trust: TrustReport,
    /// Recorded artifact digests not present among supplied files.
    pub missing_artifacts: Vec<String>,
    /// Content digests computed for supplied files.
    pub supplied_artifacts: Vec<String>,
    /// Non-failing privacy or compatibility notes.
    pub warnings: Vec<String>,
}

impl ReceiptInspection {
    /// Stable JSON report.
    pub fn to_json(&self) -> String {
        let strings = |values: &[String]| {
            format!(
                "[{}]",
                values
                    .iter()
                    .map(|value| format!("\"{}\"", escape(value)))
                    .collect::<Vec<_>>()
                    .join(",")
            )
        };
        format!(
            "{{\"schema\":1,\"trust\":{},\"missing_artifacts\":{},\"supplied_artifacts\":{},\"warnings\":{}}}",
            self.trust.to_json(),
            strings(&self.missing_artifacts),
            strings(&self.supplied_artifacts),
            strings(&self.warnings),
        )
    }
}

/// Inspect one receipt, content-match supplied artifacts, and validate an optional
/// reproduction certificate against the receipt's claim ID.
pub fn inspect_receipt(
    manifest_path: &Path,
    artifacts: &[PathBuf],
    certificate_path: Option<&Path>,
) -> Result<ReceiptInspection> {
    let text = std::fs::read_to_string(manifest_path)
        .with_context(|| format!("failed to read receipt {}", manifest_path.display()))?;
    let manifest = RunManifest::from_canonical_json(&text).map_err(|error| {
        anyhow!(
            "failed to parse receipt {}: {error}",
            manifest_path.display()
        )
    })?;

    let mut supplied = BTreeSet::new();
    for artifact in artifacts {
        supplied.insert(
            blake3_file(artifact)
                .with_context(|| format!("failed to hash artifact {}", artifact.display()))?,
        );
    }
    let expected = manifest
        .inputs
        .iter()
        .chain(&manifest.outputs)
        .map(|file| file.blake3.clone())
        .collect::<BTreeSet<_>>();
    let missing_artifacts = expected.difference(&supplied).cloned().collect::<Vec<_>>();
    let artifact_evidence = if artifacts.is_empty() {
        ArtifactEvidence::NotChecked
    } else if missing_artifacts.is_empty() {
        ArtifactEvidence::Complete
    } else {
        ArtifactEvidence::Incomplete
    };

    let certificate = certificate_path.map(|path| {
        std::fs::read_to_string(path)
            .with_context(|| format!("failed to read certificate {}", path.display()))
            .and_then(|text| {
                ReproReceipt::from_canonical_json(&text).map_err(|error| {
                    anyhow!("failed to parse certificate {}: {error}", path.display())
                })
            })
    });
    let certificate_evidence = match certificate.as_ref() {
        None => CertificateEvidence::NotSupplied,
        Some(Ok(certificate)) => CertificateEvidence::Parsed(certificate),
        Some(Err(error)) => CertificateEvidence::Invalid(error.to_string()),
    };
    let trust = TrustReport::evaluate(&manifest, artifact_evidence, certificate_evidence);

    Ok(ReceiptInspection {
        trust,
        missing_artifacts,
        supplied_artifacts: supplied.into_iter().collect(),
        warnings: vec![
            "Recorded paths are relocatable metadata and are not part of schema-3+ claim hashes."
                .to_string(),
        ],
    })
}

/// Replace schema-3+ recorded paths with stable role labels while proving the claim
/// hash remains unchanged. Claim-bearing extension parameters are intentionally kept.
pub fn sanitize_receipt(manifest_path: &Path, output: &Path) -> Result<String> {
    if output.exists() {
        bail!(
            "sanitized receipt already exists: {} (refusing to overwrite)",
            output.display()
        );
    }
    let text = std::fs::read_to_string(manifest_path)
        .with_context(|| format!("failed to read receipt {}", manifest_path.display()))?;
    let mut manifest = RunManifest::from_canonical_json(&text).map_err(|error| {
        anyhow!(
            "failed to parse receipt {}: {error}",
            manifest_path.display()
        )
    })?;
    let schema = manifest
        .params
        .get("schema_version")
        .and_then(|value| value.parse::<u32>().ok())
        .ok_or_else(|| anyhow!("receipt has no valid schema_version"))?;
    if !(3..=5).contains(&schema) {
        bail!("receipt sanitize supports schemas 3 through 5, found schema {schema}");
    }
    if manifest.self_hash_ok() != Some(true) || manifest.measurement_hash_ok() == Some(false) {
        bail!("refusing to sanitize a receipt whose protected hashes are not intact");
    }
    let before = manifest.content_hash();
    let input_roles = (0..manifest.inputs.len())
        .map(|index| artifact_role(&manifest, "input", index))
        .collect::<Vec<_>>();
    let output_roles = (0..manifest.outputs.len())
        .map(|index| artifact_role(&manifest, "output", index))
        .collect::<Vec<_>>();
    for (file, role) in manifest.inputs.iter_mut().zip(input_roles) {
        file.path = format!("role://input/{role}");
    }
    for (file, role) in manifest.outputs.iter_mut().zip(output_roles) {
        file.path = format!("role://output/{role}");
    }
    let after = manifest.content_hash();
    if before != after {
        bail!("sanitization changed the claim hash; no output was written");
    }
    write_atomic(output, manifest.to_canonical_json().as_bytes(), false)
        .with_context(|| format!("failed to write sanitized receipt {}", output.display()))?;
    Ok(before)
}

/// Export an intact native schema-5 run receipt as an unsigned in-toto Statement v1.
pub fn export_intoto(manifest_path: &Path, output: &Path) -> Result<String> {
    if output.exists() {
        bail!("in-toto output already exists: {}", output.display());
    }
    let text = std::fs::read_to_string(manifest_path)
        .with_context(|| format!("failed to read receipt {}", manifest_path.display()))?;
    let manifest = RunManifest::from_canonical_json(&text).map_err(|error| {
        anyhow!(
            "failed to parse receipt {}: {error}",
            manifest_path.display()
        )
    })?;
    if manifest.self_hash_ok() != Some(true) || manifest.measurement_hash_ok() == Some(false) {
        bail!("refusing to export a receipt whose protected hashes are not intact");
    }
    let schema = manifest
        .params
        .get("schema_version")
        .and_then(|value| value.parse::<u32>().ok())
        .ok_or_else(|| anyhow!("receipt has no valid schema_version"))?;
    if schema != 5 {
        bail!("in-toto export currently requires native receipt schema 5");
    }
    let argv = command_template_tokens(&manifest)
        .map_err(|error| anyhow!("receipt has no valid tokenized invocation: {error}"))?;
    let argv_json = json_strings(&argv);
    let subjects = manifest
        .outputs
        .iter()
        .enumerate()
        .map(|(index, file)| {
            format!(
                "{{\"name\":\"{}\",\"digest\":{{\"blake3\":\"{}\"}}}}",
                escape(&artifact_role(&manifest, "output", index)),
                escape(&file.blake3)
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    let materials = manifest
        .inputs
        .iter()
        .enumerate()
        .map(|(index, file)| {
            format!(
                "{{\"uri\":\"{}\",\"digest\":{{\"blake3\":\"{}\"}}}}",
                escape(&artifact_role(&manifest, "input", index)),
                escape(&file.blake3)
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    let identity = selected_params(
        &manifest,
        &[
            "producer.name",
            "producer.version",
            "producer.repository",
            "producer.binary",
            "analyzer.id",
            "analyzer.version",
            "code_git_sha",
            "code_dirty",
            "rustc_version",
            "target_triple",
            "deps_lock_blake3",
        ],
    );
    let measurements = string_map_json(&manifest.measurements);
    let roles = manifest
        .inputs
        .iter()
        .enumerate()
        .map(|(index, _)| {
            format!(
                "{{\"direction\":\"input\",\"index\":{index},\"role\":\"{}\"}}",
                escape(&artifact_role(&manifest, "input", index))
            )
        })
        .chain(manifest.outputs.iter().enumerate().map(|(index, _)| {
            format!(
                "{{\"direction\":\"output\",\"index\":{index},\"role\":\"{}\"}}",
                escape(&artifact_role(&manifest, "output", index))
            )
        }))
        .collect::<Vec<_>>()
        .join(",");
    let statement = format!(
        "{{\"_type\":\"https://in-toto.io/Statement/v1\",\"subject\":[{subjects}],\"predicateType\":\"https://logannye.github.io/rosalind/schema/run-attestation-v1\",\"predicate\":{{\"nativeClaim\":\"{}\",\"nativeReceiptSchema\":{schema},\"subcommand\":\"{}\",\"identity\":{identity},\"invocation\":{{\"argv\":{argv_json}}},\"artifactRoles\":[{roles}],\"materials\":[{materials}],\"resourceMeasurements\":{measurements},\"enforcementAssurance\":{}}}}}",
        manifest.content_hash(),
        escape(&manifest.subcommand),
        optional_json_string(manifest.params.get("contract.assurance").map(String::as_str)),
    );
    write_atomic(output, statement.as_bytes(), false)
        .with_context(|| format!("failed to write in-toto statement {}", output.display()))?;
    Ok(manifest.content_hash())
}

fn artifact_role(manifest: &RunManifest, direction: &str, index: usize) -> String {
    manifest
        .params
        .get(&format!("artifact.{direction}.{index}.role"))
        .cloned()
        .unwrap_or_else(|| format!("{direction}-{index}"))
}

fn selected_params(manifest: &RunManifest, keys: &[&str]) -> String {
    let entries = keys
        .iter()
        .filter_map(|key| {
            manifest
                .params
                .get(*key)
                .map(|value| format!("\"{}\":\"{}\"", escape(key), escape(value)))
        })
        .collect::<Vec<_>>()
        .join(",");
    format!("{{{entries}}}")
}

fn string_map_json(values: &std::collections::BTreeMap<String, String>) -> String {
    let entries = values
        .iter()
        .map(|(key, value)| format!("\"{}\":\"{}\"", escape(key), escape(value)))
        .collect::<Vec<_>>()
        .join(",");
    format!("{{{entries}}}")
}

fn json_strings(values: &[String]) -> String {
    format!(
        "[{}]",
        values
            .iter()
            .map(|value| format!("\"{}\"", escape(value)))
            .collect::<Vec<_>>()
            .join(",")
    )
}

fn optional_json_string(value: Option<&str>) -> String {
    value
        .map(|value| format!("\"{}\"", escape(value)))
        .unwrap_or_else(|| "null".to_string())
}

fn escape(value: &str) -> String {
    let mut output = String::with_capacity(value.len());
    for character in value.chars() {
        match character {
            '"' => output.push_str("\\\""),
            '\\' => output.push_str("\\\\"),
            '\n' => output.push_str("\\n"),
            '\r' => output.push_str("\\r"),
            '\t' => output.push_str("\\t"),
            value if (value as u32) < 0x20 => {
                output.push_str(&format!("\\u{:04x}", value as u32));
            }
            value => output.push(value),
        }
    }
    output
}
