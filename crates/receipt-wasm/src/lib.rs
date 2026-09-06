//! Client-side bindings for Rosalind Receipt Studio.

use rosalind_receipt::{
    diff_receipts as diff_core, verify_manifest_str, walk_chain, ArtifactEvidence,
    CertificateEvidence, EdgeStatus, ReceiptVerdict, ReproReceipt, RunManifest, TrustReport,
};
use wasm_bindgen::prelude::*;

/// Check receipt self-integrity. This intentionally does not claim that artifact
/// files were re-hashed or that the computation was reproduced.
#[wasm_bindgen]
pub fn verify(json: &str) -> String {
    let check = verify_manifest_str(json);
    let verdict = match check.verdict {
        ReceiptVerdict::Verified => "intact",
        ReceiptVerdict::Tampered => "tampered",
        ReceiptVerdict::Unverifiable => "unverifiable",
        ReceiptVerdict::Unparseable => "unparseable",
    };
    format!(
        "{{\"verdict\":\"{}\",\"detail\":\"{}\",\"self_hash_ok\":{},\"measurement_hash_ok\":{},\"schema_version\":{},\"subcommand\":{}}}",
        verdict,
        json_escape(&check.detail),
        tri_bool(check.self_hash_ok),
        tri_bool(check.measurement_hash_ok),
        check
            .schema_version
            .map(|v| v.to_string())
            .unwrap_or_else(|| "null".to_string()),
        optional_string(check.subcommand.as_deref()),
    )
}

/// Inspect claim identity and resource/reproduction metadata for one receipt.
#[wasm_bindgen]
pub fn inspect_receipt(json: &str) -> String {
    let check = verify_manifest_str(json);
    let Ok(manifest) = RunManifest::from_canonical_json(json) else {
        return format!(
            "{{\"ok\":false,\"integrity\":{},\"claim\":null}}",
            verify(json)
        );
    };
    let budget = manifest
        .get_recorded("memory_budget_mb")
        .and_then(|v| v.parse::<u64>().ok());
    let peak = manifest
        .get_recorded("peak_rss_bytes")
        .and_then(|v| v.parse::<u64>().ok());
    let budget_bytes = manifest.memory_budget_bytes();
    let resource = match (&budget_bytes, peak) {
        (Err(_), _) => "invalid",
        (Ok(Some(limit)), Some(bytes)) if bytes <= *limit => "within",
        (Ok(Some(_)), Some(_)) => "over",
        _ => "unknown",
    };
    let integrity = match check.verdict {
        ReceiptVerdict::Verified => "intact",
        ReceiptVerdict::Tampered => "tampered",
        ReceiptVerdict::Unverifiable => "unverifiable",
        ReceiptVerdict::Unparseable => "unparseable",
    };
    let trust = TrustReport::evaluate(
        &manifest,
        ArtifactEvidence::NotChecked,
        CertificateEvidence::NotSupplied,
    );
    format!(
        "{{\"ok\":true,\"integrity\":\"{}\",\"detail\":\"{}\",\"claim\":\"{}\",\"subcommand\":\"{}\",\"tool_version\":\"{}\",\"producer_name\":{},\"producer_version\":{},\"analyzer_id\":{},\"analyzer_version\":{},\"budget_mb\":{},\"budget_bytes\":{},\"peak_rss_bytes\":{},\"resource\":\"{}\",\"parent_claim\":{},\"reproduction_verdict\":{},\"inputs\":{},\"outputs\":{},\"input_files\":{},\"output_files\":{},\"parameters\":{},\"measurements\":{},\"trust\":{}}}",
        integrity,
        json_escape(&check.detail),
        manifest.content_hash(),
        json_escape(&manifest.subcommand),
        json_escape(&manifest.tool_version),
        optional_param(&manifest, "producer.name"),
        optional_param(&manifest, "producer.version"),
        optional_param(&manifest, "analyzer.id"),
        optional_param(&manifest, "analyzer.version"),
        budget.map_or_else(|| "null".to_string(), |v| v.to_string()),
        budget_bytes.ok().flatten().map_or_else(|| "null".to_string(), |v| v.to_string()),
        peak.map_or_else(|| "null".to_string(), |v| v.to_string()),
        resource,
        optional_param(&manifest, "parent_claim"),
        optional_param(&manifest, "verdict"),
        manifest.inputs.len(),
        manifest.outputs.len(),
        files_json(&manifest.inputs),
        files_json(&manifest.outputs),
        string_map_json(&manifest.params),
        string_map_json(&manifest.measurements),
        trust.to_json(),
    )
}

/// Evaluate the shared trust model after JavaScript content-matches artifacts and
/// optionally supplies a reproduction certificate. `artifact_status` is one of
/// `not-checked`, `complete`, or `incomplete`.
#[wasm_bindgen]
pub fn evaluate_trust(
    receipt_json: &str,
    artifact_status: &str,
    certificate_json: &str,
) -> String {
    let manifest = match RunManifest::from_canonical_json(receipt_json) {
        Ok(manifest) => manifest,
        Err(error) => return TrustReport::unparseable(error.to_string()).to_json(),
    };
    let artifacts = match artifact_status {
        "complete" => ArtifactEvidence::Complete,
        "incomplete" => ArtifactEvidence::Incomplete,
        _ => ArtifactEvidence::NotChecked,
    };
    let certificate = if certificate_json.trim().is_empty() {
        None
    } else {
        Some(ReproReceipt::from_canonical_json(certificate_json))
    };
    let evidence = match certificate.as_ref() {
        None => CertificateEvidence::NotSupplied,
        Some(Ok(certificate)) => CertificateEvidence::Parsed(certificate),
        Some(Err(error)) => CertificateEvidence::Invalid(error.to_string()),
    };
    TrustReport::evaluate(&manifest, artifacts, evidence).to_json()
}

/// Diff two receipts using the same cause/effect/noise localizer as the CLI.
#[wasm_bindgen]
pub fn diff_receipts(a: &str, b: &str) -> String {
    let Ok(a) = RunManifest::from_canonical_json(a) else {
        return "{\"ok\":false,\"error\":\"first receipt is not parseable\"}".to_string();
    };
    let Ok(b) = RunManifest::from_canonical_json(b) else {
        return "{\"ok\":false,\"error\":\"second receipt is not parseable\"}".to_string();
    };
    let diff = diff_core(&a, &b);
    format!(
        "{{\"ok\":true,\"claims_identical\":{},\"verdict\":\"{}\",\"inputs\":{},\"outputs\":{},\"code_identity\":{},\"science_params\":{},\"measurements\":{}}}",
        diff.claims_identical,
        json_escape(&diff.verdict()),
        operand_changes_json(&diff.inputs),
        operand_changes_json(&diff.outputs),
        field_changes_json(&diff.code_identity),
        field_changes_json(&diff.science_params),
        field_changes_json(&diff.measurements),
    )
}

/// Walk newline-delimited canonical receipts as a provenance DAG. Canonical
/// receipts contain no literal newlines, making JSONL a dependency-free bridge.
#[wasm_bindgen]
pub fn walk_receipt_chain(receipts_jsonl: &str) -> String {
    let mut receipts = Vec::new();
    for (line, text) in receipts_jsonl.lines().filter(|line| !line.trim().is_empty()).enumerate() {
        match RunManifest::from_canonical_json(text) {
            Ok(receipt) => receipts.push(receipt),
            Err(_) => {
                return format!(
                    "{{\"ok\":false,\"error\":\"receipt {} is not parseable\"}}",
                    line + 1
                )
            }
        }
    }
    let report = walk_chain(&receipts);
    let nodes = report
        .nodes
        .iter()
        .map(|node| {
            format!(
                "{{\"id\":\"{}\",\"subcommand\":\"{}\",\"self_hash_ok\":{}}}",
                node.id,
                json_escape(&node.subcommand),
                tri_bool(node.self_hash)
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    let mut chain_intact = report.intact;
    let mut edge_values = report
        .edges
        .iter()
        .map(|edge| {
            let (status, parent) = match &edge.status {
                EdgeStatus::Resolved { parent_id } => {
                    ("resolved", format!("\"{}\"", parent_id))
                }
                EdgeStatus::External => ("external", "null".to_string()),
                EdgeStatus::Broken => ("broken", "null".to_string()),
            };
            format!(
                "{{\"child\":\"{}\",\"subcommand\":\"{}\",\"flag\":\"{}\",\"input_blake3\":\"{}\",\"status\":\"{}\",\"parent\":{}}}",
                edge.child_id,
                json_escape(&edge.child_subcommand),
                json_escape(&edge.flag),
                edge.input_blake3,
                status,
                parent
            )
        })
        .collect::<Vec<_>>();
    for certificate in receipts
        .iter()
        .filter(|receipt| receipt.subcommand == "reproduce")
    {
        let parent_claim = certificate.params.get("parent_claim");
        let certificate_ok = ReproReceipt::from_canonical_json(&certificate.to_canonical_json())
            .map(|certificate| {
                certificate.integrity_ok()
                    && certificate.verdict() == Some("REPRODUCED")
                    && certificate.outputs_match()
                    && parent_claim.is_some()
            })
            .unwrap_or(false);
        let resolved_parent = parent_claim.and_then(|parent| {
            receipts
                .iter()
                .find(|candidate| candidate.content_hash() == parent.as_str())
                .map(RunManifest::content_hash)
        });
        let status = if !certificate_ok {
            chain_intact = false;
            "broken"
        } else if resolved_parent.is_some() {
            "resolved"
        } else {
            "external"
        };
        edge_values.push(format!(
            "{{\"child\":\"{}\",\"subcommand\":\"reproduce\",\"flag\":\"parent_claim\",\"input_blake3\":{},\"status\":\"{}\",\"parent\":{}}}",
            certificate.content_hash(),
            optional_string(parent_claim.map(String::as_str)),
            status,
            optional_string(resolved_parent.as_deref()),
        ));
    }
    let edges = edge_values.join(",");
    format!(
        "{{\"ok\":true,\"intact\":{},\"nodes\":[{}],\"edges\":[{}]}}",
        chain_intact, nodes, edges
    )
}

/// Incremental BLAKE3 hashing for large dropped artifacts. JavaScript should
/// feed fixed-size `File.slice()` chunks rather than materializing the whole file.
#[wasm_bindgen]
pub struct Blake3Hasher {
    inner: blake3::Hasher,
}

#[wasm_bindgen]
impl Blake3Hasher {
    /// Create an empty hasher.
    #[wasm_bindgen(constructor)]
    pub fn new() -> Self {
        Self {
            inner: blake3::Hasher::new(),
        }
    }

    /// Add one artifact chunk.
    pub fn update(&mut self, bytes: &[u8]) {
        self.inner.update(bytes);
    }

    /// Return the current digest without consuming the hasher.
    pub fn finalize_hex(&self) -> String {
        self.inner.clone().finalize().to_hex().to_string()
    }
}

impl Default for Blake3Hasher {
    fn default() -> Self {
        Self::new()
    }
}

fn optional_param(manifest: &RunManifest, key: &str) -> String {
    optional_string(manifest.params.get(key).map(String::as_str))
}

fn files_json(files: &[rosalind_receipt::FileHash]) -> String {
    let body = files
        .iter()
        .map(|file| {
            format!(
                "{{\"path\":\"{}\",\"blake3\":\"{}\"}}",
                json_escape(&file.path),
                json_escape(&file.blake3)
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    format!("[{body}]")
}

fn string_map_json(map: &std::collections::BTreeMap<String, String>) -> String {
    let body = map
        .iter()
        .map(|(key, value)| {
            format!(
                "\"{}\":\"{}\"",
                json_escape(key),
                json_escape(value)
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    format!("{{{body}}}")
}

fn optional_string(value: Option<&str>) -> String {
    value
        .map(|value| format!("\"{}\"", json_escape(value)))
        .unwrap_or_else(|| "null".to_string())
}

fn tri_bool(value: Option<bool>) -> &'static str {
    match value {
        Some(true) => "true",
        Some(false) => "false",
        None => "null",
    }
}

fn field_changes_json(changes: &[rosalind_receipt::FieldChange]) -> String {
    let body = changes
        .iter()
        .map(|change| {
            format!(
                "{{\"key\":\"{}\",\"a\":{},\"b\":{}}}",
                json_escape(&change.key),
                optional_string(change.a.as_deref()),
                optional_string(change.b.as_deref())
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    format!("[{body}]")
}

fn operand_changes_json(changes: &[rosalind_receipt::OperandChange]) -> String {
    let body = changes
        .iter()
        .map(|change| {
            format!(
                "{{\"flag\":\"{}\",\"a\":{},\"b\":{}}}",
                json_escape(&change.flag),
                optional_string(change.a.as_deref()),
                optional_string(change.b.as_deref())
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    format!("[{body}]")
}

fn json_escape(value: &str) -> String {
    let mut out = String::with_capacity(value.len());
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

#[cfg(test)]
mod tests {
    use super::*;
    use rosalind_receipt::{FileHash, ReproOutput, ReproReceipt};

    fn run_receipt(depth: &str) -> RunManifest {
        let mut receipt = RunManifest::new("features");
        receipt.inputs.push(FileHash {
            path: "/tmp/input.bam".to_string(),
            blake3: "1111111111111111111111111111111111111111111111111111111111111111"
                .to_string(),
        });
        receipt.outputs.push(FileHash {
            path: "/tmp/output.tsv".to_string(),
            blake3: "2222222222222222222222222222222222222222222222222222222222222222"
                .to_string(),
        });
        receipt
            .params
            .insert("max_depth".to_string(), depth.to_string());
        receipt.record_measurement("peak_rss_bytes", "1048576");
        receipt.finalize();
        receipt
    }

    #[test]
    fn inspection_exposes_the_full_receipt_payload() {
        let receipt = run_receipt("1000");
        let inspected = inspect_receipt(&receipt.to_canonical_json());
        for expected in [
            "\"integrity\":\"intact\"",
            "\"input_files\":[",
            "\"output_files\":[",
            "\"parameters\":{",
            "\"measurements\":{",
            "\"max_depth\":\"1000\"",
        ] {
            assert!(inspected.contains(expected), "missing {expected}: {inspected}");
        }
    }

    #[test]
    fn diff_and_reproduction_chain_are_exposed() {
        let first = run_receipt("1000");
        let second = run_receipt("500");
        let diff = diff_receipts(&first.to_canonical_json(), &second.to_canonical_json());
        assert!(diff.contains("\"claims_identical\":false"));
        assert!(diff.contains("\"key\":\"max_depth\""));

        let certificate = ReproReceipt::build(
            &first.content_hash(),
            "features",
            "REPRODUCED",
            1,
            &[ReproOutput {
                role: "output[0]".to_string(),
                recorded_blake3: "same".to_string(),
                observed_blake3: "same".to_string(),
                matched: true,
            }],
            None,
            None,
        );
        let chain = walk_receipt_chain(&format!(
            "{}\n{}",
            first.to_canonical_json(),
            certificate.to_canonical_json()
        ));
        assert!(chain.contains("\"intact\":true"), "{chain}");
        assert!(chain.contains("\"flag\":\"parent_claim\""), "{chain}");
        assert!(chain.contains("\"status\":\"resolved\""), "{chain}");
    }

    #[test]
    fn trust_evaluation_links_certificate_and_artifact_evidence() {
        let receipt = run_receipt("1000");
        let certificate = ReproReceipt::build(
            &receipt.content_hash(),
            "features",
            "REPRODUCED",
            1,
            &[ReproOutput {
                role: "output[0]".to_string(),
                recorded_blake3: "same".to_string(),
                observed_blake3: "same".to_string(),
                matched: true,
            }],
            None,
            None,
        );
        let trust = evaluate_trust(
            &receipt.to_canonical_json(),
            "complete",
            &certificate.to_canonical_json(),
        );
        assert!(trust.contains("\"status\":\"complete\""), "{trust}");
        assert!(trust.contains("\"status\":\"reproduced\""), "{trust}");
    }

    #[test]
    fn streaming_hasher_matches_one_shot_blake3() {
        let mut hasher = Blake3Hasher::new();
        hasher.update(b"large ");
        hasher.update(b"artifact");
        assert_eq!(
            hasher.finalize_hex(),
            blake3::hash(b"large artifact").to_hex().to_string()
        );
    }
}
