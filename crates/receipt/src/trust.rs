//! One evidence model shared by command-line, browser, and badge consumers.

use super::{json_escape, ReproReceipt, RunManifest};

/// Coarse state of one independent trust dimension.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TrustState {
    /// The supplied evidence satisfies this dimension.
    Satisfied,
    /// Evidence was not supplied or the receipt predates the capability.
    Missing,
    /// Supplied evidence is invalid, contradictory, or records a failure.
    Failed,
    /// The dimension does not apply to this receipt.
    NotApplicable,
}

impl TrustState {
    /// Stable machine-readable value.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Satisfied => "satisfied",
            Self::Missing => "missing",
            Self::Failed => "failed",
            Self::NotApplicable => "not-applicable",
        }
    }
}

/// Status and plain-language explanation for one trust dimension.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TrustFacet {
    /// Stable, dimension-specific status such as `intact` or `not-checked`.
    pub status: String,
    /// Whether that status satisfies, lacks, or contradicts evidence.
    pub state: TrustState,
    /// Human-readable explanation suitable for a "why this status?" disclosure.
    pub explanation: String,
}

impl TrustFacet {
    fn new(status: &str, state: TrustState, explanation: impl Into<String>) -> Self {
        Self {
            status: status.to_string(),
            state,
            explanation: explanation.into(),
        }
    }

    fn to_json(&self) -> String {
        format!(
            "{{\"status\":\"{}\",\"state\":\"{}\",\"explanation\":\"{}\"}}",
            json_escape(&self.status),
            self.state.as_str(),
            json_escape(&self.explanation)
        )
    }
}

/// Result of matching the receipt's inputs and outputs to supplied artifacts.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ArtifactEvidence {
    /// Artifact bytes were not supplied to this verifier.
    NotChecked,
    /// Every recorded input and output digest matched supplied bytes.
    Complete,
    /// At least one recorded artifact was absent or had a different digest.
    Incomplete,
}

/// Optional reproduction-certificate evidence.
#[derive(Debug)]
pub enum CertificateEvidence<'a> {
    /// No certificate was supplied. This is neutral, not a negative finding.
    NotSupplied,
    /// A certificate parsed and can be linked to the run receipt.
    Parsed(&'a ReproReceipt),
    /// A supplied certificate could not be parsed.
    Invalid(String),
}

/// Independent trust dimensions for one run receipt.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TrustReport {
    /// Receipt's deterministic claim ID, when derivable.
    pub claim_id: Option<String>,
    /// Claim and measurement self-hash status.
    pub receipt_integrity: TrustFacet,
    /// Whether all recorded artifact bytes were supplied and matched.
    pub artifact_completeness: TrustFacet,
    /// Budget result and enforcement evidence.
    pub resource_contract: TrustFacet,
    /// Linked byte-reproduction evidence.
    pub reproduction_evidence: TrustFacet,
    /// Signature evidence. Signing is not implemented in native schema 5.
    pub signature: TrustFacet,
    /// Recorded enforcement assurance, when present.
    pub enforcement_assurance: Option<String>,
    /// Recorded memory budget, when present and valid.
    pub budget_mb: Option<u64>,
}

impl TrustReport {
    /// Construct a report for a parse failure.
    pub fn unparseable(detail: impl Into<String>) -> Self {
        Self {
            claim_id: None,
            receipt_integrity: TrustFacet::new("unparseable", TrustState::Failed, detail),
            artifact_completeness: TrustFacet::new(
                "not-checked",
                TrustState::Missing,
                "Artifacts cannot be checked until the receipt parses.",
            ),
            resource_contract: TrustFacet::new(
                "unavailable",
                TrustState::Missing,
                "Resource evidence cannot be interpreted until the receipt parses.",
            ),
            reproduction_evidence: TrustFacet::new(
                "not-supplied",
                TrustState::Missing,
                "No valid run receipt is available to link reproduction evidence to.",
            ),
            signature: TrustFacet::new(
                "not-supplied",
                TrustState::Missing,
                "No signature evidence was supplied.",
            ),
            enforcement_assurance: None,
            budget_mb: None,
        }
    }

    /// Build a report from one run receipt and optional external evidence.
    pub fn evaluate(
        manifest: &RunManifest,
        artifacts: ArtifactEvidence,
        certificate: CertificateEvidence<'_>,
    ) -> Self {
        let measurement_required = manifest.claims_measurements();
        let receipt_integrity = match (
            manifest.self_hash_ok(),
            manifest.measurement_hash_ok(),
            measurement_required,
        ) {
            (Some(true), Some(false), _) | (Some(false), _, _) => TrustFacet::new(
                "tampered",
                TrustState::Failed,
                "A claim or measurement self-hash does not re-derive; the receipt was altered after sealing.",
            ),
            (Some(true), None, true) => TrustFacet::new(
                "tampered",
                TrustState::Failed,
                "The claim says measurements exist, but their protected block is missing.",
            ),
            (Some(true), _, _) => TrustFacet::new(
                "intact",
                TrustState::Satisfied,
                "The claim hash and every present measurement hash re-derive and match.",
            ),
            (None, _, _) => TrustFacet::new(
                "unverifiable",
                TrustState::Missing,
                "This legacy receipt has no claim self-hash, so its integrity cannot be established.",
            ),
        };

        let artifact_completeness = match artifacts {
            ArtifactEvidence::NotChecked => TrustFacet::new(
                "not-checked",
                TrustState::Missing,
                "Artifact bytes were not supplied; the receipt can be intact without proving that local files still match it.",
            ),
            ArtifactEvidence::Complete => TrustFacet::new(
                "complete",
                TrustState::Satisfied,
                "Every recorded input and output was supplied and matched by content hash.",
            ),
            ArtifactEvidence::Incomplete => TrustFacet::new(
                "incomplete-or-mismatched",
                TrustState::Failed,
                "At least one recorded artifact is missing or does not match its recorded content hash.",
            ),
        };

        let budget_mb = manifest
            .get_recorded("memory_budget_mb")
            .and_then(|value| value.parse::<u64>().ok());
        let peak = manifest
            .get_recorded("peak_rss_bytes")
            .and_then(|value| value.parse::<u64>().ok());
        let enforcement_assurance =
            manifest
                .params
                .get("contract.assurance")
                .cloned()
                .or_else(|| {
                    (manifest.params.get("enforce").map(String::as_str) == Some("true"))
                        .then(|| "declared-bound-cooperative".to_string())
                });
        let resource_contract = match (budget_mb, peak) {
            (Some(budget), Some(realized)) if realized <= budget.saturating_mul(1024 * 1024) => {
                let assurance = enforcement_assurance.as_deref().unwrap_or("observed-only");
                TrustFacet::new(
                    "within-budget",
                    TrustState::Satisfied,
                    format!(
                        "The recorded peak fits the {budget} MiB budget; enforcement assurance is {assurance}."
                    ),
                )
            }
            (Some(budget), Some(realized)) => TrustFacet::new(
                "over-budget",
                TrustState::Failed,
                format!(
                    "The recorded peak ({} MiB) exceeds the {budget} MiB budget.",
                    realized / (1 << 20)
                ),
            ),
            (None, _) => TrustFacet::new(
                "not-declared",
                TrustState::NotApplicable,
                "This run did not declare a memory budget; measurements may still be present.",
            ),
            (Some(_), None) => TrustFacet::new(
                "measurement-missing",
                TrustState::Failed,
                "A memory budget was declared, but no realized peak measurement is available.",
            ),
        };

        let reproduction_evidence = match certificate {
            CertificateEvidence::NotSupplied => TrustFacet::new(
                "not-supplied",
                TrustState::Missing,
                "No reproduction certificate was supplied; this does not mean the run is unreproducible.",
            ),
            CertificateEvidence::Invalid(detail) => TrustFacet::new(
                "invalid",
                TrustState::Failed,
                format!("The supplied reproduction certificate is invalid: {detail}"),
            ),
            CertificateEvidence::Parsed(cert) if !cert.integrity_ok() => TrustFacet::new(
                "invalid",
                TrustState::Failed,
                "The reproduction certificate's protected claim or measurements were altered.",
            ),
            CertificateEvidence::Parsed(cert)
                if cert.parent_claim() != Some(manifest.content_hash().as_str()) =>
            {
                TrustFacet::new(
                    "wrong-parent",
                    TrustState::Failed,
                    "The reproduction certificate names a different parent claim.",
                )
            }
            CertificateEvidence::Parsed(cert)
                if cert.verdict() == Some("REPRODUCED") && cert.outputs_match() =>
            {
                TrustFacet::new(
                    "reproduced",
                    TrustState::Satisfied,
                    "A valid linked certificate records byte-identical reproduced outputs.",
                )
            }
            CertificateEvidence::Parsed(cert) if cert.verdict() == Some("DIVERGED") => {
                TrustFacet::new(
                    "diverged",
                    TrustState::Failed,
                    "A valid linked certificate records output bytes that diverged.",
                )
            }
            CertificateEvidence::Parsed(_) => TrustFacet::new(
                "invalid",
                TrustState::Failed,
                "The certificate does not contain a supported REPRODUCED result with matching outputs.",
            ),
        };

        Self {
            claim_id: Some(manifest.content_hash()),
            receipt_integrity,
            artifact_completeness,
            resource_contract,
            reproduction_evidence,
            signature: TrustFacet::new(
                "not-supplied",
                TrustState::Missing,
                "Native schema 5 carries no signature; receipt integrity is hash-based, not identity-authenticated.",
            ),
            enforcement_assurance,
            budget_mb,
        }
    }

    /// Stable JSON object used inside schema-2 verification and inspection reports.
    pub fn to_json(&self) -> String {
        let claim = self
            .claim_id
            .as_ref()
            .map(|value| format!("\"{}\"", json_escape(value)))
            .unwrap_or_else(|| "null".to_string());
        let assurance = self
            .enforcement_assurance
            .as_ref()
            .map(|value| format!("\"{}\"", json_escape(value)))
            .unwrap_or_else(|| "null".to_string());
        let budget = self
            .budget_mb
            .map(|value| value.to_string())
            .unwrap_or_else(|| "null".to_string());
        format!(
            "{{\"claim_id\":{claim},\"receipt_integrity\":{},\"artifact_completeness\":{},\"resource_contract\":{},\"reproduction_evidence\":{},\"signature\":{},\"enforcement_assurance\":{assurance},\"budget_mb\":{budget}}}",
            self.receipt_integrity.to_json(),
            self.artifact_completeness.to_json(),
            self.resource_contract.to_json(),
            self.reproduction_evidence.to_json(),
            self.signature.to_json(),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ReproOutput;

    fn run() -> RunManifest {
        let mut manifest = RunManifest::new("variants");
        manifest
            .params
            .insert("memory_budget_mb".to_string(), "128".to_string());
        manifest.record_measurement("peak_rss_bytes", (32_u64 << 20).to_string());
        manifest.finalize();
        manifest
    }

    #[test]
    fn intact_receipt_without_certificate_is_neutral_about_reproduction() {
        let report = TrustReport::evaluate(
            &run(),
            ArtifactEvidence::NotChecked,
            CertificateEvidence::NotSupplied,
        );
        assert_eq!(report.receipt_integrity.status, "intact");
        assert_eq!(report.resource_contract.status, "within-budget");
        assert_eq!(report.reproduction_evidence.status, "not-supplied");
        assert_eq!(report.reproduction_evidence.state, TrustState::Missing);
    }

    #[test]
    fn valid_linked_certificate_is_reproduced() {
        let run = run();
        let output = ReproOutput {
            role: "output[0]".to_string(),
            recorded_blake3: "same".to_string(),
            observed_blake3: "same".to_string(),
            matched: true,
        };
        let certificate = ReproReceipt::build(
            &run.content_hash(),
            "variants",
            "REPRODUCED",
            1,
            &[output],
            None,
            Some(128),
        );
        let report = TrustReport::evaluate(
            &run,
            ArtifactEvidence::Complete,
            CertificateEvidence::Parsed(&certificate),
        );
        assert_eq!(report.reproduction_evidence.status, "reproduced");
        assert_eq!(report.reproduction_evidence.state, TrustState::Satisfied);
    }

    #[test]
    fn wrong_parent_certificate_fails_independently() {
        let run = run();
        let certificate = ReproReceipt::build(
            "another-claim",
            "variants",
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
        let report = TrustReport::evaluate(
            &run,
            ArtifactEvidence::NotChecked,
            CertificateEvidence::Parsed(&certificate),
        );
        assert_eq!(report.receipt_integrity.status, "intact");
        assert_eq!(report.reproduction_evidence.status, "wrong-parent");
    }
}
