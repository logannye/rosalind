//! `ReproReceipt` — the reproduction certificate written by `rosalind reproduce`.
//!
//! A small, content-addressed attestation that an independent party re-derived a
//! recorded result. It **chains to the original** by recording the parent receipt's
//! `content_hash` (`parent_claim`): N certificates that name the same `parent_claim` are
//! N independent confirmations — a serverless "reproducibility web".
//!
//! It is backed by a [`RunManifest`] (`subcommand = "reproduce"`) so it inherits the
//! proven canonical-JSON + self-hash + build-identity machinery: the certificate is
//! tamper-evident (a `manifest_blake3` self-hash), records the reproducing binary's
//! build-identity for free (via `finalize`), is signing-ready for the later Ed25519
//! track, and is itself checkable with `rosalind verify`.

use super::{ManifestError, RunManifest};

/// One output's recorded-vs-observed comparison, as carried by a certificate.
#[derive(Debug, Clone)]
pub struct ReproOutput {
    /// Output role label (e.g. `output[0]`).
    pub role: String,
    /// The hash the original receipt recorded.
    pub recorded_blake3: String,
    /// The hash this reproduction observed.
    pub observed_blake3: String,
    /// Whether they matched.
    pub matched: bool,
}

/// A reproduction certificate. Construct with [`ReproReceipt::build`]; serialize with
/// [`ReproReceipt::to_canonical_json`]; reload with [`ReproReceipt::from_canonical_json`].
#[derive(Debug, Clone)]
pub struct ReproReceipt {
    inner: RunManifest,
}

impl ReproReceipt {
    /// Build (and seal) a certificate for one reproduction.
    ///
    /// `parent_claim` is the original receipt's `content_hash()`. `peak_rss_bytes` is the
    /// reproducing machine's realized peak (a machine-dependent *measurement*, excluded
    /// from the certificate's content-address); `declared_budget_mb` is the original's
    /// declared budget. The reproducing binary's build-identity is stamped automatically.
    pub fn build(
        parent_claim: &str,
        parent_subcommand: &str,
        verdict: &str,
        chain_depth: u32,
        outputs: &[ReproOutput],
        peak_rss_bytes: Option<u64>,
        declared_budget_mb: Option<u64>,
    ) -> Self {
        let mut m = RunManifest::new("reproduce");
        m.params
            .insert("parent_claim".to_string(), parent_claim.to_string());
        m.params.insert(
            "parent_subcommand".to_string(),
            parent_subcommand.to_string(),
        );
        m.params.insert("verdict".to_string(), verdict.to_string());
        m.params
            .insert("chain_depth".to_string(), chain_depth.to_string());
        m.params
            .insert("out_count".to_string(), outputs.len().to_string());
        for (i, o) in outputs.iter().enumerate() {
            m.params.insert(format!("out{i}_role"), o.role.clone());
            m.params
                .insert(format!("out{i}_recorded"), o.recorded_blake3.clone());
            m.params
                .insert(format!("out{i}_observed"), o.observed_blake3.clone());
            m.params
                .insert(format!("out{i}_matched"), o.matched.to_string());
        }
        if let Some(mb) = declared_budget_mb {
            m.params
                .insert("declared_budget_mb".to_string(), mb.to_string());
        }
        // The reproducing machine's peak is a measurement (excluded from the claim hash),
        // so the certificate's content-address is cross-machine stable.
        if let Some(peak) = peak_rss_bytes {
            m.record_measurement("peak_rss_bytes", peak.to_string());
        }
        m.finalize(); // stamps reproducer build-identity, schema_version, the self-hash
        Self { inner: m }
    }

    /// The certificate's canonical JSON (what `<manifest>.repro.json` holds).
    pub fn to_canonical_json(&self) -> String {
        self.inner.to_canonical_json()
    }

    /// Parse a certificate from its canonical JSON.
    pub fn from_canonical_json(s: &str) -> Result<Self, ManifestError> {
        Ok(Self {
            inner: RunManifest::from_canonical_json(s)?,
        })
    }

    /// Whether the certificate's self-hash matches (tamper-evidence).
    pub fn self_hash_ok(&self) -> bool {
        self.inner.self_hash_ok() == Some(true)
    }

    /// The parent receipt's `content_hash` this certificate chains to.
    pub fn parent_claim(&self) -> Option<&str> {
        self.inner.params.get("parent_claim").map(String::as_str)
    }

    /// The recorded verdict (`REPRODUCED` / `DIVERGED`).
    pub fn verdict(&self) -> Option<&str> {
        self.inner.params.get("verdict").map(String::as_str)
    }

    /// The chain depth (1 for a reproduction of an original run receipt).
    pub fn chain_depth(&self) -> u32 {
        self.inner
            .params
            .get("chain_depth")
            .and_then(|v| v.parse().ok())
            .unwrap_or(0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_outputs() -> Vec<ReproOutput> {
        vec![ReproOutput {
            role: "output[0]".to_string(),
            recorded_blake3: "h".to_string(),
            observed_blake3: "h".to_string(),
            matched: true,
        }]
    }

    #[test]
    fn certificate_self_hashes_and_roundtrips() {
        let c = ReproReceipt::build(
            "a1b2",
            "variants",
            "REPRODUCED",
            1,
            &sample_outputs(),
            Some(22 * 1024 * 1024),
            Some(256),
        );
        let json = c.to_canonical_json();
        let back = ReproReceipt::from_canonical_json(&json).unwrap();
        assert_eq!(back.parent_claim(), Some("a1b2"));
        assert_eq!(back.verdict(), Some("REPRODUCED"));
        assert_eq!(back.chain_depth(), 1);
        assert!(back.self_hash_ok());
    }

    #[test]
    fn tamper_breaks_the_self_hash() {
        let c = ReproReceipt::build("a1b2", "variants", "REPRODUCED", 1, &[], None, None);
        // Edit a claim field without re-sealing — the self-hash must catch it.
        let json = c.to_canonical_json().replace("REPRODUCED", "DIVERGED");
        let back = ReproReceipt::from_canonical_json(&json).unwrap();
        assert!(!back.self_hash_ok());
    }
}
