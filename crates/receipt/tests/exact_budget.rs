use rosalind_receipt::{verify_receipt, RunManifest, TrustState, VerifyOpts};

fn report(bytes: Option<&str>, mib: Option<&str>, peak: u64) -> rosalind_receipt::VerifyReport {
    let mut receipt = RunManifest::new("budget-test");
    for (key, value) in [("memory_budget_bytes", bytes), ("memory_budget_mb", mib)] {
        if let Some(value) = value {
            receipt.params.insert(key.into(), value.into());
        }
    }
    receipt.record_measurement("peak_rss_bytes", peak.to_string());
    receipt.finalize();
    verify_receipt(&receipt.to_canonical_json(), &VerifyOpts::default())
}

#[test]
fn exact_byte_boundary_is_shared_by_verifier_and_trust() {
    let within = report(Some("1048577"), None, 1048577);
    assert!(within.ok, "{:?}", within.problems);
    assert_eq!(within.trust.resource_contract.state, TrustState::Satisfied);
    assert_eq!(within.trust.budget_bytes, Some(1048577));
    assert_eq!(within.trust.budget_mb, None);
    let over = report(Some("1048577"), None, 1048578);
    assert!(!over.ok);
    assert_eq!(over.trust.resource_contract.state, TrustState::Failed);
}

#[test]
fn legacy_mib_and_matching_dual_claims_remain_valid() {
    for bytes in [None, Some("1048576")] {
        let report = report(bytes, Some("1"), 1048576);
        assert!(report.ok, "{:?}", report.problems);
        assert_eq!(report.trust.budget_bytes, Some(1048576));
        assert_eq!(report.trust.budget_mb, Some(1));
    }
}

#[test]
fn conflicting_malformed_and_overflowing_budgets_fail() {
    for (bytes, mib) in [
        (Some("1048577"), Some("1")),
        (Some("not-a-budget"), None),
        (Some("1048576"), Some("invalid")),
        (None, Some("18446744073709551615")),
    ] {
        let report = report(bytes, mib, 1);
        assert!(!report.ok, "{bytes:?} / {mib:?}");
        assert_eq!(report.trust.resource_contract.state, TrustState::Failed);
        assert!(report.manifest.unwrap().memory_budget_bytes().is_err());
    }
}
