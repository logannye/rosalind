//! `rosalind diff` — a pure, claim-level divergence localizer over two run receipts.
//! It names WHICH hashed field differs, bucketed by causal role: inputs / code-identity /
//! science-params (causes), outputs (effect), measurements (machine-dependent noise).
//! No I/O, no htslib — wasm-portable. Compares two receipts to each other (distinct from
//! `reproduce`, which compares a receipt's recorded outputs against freshly-produced ones).

use std::collections::BTreeMap;

use crate::{command::manifest_operands, RunManifest, BUILD_IDENTITY_KEYS};

/// Params that are derived or redundant, excluded from the science-params bucket:
/// `manifest_blake3` is the self-hash (a function of everything else); `command` is the
/// recipe whose operands/opts are already surfaced by the input/output/param buckets.
const SKIP_PARAMS: &[&str] = &[
    "manifest_blake3",
    "command",
    "command_argv",
    "replay_schema",
];

/// One differing scalar claim/measurement field. `None` = absent on that side.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FieldChange {
    pub key: String,
    pub a: Option<String>,
    pub b: Option<String>,
}

/// One differing input/output operand, labeled by its CLI flag (recovered from `command`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OperandChange {
    pub flag: String,
    pub a: Option<String>,
    pub b: Option<String>,
}

/// The bucketed difference between two receipts' claims.
#[derive(Debug, Clone)]
pub struct ReceiptDiff {
    pub subcommand: Option<(String, String)>,
    pub inputs: Vec<OperandChange>,
    pub outputs: Vec<OperandChange>,
    pub code_identity: Vec<FieldChange>,
    pub science_params: Vec<FieldChange>,
    pub measurements: Vec<FieldChange>,
    /// `content_hash(a) == content_hash(b)` — the cross-machine claim addresses match.
    pub claims_identical: bool,
}

fn operand_changes(a: &RunManifest, b: &RunManifest, marker: &str) -> Vec<OperandChange> {
    let ma: BTreeMap<String, String> = manifest_operands(a, marker).into_iter().collect();
    let mb: BTreeMap<String, String> = manifest_operands(b, marker).into_iter().collect();
    let mut flags: Vec<String> = ma.keys().chain(mb.keys()).cloned().collect();
    flags.sort();
    flags.dedup();
    flags
        .into_iter()
        .filter(|f| ma.get(f) != mb.get(f))
        .map(|f| OperandChange {
            a: ma.get(&f).cloned(),
            b: mb.get(&f).cloned(),
            flag: f,
        })
        .collect()
}

fn map_changes(
    a: &BTreeMap<String, String>,
    b: &BTreeMap<String, String>,
    keep: impl Fn(&str) -> bool,
) -> Vec<FieldChange> {
    let mut keys: Vec<String> = a.keys().chain(b.keys()).cloned().collect();
    keys.sort();
    keys.dedup();
    keys.into_iter()
        .filter(|k| keep(k))
        .filter(|k| a.get(k) != b.get(k))
        .map(|k| FieldChange {
            a: a.get(&k).cloned(),
            b: b.get(&k).cloned(),
            key: k,
        })
        .collect()
}

/// Bucket the difference between two receipts' claims by causal role. Pure.
pub fn diff_receipts(a: &RunManifest, b: &RunManifest) -> ReceiptDiff {
    let subcommand = if a.subcommand != b.subcommand {
        Some((a.subcommand.clone(), b.subcommand.clone()))
    } else {
        None
    };
    ReceiptDiff {
        subcommand,
        inputs: operand_changes(a, b, "@in:"),
        outputs: operand_changes(a, b, "@out:"),
        code_identity: map_changes(&a.params, &b.params, |k| BUILD_IDENTITY_KEYS.contains(&k)),
        science_params: map_changes(&a.params, &b.params, |k| {
            !BUILD_IDENTITY_KEYS.contains(&k) && !SKIP_PARAMS.contains(&k)
        }),
        measurements: map_changes(&a.measurements, &b.measurements, |k| {
            k != "measurement_blake3"
        }),
        claims_identical: a.content_hash() == b.content_hash(),
    }
}

impl ReceiptDiff {
    /// `0` if the claims are identical, else `1` (read/parse errors are `2`, in the CLI).
    pub fn exit_code(&self) -> i32 {
        if self.claims_identical {
            0
        } else {
            1
        }
    }

    /// Whether an upstream CAUSE differs (vs only the output effect / measurement noise).
    pub fn has_cause(&self) -> bool {
        self.subcommand.is_some()
            || !self.inputs.is_empty()
            || !self.science_params.is_empty()
            || !self.code_identity.is_empty()
    }

    /// A one-line causal localization.
    pub fn verdict(&self) -> String {
        if self.claims_identical {
            if self.measurements.is_empty() {
                return "IDENTICAL claims".to_string();
            }
            let keys: Vec<&str> = self.measurements.iter().map(|c| c.key.as_str()).collect();
            return format!(
                "IDENTICAL claims — only machine-dependent measurements differ ({})",
                keys.join(", ")
            );
        }
        if self.has_cause() {
            let mut causes: Vec<String> = Vec::new();
            if self.subcommand.is_some() {
                causes.push("subcommand".to_string());
            }
            if !self.code_identity.is_empty() {
                causes.push(format!(
                    "code-identity ({})",
                    self.code_identity
                        .iter()
                        .map(|c| c.key.as_str())
                        .collect::<Vec<_>>()
                        .join(", ")
                ));
            }
            if !self.inputs.is_empty() {
                causes.push(format!(
                    "{} input(s) ({})",
                    self.inputs.len(),
                    self.inputs
                        .iter()
                        .map(|c| c.flag.as_str())
                        .collect::<Vec<_>>()
                        .join(", ")
                ));
            }
            if !self.science_params.is_empty() {
                causes.push(format!(
                    "{} science param(s) ({})",
                    self.science_params.len(),
                    self.science_params
                        .iter()
                        .map(|c| c.key.as_str())
                        .collect::<Vec<_>>()
                        .join(", ")
                ));
            }
            let effect = if self.outputs.is_empty() {
                String::new()
            } else {
                format!("; effect: {} output(s)", self.outputs.len())
            };
            return format!("claims DIFFER — cause: {}{}", causes.join(", "), effect);
        }
        if !self.outputs.is_empty() {
            return "claims DIFFER — outputs differ with identical inputs/params/code → nondeterminism or corruption".to_string();
        }
        "claims DIFFER".to_string()
    }

    /// A compact, dependency-free JSON summary for `--json`.
    pub fn to_json(&self) -> String {
        format!(
            "{{\"claims_identical\":{},\"inputs\":{},\"outputs\":{},\"code_identity\":{},\"science_params\":{},\"measurements\":{},\"exit_code\":{}}}",
            self.claims_identical,
            self.inputs.len(),
            self.outputs.len(),
            self.code_identity.len(),
            self.science_params.len(),
            self.measurements.len(),
            self.exit_code()
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build + finalize a manifest from a command, extra params, and measurements.
    /// (Empty inputs[]/outputs[] — the operand diff reads the `command` param, and
    /// `command` is itself a claim param so content_hash reflects it.)
    fn mk(
        sub: &str,
        command: &str,
        params: &[(&str, &str)],
        measurements: &[(&str, &str)],
    ) -> RunManifest {
        let mut m = RunManifest::new(sub);
        m.params.insert("command".to_string(), command.to_string());
        for (k, v) in params {
            m.params.insert(k.to_string(), v.to_string());
        }
        for (k, v) in measurements {
            m.record_measurement(*k, *v);
        }
        m.finalize();
        m
    }

    #[test]
    fn input_change_is_a_cause_with_the_output_as_effect() {
        let a = mk(
            "variants",
            "variants --index @in:IDX --alignments @in:B1 -o @out:O1",
            &[],
            &[],
        );
        let b = mk(
            "variants",
            "variants --index @in:IDX --alignments @in:B2 -o @out:O2",
            &[],
            &[],
        );
        let d = diff_receipts(&a, &b);
        assert!(!d.claims_identical);
        assert_eq!(d.exit_code(), 1);
        assert_eq!(d.inputs.len(), 1);
        assert_eq!(d.inputs[0].flag, "--alignments");
        assert_eq!(d.inputs[0].a.as_deref(), Some("B1"));
        assert_eq!(d.inputs[0].b.as_deref(), Some("B2"));
        assert_eq!(d.outputs.len(), 1, "the -o change is the effect");
        assert!(d.has_cause());
        assert!(d.verdict().contains("--alignments"), "{}", d.verdict());
    }

    #[test]
    fn science_param_change_is_isolated_from_code_identity() {
        let a = mk(
            "features",
            "features --index @in:I -o @out:O",
            &[("max_depth", "1000")],
            &[],
        );
        let b = mk(
            "features",
            "features --index @in:I -o @out:O",
            &[("max_depth", "500")],
            &[],
        );
        let d = diff_receipts(&a, &b);
        assert_eq!(d.science_params.len(), 1);
        assert_eq!(d.science_params[0].key, "max_depth");
        assert!(d.code_identity.is_empty());
        assert!(d.inputs.is_empty());
        assert!(d.outputs.is_empty(), "same @out:O");
        assert_eq!(d.exit_code(), 1);
    }

    #[test]
    fn code_identity_drift_is_its_own_bucket() {
        // Set code_git_sha AFTER finalize (finalize overwrites build-identity keys).
        let mut a = mk("variants", "variants --index @in:I -o @out:O", &[], &[]);
        a.params
            .insert("code_git_sha".to_string(), "aaaa111".to_string());
        let mut b = mk("variants", "variants --index @in:I -o @out:O", &[], &[]);
        b.params
            .insert("code_git_sha".to_string(), "bbbb222".to_string());
        let d = diff_receipts(&a, &b);
        assert_eq!(d.code_identity.len(), 1);
        assert_eq!(d.code_identity[0].key, "code_git_sha");
        assert!(
            d.science_params.is_empty(),
            "code-identity is segregated from science params"
        );
        assert!(d.verdict().contains("code-identity"), "{}", d.verdict());
    }

    #[test]
    fn outputs_differ_with_identical_cause_flags_nondeterminism() {
        let a = mk("variants", "variants --index @in:I -o @out:O1", &[], &[]);
        let b = mk("variants", "variants --index @in:I -o @out:O2", &[], &[]);
        let d = diff_receipts(&a, &b);
        assert!(!d.claims_identical);
        assert!(d.inputs.is_empty());
        assert!(d.science_params.is_empty());
        assert!(d.code_identity.is_empty());
        assert!(!d.has_cause());
        assert_eq!(d.outputs.len(), 1);
        assert!(d.verdict().contains("nondeterminism"), "{}", d.verdict());
    }

    #[test]
    fn identical_claims_with_only_a_measurement_diff_exit_zero() {
        let a = mk(
            "features",
            "features --index @in:I -o @out:O",
            &[],
            &[("peak_rss_bytes", "100")],
        );
        let b = mk(
            "features",
            "features --index @in:I -o @out:O",
            &[],
            &[("peak_rss_bytes", "200")],
        );
        let d = diff_receipts(&a, &b);
        assert!(
            d.claims_identical,
            "measurements are excluded from the claim hash"
        );
        assert_eq!(d.exit_code(), 0);
        assert_eq!(d.measurements.len(), 1);
        assert_eq!(d.measurements[0].key, "peak_rss_bytes");
        assert!(
            d.verdict().contains("measurements differ"),
            "{}",
            d.verdict()
        );
        assert_eq!(
            d.to_json(),
            "{\"claims_identical\":true,\"inputs\":0,\"outputs\":0,\"code_identity\":0,\"science_params\":0,\"measurements\":1,\"exit_code\":0}"
        );
    }
}
