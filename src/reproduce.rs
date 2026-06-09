//! `reproduce` — re-derive a recorded result and compare it byte-for-byte. Standalone:
//! filesystem + blake3 + subprocess only (no htslib). The verdict is over OUTPUT bytes;
//! code / inputs / resource are diagnostic context. v1 supports the deterministic text
//! outputs Rosalind emits (VCF, TSV); BAM/bgzf is reported INCONCLUSIVE, never a false
//! DIVERGED.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use anyhow::{anyhow, Context, Result};

use crate::provenance::{blake3_file, CommandCapture, FileHash, RunManifest};

/// The output-byte verdict. Exit codes are stable: a stranger's CI can branch on them.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Verdict {
    /// Every recorded output re-derived byte-identically.
    Reproduced,
    /// At least one output differs from what the receipt recorded.
    Diverged,
    /// Could not run the comparison (pre-schema-5, input not located, unsupported output).
    Inconclusive,
}

impl Verdict {
    /// The process exit code for this verdict (0 / 6 / 7).
    pub fn exit_code(self) -> i32 {
        match self {
            Verdict::Reproduced => 0,
            Verdict::Diverged => 6,
            Verdict::Inconclusive => 7,
        }
    }
    /// The human-readable label.
    pub fn label(self) -> &'static str {
        match self {
            Verdict::Reproduced => "REPRODUCED",
            Verdict::Diverged => "DIVERGED",
            Verdict::Inconclusive => "INCONCLUSIVE",
        }
    }
}

/// The result of comparing produced outputs against the recorded ones.
#[derive(Debug)]
pub struct Outcome {
    /// REPRODUCED iff every recorded output matched; DIVERGED otherwise.
    pub verdict: Verdict,
    /// One human-readable line per mismatch (empty when REPRODUCED).
    pub diffs: Vec<String>,
}

/// Pure: compare produced `(role, hash)` pairs against the recorded outputs. REPRODUCED
/// iff every recorded output has a positional produced match; DIVERGED otherwise, naming
/// each mismatch.
pub fn classify_outputs(recorded: &[FileHash], produced: &[(String, String)]) -> Outcome {
    let mut diffs = Vec::new();
    for (i, rec) in recorded.iter().enumerate() {
        match produced.get(i) {
            Some((_, got)) if *got == rec.blake3 => {}
            Some((role, got)) => {
                diffs.push(format!("output {role}: recorded {} got {got}", rec.blake3))
            }
            None => diffs.push(format!(
                "output {i}: recorded {} but not produced",
                rec.blake3
            )),
        }
    }
    let verdict = if diffs.is_empty() {
        Verdict::Reproduced
    } else {
        Verdict::Diverged
    };
    Outcome { verdict, diffs }
}

/// Per-output comparison detail (carried for the reproduction certificate).
#[derive(Debug, Clone)]
pub struct OutputCmp {
    /// Output role label (e.g. `output[0]`).
    pub role: String,
    /// The hash the receipt recorded.
    pub recorded: String,
    /// The hash the re-derivation produced.
    pub observed: String,
    /// Whether they match.
    pub matched: bool,
}

/// Machine-local resource observation of the re-run (explicitly NOT cross-machine).
#[derive(Debug, Clone, Default)]
pub struct ResourceHere {
    /// Realized peak RSS of the re-run on this machine, if its receipt was readable.
    pub peak_rss_bytes: Option<u64>,
    /// The budget the original declared (MiB), if any.
    pub declared_budget_mb: Option<u64>,
    /// Whether the re-run's peak fit the declared budget here.
    pub fit: Option<bool>,
}

/// The full result of a `reproduce` run: the printable report, the exit code, and the
/// structured data the reproduction certificate is minted from.
#[derive(Debug)]
pub struct ReproReport {
    /// `REPRODUCED` | `DIVERGED` | `INCONCLUSIVE` | `TAMPERED`.
    pub verdict_label: String,
    /// Process exit code: 0 | 6 | 7 | 5 (tamper).
    pub exit_code: i32,
    /// Aligned, human-readable report lines.
    pub lines: Vec<String>,
    /// The parent receipt's `content_hash()` — the cross-machine id the certificate chains to.
    pub parent_claim: String,
    /// The parent receipt's subcommand.
    pub parent_subcommand: String,
    /// Per-output comparison (populated for REPRODUCED / DIVERGED).
    pub outputs: Vec<OutputCmp>,
    /// Machine-local resource observation of the re-run.
    pub resource_here: ResourceHere,
    /// The build-identity of the binary that performed this reproduction.
    pub reproducer_code: BTreeMap<String, String>,
    /// True only for REPRODUCED / DIVERGED (an actual byte comparison happened).
    pub compared: bool,
}

/// Output types `reproduce` can byte-compare in v1. BAM/bgzf is out of scope (a C zlib
/// not captured by `deps_lock_blake3`).
fn output_is_text(path: &str) -> bool {
    let p = path.to_ascii_lowercase();
    p.ends_with(".vcf") || p.ends_with(".tsv") || p.ends_with(".txt") || p.ends_with(".gvcf")
}

/// Build a content-hash → path index of every file directly under `inputs_dir`.
fn index_inputs(inputs_dir: &Path) -> Result<BTreeMap<String, PathBuf>> {
    let mut idx = BTreeMap::new();
    let rd = std::fs::read_dir(inputs_dir)
        .with_context(|| format!("failed to read --inputs dir {}", inputs_dir.display()))?;
    for entry in rd {
        let p = entry?.path();
        if p.is_file() {
            if let Ok(h) = blake3_file(&p) {
                idx.entry(h).or_insert(p);
            }
        }
    }
    Ok(idx)
}

/// A fresh, unique temp directory for the re-run's outputs.
fn make_temp_dir() -> Result<PathBuf> {
    static COUNTER: AtomicU64 = AtomicU64::new(0);
    let n = COUNTER.fetch_add(1, Ordering::Relaxed);
    let d = std::env::temp_dir().join(format!("rosalind-reproduce-{}-{n}", std::process::id()));
    std::fs::create_dir_all(&d)
        .with_context(|| format!("failed to create temp dir {}", d.display()))?;
    Ok(d)
}

/// The build-identity of the *current* binary (the one performing the reproduction),
/// read by finalizing a throwaway manifest (which stamps build-identity from `build.rs`).
fn current_build_identity() -> BTreeMap<String, String> {
    let mut m = RunManifest::new("reproduce-probe");
    m.finalize();
    let mut out = BTreeMap::new();
    for k in [
        "code_git_sha",
        "code_dirty",
        "rustc_version",
        "target_triple",
        "deps_lock_blake3",
    ] {
        if let Some(v) = m.params.get(k) {
            out.insert(k.to_string(), v.clone());
        }
    }
    out
}

/// Read `peak_rss_bytes` from the receipt the re-run wrote next to `output_temp`.
fn read_rerun_peak(output_temp: &Path) -> Option<u64> {
    let mut name = output_temp.as_os_str().to_os_string();
    name.push(".manifest.json");
    let text = std::fs::read_to_string(PathBuf::from(name)).ok()?;
    let m = RunManifest::from_canonical_json(&text).ok()?;
    m.get_recorded("peak_rss_bytes")
        .and_then(|v| v.parse().ok())
}

/// Re-derive the result recorded in `manifest_path` from inputs content-located under
/// `inputs_dir`, and compare byte-for-byte. See the module docs for scope.
pub fn reproduce(manifest_path: &Path, inputs_dir: &Path) -> Result<ReproReport> {
    let text = std::fs::read_to_string(manifest_path)
        .with_context(|| format!("failed to read receipt {}", manifest_path.display()))?;
    let manifest = RunManifest::from_canonical_json(&text)
        .map_err(|e| anyhow!("malformed receipt {}: {e}", manifest_path.display()))?;

    let parent_claim = manifest.content_hash();
    let parent_subcommand = manifest.subcommand.clone();
    let reproducer_code = current_build_identity();

    let report =
        |verdict_label: &str, exit_code: i32, lines: Vec<String>, compared: bool| ReproReport {
            verdict_label: verdict_label.to_string(),
            exit_code,
            lines,
            parent_claim: parent_claim.clone(),
            parent_subcommand: parent_subcommand.clone(),
            outputs: Vec::new(),
            resource_here: ResourceHere::default(),
            reproducer_code: reproducer_code.clone(),
            compared,
        };

    // 1. Integrity gate (tamper-evidence). Budget-fit is verify's job, NOT reproduce's —
    //    a receipt that honestly records an over-budget run is still reproducible.
    if manifest.self_hash_ok() == Some(false) || manifest.measurement_hash_ok() == Some(false) {
        return Ok(report(
            "TAMPERED",
            5,
            vec![
                "  the receipt's self-hash does not match — it was modified after it was written"
                    .to_string(),
                "  VERDICT     : TAMPERED".to_string(),
            ],
            false,
        ));
    }

    // 2. Need a recorded command (schema >= 5) to replay.
    let command = match manifest.params.get("command") {
        Some(c) => c.clone(),
        None => {
            return Ok(report(
                "INCONCLUSIVE",
                7,
                vec![
                    "  pre-schema-5 receipt: no recorded command to replay".to_string(),
                    "  VERDICT     : INCONCLUSIVE".to_string(),
                ],
                false,
            ))
        }
    };

    // 3. Supported, comparable outputs (text only in v1; at least one to compare).
    if manifest.outputs.is_empty() {
        return Ok(report(
            "INCONCLUSIVE",
            7,
            vec![
                "  receipt records no outputs to compare (e.g. a stdout run)".to_string(),
                "  VERDICT     : INCONCLUSIVE".to_string(),
            ],
            false,
        ));
    }
    if let Some(o) = manifest.outputs.iter().find(|o| !output_is_text(&o.path)) {
        return Ok(report(
            "INCONCLUSIVE",
            7,
            vec![
                format!(
                    "  output {} is not byte-comparable in v1 (BAM/bgzf rests on a C zlib outside the contract)",
                    o.path
                ),
                "  VERDICT     : INCONCLUSIVE".to_string(),
            ],
            false,
        ));
    }

    // 4. Content-locate inputs; bind outputs to fresh temp paths.
    let index = index_inputs(inputs_dir)?;
    let located = index.len();
    let locate = |h: &str| index.get(h).map(|p| p.display().to_string());

    let work = make_temp_dir()?;
    let mut out_temp: BTreeMap<String, PathBuf> = BTreeMap::new();
    for (i, o) in manifest.outputs.iter().enumerate() {
        let ext = Path::new(&o.path)
            .extension()
            .and_then(|e| e.to_str())
            .unwrap_or("out");
        out_temp.insert(o.blake3.clone(), work.join(format!("out_{i}.{ext}")));
    }
    let temp_output = |h: &str| {
        out_temp
            .get(h)
            .map(|p| p.display().to_string())
            .unwrap_or_else(|| work.join("out_unknown").display().to_string())
    };

    // 5. Reconstruct the argv; an unlocatable input is INCONCLUSIVE (not DIVERGED).
    let argv = match CommandCapture::argv_from_command(&command, &locate, &temp_output) {
        Ok(a) => a,
        Err(e) => {
            return Ok(report(
                "INCONCLUSIVE",
                7,
                vec![
                    format!("  cannot reproduce: {e}"),
                    "  (pass --inputs pointing at a directory holding the recorded files)"
                        .to_string(),
                    "  VERDICT     : INCONCLUSIVE".to_string(),
                ],
                false,
            ))
        }
    };

    // 6. Re-execute the same binary.
    let exe = std::env::current_exe().context("locating the rosalind binary to re-run")?;
    let child = std::process::Command::new(&exe)
        .args(&argv)
        .output()
        .context("re-running the recorded command")?;
    if !child.status.success() {
        return Ok(report(
            "INCONCLUSIVE",
            7,
            vec![
                format!(
                    "  the re-run did not complete cleanly (exit {:?})",
                    child.status.code()
                ),
                format!(
                    "  stderr: {}",
                    String::from_utf8_lossy(&child.stderr).trim()
                ),
                "  VERDICT     : INCONCLUSIVE".to_string(),
            ],
            false,
        ));
    }

    // 7. Hash produced outputs and classify.
    let mut produced: Vec<(String, String)> = Vec::new();
    let mut cmps: Vec<OutputCmp> = Vec::new();
    for (i, o) in manifest.outputs.iter().enumerate() {
        let temp = out_temp
            .get(&o.blake3)
            .expect("temp path per recorded output");
        let observed = blake3_file(temp).with_context(|| {
            format!(
                "the re-run did not produce the expected output {}",
                temp.display()
            )
        })?;
        produced.push((format!("output[{i}]"), observed.clone()));
        cmps.push(OutputCmp {
            role: format!("output[{i}]"),
            recorded: o.blake3.clone(),
            observed: observed.clone(),
            matched: observed == o.blake3,
        });
    }
    let outcome = classify_outputs(&manifest.outputs, &produced);

    // 8. Machine-local resource + code context (diagnostic; never fatal).
    let declared_budget_mb = manifest
        .get_recorded("memory_budget_mb")
        .and_then(|v| v.parse().ok());
    let peak_here = manifest
        .outputs
        .first()
        .and_then(|o| out_temp.get(&o.blake3))
        .and_then(|p| read_rerun_peak(p));
    let fit = match (peak_here, declared_budget_mb) {
        (Some(peak), Some(mb)) => Some(crate::core::MemoryBudget::from_mb(mb).admits(peak)),
        _ => None,
    };
    let resource_here = ResourceHere {
        peak_rss_bytes: peak_here,
        declared_budget_mb,
        fit,
    };

    // 9. Build the report.
    let recorded_sha = manifest.params.get("code_git_sha").map(String::as_str);
    let here_sha = reproducer_code.get("code_git_sha").map(String::as_str);
    let code_line = match (recorded_sha, here_sha) {
        (Some(r), Some(h)) if r == h => format!("  code        : git {} — matches", short(h)),
        (Some(r), Some(h)) => format!(
            "  code        : git {} — DIFFERS from recorded {} (byte-match still counts)",
            short(h),
            short(r)
        ),
        _ => "  code        : (build-identity unavailable)".to_string(),
    };

    let mut lines = vec![
        format!("  claim       : {} (re-derived)", short(&parent_claim)),
        code_line,
        format!("  inputs      : {located} file(s) indexed by content hash"),
    ];
    for c in &cmps {
        if c.matched {
            lines.push(format!(
                "  output      : {}  OK byte-identical (blake3 {})",
                c.role,
                short(&c.observed)
            ));
        } else {
            lines.push(format!(
                "  output      : {}  DIVERGED (recorded {} got {})",
                c.role,
                short(&c.recorded),
                short(&c.observed)
            ));
        }
    }
    if let (Some(peak), Some(mb)) = (
        resource_here.peak_rss_bytes,
        resource_here.declared_budget_mb,
    ) {
        lines.push(format!(
            "  resource    : peak {} MiB vs declared {} MiB (here: this machine)",
            peak / (1 << 20),
            mb
        ));
    }
    lines.push(format!("  VERDICT     : {}", outcome.verdict.label()));

    // Clean up the temp work dir (best-effort).
    std::fs::remove_dir_all(&work).ok();

    Ok(ReproReport {
        verdict_label: outcome.verdict.label().to_string(),
        exit_code: outcome.verdict.exit_code(),
        lines,
        parent_claim,
        parent_subcommand,
        outputs: cmps,
        resource_here,
        reproducer_code,
        compared: true,
    })
}

/// First 10 hex chars of a digest for compact display.
fn short(hex: &str) -> &str {
    &hex[..hex.len().min(10)]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::provenance::FileHash;

    fn fh(hash: &str) -> FileHash {
        FileHash {
            path: "x".into(),
            blake3: hash.into(),
        }
    }

    #[test]
    fn classify_reproduced_when_all_outputs_match() {
        let recorded = vec![fh("h1")];
        let produced = vec![("o0".to_string(), "h1".to_string())];
        let v = classify_outputs(&recorded, &produced);
        assert!(matches!(v.verdict, Verdict::Reproduced));
        assert!(v.diffs.is_empty());
    }

    #[test]
    fn classify_diverged_names_the_first_mismatch() {
        let recorded = vec![fh("h1")];
        let produced = vec![("o0".to_string(), "DIFFERENT".to_string())];
        let v = classify_outputs(&recorded, &produced);
        assert!(matches!(v.verdict, Verdict::Diverged));
        assert!(v
            .diffs
            .iter()
            .any(|d| d.contains("h1") && d.contains("DIFFERENT")));
    }

    #[test]
    fn exit_codes_are_stable() {
        assert_eq!(Verdict::Reproduced.exit_code(), 0);
        assert_eq!(Verdict::Diverged.exit_code(), 6);
        assert_eq!(Verdict::Inconclusive.exit_code(), 7);
    }
}
