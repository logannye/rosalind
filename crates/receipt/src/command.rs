//! `CommandCapture` — the single chokepoint that records a normalized, replayable
//! invocation into a receipt's claim, and reconstructs an argv from it. Recording and
//! replay share this one structure so they cannot drift (the "forgot to record flag X"
//! bug class is eliminated by construction).
//!
//! Input/output operands carry a content hash, not a path, so the recorded `command`
//! string is machine-independent and is protected by the existing claim self-hash.

use std::path::Path;

use super::{blake3_file, FileHash, RunManifest};

/// One token of a recorded invocation.
#[derive(Debug)]
enum Token {
    Flag(String),
    Opt(String, String),
    Input { flag: String, blake3: String },
    Output { flag: String, blake3: String },
}

/// Accumulates an invocation, then writes it into a `RunManifest` (claim) and/or
/// reconstructs an argv for re-execution.
#[derive(Debug)]
pub struct CommandCapture {
    subcommand: String,
    tokens: Vec<Token>,
    inputs: Vec<FileHash>,
    outputs: Vec<FileHash>,
}

impl CommandCapture {
    /// Start capturing an invocation of `subcommand` (e.g. `"variants"`).
    pub fn new(subcommand: impl Into<String>) -> Self {
        Self {
            subcommand: subcommand.into(),
            tokens: Vec::new(),
            inputs: Vec::new(),
            outputs: Vec::new(),
        }
    }

    /// A content-addressed input operand; hashes the file and records it.
    pub fn input(&mut self, flag: &str, path: &Path) -> std::io::Result<&mut Self> {
        let h = blake3_file(path)?;
        Ok(self.input_hashed(flag, &path.display().to_string(), &h))
    }

    /// A content-addressed output operand; hashes the file and records it.
    pub fn output(&mut self, flag: &str, path: &Path) -> std::io::Result<&mut Self> {
        let h = blake3_file(path)?;
        Ok(self.output_hashed(flag, &path.display().to_string(), &h))
    }

    /// Input operand with a precomputed hash (when the caller already hashed it, and in
    /// tests). `path` is recorded into `inputs[]` for humans / `verify` to re-hash.
    pub fn input_hashed(&mut self, flag: &str, path: &str, blake3: &str) -> &mut Self {
        self.inputs.push(FileHash {
            path: path.to_string(),
            blake3: blake3.to_string(),
        });
        self.tokens.push(Token::Input {
            flag: flag.to_string(),
            blake3: blake3.to_string(),
        });
        self
    }

    /// Output operand with a precomputed hash.
    pub fn output_hashed(&mut self, flag: &str, path: &str, blake3: &str) -> &mut Self {
        self.outputs.push(FileHash {
            path: path.to_string(),
            blake3: blake3.to_string(),
        });
        self.tokens.push(Token::Output {
            flag: flag.to_string(),
            blake3: blake3.to_string(),
        });
        self
    }

    /// An option with a value, e.g. `("--max-depth", 1000)`.
    pub fn opt(&mut self, flag: &str, value: impl ToString) -> &mut Self {
        self.tokens
            .push(Token::Opt(flag.to_string(), value.to_string()));
        self
    }

    /// A bare flag, e.g. `"--enforce"`.
    pub fn flag(&mut self, flag: &str) -> &mut Self {
        self.tokens.push(Token::Flag(flag.to_string()));
        self
    }

    /// Record `flag` only when `cond` (so a `false` boolean leaves no trace).
    pub fn flag_if(&mut self, cond: bool, flag: &str) -> &mut Self {
        if cond {
            self.flag(flag)
        } else {
            self
        }
    }

    /// Render the normalized, machine-independent command string. Canonical order:
    /// subcommand, then input operands (in order added), then options sorted by flag,
    /// then bare flags sorted, then output operands — so it is stable run-to-run.
    fn render_command(&self) -> String {
        let mut parts: Vec<String> = vec![self.subcommand.clone()];
        for t in &self.tokens {
            if let Token::Input { flag, blake3 } = t {
                parts.push(flag.clone());
                parts.push(format!("@in:{blake3}"));
            }
        }
        let mut opts: Vec<(&String, &String)> = self
            .tokens
            .iter()
            .filter_map(|t| match t {
                Token::Opt(f, v) => Some((f, v)),
                _ => None,
            })
            .collect();
        opts.sort_by(|a, b| a.0.cmp(b.0));
        for (f, v) in opts {
            parts.push(f.clone());
            parts.push(v.clone());
        }
        let mut flags: Vec<&String> = self
            .tokens
            .iter()
            .filter_map(|t| match t {
                Token::Flag(f) => Some(f),
                _ => None,
            })
            .collect();
        flags.sort();
        for f in flags {
            parts.push(f.clone());
        }
        for t in &self.tokens {
            if let Token::Output { flag, blake3 } = t {
                parts.push(flag.clone());
                parts.push(format!("@out:{blake3}"));
            }
        }
        parts.join(" ")
    }

    /// Write the capture into a manifest's claim: the `command` recipe, the derived
    /// `inputs[]`/`outputs[]`, the discrete params (mechanical flag→key projection), and
    /// the inferred `mode`. Call before `finalize()`. Measurement fields remain the
    /// caller's responsibility (recorded separately, relocated by `finalize`).
    pub fn record_into(self, m: &mut RunManifest) {
        m.params
            .insert("command".to_string(), self.render_command());
        for t in &self.tokens {
            match t {
                Token::Opt(f, v) => {
                    m.params.insert(flag_to_key(f), v.clone());
                }
                Token::Flag(f) => {
                    m.params.insert(flag_to_key(f), "true".to_string());
                }
                _ => {}
            }
        }
        let has = |want: &str| {
            self.tokens
                .iter()
                .any(|t| matches!(t, Token::Input { flag, .. } if flag == want))
        };
        if has("--index") {
            m.params.insert("mode".to_string(), "index".to_string());
        } else if has("--reference") {
            m.params.insert("mode".to_string(), "reference".to_string());
        }
        m.inputs = self.inputs;
        m.outputs = self.outputs;
    }

    /// Reconstruct an argv from a recorded `command`: substitute each `@in:<h>` with the
    /// located input path and each `@out:<h>` with a caller-supplied temp path. Errors
    /// (naming the hash) if an input cannot be located.
    pub fn argv_from_command(
        command: &str,
        locate_input: &dyn Fn(&str) -> Option<String>,
        temp_output: &dyn Fn(&str) -> String,
    ) -> Result<Vec<String>, String> {
        let mut argv = Vec::new();
        for tok in command.split(' ') {
            if let Some(h) = tok.strip_prefix("@in:") {
                match locate_input(h) {
                    Some(p) => argv.push(p),
                    None => return Err(format!("input not located by content hash @in:{h}")),
                }
            } else if let Some(h) = tok.strip_prefix("@out:") {
                argv.push(temp_output(h));
            } else {
                argv.push(tok.to_string());
            }
        }
        Ok(argv)
    }
}

/// Mechanical `--max-depth` → `max_depth` projection (strip leading dashes,
/// dashes → underscores). Deterministic; no per-flag config.
fn flag_to_key(flag: &str) -> String {
    flag.trim_start_matches('-').replace('-', "_")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::RunManifest;

    // Build a capture WITHOUT touching the filesystem by injecting hashes directly.
    fn sample() -> CommandCapture {
        let mut c = CommandCapture::new("variants");
        c.input_hashed("--index", "ref.idx", "h_idx");
        c.input_hashed("--alignments", "s.bam", "h_bam");
        c.opt("--mapq-threshold", 20u8);
        c.opt("--max-depth", 1000u32);
        c.flag_if(true, "--enforce");
        c.flag_if(false, "--gvcf");
        c.opt("--memory-budget-mb", 256u64);
        c.output_hashed("-o", "out.vcf", "h_out");
        c
    }

    #[test]
    fn record_into_writes_command_inputs_outputs_and_discrete_params() {
        let mut m = RunManifest::new("variants");
        sample().record_into(&mut m);

        assert_eq!(
            m.params.get("command").unwrap(),
            "variants --index @in:h_idx --alignments @in:h_bam \
             --mapq-threshold 20 --max-depth 1000 --memory-budget-mb 256 --enforce -o @out:h_out"
        );
        assert_eq!(
            m.inputs
                .iter()
                .map(|f| f.blake3.as_str())
                .collect::<Vec<_>>(),
            ["h_idx", "h_bam"]
        );
        assert_eq!(
            m.outputs
                .iter()
                .map(|f| f.blake3.as_str())
                .collect::<Vec<_>>(),
            ["h_out"]
        );
        assert_eq!(m.params.get("mapq_threshold").unwrap(), "20");
        assert_eq!(m.params.get("max_depth").unwrap(), "1000");
        assert_eq!(m.params.get("memory_budget_mb").unwrap(), "256");
        assert_eq!(m.params.get("enforce").unwrap(), "true");
        assert_eq!(m.params.get("mode").unwrap(), "index");
        assert!(!m.params.contains_key("gvcf"));
    }

    #[test]
    fn argv_roundtrips_from_the_recorded_command() {
        let mut m = RunManifest::new("variants");
        sample().record_into(&mut m);
        let command = m.params.get("command").unwrap();

        let locate = |h: &str| match h {
            "h_idx" => Some("/data/ref.idx".to_string()),
            "h_bam" => Some("/data/s.bam".to_string()),
            _ => None,
        };
        let out_temp = |_h: &str| "/tmp/out.vcf".to_string();
        let argv = CommandCapture::argv_from_command(command, &locate, &out_temp).unwrap();

        assert_eq!(
            argv,
            vec![
                "variants",
                "--index",
                "/data/ref.idx",
                "--alignments",
                "/data/s.bam",
                "--mapq-threshold",
                "20",
                "--max-depth",
                "1000",
                "--memory-budget-mb",
                "256",
                "--enforce",
                "-o",
                "/tmp/out.vcf",
            ]
        );
    }

    #[test]
    fn argv_errors_when_an_input_cannot_be_located() {
        let mut m = RunManifest::new("variants");
        sample().record_into(&mut m);
        let command = m.params.get("command").unwrap();
        let locate = |_h: &str| None;
        let out_temp = |_h: &str| "/tmp/out.vcf".to_string();
        let err = CommandCapture::argv_from_command(command, &locate, &out_temp).unwrap_err();
        assert!(
            err.contains("h_idx"),
            "error names the unresolved input hash: {err}"
        );
    }
}
