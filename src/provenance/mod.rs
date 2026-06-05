//! A minimal, deterministic reproducibility receipt for a run: tool version,
//! subcommand, BLAKE3 content hashes of inputs + outputs, and the parameters.
//! Serialized as canonical JSON (sorted keys, no timestamps) so two identical
//! runs produce a byte-identical manifest.
//!
//! The receipt is split into a deterministic **claim** (inputs, outputs, params,
//! subcommand, versions) and a machine-/run-dependent **measurement** block (peak
//! RSS, working set, verdict, …). The self-hash (`manifest_blake3`) covers the
//! claim only — so the machine-dependent measured *cost* no longer perturbs it —
//! while a second `measurement_blake3` keeps that cost locally tamper-evident.
//!
//! The claim hash is a cross-machine **content-address**: for schema-3 receipts the
//! claim hashes inputs/outputs by their sorted `blake3` digests (recorded paths are
//! dropped from the claim form, though the on-disk receipt keeps them for humans and
//! for `verify` to re-hash files), so the same data at different paths hashes
//! identically. Pre-3 receipts hashed paths into the claim; the version gate
//! reproduces their form so they still self-verify.

use std::collections::BTreeMap;
use std::io::{self, Read, Write};
use std::path::{Path, PathBuf};

mod command;
pub use command::CommandCapture;

mod repro;
pub use repro::{ReproOutput, ReproReceipt};

/// Current receipt/feature schema version. Bump on any breaking schema change.
/// v2: split into a deterministic *claim* and a machine-dependent *measurement* block;
/// the self-hash (`manifest_blake3`) covers the claim only. v3: the claim hashes
/// inputs/outputs by their sorted `blake3` digests (paths dropped), so it is a
/// cross-machine content-address — the same data at different paths hashes identically.
/// v4: the claim records build-identity (`code_git_sha`/`code_dirty`/`rustc_version`/
/// `target_triple`/`deps_lock_blake3`), committing to exactly which code, toolchain, and
/// dependencies produced the run. v5: the claim records a normalized, replayable
/// `command` recipe (via [`CommandCapture`]) plus the discrete output-affecting params and
/// `mode`, so `reproduce` can re-derive the exact invocation.
pub const MANIFEST_SCHEMA_VERSION: u32 = 5;

/// Keys whose values are machine-/run-dependent measurements, not part of the
/// deterministic claim. `finalize` relocates these out of `params` into the
/// `measurements` block so the measured cost is excluded from the claim hash. The
/// single audited source of truth for the claim/measurement split.
pub const MEASUREMENT_KEYS: &[&str] = &[
    "peak_rss_bytes",
    "max_working_set_bytes",
    "predicted_peak_rss_bytes",
    "baseline_rss_bytes",
    "rss_residual_bytes",
    "governor",
    "contract_verdict",
];

/// Failure parsing a canonical run manifest.
#[derive(Debug)]
pub struct ManifestError(pub String);

impl std::fmt::Display for ManifestError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "malformed manifest: {}", self.0)
    }
}

impl std::error::Error for ManifestError {}

/// A file referenced by a run, with its content hash.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FileHash {
    /// Path as recorded (normalized to a string by the caller).
    pub path: String,
    /// BLAKE3 hex digest of the file's contents.
    pub blake3: String,
}

/// How input/output file entries render in a canonical form. The on-disk receipt
/// keeps full `{path, blake3}`; the schema-3 claim drops the path and hashes only the
/// content digest, so the claim hash does not depend on where files live.
#[derive(Clone, Copy)]
enum FileRender {
    WithPath,
    ContentOnly,
}

/// A reproducibility receipt for a single run.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RunManifest {
    /// Rosalind version (`CARGO_PKG_VERSION`).
    pub tool_version: String,
    /// Subcommand that produced the run (e.g. `variants`, `somatic`).
    pub subcommand: String,
    /// Input files and their content hashes.
    pub inputs: Vec<FileHash>,
    /// Run parameters (sorted by key in the canonical form).
    pub params: BTreeMap<String, String>,
    /// Output files and their content hashes.
    pub outputs: Vec<FileHash>,
    /// Machine-/run-dependent measured cost (peak RSS, working set, verdict, …),
    /// excluded from the claim hash. Carries its own `measurement_blake3`.
    pub measurements: BTreeMap<String, String>,
}

impl RunManifest {
    /// A new manifest stamped with the current tool version.
    pub fn new(subcommand: impl Into<String>) -> Self {
        Self {
            tool_version: env!("CARGO_PKG_VERSION").to_string(),
            subcommand: subcommand.into(),
            inputs: Vec::new(),
            params: BTreeMap::new(),
            outputs: Vec::new(),
            measurements: BTreeMap::new(),
        }
    }

    /// Record a machine-/run-dependent measurement (excluded from the claim hash).
    pub fn record_measurement(&mut self, key: impl Into<String>, value: impl Into<String>) {
        self.measurements.insert(key.into(), value.into());
    }

    /// Serialize to canonical JSON: keys sorted, arrays sorted by path, no
    /// timestamps — so identical runs render identically. Files render as full
    /// `{path, blake3}`; includes the measurement block (when non-empty). This is the
    /// on-disk receipt — humans and `verify` read files from its recorded paths.
    pub fn to_canonical_json(&self) -> String {
        self.push_canonical(true, FileRender::WithPath)
    }

    /// The claim-only canonical JSON: never emits the measurement block. This is the
    /// portion the self-hash commits to. For schema-3 receipts, files render as their
    /// sorted `blake3` digests (paths dropped) so the claim hash is a cross-machine
    /// content-address.
    pub fn to_canonical_claim_json(&self) -> String {
        self.push_canonical(false, self.claim_file_render())
    }

    /// Schema >= 3 → content-only claim (paths dropped, cross-machine). Older receipts
    /// hashed paths into the claim; reproduce their form so they still self-verify.
    fn claim_file_render(&self) -> FileRender {
        match self
            .params
            .get("schema_version")
            .and_then(|v| v.parse::<u32>().ok())
        {
            Some(v) if v >= 3 => FileRender::ContentOnly,
            _ => FileRender::WithPath,
        }
    }

    /// Render the canonical JSON, optionally including the measurement block, with
    /// files in the requested form. Canonical key order is alphabetical, so
    /// `measurements` sits between `inputs` and `outputs`; it is emitted only when
    /// non-empty (a pre-v2 receipt and a measurement-free receipt are byte-identical).
    fn push_canonical(&self, include_measurements: bool, files: FileRender) -> String {
        let mut out = String::new();
        out.push('{');
        out.push_str("\"inputs\":");
        push_files(&mut out, &self.inputs, files);
        if include_measurements && !self.measurements.is_empty() {
            out.push_str(",\"measurements\":");
            push_string_map(&mut out, &self.measurements);
        }
        out.push_str(",\"outputs\":");
        push_files(&mut out, &self.outputs, files);
        out.push_str(",\"params\":");
        push_string_map(&mut out, &self.params);
        out.push_str(",\"subcommand\":\"");
        out.push_str(&json_escape(&self.subcommand));
        out.push_str("\",\"tool_version\":\"");
        out.push_str(&json_escape(&self.tool_version));
        out.push_str("\"}");
        out
    }

    /// Parse a manifest from its canonical JSON form (the exact shape
    /// `to_canonical_json` emits; all values are strings). A small hand-parser —
    /// no general JSON dependency. Round-trips with `to_canonical_json`.
    pub fn from_canonical_json(s: &str) -> Result<RunManifest, ManifestError> {
        let mut p = Parser {
            b: s.as_bytes(),
            i: 0,
        };
        p.parse_manifest()
    }

    /// BLAKE3 hex of the **claim** canonical JSON with the self-hash excluded — the
    /// content this manifest commits to. Excludes the machine-dependent measurement
    /// block (so the measured cost does not perturb it) and, for schema-3 receipts,
    /// recorded paths (so it is a cross-machine content-address). Deterministic;
    /// `verify` re-derives it.
    pub fn content_hash(&self) -> String {
        let mut m = self.clone();
        m.params.remove("manifest_blake3");
        blake3_hex(m.to_canonical_claim_json().as_bytes())
    }

    /// BLAKE3 hex of the measurement block with `measurement_blake3` excluded — a
    /// LOCAL attestation of the measured cost. Not cross-machine stable by design
    /// (it hashes machine-dependent numbers), so it lives inside the measurement
    /// block rather than the claim.
    pub fn measurement_hash(&self) -> String {
        let mut m = self.measurements.clone();
        m.remove("measurement_blake3");
        let mut s = String::new();
        push_string_map(&mut s, &m);
        blake3_hex(s.as_bytes())
    }

    /// Look up a recorded value by key, checking `measurements` then `params`. Lets
    /// `verify` read v2 receipts (measured fields in `measurements`) and pre-v2
    /// receipts (everything in `params`) uniformly.
    pub fn get_recorded(&self, key: &str) -> Option<&String> {
        self.measurements.get(key).or_else(|| self.params.get(key))
    }

    /// Whether the (hash-protected) claim records that a measurement block exists.
    /// `verify` uses this to detect a measurement block stripped after the run — a
    /// claim that says `has_measurements` paired with a receipt that has none.
    pub fn claims_measurements(&self) -> bool {
        self.params
            .get("has_measurements")
            .map(|v| v == "true")
            .unwrap_or(false)
    }

    /// Check the recorded `code_git_sha` against an expected commit (prefix match, so
    /// short SHAs work). Returns the problems found: a mismatch, a clean match from a
    /// DIRTY tree (not reproducible from a SHA alone), or an inability to check (no /
    /// `unknown` SHA). An empty vec means a clean, matching build.
    pub fn check_expected_code(&self, expected: &str) -> Vec<String> {
        // Reject a degenerate expected SHA up front: an empty / too-short / non-hex
        // value would match loosely (or vacuously) via `starts_with` and give false
        // confidence — a scripted `--expect-code "$MAYBE_EMPTY"` must fail, not pass.
        if expected.len() < 7 || !expected.chars().all(|c| c.is_ascii_hexdigit()) {
            return vec![format!(
                "invalid --expect-code {expected:?}: expected a hex commit SHA of at least 7 chars"
            )];
        }
        match self.params.get("code_git_sha").map(String::as_str) {
            None => vec![
                "cannot check --expect-code: the receipt records no code_git_sha (a pre-P0.3 receipt)"
                    .to_string(),
            ],
            Some("unknown") => vec![
                "cannot check --expect-code: the receipt's code_git_sha is 'unknown' (a non-git build)"
                    .to_string(),
            ],
            Some(sha) if !sha.starts_with(expected) => vec![format!(
                "code mismatch: receipt was built from {sha}, expected {expected}"
            )],
            Some(_) => {
                if self.params.get("code_dirty").map(String::as_str) == Some("true") {
                    vec![format!(
                        "code matches {expected} but the receipt was built from a DIRTY tree \
                         (uncommitted changes) — not reproducible from a commit SHA alone"
                    )]
                } else {
                    Vec::new()
                }
            }
        }
    }

    /// Seal the receipt: partition measured fields out of the claim, stamp the
    /// measurement hash, the schema version, then the claim self-hash LAST (so it
    /// covers every other claim field, including the version). Idempotent.
    pub fn finalize(&mut self) {
        // 1. Partition: relocate machine-dependent fields OUT of the claim so the
        //    measured cost no longer perturbs the claim hash.
        for key in MEASUREMENT_KEYS {
            if let Some(v) = self.params.remove(*key) {
                self.measurements.insert((*key).to_string(), v);
            }
        }
        // 2. Local measurement attestation (only when a measurement exists), plus a
        //    DETERMINISTIC marker in the claim recording that a block exists. The
        //    marker is covered by the claim self-hash and is identical for the same
        //    logical run on any machine (unlike `measurement_blake3`, which hashes
        //    machine-dependent numbers and so cannot live in the claim). It lets
        //    `verify` catch a measurement block stripped to evade the cost checks:
        //    dropping the block leaves the claim asserting one must exist.
        if !self.measurements.is_empty() {
            let mh = self.measurement_hash();
            self.measurements
                .insert("measurement_blake3".to_string(), mh);
            self.params
                .insert("has_measurements".to_string(), "true".to_string());
        }
        // 3. Build-identity (baked at compile time by build.rs) — part of the claim,
        //    so it is committed to by the self-hash and forms the reproduction key:
        //    exactly which code, toolchain, and deps produced this run.
        for (k, v) in [
            ("code_git_sha", env!("ROSALIND_GIT_SHA")),
            ("code_dirty", env!("ROSALIND_GIT_DIRTY")),
            ("rustc_version", env!("ROSALIND_RUSTC_VERSION")),
            ("target_triple", env!("ROSALIND_TARGET")),
            ("deps_lock_blake3", env!("ROSALIND_DEPS_LOCK_BLAKE3")),
        ] {
            self.params.insert(k.to_string(), v.to_string());
        }
        // 4. Stamp the version into the claim, then the claim self-hash last.
        self.params.insert(
            "schema_version".to_string(),
            MANIFEST_SCHEMA_VERSION.to_string(),
        );
        let h = self.content_hash();
        self.params.insert("manifest_blake3".to_string(), h);
    }

    /// `Some(true)`/`Some(false)` if a claim self-hash is recorded and matches /
    /// mismatches; `None` if none is recorded (a pre-1.2 receipt).
    pub fn self_hash_ok(&self) -> Option<bool> {
        self.params
            .get("manifest_blake3")
            .map(|recorded| *recorded == self.content_hash())
    }

    /// `Some(true)`/`Some(false)` if a measurement self-hash is recorded and matches /
    /// mismatches; `None` if none is recorded (no measurement block, or a pre-v2 receipt).
    pub fn measurement_hash_ok(&self) -> Option<bool> {
        self.measurements
            .get("measurement_blake3")
            .map(|recorded| *recorded == self.measurement_hash())
    }
}

/// Minimal recursive parser for the fixed canonical-manifest shape. Every value
/// is a JSON string (inputs/outputs are arrays of `{blake3, path}` objects;
/// params is an object of string→string), so the parser only needs strings,
/// arrays, and objects — no numbers/bools/null.
struct Parser<'a> {
    b: &'a [u8],
    i: usize,
}

impl Parser<'_> {
    fn err(&self, m: &str) -> ManifestError {
        ManifestError(format!("{m} at byte {}", self.i))
    }

    fn expect(&mut self, c: u8) -> Result<(), ManifestError> {
        if self.i < self.b.len() && self.b[self.i] == c {
            self.i += 1;
            Ok(())
        } else {
            Err(self.err(&format!("expected '{}'", c as char)))
        }
    }

    fn parse_string(&mut self) -> Result<String, ManifestError> {
        self.expect(b'"')?;
        let mut buf: Vec<u8> = Vec::new();
        while self.i < self.b.len() {
            let c = self.b[self.i];
            self.i += 1;
            match c {
                b'"' => {
                    return String::from_utf8(buf).map_err(|_| self.err("invalid utf-8"));
                }
                b'\\' => {
                    let e = *self
                        .b
                        .get(self.i)
                        .ok_or_else(|| self.err("trailing escape"))?;
                    self.i += 1;
                    match e {
                        b'"' => buf.push(b'"'),
                        b'\\' => buf.push(b'\\'),
                        b'n' => buf.push(b'\n'),
                        b'r' => buf.push(b'\r'),
                        b't' => buf.push(b'\t'),
                        b'u' => {
                            let hex = self
                                .b
                                .get(self.i..self.i + 4)
                                .ok_or_else(|| self.err("short \\u"))?;
                            let cp = u32::from_str_radix(
                                std::str::from_utf8(hex).map_err(|_| self.err("bad \\u"))?,
                                16,
                            )
                            .map_err(|_| self.err("bad \\u"))?;
                            let ch = char::from_u32(cp).ok_or_else(|| self.err("bad codepoint"))?;
                            let mut tmp = [0u8; 4];
                            buf.extend_from_slice(ch.encode_utf8(&mut tmp).as_bytes());
                            self.i += 4;
                        }
                        _ => return Err(self.err("bad escape")),
                    }
                }
                _ => buf.push(c),
            }
        }
        Err(self.err("unterminated string"))
    }

    fn expect_key(&mut self, key: &str) -> Result<(), ManifestError> {
        let k = self.parse_string()?;
        if k != key {
            return Err(self.err(&format!("expected key \"{key}\", got \"{k}\"")));
        }
        self.expect(b':')
    }

    fn parse_file_array(&mut self) -> Result<Vec<FileHash>, ManifestError> {
        self.expect(b'[')?;
        let mut out = Vec::new();
        if self.i < self.b.len() && self.b[self.i] == b']' {
            self.i += 1;
            return Ok(out);
        }
        loop {
            self.expect(b'{')?;
            self.expect_key("blake3")?;
            let blake3 = self.parse_string()?;
            self.expect(b',')?;
            self.expect_key("path")?;
            let path = self.parse_string()?;
            self.expect(b'}')?;
            out.push(FileHash { path, blake3 });
            match self.b.get(self.i) {
                Some(b',') => self.i += 1,
                Some(b']') => {
                    self.i += 1;
                    break;
                }
                _ => return Err(self.err("expected ',' or ']' in array")),
            }
        }
        Ok(out)
    }

    fn parse_params(&mut self) -> Result<BTreeMap<String, String>, ManifestError> {
        self.expect(b'{')?;
        let mut map = BTreeMap::new();
        if self.i < self.b.len() && self.b[self.i] == b'}' {
            self.i += 1;
            return Ok(map);
        }
        loop {
            let k = self.parse_string()?;
            self.expect(b':')?;
            let v = self.parse_string()?;
            // Reject duplicate keys: a crafted receipt with two values for one key
            // could otherwise show one value to a human reader while `verify` (and
            // `get_recorded`) act on the other.
            if map.insert(k.clone(), v).is_some() {
                return Err(self.err(&format!("duplicate key \"{k}\"")));
            }
            match self.b.get(self.i) {
                Some(b',') => self.i += 1,
                Some(b'}') => {
                    self.i += 1;
                    break;
                }
                _ => return Err(self.err("expected ',' or '}' in object")),
            }
        }
        Ok(map)
    }

    fn parse_manifest(&mut self) -> Result<RunManifest, ManifestError> {
        self.expect(b'{')?;
        self.expect_key("inputs")?;
        let inputs = self.parse_file_array()?;
        self.expect(b',')?;
        // `measurements` is optional: absent in pre-v2 receipts and measurement-free
        // runs. Read the next key and branch on whether it is the measurement block.
        let key = self.parse_string()?;
        self.expect(b':')?;
        let (measurements, outputs) = if key == "measurements" {
            let m = self.parse_params()?;
            self.expect(b',')?;
            self.expect_key("outputs")?;
            (m, self.parse_file_array()?)
        } else if key == "outputs" {
            (BTreeMap::new(), self.parse_file_array()?)
        } else {
            return Err(self.err(&format!(
                "expected \"measurements\" or \"outputs\", got \"{key}\""
            )));
        };
        self.expect(b',')?;
        self.expect_key("params")?;
        let params = self.parse_params()?;
        self.expect(b',')?;
        self.expect_key("subcommand")?;
        let subcommand = self.parse_string()?;
        self.expect(b',')?;
        self.expect_key("tool_version")?;
        let tool_version = self.parse_string()?;
        self.expect(b'}')?;
        Ok(RunManifest {
            tool_version,
            subcommand,
            inputs,
            params,
            outputs,
            measurements,
        })
    }
}

/// Render a file array in the requested form: `[{"blake3","path"}]` sorted by path
/// (on-disk), or `["<blake3>",…]` sorted by blake3 (the content-only claim).
fn push_files(out: &mut String, files: &[FileHash], render: FileRender) {
    match render {
        FileRender::WithPath => push_file_hashes(out, files),
        FileRender::ContentOnly => push_blake3_list(out, files),
    }
}

/// Render `["<blake3>",…]`, sorted by digest (a content multiset — duplicates kept).
fn push_blake3_list(out: &mut String, files: &[FileHash]) {
    let mut digests: Vec<&str> = files.iter().map(|f| f.blake3.as_str()).collect();
    digests.sort_unstable();
    out.push('[');
    for (i, d) in digests.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        out.push('"');
        out.push_str(&json_escape(d));
        out.push('"');
    }
    out.push(']');
}

/// Render a `[{"blake3":..,"path":..}, ..]` array, entries sorted by path.
fn push_file_hashes(out: &mut String, files: &[FileHash]) {
    let mut sorted: Vec<&FileHash> = files.iter().collect();
    sorted.sort_by(|a, b| a.path.cmp(&b.path));
    out.push('[');
    for (i, f) in sorted.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        out.push_str("{\"blake3\":\"");
        out.push_str(&json_escape(&f.blake3));
        out.push_str("\",\"path\":\"");
        out.push_str(&json_escape(&f.path));
        out.push_str("\"}");
    }
    out.push(']');
}

/// Render a `{"k":"v",…}` object, entries in the map's (sorted) key order. Shared by
/// the `params` and `measurements` blocks (both are `string → string`).
fn push_string_map(out: &mut String, map: &BTreeMap<String, String>) {
    out.push('{');
    for (i, (k, v)) in map.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        out.push('"');
        out.push_str(&json_escape(k));
        out.push_str("\":\"");
        out.push_str(&json_escape(v));
        out.push('"');
    }
    out.push('}');
}

/// Minimal RFC-8259 string escaping for the characters we can encounter.
fn json_escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
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

/// BLAKE3 hex digest of a byte slice.
pub fn blake3_hex(bytes: &[u8]) -> String {
    blake3::hash(bytes).to_hex().to_string()
}

/// BLAKE3 hex digest of a file's contents, streamed in fixed-size chunks
/// (bounded memory regardless of file size).
pub fn blake3_file(path: &Path) -> io::Result<String> {
    let mut hasher = blake3::Hasher::new();
    let mut file = std::fs::File::open(path)?;
    let mut buf = [0u8; 64 * 1024];
    loop {
        let n = file.read(&mut buf)?;
        if n == 0 {
            break;
        }
        hasher.update(&buf[..n]);
    }
    Ok(hasher.finalize().to_hex().to_string())
}

/// Options for [`verify_receipt`] beyond the receipt text itself.
#[derive(Debug, Default)]
pub struct VerifyOpts {
    /// Budget (MiB) to check the recorded peak against (overrides the recorded one).
    pub budget_mb: Option<u64>,
    /// Assert the receipt was built from exactly this commit SHA (prefix match).
    pub expect_code: Option<String>,
    /// Re-hash the files at the recorded input/output paths and check their digests.
    pub rehash_files: bool,
}

/// The outcome of [`verify_receipt`]: failures (empty == ok), informational notes, and
/// the parsed manifest (when parsing succeeded).
#[derive(Debug)]
pub struct VerifyReport {
    /// Whether the receipt passed every check.
    pub ok: bool,
    /// Human-readable failures (empty when `ok`).
    pub problems: Vec<String>,
    /// Informational, non-failing notes (e.g. "peak within budget").
    pub notes: Vec<String>,
    /// The parsed manifest, when parsing succeeded.
    pub manifest: Option<RunManifest>,
}

/// Check a receipt's internal integrity — self-hashes, cross-field consistency, optional
/// budget + expected-code — and, when `opts.rehash_files`, re-hash recorded files. The
/// single source of truth shared by the `verify` CLI and `reproduce` (and a future WASM
/// verifier) so they cannot drift.
pub fn verify_receipt(text: &str, opts: &VerifyOpts) -> VerifyReport {
    let manifest = match RunManifest::from_canonical_json(text) {
        Ok(m) => m,
        Err(e) => {
            return VerifyReport {
                ok: false,
                problems: vec![format!("parse error: {e}")],
                notes: Vec::new(),
                manifest: None,
            }
        }
    };
    let mut problems: Vec<String> = Vec::new();
    let mut notes: Vec<String> = Vec::new();

    // Re-hash inputs + outputs against the recorded digests (CLI sets this).
    if opts.rehash_files {
        for (kind, files) in [("input", &manifest.inputs), ("output", &manifest.outputs)] {
            for f in files {
                match blake3_file(Path::new(&f.path)) {
                    Ok(h) if h == f.blake3 => {}
                    Ok(h) => problems.push(format!(
                        "{kind} {} hash mismatch: recorded {}, now {}",
                        f.path, f.blake3, h
                    )),
                    Err(e) => problems.push(format!("{kind} {} unreadable: {e}", f.path)),
                }
            }
        }
    }

    // Strictly parse the numeric fields. A PRESENT-but-unparseable field is corruption,
    // not absence. `get_recorded` reads measurements (v2) or params (pre-v2) uniformly.
    let parse_num = |k: &str, problems: &mut Vec<String>| -> Option<u64> {
        match manifest.get_recorded(k) {
            None => None,
            Some(v) => match v.parse::<u64>() {
                Ok(n) => Some(n),
                Err(_) => {
                    problems.push(format!("malformed numeric field {k}: {v:?}"));
                    None
                }
            },
        }
    };
    let recorded_peak = parse_num("peak_rss_bytes", &mut problems);
    let recorded_ws = parse_num("max_working_set_bytes", &mut problems);
    let recorded_budget = parse_num("memory_budget_mb", &mut problems);

    // Re-check the recorded realized peak against the budget (CLI overrides manifest).
    let budget_mb = opts.budget_mb.or(recorded_budget);
    match (budget_mb, recorded_peak) {
        (Some(mb), Some(peak)) => {
            if crate::core::MemoryBudget::from_mb(mb).admits(peak) {
                notes.push(format!(
                    "peak {} MiB within budget {mb} MiB",
                    peak / (1 << 20)
                ));
            } else {
                problems.push(format!(
                    "recorded peak {} MiB exceeded budget {mb} MiB",
                    peak / (1 << 20)
                ));
            }
        }
        (None, _) => notes.push("no budget to check (none supplied or recorded)".to_string()),
        (Some(_), None) => problems.push("manifest has no recorded peak_rss_bytes".to_string()),
    }

    // Internal-consistency cross-checks: the working set cannot exceed peak RSS.
    if let (Some(ws), Some(peak)) = (recorded_ws, recorded_peak) {
        if ws > peak {
            problems.push(format!(
                "internally inconsistent: max_working_set_bytes ({ws}) exceeds peak_rss_bytes ({peak})"
            ));
        }
    }
    // A recorded verdict must agree with the recorded peak vs the recorded budget.
    if let (Some(verdict), Some(mb), Some(peak)) = (
        manifest
            .get_recorded("contract_verdict")
            .map(String::as_str),
        recorded_budget,
        recorded_peak,
    ) {
        let actually_within = crate::core::MemoryBudget::from_mb(mb).admits(peak);
        if verdict == "within" && !actually_within {
            problems.push(format!(
                "internally inconsistent: contract_verdict='within' but recorded peak {} MiB \
                 exceeds recorded budget {mb} MiB",
                peak / (1 << 20)
            ));
        }
        if verdict == "over" && actually_within {
            problems.push(format!(
                "internally inconsistent: contract_verdict='over' but recorded peak {} MiB \
                 is within recorded budget {mb} MiB",
                peak / (1 << 20)
            ));
        }
    }

    // Claim self-hash: catches any post-write edit to the claim. A pre-1.2 receipt has none.
    match manifest.self_hash_ok() {
        Some(true) => {}
        Some(false) => problems.push(
            "manifest_blake3 mismatch: the receipt was modified after it was written".to_string(),
        ),
        None => {
            notes.push("no manifest_blake3 (a pre-1.2 receipt); skipping self-hash".to_string())
        }
    }

    // Measurement self-hash: catches an edit to a measured field the claim hash can't see.
    match manifest.measurement_hash_ok() {
        Some(true) => {}
        Some(false) => problems.push(
            "measurement_blake3 mismatch: a measured field was modified after the run".to_string(),
        ),
        None => {
            if manifest.claims_measurements() {
                problems.push(
                    "measurement block missing: the claim records a measurement block but \
                     the receipt has none (stripped after the run?)"
                        .to_string(),
                );
            }
        }
    }

    // Build-identity: assert the receipt came from exactly the expected commit.
    if let Some(expected) = &opts.expect_code {
        let code_problems = manifest.check_expected_code(expected);
        if code_problems.is_empty() {
            notes.push(format!("code matches {expected} (clean build)"));
        } else {
            problems.extend(code_problems);
        }
    }

    VerifyReport {
        ok: problems.is_empty(),
        notes,
        manifest: Some(manifest),
        problems,
    }
}

/// Write `<output_path>.manifest.json` next to the output, returning its path.
pub fn write_manifest(output_path: &Path, manifest: &RunManifest) -> io::Result<PathBuf> {
    let mut name = output_path.as_os_str().to_os_string();
    name.push(".manifest.json");
    let manifest_path = PathBuf::from(name);
    let mut file = std::fs::File::create(&manifest_path)?;
    file.write_all(manifest.to_canonical_json().as_bytes())?;
    file.flush()?;
    Ok(manifest_path)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn verify_receipt_reports_a_tampered_claim() {
        let mut m = RunManifest::new("variants");
        m.inputs.push(FileHash {
            path: "a".into(),
            blake3: "aa".into(),
        });
        m.finalize();
        // Flip a byte of a recorded input digest without re-sealing — the claim
        // self-hash must catch it.
        let text = m.to_canonical_json().replace("\"aa\"", "\"ab\"");
        let report = verify_receipt(&text, &VerifyOpts::default());
        assert!(!report.ok, "tampered claim must not verify");
        assert!(
            report
                .problems
                .iter()
                .any(|p| p.contains("manifest_blake3")),
            "{:?}",
            report.problems
        );
    }

    #[test]
    fn verify_receipt_passes_a_clean_receipt() {
        let mut m = RunManifest::new("variants");
        m.finalize();
        let report = verify_receipt(&m.to_canonical_json(), &VerifyOpts::default());
        assert!(report.ok, "{:?}", report.problems);
    }

    #[test]
    fn blake3_is_deterministic_and_sensitive() {
        assert_eq!(blake3_hex(b"abc"), blake3_hex(b"abc"));
        assert_ne!(blake3_hex(b"abc"), blake3_hex(b"abd"));
        // Known BLAKE3 vector for the empty input.
        assert_eq!(
            blake3_hex(b""),
            "af1349b9f5f9a1a6a0404dea36dcc9499bcb25c9adc112b7cc9a93cae41f3262"
        );
        assert_eq!(blake3_hex(b"abc").len(), 64);
    }

    #[test]
    fn canonical_json_has_sorted_keys_and_is_exact() {
        let mut m = RunManifest::new("variants");
        m.tool_version = "0.1.0".to_string();
        m.inputs.push(FileHash {
            path: "ref.fa".to_string(),
            blake3: "aa".to_string(),
        });
        m.outputs.push(FileHash {
            path: "out.vcf".to_string(),
            blake3: "bb".to_string(),
        });
        m.params.insert("min_qual".to_string(), "30".to_string());
        m.params.insert("min_depth".to_string(), "8".to_string());

        let json = m.to_canonical_json();
        assert_eq!(
            json,
            r#"{"inputs":[{"blake3":"aa","path":"ref.fa"}],"outputs":[{"blake3":"bb","path":"out.vcf"}],"params":{"min_depth":"8","min_qual":"30"},"subcommand":"variants","tool_version":"0.1.0"}"#
        );
    }

    #[test]
    fn canonical_json_is_order_independent() {
        let mk = |order_swapped: bool| {
            let mut m = RunManifest::new("variants");
            m.tool_version = "0.1.0".to_string();
            let a = FileHash {
                path: "a.fa".to_string(),
                blake3: "1".to_string(),
            };
            let b = FileHash {
                path: "b.fa".to_string(),
                blake3: "2".to_string(),
            };
            if order_swapped {
                m.inputs.push(b);
                m.inputs.push(a);
            } else {
                m.inputs.push(a);
                m.inputs.push(b);
            }
            m.to_canonical_json()
        };
        // Inputs are sorted by path in the canonical form → order-independent.
        assert_eq!(mk(false), mk(true));
    }

    #[test]
    fn json_escapes_special_characters() {
        let mut m = RunManifest::new("variants");
        m.tool_version = "0.1.0".to_string();
        m.params.insert("note".to_string(), "a\"b\\c".to_string());
        let json = m.to_canonical_json();
        assert!(json.contains(r#""note":"a\"b\\c""#));
    }

    #[test]
    fn write_manifest_emits_sidecar_file() {
        let dir =
            std::env::temp_dir().join(format!("rosalind_manifest_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let out = dir.join("calls.vcf");
        std::fs::write(&out, b"##fileformat=VCFv4.2\n").unwrap();

        let mut m = RunManifest::new("variants");
        m.tool_version = "0.1.0".to_string();
        m.outputs.push(FileHash {
            path: out.display().to_string(),
            blake3: blake3_file(&out).unwrap(),
        });

        let manifest_path = write_manifest(&out, &m).unwrap();
        assert_eq!(manifest_path, dir.join("calls.vcf.manifest.json"));
        let written = std::fs::read_to_string(&manifest_path).unwrap();
        assert_eq!(written, m.to_canonical_json());

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn parse_round_trips_canonical_json_including_escapes() {
        let mut m = RunManifest::new("variants");
        m.tool_version = "9.9.9".to_string();
        m.inputs.push(FileHash {
            path: "weird \"path\"\twith\\escapes/和.fa".to_string(),
            blake3: "aa".to_string(),
        });
        m.inputs.push(FileHash {
            path: "a.idx".to_string(),
            blake3: "bb".to_string(),
        });
        m.outputs.push(FileHash {
            path: "out.vcf".to_string(),
            blake3: "cc".to_string(),
        });
        m.params
            .insert("contract_verdict".to_string(), "within".to_string());
        m.params
            .insert("peak_rss_bytes".to_string(), "12345".to_string());
        m.params
            .insert("note".to_string(), "line1\nline2".to_string());

        let json = m.to_canonical_json();
        let parsed = RunManifest::from_canonical_json(&json).expect("parse");
        // serialize → parse → serialize is the identity on the canonical form.
        assert_eq!(parsed.to_canonical_json(), json);
        assert_eq!(parsed.tool_version, "9.9.9");
        assert_eq!(parsed.subcommand, "variants");
        assert_eq!(parsed.params.get("note").unwrap(), "line1\nline2");
        assert_eq!(parsed.params.get("contract_verdict").unwrap(), "within");
    }

    #[test]
    fn parse_rejects_malformed() {
        assert!(RunManifest::from_canonical_json("not json").is_err());
        assert!(RunManifest::from_canonical_json("{\"inputs\":[}").is_err());
    }

    #[test]
    fn finalize_stamps_schema_version_and_a_matching_self_hash() {
        let mut m = RunManifest::new("variants");
        m.tool_version = "0.1.0".to_string();
        m.params
            .insert("peak_rss_bytes".to_string(), "123".to_string());
        m.finalize();
        assert_eq!(
            m.params.get("schema_version").map(String::as_str),
            Some(MANIFEST_SCHEMA_VERSION.to_string().as_str())
        );
        assert!(m.params.contains_key("manifest_blake3"));
        assert_eq!(m.self_hash_ok(), Some(true), "fresh finalize must verify");
    }

    #[test]
    fn tampering_any_field_breaks_the_self_hash() {
        let mut m = RunManifest::new("variants");
        m.tool_version = "0.1.0".to_string();
        m.params
            .insert("peak_rss_bytes".to_string(), "123".to_string());
        m.finalize();
        // Flip a field WITHOUT re-finalizing — the recorded hash no longer matches.
        m.params
            .insert("peak_rss_bytes".to_string(), "999".to_string());
        assert_eq!(m.self_hash_ok(), Some(false));
    }

    #[test]
    fn a_manifest_without_a_self_hash_returns_none() {
        let mut m = RunManifest::new("variants");
        m.tool_version = "0.1.0".to_string();
        assert_eq!(m.self_hash_ok(), None);
    }

    #[test]
    fn finalize_is_idempotent() {
        let mut m = RunManifest::new("variants");
        m.tool_version = "0.1.0".to_string();
        m.params
            .insert("peak_rss_bytes".to_string(), "123".to_string());
        m.finalize();
        let first = m.params.get("manifest_blake3").cloned();
        m.finalize();
        assert_eq!(m.params.get("manifest_blake3").cloned(), first);
        assert_eq!(m.self_hash_ok(), Some(true));
    }

    #[test]
    fn claim_hash_is_stable_across_machine_dependent_measurements() {
        // Same logical run on two machines: identical claim, different measured cost.
        // The claim self-hash must match; the measurement is excluded from it.
        let mk = |peak: &str, ws: &str| {
            let mut m = RunManifest::new("variants");
            m.tool_version = "0.1.0".to_string();
            m.inputs.push(FileHash {
                path: "ref.fa".to_string(),
                blake3: "aa".to_string(),
            });
            m.outputs.push(FileHash {
                path: "out.vcf".to_string(),
                blake3: "bb".to_string(),
            });
            m.params.insert("min_qual".to_string(), "30".to_string());
            m.record_measurement("peak_rss_bytes", peak);
            m.record_measurement("max_working_set_bytes", ws);
            m.finalize();
            m
        };
        let a = mk("1000000", "4096");
        let b = mk("9999999", "8192");
        assert_eq!(
            a.content_hash(),
            b.content_hash(),
            "measured cost must not change the claim hash"
        );
        assert_eq!(a.self_hash_ok(), Some(true));
        assert_eq!(b.self_hash_ok(), Some(true));
        // Differing measurements DO change the measurement hash.
        assert_ne!(
            a.measurements.get("measurement_blake3"),
            b.measurements.get("measurement_blake3")
        );
    }

    #[test]
    fn claim_excludes_but_full_form_includes_measurements() {
        let mut m = RunManifest::new("variants");
        m.tool_version = "0.1.0".to_string();
        m.record_measurement("peak_rss_bytes", "123");
        m.finalize();
        assert!(
            !m.to_canonical_claim_json().contains("peak_rss_bytes"),
            "claim form must not carry the measurement"
        );
        assert!(
            m.to_canonical_json().contains("peak_rss_bytes"),
            "full form must record the measurement"
        );
    }

    #[test]
    fn editing_a_measurement_breaks_only_the_measurement_hash() {
        let mut m = RunManifest::new("variants");
        m.tool_version = "0.1.0".to_string();
        m.record_measurement("peak_rss_bytes", "123");
        m.finalize();
        assert_eq!(m.self_hash_ok(), Some(true));
        assert_eq!(m.measurement_hash_ok(), Some(true));
        // Lower the recorded peak WITHOUT re-finalizing (a tampered receipt).
        m.measurements
            .insert("peak_rss_bytes".to_string(), "1".to_string());
        assert_eq!(
            m.self_hash_ok(),
            Some(true),
            "claim hash is unaffected by the measurement edit"
        );
        assert_eq!(
            m.measurement_hash_ok(),
            Some(false),
            "measurement hash must catch the edit"
        );
    }

    #[test]
    fn finalize_relocates_measured_keys_out_of_the_claim() {
        let mut m = RunManifest::new("variants");
        m.tool_version = "0.1.0".to_string();
        // Insert measured fields the legacy way (into params); finalize must relocate.
        m.params
            .insert("peak_rss_bytes".to_string(), "555".to_string());
        m.params
            .insert("contract_verdict".to_string(), "within".to_string());
        m.params.insert("min_qual".to_string(), "30".to_string());
        m.finalize();
        for k in ["peak_rss_bytes", "contract_verdict"] {
            assert!(!m.params.contains_key(k), "{k} must leave the claim");
            assert!(
                m.measurements.contains_key(k),
                "{k} must enter measurements"
            );
        }
        assert!(m.params.contains_key("min_qual"), "claim params stay put");
    }

    #[test]
    fn pre_v2_receipt_with_measurements_in_params_still_verifies() {
        // A pre-v2 receipt: peak in params, no measurements key, schema 1, self-hash
        // over the params-inclusive claim. It must still self-verify (graceful degrade).
        let mut m = RunManifest::new("variants");
        m.tool_version = "0.1.0".to_string();
        m.params
            .insert("peak_rss_bytes".to_string(), "123".to_string());
        m.params
            .insert("schema_version".to_string(), "1".to_string());
        let h = m.content_hash();
        m.params.insert("manifest_blake3".to_string(), h);

        assert_eq!(m.self_hash_ok(), Some(true));
        assert_eq!(m.measurement_hash_ok(), None, "no measurement block in v1");
        let json = m.to_canonical_json();
        assert!(
            !json.contains("\"measurements\""),
            "v1 emits no measurements key"
        );
        let parsed = RunManifest::from_canonical_json(&json).expect("parse v1");
        assert!(parsed.measurements.is_empty());
        assert_eq!(parsed.self_hash_ok(), Some(true));
    }

    #[test]
    fn v2_receipt_round_trips_through_the_parser() {
        let mut m = RunManifest::new("features");
        m.tool_version = "9.9.9".to_string();
        m.inputs.push(FileHash {
            path: "a.idx".to_string(),
            blake3: "aa".to_string(),
        });
        m.outputs.push(FileHash {
            path: "out.tsv".to_string(),
            blake3: "bb".to_string(),
        });
        m.params
            .insert("feature_rows".to_string(), "42".to_string());
        m.record_measurement("peak_rss_bytes", "1000");
        m.record_measurement("governor", "enforced");
        m.finalize();

        let json = m.to_canonical_json();
        let parsed = RunManifest::from_canonical_json(&json).expect("parse v2");
        assert_eq!(
            parsed.to_canonical_json(),
            json,
            "round-trip is the identity"
        );
        assert_eq!(
            parsed
                .measurements
                .get("peak_rss_bytes")
                .map(String::as_str),
            Some("1000")
        );
        assert_eq!(
            parsed.params.get("feature_rows").map(String::as_str),
            Some("42")
        );
        assert_eq!(parsed.self_hash_ok(), Some(true));
        assert_eq!(parsed.measurement_hash_ok(), Some(true));
    }

    #[test]
    fn finalize_records_a_claim_marker_only_when_a_measurement_exists() {
        // With a measurement: the claim records has_measurements (so stripping the
        // block is detectable), and the marker is covered by the claim self-hash.
        let mut with = RunManifest::new("variants");
        with.tool_version = "0.1.0".to_string();
        with.record_measurement("peak_rss_bytes", "123");
        with.finalize();
        assert!(with.claims_measurements());
        assert_eq!(
            with.params.get("has_measurements").map(String::as_str),
            Some("true")
        );
        assert_eq!(with.self_hash_ok(), Some(true));

        // Without a measurement (e.g. a somatic run): no marker, nothing to strip.
        let mut without = RunManifest::new("somatic");
        without.tool_version = "0.1.0".to_string();
        without.finalize();
        assert!(!without.claims_measurements());
        assert!(!without.params.contains_key("has_measurements"));
    }

    #[test]
    fn stripping_the_measurement_block_leaves_the_claim_asserting_one_exists() {
        // The deletion-bypass guard: clearing the measurement block keeps the claim
        // self-hash valid (the claim never carried the block), but the claim still
        // records has_measurements while measurement_hash_ok drops to None — the
        // exact signal `verify` keys on.
        let mut m = RunManifest::new("variants");
        m.tool_version = "0.1.0".to_string();
        m.record_measurement("peak_rss_bytes", "900000000");
        m.finalize();
        m.measurements.clear();
        assert_eq!(
            m.self_hash_ok(),
            Some(true),
            "claim is intact after stripping"
        );
        assert_eq!(m.measurement_hash_ok(), None, "no block to hash");
        assert!(
            m.claims_measurements(),
            "the claim still asserts a measurement block must exist"
        );
    }

    #[test]
    fn parser_rejects_a_duplicate_key() {
        // A crafted receipt cannot carry two values for one key (reader/verifier
        // shadowing). `params` and `measurements` both go through `parse_params`.
        let dup = r#"{"inputs":[],"outputs":[],"params":{"k":"1","k":"2"},"subcommand":"x","tool_version":"0.1.0"}"#;
        assert!(RunManifest::from_canonical_json(dup).is_err());
    }

    #[test]
    fn claim_hash_is_stable_across_machine_dependent_paths() {
        // Same content (blake3) at DIFFERENT paths on two machines → SAME claim hash.
        // The keystone of P0.2b: paths are dropped from the claim form.
        let mk = |idx_path: &str, out_path: &str| {
            let mut m = RunManifest::new("variants");
            m.tool_version = "0.1.0".to_string();
            m.inputs.push(FileHash {
                path: idx_path.to_string(),
                blake3: "aa".to_string(),
            });
            m.outputs.push(FileHash {
                path: out_path.to_string(),
                blake3: "bb".to_string(),
            });
            m.params.insert("min_qual".to_string(), "30".to_string());
            m.finalize();
            m
        };
        let a = mk("/home/alice/ref.idx", "/tmp/run-1/out.vcf");
        let b = mk("/data/ref.idx", "out.vcf");
        assert_eq!(
            a.content_hash(),
            b.content_hash(),
            "the claim hash must not depend on recorded paths"
        );
        assert_eq!(a.self_hash_ok(), Some(true));
        assert_eq!(b.self_hash_ok(), Some(true));
    }

    #[test]
    fn claim_hash_still_tracks_content() {
        // Different content (blake3) → different claim hash (the address is the content).
        let mk = |digest: &str| {
            let mut m = RunManifest::new("variants");
            m.tool_version = "0.1.0".to_string();
            m.inputs.push(FileHash {
                path: "ref.idx".to_string(),
                blake3: digest.to_string(),
            });
            m.finalize();
            m
        };
        assert_ne!(mk("aa").content_hash(), mk("bb").content_hash());
    }

    #[test]
    fn claim_drops_paths_but_the_on_disk_form_keeps_them() {
        let mut m = RunManifest::new("variants");
        m.tool_version = "0.1.0".to_string();
        m.inputs.push(FileHash {
            path: "/home/alice/secret/ref.idx".to_string(),
            blake3: "aa".to_string(),
        });
        m.finalize();
        assert!(
            !m.to_canonical_claim_json().contains("/home/alice"),
            "the claim form must not carry the recorded path"
        );
        assert!(
            m.to_canonical_json().contains("/home/alice"),
            "the on-disk form must keep the recorded path"
        );
        // The content digest is present in BOTH.
        assert!(m.to_canonical_claim_json().contains("aa"));
    }

    #[test]
    fn the_schema_gate_keeps_v2_verifying_and_makes_only_v3_path_independent() {
        // Seal a receipt at a given schema with a given input path.
        let seal = |schema: &str, path: &str| {
            let mut m = RunManifest::new("variants");
            m.tool_version = "0.1.0".to_string();
            m.inputs.push(FileHash {
                path: path.to_string(),
                blake3: "aa".to_string(),
            });
            m.params
                .insert("schema_version".to_string(), schema.to_string());
            let h = m.content_hash();
            m.params.insert("manifest_blake3".to_string(), h);
            m
        };
        // Back-compat: the gate reproduces each schema's own claim form, so both
        // self-verify — a pre-P0.2b receipt does not break.
        assert_eq!(
            seal("2", "ref.idx").self_hash_ok(),
            Some(true),
            "schema-2 must still self-verify"
        );
        assert_eq!(
            seal("3", "ref.idx").self_hash_ok(),
            Some(true),
            "schema-3 must self-verify"
        );
        // Isolate the file-render gate: hold the schema string fixed, vary ONLY the
        // path. schema 2 hashed paths into the claim → path-SENSITIVE.
        assert_ne!(
            seal("2", "/a/ref.idx").content_hash(),
            seal("2", "/b/ref.idx").content_hash(),
            "the pre-P0.2b claim form is path-inclusive"
        );
        // schema 3 drops paths → path-INDEPENDENT: the gate fired for the right reason
        // (the render change, not the schema string sitting inside the hashed claim).
        assert_eq!(
            seal("3", "/a/ref.idx").content_hash(),
            seal("3", "/b/ref.idx").content_hash(),
            "the schema-3 claim form is content-only"
        );
    }

    #[test]
    fn finalize_stamps_build_identity_into_the_claim() {
        let mut m = RunManifest::new("variants");
        m.finalize();
        for k in [
            "code_git_sha",
            "code_dirty",
            "rustc_version",
            "target_triple",
            "deps_lock_blake3",
        ] {
            assert!(
                m.params.get(k).is_some_and(|v| !v.is_empty()),
                "finalize must stamp {k}"
            );
        }
        // Build-identity is in the claim → tampering it breaks the claim self-hash.
        assert_eq!(m.self_hash_ok(), Some(true));
        m.params
            .insert("code_git_sha".to_string(), "tampered".to_string());
        assert_eq!(m.self_hash_ok(), Some(false));
    }

    #[test]
    fn check_expected_code_matches_mismatches_and_flags_dirty() {
        let mk = |sha: &str, dirty: &str| {
            let mut m = RunManifest::new("variants");
            m.params.insert("code_git_sha".to_string(), sha.to_string());
            m.params.insert("code_dirty".to_string(), dirty.to_string());
            m
        };
        // Exact + (>=7-char) prefix match on a clean build → no problems.
        assert!(mk("abc123def456", "false")
            .check_expected_code("abc123def456")
            .is_empty());
        assert!(mk("abc123def456", "false")
            .check_expected_code("abc123d")
            .is_empty());
        // Mismatch → one problem mentioning "mismatch".
        let mm = mk("abc123def456", "false").check_expected_code("deadbeef");
        assert_eq!(mm.len(), 1);
        assert!(mm[0].contains("mismatch"));
        // Match but dirty → one problem mentioning the DIRTY tree.
        let dirty = mk("abc123def456", "true").check_expected_code("abc123def456");
        assert_eq!(dirty.len(), 1);
        assert!(dirty[0].contains("DIRTY"));
        // Absent / unknown → cannot check (one problem each).
        assert_eq!(
            RunManifest::new("variants")
                .check_expected_code("abc1234")
                .len(),
            1
        );
        assert_eq!(
            mk("unknown", "false").check_expected_code("abc1234").len(),
            1
        );
    }

    #[test]
    fn check_expected_code_rejects_a_degenerate_expected_sha() {
        // An empty / too-short / non-hex expected must FAIL (not vacuously pass via
        // starts_with) — otherwise `--expect-code "$MAYBE_EMPTY"` is a false-confidence
        // footgun.
        let m = {
            let mut m = RunManifest::new("variants");
            m.params
                .insert("code_git_sha".to_string(), "abc123def456".to_string());
            m.params
                .insert("code_dirty".to_string(), "false".to_string());
            m
        };
        for bad in ["", "9", "abc", "zzzzzzz"] {
            let problems = m.check_expected_code(bad);
            assert_eq!(problems.len(), 1, "{bad:?} must be rejected");
            assert!(
                problems[0].contains("invalid --expect-code"),
                "{bad:?}: {}",
                problems[0]
            );
        }
    }
}
