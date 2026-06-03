//! A minimal, deterministic reproducibility receipt for a run: tool version,
//! subcommand, BLAKE3 content hashes of inputs + outputs, and the parameters.
//! Serialized as canonical JSON (sorted keys, no timestamps) so two identical
//! runs produce a byte-identical manifest. Full `rosalind verify` is a later
//! phase; this phase emits the receipt and proves it is deterministic.

use std::collections::BTreeMap;
use std::io::{self, Read, Write};
use std::path::{Path, PathBuf};

/// Current receipt/feature schema version. Bump on any breaking schema change.
pub const MANIFEST_SCHEMA_VERSION: u32 = 1;

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
        }
    }

    /// Serialize to canonical JSON: keys sorted, `inputs`/`outputs` sorted by
    /// path, no timestamps — so identical runs hash and render identically.
    pub fn to_canonical_json(&self) -> String {
        let mut out = String::new();
        out.push('{');

        out.push_str("\"inputs\":");
        push_file_hashes(&mut out, &self.inputs);

        out.push_str(",\"outputs\":");
        push_file_hashes(&mut out, &self.outputs);

        out.push_str(",\"params\":{");
        for (i, (k, v)) in self.params.iter().enumerate() {
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

    /// BLAKE3 hex of the canonical JSON with the self-hash field excluded — the
    /// content this manifest commits to. Deterministic; `verify` re-derives it.
    pub fn content_hash(&self) -> String {
        let mut m = self.clone();
        m.params.remove("manifest_blake3");
        blake3_hex(m.to_canonical_json().as_bytes())
    }

    /// Stamp the schema version + the self-hash. Call LAST, immediately before
    /// serialization, so the hash covers every other field (including the version).
    pub fn finalize(&mut self) {
        self.params.insert(
            "schema_version".to_string(),
            MANIFEST_SCHEMA_VERSION.to_string(),
        );
        let h = self.content_hash();
        self.params.insert("manifest_blake3".to_string(), h);
    }

    /// `Some(true)`/`Some(false)` if a self-hash is recorded and matches / mismatches;
    /// `None` if none is recorded (a pre-1.2 receipt).
    pub fn self_hash_ok(&self) -> Option<bool> {
        self.params
            .get("manifest_blake3")
            .map(|recorded| *recorded == self.content_hash())
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
            map.insert(k, v);
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
        self.expect_key("outputs")?;
        let outputs = self.parse_file_array()?;
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
        })
    }
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
}
