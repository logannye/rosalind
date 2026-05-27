//! A minimal, deterministic reproducibility receipt for a run: tool version,
//! subcommand, BLAKE3 content hashes of inputs + outputs, and the parameters.
//! Serialized as canonical JSON (sorted keys, no timestamps) so two identical
//! runs produce a byte-identical manifest. Full `rosalind verify` is a later
//! phase; this phase emits the receipt and proves it is deterministic.

use std::collections::BTreeMap;
use std::io::{self, Read, Write};
use std::path::{Path, PathBuf};

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
}
