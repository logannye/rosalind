//! Browser (wasm) bindings for the Rosalind receipt verifier.
//!
//! [`verify`] runs the SAME canonical-JSON + BLAKE3 self-hash check that ships in the
//! `rosalind verify` CLI (via the `rosalind-receipt` crate) — entirely client-side, no
//! upload, no server. It returns a small JSON string the page renders.

use rosalind_receipt::{verify_manifest_str, ReceiptVerdict};
use wasm_bindgen::prelude::*;

/// Verify a receipt's JSON text. Returns a JSON object string:
/// `{"verdict":"verified|tampered|unverifiable|unparseable","detail":"…",
/// "self_hash_ok":true|false|null,"measurement_hash_ok":…,"schema_version":N|null,
/// "subcommand":"…"|null}` — the page does `JSON.parse(verify(text))`.
#[wasm_bindgen]
pub fn verify(json: &str) -> String {
    let c = verify_manifest_str(json);
    let verdict = match c.verdict {
        ReceiptVerdict::Verified => "verified",
        ReceiptVerdict::Tampered => "tampered",
        ReceiptVerdict::Unverifiable => "unverifiable",
        ReceiptVerdict::Unparseable => "unparseable",
    };
    let tribool = |b: Option<bool>| match b {
        Some(true) => "true",
        Some(false) => "false",
        None => "null",
    };
    let opt_u32 = |n: Option<u32>| n.map(|v| v.to_string()).unwrap_or_else(|| "null".to_string());
    let opt_str = |s: Option<String>| match s {
        Some(v) => format!("\"{}\"", json_escape(&v)),
        None => "null".to_string(),
    };
    format!(
        "{{\"verdict\":\"{}\",\"detail\":\"{}\",\"self_hash_ok\":{},\"measurement_hash_ok\":{},\"schema_version\":{},\"subcommand\":{}}}",
        verdict,
        json_escape(&c.detail),
        tribool(c.self_hash_ok),
        tribool(c.measurement_hash_ok),
        opt_u32(c.schema_version),
        opt_str(c.subcommand),
    )
}

/// Minimal JSON string escaping for the small set of characters in our detail/subcommand.
fn json_escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for ch in s.chars() {
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
