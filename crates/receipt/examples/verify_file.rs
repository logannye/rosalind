//! Check a receipt's tamper-evident integrity from its `*.manifest.json` alone — the
//! same self-hash check the in-browser verifier (and `rosalind verify`) runs, without
//! re-hashing the input/output files.
//!
//!   cargo run -p rosalind-receipt --example verify_file -- path/to/sample.manifest.json
//!
//! Exit code: 0 = Verified, 1 = anything else (Tampered / Unverifiable / Unparseable).

use std::process::ExitCode;

use rosalind_receipt::{verify_manifest_str, ReceiptVerdict};

fn main() -> ExitCode {
    let Some(path) = std::env::args().nth(1) else {
        eprintln!("usage: verify_file <path/to/*.manifest.json>");
        return ExitCode::from(2);
    };
    let json = match std::fs::read_to_string(&path) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("cannot read {path}: {e}");
            return ExitCode::from(2);
        }
    };
    let c = verify_manifest_str(&json);
    println!(
        "{:?}  (self_hash={:?}, measurement_hash={:?}, schema={:?}, subcommand={:?})\n  {}",
        c.verdict, c.self_hash_ok, c.measurement_hash_ok, c.schema_version, c.subcommand, c.detail
    );
    match c.verdict {
        ReceiptVerdict::Verified => ExitCode::SUCCESS,
        _ => ExitCode::FAILURE,
    }
}
