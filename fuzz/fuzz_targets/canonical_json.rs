#![no_main]

use libfuzzer_sys::fuzz_target;
use rosalind_receipt::RunManifest;

fuzz_target!(|data: &[u8]| {
    if let Ok(text) = std::str::from_utf8(data) {
        let _ = RunManifest::from_canonical_json(text);
    }
});
