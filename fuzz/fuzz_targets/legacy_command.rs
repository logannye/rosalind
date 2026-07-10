#![no_main]

use libfuzzer_sys::fuzz_target;
use rosalind_receipt::{command_template_tokens, RunManifest};

fuzz_target!(|data: &[u8]| {
    if let Ok(command) = std::str::from_utf8(data) {
        let mut receipt = RunManifest::new("variants");
        receipt
            .params
            .insert("schema_version".to_string(), "5".to_string());
        receipt
            .params
            .insert("command".to_string(), command.to_string());
        let _ = command_template_tokens(&receipt);
    }
});
