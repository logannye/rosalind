use rosalind_receipt::{FileHash, RunManifest};

fn base(schema: u32) -> RunManifest {
    let mut manifest = RunManifest::new("variants");
    manifest.inputs = vec![FileHash {
        path: "/legacy/input.bam".to_string(),
        blake3: "1111111111111111111111111111111111111111111111111111111111111111".to_string(),
    }];
    manifest.outputs = vec![FileHash {
        path: "/legacy/output.vcf".to_string(),
        blake3: "2222222222222222222222222222222222222222222222222222222222222222".to_string(),
    }];
    manifest
        .params
        .insert("memory_budget_mb".to_string(), "128".to_string());
    manifest
        .params
        .insert("schema_version".to_string(), schema.to_string());
    manifest
}

fn seal_claim(manifest: &mut RunManifest) {
    let hash = manifest.content_hash();
    manifest.params.insert("manifest_blake3".to_string(), hash);
}

fn add_measurements(manifest: &mut RunManifest) {
    manifest
        .measurements
        .insert("contract_verdict".to_string(), "within".to_string());
    manifest
        .measurements
        .insert("peak_rss_bytes".to_string(), "33554432".to_string());
    let hash = manifest.measurement_hash();
    manifest
        .measurements
        .insert("measurement_blake3".to_string(), hash);
    manifest
        .params
        .insert("has_measurements".to_string(), "true".to_string());
}

fn main() {
    for schema in 1..=5 {
        let mut manifest = base(schema);
        if schema == 1 {
            manifest
                .params
                .insert("contract_verdict".to_string(), "within".to_string());
            manifest
                .params
                .insert("peak_rss_bytes".to_string(), "33554432".to_string());
        } else {
            add_measurements(&mut manifest);
        }
        if schema >= 4 {
            for (key, value) in [
                ("code_git_sha", "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"),
                ("code_dirty", "false"),
                ("rustc_version", "rustc fixture"),
                ("target_triple", "x86_64-unknown-linux-gnu"),
                (
                    "deps_lock_blake3",
                    "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb",
                ),
            ] {
                manifest.params.insert(key.to_string(), value.to_string());
            }
        }
        if schema >= 5 {
            manifest.params.insert(
                "command".to_string(),
                "variants --alignments @in:1111111111111111111111111111111111111111111111111111111111111111 -o @out:2222222222222222222222222222222222222222222222222222222222222222"
                    .to_string(),
            );
        }
        seal_claim(&mut manifest);
        println!("SCHEMA {schema}\n{}", manifest.to_canonical_json());
    }
}
