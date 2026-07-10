use rosalind_receipt::{command_template_tokens, RunManifest};

#[test]
fn schemas_one_through_five_remain_parseable_and_historically_verifiable() {
    let fixtures = [
        include_str!("../fixtures/receipt-schema-1.json"),
        include_str!("../fixtures/receipt-schema-2.json"),
        include_str!("../fixtures/receipt-schema-3.json"),
        include_str!("../fixtures/receipt-schema-4.json"),
        include_str!("../fixtures/receipt-schema-5.json"),
    ];

    for (offset, fixture) in fixtures.iter().enumerate() {
        let schema = offset + 1;
        let manifest = RunManifest::from_canonical_json(fixture.trim()).unwrap();
        assert_eq!(
            manifest
                .params
                .get("schema_version")
                .and_then(|value| value.parse::<usize>().ok()),
            Some(schema)
        );
        assert_eq!(manifest.self_hash_ok(), Some(true), "schema {schema}");
        if schema == 1 {
            assert_eq!(manifest.measurement_hash_ok(), None);
            assert!(manifest.measurements.is_empty());
        } else {
            assert_eq!(
                manifest.measurement_hash_ok(),
                Some(true),
                "schema {schema}"
            );
        }
        if schema < 5 {
            assert!(command_template_tokens(&manifest).is_err());
        } else {
            assert!(command_template_tokens(&manifest).is_ok());
        }
    }
}

#[test]
fn schema_three_and_newer_paths_are_portable_metadata() {
    for fixture in [
        include_str!("../fixtures/receipt-schema-3.json"),
        include_str!("../fixtures/receipt-schema-4.json"),
        include_str!("../fixtures/receipt-schema-5.json"),
    ] {
        let manifest = RunManifest::from_canonical_json(fixture.trim()).unwrap();
        let claim = manifest.content_hash();
        let relocated = fixture.replace("/legacy/", "/relocated/data/");
        let relocated = RunManifest::from_canonical_json(relocated.trim()).unwrap();
        assert_eq!(relocated.content_hash(), claim);
        assert_eq!(relocated.self_hash_ok(), Some(true));
    }
}
