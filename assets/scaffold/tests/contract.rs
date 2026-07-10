use rosalind::contract::testkit::assert_paths_are_portable;
use rosalind::provenance::{FileHash, RunManifest};

#[test]
fn receipt_claim_is_path_portable() {
    let mut manifest = RunManifest::new("run");
    manifest.inputs.push(FileHash {
        path: "/one/input.bam".into(),
        blake3: "abc".into(),
    });
    manifest.outputs.push(FileHash {
        path: "/one/output.tsv".into(),
        blake3: "def".into(),
    });
    manifest.finalize();
    assert_paths_are_portable(&manifest);
}
