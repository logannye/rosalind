use std::fs;
use std::path::{Path, PathBuf};

fn snapshot_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("snapshots")
}

pub fn assert_snapshot(name: &str, actual: &str) {
    let path = snapshot_root().join(name);
    if std::env::var("ROSALIND_UPDATE_SNAPSHOTS").is_ok() {
        if let Some(parent) = path.parent() {
            fs::create_dir_all(parent).expect("create snapshot directory");
        }
        fs::write(&path, actual).expect("write snapshot");
        return;
    }

    let expected =
        fs::read_to_string(&path).unwrap_or_else(|_| panic!("snapshot {:?} not found", path));
    if normalize(&expected) != normalize(actual) {
        panic!(
            "Snapshot mismatch for {:?}. Set ROSALIND_UPDATE_SNAPSHOTS=1 to regenerate.\nExpected:\n{}\nActual:\n{}",
            path,
            expected,
            actual
        );
    }
}

fn normalize(input: &str) -> String {
    input.replace("\r\n", "\n")
}
