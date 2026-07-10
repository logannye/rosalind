use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicU64, Ordering};

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_rosalind")
}

fn run(args: &[&str]) -> Output {
    Command::new(bin()).args(args).output().unwrap()
}

fn ok(args: &[&str]) {
    let output = run(args);
    assert!(
        output.status.success(),
        "stdout: {}\nstderr: {}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
}

fn dir() -> PathBuf {
    static COUNTER: AtomicU64 = AtomicU64::new(0);
    let path = std::env::temp_dir().join(format!(
        "rosalind-atomic-cli-{}-{}",
        std::process::id(),
        COUNTER.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir_all(&path).unwrap();
    path
}

fn write_inputs(dir: &Path) -> (PathBuf, PathBuf) {
    let reference = dir.join("reference.fa");
    let reads = dir.join("reads.fastq");
    std::fs::write(&reference, ">chr1\nACGTACGTACGTACGTACGTACGTACGT\n").unwrap();
    std::fs::write(&reads, "@read\nACGTACGT\n+\nIIIIIIII\n").unwrap();
    (reference, reads)
}

#[test]
fn index_align_sort_and_variants_refuse_collisions_and_force_replaces() {
    let dir = dir();
    let (reference, reads) = write_inputs(&dir);
    let index = dir.join("reference.idx");
    let raw = dir.join("raw.bam");
    let sorted = dir.join("sorted.bam");
    let calls = dir.join("calls.vcf");

    let index_args = [
        "index",
        "--reference",
        reference.to_str().unwrap(),
        "--output",
        index.to_str().unwrap(),
    ];
    ok(&index_args);
    let original_index = std::fs::read(&index).unwrap();
    assert_eq!(run(&index_args).status.code(), Some(2));
    assert_eq!(std::fs::read(&index).unwrap(), original_index);
    let mut forced_index = index_args.to_vec();
    forced_index.push("--force");
    ok(&forced_index);
    assert_eq!(std::fs::read(&index).unwrap(), original_index);

    let align_args = [
        "align",
        "--reference",
        reference.to_str().unwrap(),
        "--reads",
        reads.to_str().unwrap(),
        "--format",
        "bam",
        "--output",
        raw.to_str().unwrap(),
    ];
    ok(&align_args);
    assert_eq!(run(&align_args).status.code(), Some(2));

    let sort_args = [
        "sort",
        "--input",
        raw.to_str().unwrap(),
        "--output",
        sorted.to_str().unwrap(),
    ];
    ok(&sort_args);
    assert_eq!(run(&sort_args).status.code(), Some(2));

    let variants_args = [
        "variants",
        "--index",
        index.to_str().unwrap(),
        "--alignments",
        sorted.to_str().unwrap(),
        "--output",
        calls.to_str().unwrap(),
    ];
    ok(&variants_args);
    let original_calls = std::fs::read(&calls).unwrap();
    assert_eq!(run(&variants_args).status.code(), Some(2));
    assert_eq!(std::fs::read(&calls).unwrap(), original_calls);
    let mut forced_variants = variants_args.to_vec();
    forced_variants.push("--force");
    ok(&forced_variants);
    assert_eq!(std::fs::read(&calls).unwrap(), original_calls);

    std::fs::remove_dir_all(dir).ok();
}

#[test]
fn ordinary_input_failure_leaves_no_output_or_temporary_file() {
    let dir = dir();
    let output = dir.join("never.bam");
    let result = run(&[
        "align",
        "--reference",
        dir.join("missing.fa").to_str().unwrap(),
        "--reads",
        dir.join("missing.fastq").to_str().unwrap(),
        "--format",
        "bam",
        "--output",
        output.to_str().unwrap(),
    ]);
    assert!(!result.status.success());
    assert!(!output.exists());
    let names = std::fs::read_dir(&dir)
        .unwrap()
        .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
        .collect::<Vec<_>>();
    assert!(
        names.iter().all(|name| !name.contains(".rosalind-")),
        "temporary files survived: {names:?}"
    );
    std::fs::remove_dir_all(dir).ok();
}
