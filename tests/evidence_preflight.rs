use rosalind::provenance::RunManifest;
use rust_htslib::bam::{self, header::HeaderRecord, Header};
use std::path::PathBuf;
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};
static FIXTURE_ID: AtomicU64 = AtomicU64::new(0);

struct Fixture {
    dir: PathBuf,
    reference: PathBuf,
    bam: PathBuf,
    bed: PathBuf,
}
impl Fixture {
    fn new() -> Self {
        let dir = std::env::temp_dir().join(format!(
            "rosalind-preflight-{}-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos(),
            FIXTURE_ID.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir(&dir).unwrap();
        let reference = dir.join("ref.fa");
        std::fs::write(&reference, ">chr1\nAAAAAAAAAA\n").unwrap();
        std::fs::write(dir.join("ref.fa.fai"), "chr1\t10\t6\t10\t11\n").unwrap();
        let bam = dir.join("reads.bam");
        let mut header = Header::new();
        header.push_record(
            HeaderRecord::new(b"HD")
                .push_tag(b"VN", "1.6")
                .push_tag(b"SO", "coordinate"),
        );
        header.push_record(
            HeaderRecord::new(b"SQ")
                .push_tag(b"SN", "chr1")
                .push_tag(b"LN", 10),
        );
        drop(bam::Writer::from_path(&bam, &header, bam::Format::Bam).unwrap());
        bam::index::build(&bam, None::<&PathBuf>, bam::index::Type::Bai, 1).unwrap();
        let bed = dir.join("targets.bed");
        std::fs::write(&bed, "chr1\t0\t10\tfirst\n").unwrap();
        Self {
            dir,
            reference,
            bam,
            bed,
        }
    }
    fn command(&self, kind: &str) -> Command {
        let mut command = Command::new(env!("CARGO_BIN_EXE_rosalind"));
        command
            .args(["analyze", kind, "--alignments"])
            .arg(&self.bam)
            .arg("--reference")
            .arg(&self.reference)
            .arg("--regions")
            .arg(&self.bed);
        command
    }
}
impl Drop for Fixture {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.dir);
    }
}
fn refused(command: &mut Command, fragment: &str) {
    let result = command.output().unwrap();
    assert_eq!(
        result.status.code(),
        Some(2),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert!(
        String::from_utf8_lossy(&result.stderr).contains(fragment),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
}

#[test]
fn preflight_rejects_output_receipt_partial_and_input_aliases_before_writing() {
    let f = Fixture::new();
    let output = f.dir.join("out.arrow");
    let original = std::fs::read(&f.bam).unwrap();
    refused(
        f.command("evidence").arg("--force").arg("-o").arg(&f.bam),
        "must not overwrite an input",
    );
    assert_eq!(std::fs::read(&f.bam).unwrap(), original);
    refused(
        f.command("panel-qc")
            .arg("-o")
            .arg(&output)
            .arg("--position-output")
            .arg(f.dir.join("./out.arrow")),
        "must use different paths",
    );
    refused(
        f.command("evidence")
            .arg("-o")
            .arg(&output)
            .arg("--manifest")
            .arg(f.dir.join("out.arrow.partial")),
        "must use different paths",
    );
    assert!(!output.exists());
    std::fs::write(f.dir.join("out.arrow.partial"), b"old diagnostic").unwrap();
    refused(
        f.command("evidence").arg("-o").arg(&output),
        "destination already exists",
    );
    assert_eq!(
        std::fs::read(f.dir.join("out.arrow.partial")).unwrap(),
        b"old diagnostic"
    );
}

#[test]
fn irrelevant_reference_and_panel_options_are_rejected() {
    let f = Fixture::new();
    refused(
        f.command("evidence")
            .arg("--cram-reference")
            .arg(&f.reference),
        "require CRAM input",
    );
    refused(
        f.command("evidence").args(["--min-callable-depth", "17"]),
        "panel-qc",
    );
    let mut command = Command::new(env!("CARGO_BIN_EXE_rosalind"));
    command
        .args(["analyze", "panel-qc", "--alignments"])
        .arg(&f.bam)
        .arg("--regions")
        .arg(&f.bed)
        .arg("--reference-fai")
        .arg(f.dir.join("ref.fa.fai"));
    refused(&mut command, "requires an indexed FASTA analysis reference");
}

#[test]
fn scientific_identity_uses_normalized_selection_and_captures_index_content() {
    let f = Fixture::new();
    let one = f.dir.join("one.tsv");
    let two = f.dir.join("two.tsv");
    for (output, bed) in [
        (&one, "chr1\t0\t10\tfirst\n"),
        (
            &two,
            "chr1\t5\t10\tsecond\nchr1\t0\t8\tfirst\nchr1\t1\t3\tduplicate\n",
        ),
    ] {
        std::fs::write(&f.bed, bed).unwrap();
        let result = f
            .command("evidence")
            .arg("-o")
            .arg(output)
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
    }
    let receipt = |p: &PathBuf| {
        RunManifest::from_canonical_json(
            &std::fs::read_to_string(format!("{}.manifest.json", p.display())).unwrap(),
        )
        .unwrap()
    };
    let a = receipt(&one);
    let b = receipt(&two);
    assert_eq!(
        a.params["evidence.science_blake3"],
        b.params["evidence.science_blake3"]
    );
    assert_eq!(std::fs::read(one).unwrap(), std::fs::read(two).unwrap());
    assert!(a.inputs.iter().any(|i| i.path.ends_with(".bai")));
}
