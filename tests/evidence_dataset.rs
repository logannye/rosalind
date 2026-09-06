use rosalind::dataset::{plan_dataset, run_dataset, DatasetError, DatasetOptions};
use rosalind::evidence::*;
use rosalind::provenance::RunManifest;
use rosalind::selection::GenomicInterval;
use rust_htslib::bam::record::{Cigar, CigarString, Record};
use rust_htslib::bam::{self};
use std::path::PathBuf;
use std::time::{SystemTime, UNIX_EPOCH};

struct Fixture {
    root: PathBuf,
    fasta: PathBuf,
    bam: PathBuf,
}
impl Fixture {
    fn new() -> Self {
        static COUNTER: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
        let root = std::env::temp_dir().join(format!(
            "rosalind-dataset-{}-{}-{}",
            std::process::id(),
            COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed),
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&root).unwrap();
        let fasta = root.join("reference.fa");
        let bam = root.join("reads.bam");
        let length = 32780;
        std::fs::write(&fasta, format!(">chr1\n{}\n", "A".repeat(length))).unwrap();
        std::fs::write(
            root.join("reference.fa.fai"),
            format!("chr1\t{length}\t6\t{length}\t{}\n", length + 1),
        )
        .unwrap();
        let mut header = bam::Header::new();
        let mut hd = bam::header::HeaderRecord::new(b"HD");
        hd.push_tag(b"VN", "1.6");
        hd.push_tag(b"SO", "coordinate");
        header.push_record(&hd);
        let mut sq = bam::header::HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr1");
        sq.push_tag(b"LN", length);
        header.push_record(&sq);
        let mut writer = bam::Writer::from_path(&bam, &header, bam::Format::Bam).unwrap();
        for (index, position) in [16378, 16381, 32764].into_iter().enumerate() {
            let mut record = Record::new();
            record.set(
                format!("r{index}").as_bytes(),
                Some(&CigarString(vec![Cigar::Match(12)])),
                b"AAACAAAAGAAA",
                &[30; 12],
            );
            record.set_tid(0);
            record.set_pos(position);
            record.set_mapq(60);
            writer.write(&record).unwrap();
        }
        drop(writer);
        bam::index::build(&bam, None::<&PathBuf>, bam::index::Type::Bai, 1).unwrap();
        Self { root, fasta, bam }
    }
    fn request(&self, width: u32) -> EvidenceRequest {
        let mut request = EvidenceRequest::new(&self.bam, &self.fasta);
        request.execution.max_microtile_bases = width;
        request.selection = EvidenceSelection::Intervals(vec![
            GenomicInterval {
                contig: 0,
                start: 16380,
                end: 16390,
            },
            GenomicInterval {
                contig: 0,
                start: 32766,
                end: 32770,
            },
        ]);
        request
    }
    fn options(&self, workers: usize, resume: bool) -> DatasetOptions {
        DatasetOptions {
            cache_dir: self.root.join("cache"),
            workers,
            resume,
        }
    }
    fn science(&self) -> String {
        blake3::hash(b"verified fixture science identity")
            .to_hex()
            .to_string()
    }
    fn run(
        &self,
        request: EvidenceRequest,
        options: DatasetOptions,
    ) -> (rosalind::dataset::DatasetOutcome, Vec<u8>) {
        let mut engine = EvidenceEngine::open(request).unwrap();
        let mut writer = EvidenceArrowWriter::new(Vec::new());
        let outcome = run_dataset(&mut engine, &self.science(), &options, &mut writer).unwrap();
        (outcome, writer.into_inner().unwrap())
    }
}
impl Drop for Fixture {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.root);
    }
}

#[test]
fn canonical_cache_workers_and_microtiles_produce_identical_evidence() {
    let fixture = Fixture::new();
    let mut direct = EvidenceEngine::open(fixture.request(7)).unwrap();
    let mut writer = EvidenceArrowWriter::new(Vec::new());
    direct.run(&mut writer).unwrap();
    let expected = writer.into_inner().unwrap();
    let (first, one) = fixture.run(fixture.request(1), fixture.options(1, false));
    assert_eq!(first.plan.partition_count, 3);
    assert_eq!(first.computed_partitions, 3);
    assert_eq!(first.stats.emitted_loci, 14);
    assert_eq!(one, expected);
    let mut parallel_options = fixture.options(2, false);
    parallel_options.cache_dir = fixture.root.join("parallel");
    let (parallel, two) = fixture.run(fixture.request(16384), parallel_options);
    assert_eq!(parallel.plan.worker_count, 2);
    assert_eq!(two, expected);
    let (reused, three) = fixture.run(fixture.request(3), fixture.options(2, false));
    assert_eq!(reused.computed_partitions, 0);
    assert_eq!(reused.reused_partitions, 3);
    assert_eq!(reused.stats.record_visits, 0);
    assert_eq!(three, expected);
}

#[test]
fn incomplete_cache_requires_resume_and_recomputes_only_missing_partition() {
    let fixture = Fixture::new();
    let (complete, expected) = fixture.run(fixture.request(2), fixture.options(2, false));
    let root = complete.dataset_manifest.parent().unwrap();
    std::fs::remove_file(&complete.dataset_manifest).unwrap();
    std::fs::remove_dir_all(root.join("c00000000-p0000016384")).unwrap();
    let mut engine = EvidenceEngine::open(fixture.request(11)).unwrap();
    let mut sink = EvidenceArrowWriter::new(Vec::new());
    let error = run_dataset(
        &mut engine,
        &fixture.science(),
        &fixture.options(2, false),
        &mut sink,
    )
    .unwrap_err();
    assert!(matches!(error, DatasetError::Incompatible(_)));
    let (resumed, actual) = fixture.run(fixture.request(5), fixture.options(2, true));
    assert_eq!(resumed.reused_partitions, 2);
    assert_eq!(resumed.computed_partitions, 1);
    assert_eq!(actual, expected);
}

#[test]
fn resumed_cache_hash_corruption_is_refused_before_analyzer_receives_rows() {
    let fixture = Fixture::new();
    let (complete, _) = fixture.run(fixture.request(2), fixture.options(2, false));
    let root = complete.dataset_manifest.parent().unwrap();
    let artifact = root.join("c00000000-p0000032768/evidence.arrow");
    let mut bytes = std::fs::read(&artifact).unwrap();
    bytes[100] ^= 1;
    std::fs::write(artifact, bytes).unwrap();
    #[derive(Default)]
    struct Count(usize);
    impl EvidenceAnalyzer for Count {
        fn on_batch(&mut self, batch: &EvidenceBatch) -> Result<(), EvidenceError> {
            self.0 += batch.len();
            Ok(())
        }
        fn additional_memory_bytes(&self) -> Option<u64> {
            Some(0)
        }
    }
    let mut analyzer = Count::default();
    let mut engine = EvidenceEngine::open(fixture.request(4)).unwrap();
    let error = run_dataset(
        &mut engine,
        &fixture.science(),
        &fixture.options(2, true),
        &mut analyzer,
    )
    .unwrap_err();
    assert!(matches!(error, DatasetError::Corrupt(_)));
    assert_eq!(error.exit_code(), 5);
    assert_eq!(analyzer.0, 0);
}

#[test]
fn completed_dataset_binds_partition_receipts_and_missing_partitions() {
    let fixture = Fixture::new();
    let (complete, _) = fixture.run(fixture.request(2), fixture.options(2, false));
    let root = complete.dataset_manifest.parent().unwrap();
    let receipt = root.join("c00000000-p0000000000/manifest.json");
    let mut manifest =
        RunManifest::from_canonical_json(&std::fs::read_to_string(&receipt).unwrap()).unwrap();
    manifest
        .params
        .insert("unbound.change".into(), "tampered".into());
    manifest.finalize();
    std::fs::write(receipt, manifest.to_canonical_json()).unwrap();
    let mut engine = EvidenceEngine::open(fixture.request(4)).unwrap();
    let error = run_dataset(
        &mut engine,
        &fixture.science(),
        &fixture.options(1, true),
        &mut EvidenceArrowWriter::new(std::io::sink()),
    )
    .unwrap_err();
    assert!(matches!(error, DatasetError::Corrupt(_)));
    std::fs::remove_dir_all(root.join("c00000000-p0000000000")).unwrap();
    let error = run_dataset(
        &mut engine,
        &fixture.science(),
        &fixture.options(1, true),
        &mut EvidenceArrowWriter::new(std::io::sink()),
    )
    .unwrap_err();
    assert!(matches!(error, DatasetError::Corrupt(_)));
}

#[test]
fn snv_site_annotations_survive_partition_cache_and_execution_changes() {
    let fixture = Fixture::new();
    let mut request = fixture.request(1);
    request.selection = EvidenceSelection::Sites(vec![
        SnvSite {
            contig: 0,
            position: 16383,
            reference: b'A',
            alternates: vec![b'C', b'G'],
        },
        SnvSite {
            contig: 0,
            position: 16384,
            reference: b'A',
            alternates: vec![b'T'],
        },
        SnvSite {
            contig: 0,
            position: 32768,
            reference: b'A',
            alternates: vec![b'G'],
        },
    ]);
    let mut direct = EvidenceEngine::open(request.clone()).unwrap();
    let mut writer = EvidenceArrowWriter::new(Vec::new());
    direct.run(&mut writer).unwrap();
    let expected = writer.into_inner().unwrap();
    let (_, actual) = fixture.run(request.clone(), fixture.options(3, false));
    assert_eq!(actual, expected);
    request.execution.max_microtile_bases = 16384;
    let (cached, actual) = fixture.run(request, fixture.options(1, false));
    assert_eq!(cached.reused_partitions, 3);
    assert_eq!(actual, expected);
}

#[test]
fn worker_admission_counts_single_baseline_and_refuses_before_cache_creation() {
    // The admission model deliberately includes process RSS. Isolate this
    // fixed-budget case from concurrent tests' Arrow buffers and allocator
    // high-water marks, just as a real CLI invocation has its own process.
    const ISOLATED: &str = "ROSALIND_TEST_WORKER_ADMISSION_ISOLATED";
    if std::env::var_os(ISOLATED).is_none() {
        let output = std::process::Command::new(std::env::current_exe().unwrap())
            .args([
                "--exact",
                "worker_admission_counts_single_baseline_and_refuses_before_cache_creation",
                "--nocapture",
            ])
            .env(ISOLATED, "1")
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "isolated admission regression failed:\n{}\n{}",
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        );
        return;
    }
    let fixture = Fixture::new();
    let mut request = fixture.request(16384);
    request.execution.memory_budget_bytes = Some(100 << 20);
    let mut engine = EvidenceEngine::open(request).unwrap();
    let writer = EvidenceArrowWriter::new(std::io::sink());
    for workers in [0, 65] {
        assert!(matches!(
            plan_dataset(&mut engine, &fixture.options(workers, false), &writer),
            Err(DatasetError::Incompatible(_))
        ));
    }
    let error = run_dataset(
        &mut engine,
        &fixture.science(),
        &fixture.options(3, false),
        &mut EvidenceArrowWriter::new(std::io::sink()),
    )
    .unwrap_err();
    assert_eq!(error.exit_code(), 3);
    assert!(!fixture.root.join("cache").exists());
    let mut engine = EvidenceEngine::open(fixture.request(4)).unwrap();
    let plan = plan_dataset(&mut engine, &fixture.options(3, false), &writer).unwrap();
    assert_eq!(
        plan.predicted_peak_rss_bytes,
        plan.baseline_rss_bytes
            + plan.worker_bytes * 3
            + plan.reducer_bytes
            + plan.metadata_queue_bytes
    );
}

fn cli_evidence(fixture: &Fixture, name: &str, format: &str, region: &str, baseq: &str) -> PathBuf {
    let output = fixture.root.join(name);
    let bed = fixture.root.join(format!("{name}.bed"));
    let (contig, span) = region.split_once(':').unwrap();
    let (start, end) = span.split_once('-').unwrap();
    std::fs::write(
        &bed,
        format!("{contig}\t{}\t{end}\n", start.parse::<u32>().unwrap() - 1),
    )
    .unwrap();
    let result = std::process::Command::new(env!("CARGO_BIN_EXE_rosalind"))
        .args(["analyze", "evidence", "--reference"])
        .arg(&fixture.fasta)
        .arg("--alignments")
        .arg(&fixture.bam)
        .arg("--regions")
        .arg(&bed)
        .args(["--format", format, "--base-quality-threshold", baseq, "-o"])
        .arg(&output)
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    PathBuf::from(format!("{}.manifest.json", output.display()))
}

#[test]
fn dataset_diff_is_empty_across_arrow_and_tsv_and_reports_profile_reference_and_locus_changes() {
    use rosalind::dataset_diff::diff_evidence_datasets;
    let first = Fixture::new();
    let arrow = cli_evidence(&first, "same.arrow", "arrow-ipc", "chr1:16381-16390", "20");
    let tsv = cli_evidence(&first, "same.tsv", "tsv", "chr1:16381-16390", "20");
    let changes = first.root.join("same-delta.tsv");
    let summary = diff_evidence_datasets(&arrow, &tsv, &changes).unwrap();
    assert_eq!(
        summary,
        rosalind::dataset_diff::DatasetDiffSummary::default()
    );
    assert_eq!(std::fs::read_to_string(changes).unwrap().lines().count(), 1);
    let report = std::process::Command::new(env!("CARGO_BIN_EXE_rosalind"))
        .arg("diff")
        .arg(&arrow)
        .arg(&tsv)
        .arg("--json")
        .output()
        .unwrap();
    assert_eq!(report.status.code(), Some(1));
    let report = String::from_utf8(report.stdout).unwrap();
    assert!(report.contains("\"science_params\":0"), "{report}");
    let execution_count = report
        .split("\"execution_params\":")
        .nth(1)
        .unwrap()
        .split(',')
        .next()
        .unwrap()
        .parse::<u64>()
        .unwrap();
    assert!(execution_count >= 2);
    let second = Fixture::new();
    let mut reference = std::fs::read(&second.fasta).unwrap();
    reference[6 + 16383] = b'G';
    std::fs::write(&second.fasta, reference).unwrap();
    let changed = cli_evidence(
        &second,
        "changed.arrow",
        "arrow-ipc",
        "chr1:16384-16392",
        "31",
    );
    let changes = first.root.join("changed-delta.tsv");
    let result = std::process::Command::new(env!("CARGO_BIN_EXE_rosalind"))
        .arg("diff")
        .arg(&tsv)
        .arg(&changed)
        .arg("--loci-output")
        .arg(&changes)
        .output()
        .unwrap();
    assert_eq!(
        result.status.code(),
        Some(1),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let causes = String::from_utf8_lossy(&result.stdout);
    assert!(
        causes.contains("CAUSE") && causes.contains("input") && causes.contains("base_quality"),
        "{causes}"
    );
    let deltas = std::fs::read_to_string(&changes).unwrap();
    assert!(
        deltas.contains("chr1\t16384\tchanged\tref\tA\tG"),
        "{deltas}"
    );
    assert!(deltas.contains("changed\tcallable_depth\t2\t0"), "{deltas}");
    assert!(
        deltas.contains("chr1\t16381\tremoved\tref\tA\t."),
        "{deltas}"
    );
    assert!(deltas.contains("chr1\t16391\tadded\tref\t.\tA"), "{deltas}");
}

#[test]
fn dataset_diff_refuses_tampered_artifacts_before_publishing_delta_output() {
    use rosalind::dataset_diff::{diff_evidence_datasets, DatasetDiffError};
    let fixture = Fixture::new();
    let receipt = cli_evidence(&fixture, "exact.tsv", "tsv", "chr1:16381-16390", "20");
    let artifact = fixture.root.join("exact.tsv");
    let mut bytes = std::fs::read(&artifact).unwrap();
    bytes[0] ^= 1;
    std::fs::write(artifact, bytes).unwrap();
    let changes = fixture.root.join("tampered-delta.tsv");
    let error = diff_evidence_datasets(&receipt, &receipt, &changes).unwrap_err();
    assert!(matches!(error, DatasetDiffError::Integrity(_)));
    assert_eq!(error.exit_code(), 5);
    assert!(!changes.exists());
}

#[test]
fn dataset_diff_compares_panel_position_artifact_and_refuses_missing_positions() {
    use rosalind::dataset_diff::{diff_evidence_datasets, DatasetDiffError};
    let fixture = Fixture::new();
    let bed = fixture.root.join("targets.bed");
    std::fs::write(&bed, "chr1\t16380\t16390\ttarget\n").unwrap();
    let exact = cli_evidence(
        &fixture,
        "exact.arrow",
        "arrow-ipc",
        "chr1:16381-16390",
        "20",
    );
    let mut receipts = Vec::new();
    for positions in [true, false] {
        let output = fixture
            .root
            .join(if positions { "with.tsv" } else { "without.tsv" });
        let mut command = std::process::Command::new(env!("CARGO_BIN_EXE_rosalind"));
        command
            .args(["analyze", "panel-qc", "--reference"])
            .arg(&fixture.fasta)
            .arg("--alignments")
            .arg(&fixture.bam)
            .arg("--regions")
            .arg(&bed)
            .arg("-o")
            .arg(&output);
        if positions {
            command
                .arg("--position-output")
                .arg(fixture.root.join("positions.arrow"));
        }
        let result = command.output().unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        receipts.push(PathBuf::from(format!("{}.manifest.json", output.display())));
    }
    let changes = fixture.root.join("panel-delta.tsv");
    let summary = diff_evidence_datasets(&exact, &receipts[0], &changes).unwrap();
    assert_eq!(summary.metric_changes, 0);
    let changes = fixture.root.join("missing-delta.tsv");
    let error = diff_evidence_datasets(&exact, &receipts[1], &changes).unwrap_err();
    assert!(matches!(error, DatasetDiffError::Incompatible(_)));
    assert!(!changes.exists());
}

#[test]
fn dataset_diff_refuses_different_reference_dictionaries() {
    use rosalind::dataset_diff::{diff_evidence_datasets, DatasetDiffError};
    let first = Fixture::new();
    let second = Fixture::new();
    let length = 32781;
    std::fs::write(&second.fasta, format!(">chr1\n{}\n", "A".repeat(length))).unwrap();
    std::fs::write(
        second.root.join("reference.fa.fai"),
        format!("chr1\t{length}\t6\t{length}\t{}\n", length + 1),
    )
    .unwrap();
    let mut header = bam::Header::new();
    let mut hd = bam::header::HeaderRecord::new(b"HD");
    hd.push_tag(b"SO", "coordinate");
    header.push_record(&hd);
    let mut sq = bam::header::HeaderRecord::new(b"SQ");
    sq.push_tag(b"SN", "chr1");
    sq.push_tag(b"LN", length);
    header.push_record(&sq);
    drop(bam::Writer::from_path(&second.bam, &header, bam::Format::Bam).unwrap());
    std::fs::remove_file(second.root.join("reads.bam.bai")).unwrap();
    bam::index::build(&second.bam, None::<&PathBuf>, bam::index::Type::Bai, 1).unwrap();
    let a = cli_evidence(&first, "a.tsv", "tsv", "chr1:16381-16390", "20");
    let b = cli_evidence(&second, "b.tsv", "tsv", "chr1:16381-16390", "20");
    let output = first.root.join("dictionary-delta.tsv");
    let error = diff_evidence_datasets(&a, &b, &output).unwrap_err();
    assert!(matches!(error, DatasetDiffError::Incompatible(_)));
    assert!(error.to_string().contains("dictionaries differ"));
    assert!(!output.exists());
}

#[test]
fn failed_worker_disconnects_a_full_job_queue_and_leaves_no_complete_dataset() {
    let fixture = Fixture::new();
    let cache = fixture.root.join("cache").join(fixture.science());
    let mut request = fixture.request(1);
    request.execution.max_read_len = 1;
    let options = fixture.options(1, false);
    let science = fixture.science();
    let (sender, receiver) = std::sync::mpsc::channel();
    let worker = std::thread::spawn(move || {
        let mut engine = EvidenceEngine::open(request).unwrap();
        let result = run_dataset(
            &mut engine,
            &science,
            &options,
            &mut EvidenceArrowWriter::new(std::io::sink()),
        );
        sender.send(result).unwrap();
    });
    let result = receiver
        .recv_timeout(std::time::Duration::from_secs(5))
        .expect("failed worker must disconnect the bounded queue without deadlock");
    assert!(matches!(
        result,
        Err(DatasetError::Evidence(EvidenceError::RecordLimit(_)))
    ));
    worker.join().unwrap();
    assert!(!cache.join("dataset.manifest.json").exists());
    assert!(!cache.join("c00000000-p0000000000").exists());
}

#[test]
fn reused_caller_cache_key_cannot_override_request_scientific_semantics() {
    let fixture = Fixture::new();
    let (complete, _) = fixture.run(fixture.request(3), fixture.options(1, false));
    for change in 0..3 {
        let mut request = fixture.request(7);
        match change {
            0 => request.profile.min_mapq = 61,
            1 => request.profile.exclude_duplicates = false,
            _ => {
                request.selection = EvidenceSelection::Intervals(vec![
                    GenomicInterval {
                        contig: 0,
                        start: 16381,
                        end: 16390,
                    },
                    GenomicInterval {
                        contig: 0,
                        start: 32766,
                        end: 32770,
                    },
                ])
            }
        }
        let mut engine = EvidenceEngine::open(request).unwrap();
        let mut called = false;
        let mut analyzer = EvidenceCallback::new(
            |_: &EvidenceBatch| {
                called = true;
                Ok(())
            },
            0,
        );
        let error = run_dataset(
            &mut engine,
            &fixture.science(),
            &fixture.options(1, false),
            &mut analyzer,
        )
        .unwrap_err();
        assert!(matches!(error, DatasetError::Incompatible(_)), "{error}");
        assert!(!called);
    }
    // The same protection applies to a resumable cache with no final receipt.
    std::fs::remove_file(&complete.dataset_manifest).unwrap();
    let mut request = fixture.request(1);
    request.profile.min_base_quality = 31;
    let mut engine = EvidenceEngine::open(request).unwrap();
    let mut analyzer = EvidenceCallback::new(|_: &EvidenceBatch| Ok(()), 0);
    let error = run_dataset(
        &mut engine,
        &fixture.science(),
        &fixture.options(1, true),
        &mut analyzer,
    )
    .unwrap_err();
    assert!(matches!(error, DatasetError::Incompatible(_)), "{error}");
}

#[test]
fn input_snapshot_refuses_replacement_and_mutation_before_cache_publication() {
    use rosalind::dataset::{run_dataset_with_snapshot, InputSnapshot};
    let fixture = Fixture::new();
    let snapshot = InputSnapshot::capture([fixture.fasta.clone(), fixture.bam.clone()]).unwrap();
    let mut engine = EvidenceEngine::open(fixture.request(1)).unwrap();
    let replacement = fixture.root.join("replacement.fa");
    std::fs::copy(&fixture.fasta, &replacement).unwrap();
    std::fs::rename(replacement, &fixture.fasta).unwrap();
    let mut analyzer = EvidenceCallback::new(|_: &EvidenceBatch| Ok(()), 0);
    let error = run_dataset_with_snapshot(
        &mut engine,
        &fixture.science(),
        &fixture.options(1, false),
        &mut analyzer,
        &snapshot,
    )
    .unwrap_err();
    assert!(
        error.to_string().contains("input changed during analysis"),
        "{error}"
    );
    assert!(!fixture.options(1, false).cache_dir.exists());
}

#[test]
fn worker_rejects_source_mutation_before_publishing_its_partition() {
    let fixture = Fixture::new();
    let mut request = fixture.request(1);
    request.selection = EvidenceSelection::Intervals(vec![GenomicInterval {
        contig: 0,
        start: 0,
        end: CANONICAL_TILE_BASES,
    }]);
    let mut engine = EvidenceEngine::open(request).unwrap();
    let options = fixture.options(1, false);
    let cache = options.cache_dir.join(fixture.science());
    let watch_cache = cache.clone();
    let fasta = fixture.fasta.clone();
    let mutator = std::thread::spawn(move || {
        let deadline = std::time::Instant::now() + std::time::Duration::from_secs(5);
        while std::time::Instant::now() < deadline {
            if std::fs::read_dir(&watch_cache).is_ok_and(|entries| {
                entries
                    .flatten()
                    .any(|entry| entry.file_name().to_string_lossy().starts_with(".c"))
            }) {
                use std::io::Write;
                std::fs::OpenOptions::new()
                    .append(true)
                    .open(&fasta)
                    .unwrap()
                    .write_all(b"\n")
                    .unwrap();
                return true;
            }
            std::thread::sleep(std::time::Duration::from_millis(1));
        }
        false
    });
    let mut analyzer = EvidenceCallback::new(|_: &EvidenceBatch| Ok(()), 0);
    let result = run_dataset(&mut engine, &fixture.science(), &options, &mut analyzer);
    assert!(
        mutator.join().unwrap(),
        "worker staging directory was never observed"
    );
    let error = result.unwrap_err();
    assert!(
        error.to_string().contains("input changed during analysis"),
        "{error}"
    );
    assert!(!cache.join("dataset.manifest.json").exists());
    assert!(!cache.join("c00000000-p0000000000").exists());
}
