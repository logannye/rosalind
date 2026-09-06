use rosalind::evidence::*;
use rosalind::selection::GenomicInterval;
use rust_htslib::bam::record::{Aux, Cigar, CigarString, Record};
use rust_htslib::bam::{self, header::HeaderRecord};
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

static NEXT_FIXTURE: AtomicU64 = AtomicU64::new(0);
const LENGTH: u32 = CANONICAL_TILE_BASES + 32;

struct Fixture {
    root: PathBuf,
    bam: PathBuf,
    fasta: PathBuf,
}

impl Fixture {
    fn new(groups: &[(&str, Option<&str>)], records: &[Record]) -> Self {
        let root = std::env::temp_dir().join(format!(
            "rosalind-evidence-samples-{}-{}",
            std::process::id(),
            NEXT_FIXTURE.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir(&root).unwrap();
        let bam = root.join("reads.bam");
        let fasta = root.join("reference.fa");
        std::fs::write(&fasta, format!(">chr1\n{}\n", "A".repeat(LENGTH as usize))).unwrap();
        std::fs::write(
            root.join("reference.fa.fai"),
            format!("chr1\t{LENGTH}\t6\t{LENGTH}\t{}\n", LENGTH + 1),
        )
        .unwrap();
        let mut header = bam::Header::new();
        header.push_record(
            HeaderRecord::new(b"HD")
                .push_tag(b"VN", "1.6")
                .push_tag(b"SO", "coordinate"),
        );
        header.push_record(
            HeaderRecord::new(b"SQ")
                .push_tag(b"SN", "chr1")
                .push_tag(b"LN", LENGTH),
        );
        for (id, sample) in groups {
            let mut group = HeaderRecord::new(b"RG");
            group.push_tag(b"ID", id);
            if let Some(sample) = sample {
                group.push_tag(b"SM", sample);
            }
            header.push_record(&group);
        }
        let mut writer = bam::Writer::from_path(&bam, &header, bam::Format::Bam).unwrap();
        for record in records {
            writer.write(record).unwrap();
        }
        drop(writer);
        bam::index::build(&bam, None::<&PathBuf>, bam::index::Type::Bai, 1).unwrap();
        Self { root, bam, fasta }
    }

    fn request(&self, sample_selection: EvidenceSampleSelection) -> EvidenceRequest {
        let mut request = EvidenceRequest::new(&self.bam, &self.fasta);
        request.sample_selection = sample_selection;
        request.selection = interval(0, 4);
        request
    }
}

impl Drop for Fixture {
    fn drop(&mut self) {
        std::fs::remove_dir_all(&self.root).ok();
    }
}

fn interval(start: u32, end: u32) -> EvidenceSelection {
    EvidenceSelection::Intervals(vec![GenomicInterval {
        contig: 0,
        start,
        end,
    }])
}

fn read(name: &str, group: Option<&str>, base: u8, start: u32, length: usize) -> Record {
    let mut record = Record::new();
    record.set(
        name.as_bytes(),
        Some(&CigarString(vec![Cigar::Match(length as u32)])),
        &vec![base; length],
        &vec![30; length],
    );
    record.set_tid(0);
    record.set_pos(start as i64);
    record.set_flags(0);
    record.set_mapq(60);
    if let Some(group) = group {
        record.push_aux(b"RG", Aux::String(group)).unwrap();
    }
    record
}

#[derive(Default)]
struct Capture(Vec<EvidenceRow>);
impl EvidenceAnalyzer for Capture {
    fn on_batch(&mut self, batch: &EvidenceBatch) -> Result<(), EvidenceError> {
        for row in batch.rows() {
            self.0.push(row.try_to_full_row()?);
        }
        Ok(())
    }
    fn additional_memory_bytes(&self) -> Option<u64> {
        // Fixtures select at most 12 loci, retaining their full rows.
        Some(64 * std::mem::size_of::<EvidenceRow>() as u64)
    }
}

fn collect(request: EvidenceRequest) -> (Vec<EvidenceRow>, EvidenceRunStats, EvidenceSampleScope) {
    let mut engine = EvidenceEngine::open(request).unwrap();
    let mut capture = Capture::default();
    let stats = engine.run(&mut capture).unwrap();
    (capture.0, stats, engine.sample_scope().clone())
}

#[test]
fn auto_unifies_read_groups_of_one_sample_and_matches_explicit_selection() {
    let fixture = Fixture::new(
        &[("lane2", Some("sample-a")), ("lane1", Some("sample-a"))],
        &[
            read("one", Some("lane1"), b'A', 0, 4),
            read("two", Some("lane2"), b'C', 0, 4),
        ],
    );
    let (automatic, _, scope) = collect(fixture.request(EvidenceSampleSelection::Auto));
    let (explicit, _, selected_scope) =
        collect(fixture.request(EvidenceSampleSelection::Named("sample-a".into())));
    assert_eq!(automatic, explicit);
    assert_eq!(scope, selected_scope);
    assert_eq!(scope.mode, EvidenceSampleMode::Named);
    assert_eq!(scope.selected_sample.as_deref(), Some("sample-a"));
    assert_eq!(scope.declared_samples, ["sample-a"]);
    assert_eq!(scope.read_groups[0].id, "lane1");
    assert!(!scope.allows_unassigned);
    assert_eq!(automatic[0].allele_counts, [1, 1, 0, 0]);
    assert_eq!(scope.canonical_json(), "{\"version\":1,\"mode\":\"named\",\"selected_sample\":\"sample-a\",\"declared_samples\":[\"sample-a\"],\"read_groups\":[{\"id\":\"lane1\",\"sample\":\"sample-a\"},{\"id\":\"lane2\",\"sample\":\"sample-a\"}],\"allows_unassigned\":false}");
}

#[test]
fn auto_refuses_multiple_samples_and_mixed_named_unnamed_groups() {
    for groups in [
        vec![("one", Some("sample-a")), ("two", Some("sample-b"))],
        vec![("one", Some("sample-a")), ("two", None)],
    ] {
        let fixture = Fixture::new(&groups, &[]);
        let error =
            EvidenceEngine::open(fixture.request(EvidenceSampleSelection::Auto)).unwrap_err();
        assert!(matches!(error, EvidenceError::InvalidRequest(_)));
        let message = error.to_string();
        assert!(message.contains("--sample"), "{message}");
        assert!(message.contains("--pool-samples"), "{message}");
    }
}

#[test]
fn named_selection_precedes_every_locus_depth_and_filter_counter() {
    let mut other_sample = read("other", Some("two"), b'C', 0, 4);
    other_sample.set_mapq(0);
    other_sample.set_flags(0x400);
    let fixture = Fixture::new(
        &[("one", Some("sample-a")), ("two", Some("sample-b"))],
        &[read("selected", Some("one"), b'A', 0, 4), other_sample],
    );
    let (rows, stats, scope) =
        collect(fixture.request(EvidenceSampleSelection::Named("sample-a".into())));
    assert_eq!(stats.record_visits, 2);
    assert_eq!(stats.sample_filtered_record_visits, 1);
    assert_eq!(stats.filtered_record_visits, 0);
    assert_eq!(scope.declared_samples, ["sample-a", "sample-b"]);
    for row in rows {
        assert_eq!(row.prefilter_depth, 1);
        assert_eq!(row.aligned_depth, 1);
        assert_eq!(row.callable_depth, 1);
        assert_eq!(row.allele_counts, [1, 0, 0, 0]);
        assert_eq!(row.filters, EvidenceFilterCounts::default());
    }
}

#[test]
fn named_runs_refuse_missing_undeclared_unnamed_or_nonstring_rg() {
    let mut wrong_type = read("wrong-type", None, b'A', 0, 4);
    wrong_type.push_aux(b"RG", Aux::I32(7)).unwrap();
    let cases = [
        read("missing", None, b'A', 0, 4),
        read("unknown", Some("absent"), b'A', 0, 4),
        read("unnamed", Some("unnamed"), b'A', 0, 4),
        wrong_type,
    ];
    for record in cases {
        let fixture = Fixture::new(&[("named", Some("sample-a")), ("unnamed", None)], &[record]);
        let mut engine = EvidenceEngine::open(
            fixture.request(EvidenceSampleSelection::Named("sample-a".into())),
        )
        .unwrap();
        let mut capture = Capture::default();
        let error = engine.run(&mut capture).unwrap_err();
        assert!(matches!(error, EvidenceError::InvalidInput(_)));
        assert!(error.to_string().contains("cannot be assigned"));
        assert!(
            capture.0.is_empty(),
            "invalid tile must not reach consumers"
        );
    }
    let fixture = Fixture::new(
        &[("named", Some("sample-a"))],
        &[read("missing", None, b'A', 0, 4)],
    );
    let mut engine = EvidenceEngine::open(fixture.request(EvidenceSampleSelection::Auto)).unwrap();
    assert!(engine.run(&mut Capture::default()).is_err());
}

#[test]
fn named_selection_must_exist_in_the_header() {
    for groups in [vec![], vec![("one", Some("sample-a"))]] {
        let fixture = Fixture::new(&groups, &[]);
        let error =
            EvidenceEngine::open(fixture.request(EvidenceSampleSelection::Named("absent".into())))
                .unwrap_err();
        assert!(error.to_string().contains("not declared"));
    }
}

#[test]
fn deliberate_pooling_includes_named_and_unassigned_records() {
    let fixture = Fixture::new(
        &[
            ("one", Some("sample-a")),
            ("two", Some("sample-b")),
            ("unnamed", None),
        ],
        &[
            read("one", Some("one"), b'A', 0, 4),
            read("two", Some("two"), b'C', 0, 4),
            read("unnamed", Some("unnamed"), b'G', 0, 4),
            read("missing", None, b'T', 0, 4),
            read("unknown", Some("absent"), b'T', 0, 4),
        ],
    );
    let (rows, stats, scope) = collect(fixture.request(EvidenceSampleSelection::Pool));
    assert_eq!(scope.mode, EvidenceSampleMode::Pooled);
    assert_eq!(scope.selected_sample, None);
    assert_eq!(scope.declared_samples, ["sample-a", "sample-b"]);
    assert!(scope.allows_unassigned);
    assert_eq!(stats.sample_filtered_record_visits, 0);
    assert_eq!(rows[0].allele_counts, [1, 1, 1, 2]);
    assert_eq!(rows[0].prefilter_depth, 5);
}

#[test]
fn unlabeled_inputs_keep_explicit_unknown_identity_and_existing_counts() {
    for groups in [vec![], vec![("unnamed", None)]] {
        let fixture = Fixture::new(&groups, &[read("missing", None, b'A', 0, 4)]);
        let (rows, _, scope) = collect(fixture.request(EvidenceSampleSelection::Auto));
        assert_eq!(scope.mode, EvidenceSampleMode::Unknown);
        assert!(scope.allows_unassigned);
        assert!(scope.declared_samples.is_empty());
        assert_eq!(rows[0].callable_depth, 1);
        let (_, _, pooled) = collect(fixture.request(EvidenceSampleSelection::Pool));
        assert_ne!(scope.canonical_json(), pooled.canonical_json());
    }
}

#[test]
fn sample_scope_and_rows_are_invariant_across_microtiles_workers_and_boundaries() {
    let start = CANONICAL_TILE_BASES - 4;
    let end = CANONICAL_TILE_BASES + 8;
    let fixture = Fixture::new(
        &[("one", Some("sample-a")), ("two", Some("sample-b"))],
        &[
            read("selected", Some("one"), b'A', start, 12),
            read("other", Some("two"), b'C', start, 12),
        ],
    );
    let mut request = fixture.request(EvidenceSampleSelection::Named("sample-a".into()));
    request.selection = interval(start, end);
    request.execution.memory_budget_bytes = Some(512 << 20);
    let (expected, _, scope) = collect(request.clone());
    for width in [1, 3, CANONICAL_TILE_BASES] {
        request.execution.max_microtile_bases = width;
        let (actual, _, actual_scope) = collect(request.clone());
        assert_eq!(actual, expected);
        assert_eq!(actual_scope.canonical_json(), scope.canonical_json());
    }
    let engine = EvidenceEngine::open(request).unwrap();
    let factory = engine.worker_factory();
    let handles: Vec<_> = [(start, CANONICAL_TILE_BASES), (CANONICAL_TILE_BASES, end)]
        .into_iter()
        .map(|(begin, finish)| {
            let factory = factory.clone();
            let scope = scope.clone();
            std::thread::spawn(move || {
                let mut worker = factory.open(interval(begin, finish)).unwrap();
                assert_eq!(worker.sample_scope(), &scope);
                assert!(worker.plan().sample_scope_bytes > 0);
                let mut capture = Capture::default();
                worker.run(&mut capture).unwrap();
                capture.0
            })
        })
        .collect();
    let actual: Vec<_> = handles
        .into_iter()
        .flat_map(|handle| handle.join().unwrap())
        .collect();
    assert_eq!(actual, expected);
}

#[test]
fn canonical_scope_ignores_header_order_and_escapes_identifiers() {
    let groups = [("z", Some("a\"b\\c")), ("a", Some("a\"b\\c"))];
    let first = Fixture::new(&groups, &[]);
    let second = Fixture::new(&[groups[1], groups[0]], &[]);
    let a = EvidenceEngine::open(first.request(EvidenceSampleSelection::Auto)).unwrap();
    let b = EvidenceEngine::open(second.request(EvidenceSampleSelection::Auto)).unwrap();
    assert_eq!(
        a.sample_scope().canonical_json(),
        b.sample_scope().canonical_json()
    );
    assert!(a.sample_scope().canonical_json().contains("a\\\"b\\\\c"));
}

#[test]
fn sample_header_state_is_reserved_through_replanning() {
    let names: Vec<_> = (0..200).map(|n| format!("read-group-{n:04}")).collect();
    let groups: Vec<_> = names
        .iter()
        .map(|name| (name.as_str(), Some("sample-a")))
        .collect();
    let fixture = Fixture::new(&groups, &[]);
    let mut engine = EvidenceEngine::open(fixture.request(EvidenceSampleSelection::Auto)).unwrap();
    let reserved = engine.plan().sample_scope_bytes;
    assert!(reserved > 200 * 256);
    assert!(engine.plan().fixed_bytes >= reserved);
    engine.set_selection(interval(1, 3)).unwrap();
    assert_eq!(engine.plan().sample_scope_bytes, reserved);
    engine.plan_for_analyzer(&Capture::default()).unwrap();
    assert_eq!(engine.plan().sample_scope_bytes, reserved);
}

#[test]
fn duplicate_read_group_ids_are_rejected_even_for_explicit_pooling() {
    let fixture = Fixture::new(
        &[("same", Some("sample-a")), ("same", Some("sample-b"))],
        &[],
    );
    let error = EvidenceEngine::open(fixture.request(EvidenceSampleSelection::Pool)).unwrap_err();
    assert!(
        error.to_string().contains("declared more than once"),
        "{error}"
    );
}
