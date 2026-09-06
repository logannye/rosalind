use rosalind::evidence::*;
use rosalind::selection::GenomicInterval;
use rust_htslib::bam::{
    self,
    record::{Aux, Cigar, CigarString},
    Read,
};
use std::{
    path::PathBuf,
    sync::atomic::{AtomicU64, Ordering},
};
static NEXT: AtomicU64 = AtomicU64::new(0);
struct Fixture {
    root: PathBuf,
    fasta: PathBuf,
    bam: PathBuf,
    cram: PathBuf,
}
impl Fixture {
    fn new(records: usize, oversized_aux: bool) -> Self {
        let root = std::env::temp_dir().join(format!(
            "rosalind-cram-envelope-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir_all(&root).unwrap();
        let fasta = root.join("ref.fa");
        std::fs::write(&fasta, format!(">chr1\n{}\n", "A".repeat(2000))).unwrap();
        std::fs::write(root.join("ref.fa.fai"), "chr1\t2000\t6\t2000\t2001\n").unwrap();
        let bam = root.join("reads.bam");
        let cram = root.join("reads.cram");
        let mut header = bam::Header::new();
        header.push_record(
            bam::header::HeaderRecord::new(b"HD")
                .push_tag(b"VN", "1.6")
                .push_tag(b"SO", "coordinate"),
        );
        header.push_record(
            bam::header::HeaderRecord::new(b"SQ")
                .push_tag(b"SN", "chr1")
                .push_tag(b"LN", 2000),
        );
        let mut writer = bam::Writer::from_path(&bam, &header, bam::Format::Bam).unwrap();
        for i in 0..records {
            let mut record = bam::Record::new();
            record.set(
                format!("r{i}").as_bytes(),
                Some(&CigarString(vec![Cigar::Match(100)])),
                &[b'A'; 100],
                &[30; 100],
            );
            record.set_tid(0);
            record.set_pos(500);
            record.set_mapq(60);
            record.set_flags(0);
            if oversized_aux && i + 1 == records {
                record
                    .push_aux(b"ZZ", Aux::String(&"x".repeat(4096)))
                    .unwrap();
            }
            writer.write(&record).unwrap();
        }
        drop(writer);
        bam::index::build(&bam, None::<&PathBuf>, bam::index::Type::Bai, 1).unwrap();
        let mut writer = bam::Writer::from_path(&cram, &header, bam::Format::Cram).unwrap();
        writer.set_reference(&fasta).unwrap();
        let mut reader = bam::Reader::from_path(&bam).unwrap();
        for record in reader.records() {
            writer.write(&record.unwrap()).unwrap();
        }
        drop(writer);
        bam::index::build(&cram, None::<&PathBuf>, bam::index::Type::Bai, 1).unwrap();
        Self {
            root,
            fasta,
            bam,
            cram,
        }
    }
    fn request(&self) -> EvidenceRequest {
        let mut request = EvidenceRequest::new(&self.cram, &self.fasta);
        request.selection = EvidenceSelection::Intervals(vec![GenomicInterval {
            contig: 0,
            start: 500,
            end: 510,
        }]);
        request
    }
}
impl Drop for Fixture {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.root);
    }
}
fn tsv(request: EvidenceRequest) -> Vec<u8> {
    let mut engine = EvidenceEngine::open(request).unwrap();
    let mut bytes = Vec::new();
    engine.run(&mut EvidenceTsvWriter::new(&mut bytes)).unwrap();
    bytes
}
#[test]
fn whole_file_maxima_are_measured_once_and_scientific_bytes_match_bam() {
    let fixture = Fixture::new(200, false);
    let request = fixture.request();
    let engine = EvidenceEngine::open(request.clone()).unwrap();
    let envelope = engine.plan().cram_envelope.as_ref().unwrap();
    assert_eq!(envelope.validated_records, 200);
    assert_eq!(envelope.validated_bases, 20_000);
    assert_eq!(envelope.validated_max_read_length, 100);
    assert!(envelope.validated_max_record_bytes < 512);
    assert!(engine.plan().decoder_bytes < 32 << 20);
    assert_eq!(
        engine.plan().decoder_measurements()["execution.cram.validated_records"],
        "200"
    );
    let worker = engine
        .worker_factory()
        .open(request.selection.clone())
        .unwrap();
    assert_eq!(
        worker
            .plan()
            .cram_envelope
            .as_ref()
            .unwrap()
            .validation_wall_micros,
        envelope.validation_wall_micros
    );
    let mut bam_request = request.clone();
    bam_request.alignments = fixture.bam.clone();
    let expected = tsv(bam_request);
    for width in [1, 7, 1024] {
        let mut request = request.clone();
        request.execution.max_microtile_bases = width;
        assert_eq!(tsv(request), expected);
    }
    let mut projected = request.clone();
    projected.fields = EvidenceFields::DEPTHS;
    assert_eq!(
        EvidenceEngine::open(projected)
            .unwrap()
            .plan()
            .decoder_bytes,
        engine.plan().decoder_bytes
    );
}
#[test]
fn off_target_payload_is_checked_before_successful_indexed_analysis() {
    let fixture = Fixture::new(2, true);
    let mut request = fixture.request();
    request.execution.max_record_bytes = 1024;
    request.selection = EvidenceSelection::Intervals(vec![GenomicInterval {
        contig: 0,
        start: 0,
        end: 10,
    }]);
    let mut bam_request = request.clone();
    bam_request.alignments = fixture.bam.clone();
    assert!(!tsv(bam_request).is_empty());
    let error = EvidenceEngine::open(request).unwrap_err();
    assert!(matches!(error, EvidenceError::RecordLimit(_)), "{error}");
    assert!(error.to_string().contains("off-target"));
}
#[test]
fn worker_cannot_tighten_limits_below_validated_records_or_ignore_index_mutation() {
    let fixture = Fixture::new(3, false);
    let request = fixture.request();
    let engine = EvidenceEngine::open(request.clone()).unwrap();
    let factory = engine.worker_factory();
    let mut execution = request.execution.clone();
    execution.max_record_bytes = 1;
    assert!(matches!(
        factory.open_with_execution(request.selection.clone(), execution),
        Err(EvidenceError::RecordLimit(_))
    ));
    let mut execution = request.execution.clone();
    execution.memory_budget_bytes = Some(1);
    assert!(matches!(
        factory.open_with_execution(request.selection.clone(), execution),
        Err(EvidenceError::Refused { .. })
    ));
    let index = PathBuf::from(format!("{}.crai", fixture.cram.display()));
    use std::io::Write;
    std::fs::OpenOptions::new()
        .append(true)
        .open(index)
        .unwrap()
        .write_all(b"changed")
        .unwrap();
    assert!(factory.open(request.selection).is_err());
}

#[test]
fn omitted_index_slice_is_rejected_instead_of_becoming_false_zero_evidence() {
    let fixture = Fixture::new(10_001, false);
    let path = PathBuf::from(format!("{}.crai", fixture.cram.display()));
    use std::io::{Read as _, Write};
    let mut text = String::new();
    flate2::read::MultiGzDecoder::new(std::fs::File::open(&path).unwrap())
        .read_to_string(&mut text)
        .unwrap();
    assert!(
        text.lines().count() >= 2,
        "fixture must contain multiple indexed slices"
    );
    let mut writer = flate2::write::GzEncoder::new(
        std::fs::File::create(path).unwrap(),
        flate2::Compression::default(),
    );
    writeln!(writer, "{}", text.lines().next().unwrap()).unwrap();
    writer.finish().unwrap();
    let error = EvidenceEngine::open(fixture.request()).unwrap_err();
    assert!(
        error.to_string().contains("CRAI entry disagrees"),
        "{error}"
    );
}
