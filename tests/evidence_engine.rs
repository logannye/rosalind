use rosalind::evidence::*;
use rosalind::selection::GenomicInterval;
use rust_htslib::bam::record::{Cigar, CigarString, Record};
use rust_htslib::bam::{self, Read};
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{SystemTime, UNIX_EPOCH};
static FIXTURE_ID: AtomicU64 = AtomicU64::new(0);

struct Fixture {
    root: PathBuf,
    fasta: PathBuf,
    bam: PathBuf,
    length: u32,
}
impl Fixture {
    fn new(length: u32) -> Self {
        let root = std::env::temp_dir().join(format!(
            "rosalind-exact-evidence-{}-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos(),
            FIXTURE_ID.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir_all(&root).unwrap();
        let fasta = root.join("reference.fa");
        let bam = root.join("reads.bam");
        std::fs::write(&fasta, format!(">chr1\n{}\n", "A".repeat(length as usize))).unwrap();
        std::fs::write(
            root.join("reference.fa.fai"),
            format!("chr1\t{length}\t6\t{length}\t{}\n", length + 1),
        )
        .unwrap();
        Self {
            root,
            fasta,
            bam,
            length,
        }
    }
    fn header(&self) -> bam::Header {
        let mut header = bam::Header::new();
        let mut hd = bam::header::HeaderRecord::new(b"HD");
        hd.push_tag(b"VN", "1.6");
        hd.push_tag(b"SO", "coordinate");
        header.push_record(&hd);
        let mut sq = bam::header::HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr1");
        sq.push_tag(b"LN", self.length);
        header.push_record(&sq);
        header
    }
    fn write(&self, records: Vec<Record>, index: bam::index::Type) {
        let mut writer =
            bam::Writer::from_path(&self.bam, &self.header(), bam::Format::Bam).unwrap();
        for record in records {
            writer.write(&record).unwrap();
        }
        drop(writer);
        bam::index::build(&self.bam, None::<&PathBuf>, index, 1).unwrap();
    }
    fn request(&self, width: u32) -> EvidenceRequest {
        let mut request = EvidenceRequest::new(&self.bam, &self.fasta);
        request.execution.max_microtile_bases = width;
        request
    }
}
impl Drop for Fixture {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.root);
    }
}
fn record(
    name: &str,
    pos: i64,
    seq: &[u8],
    qual: &[u8],
    mapq: u8,
    flags: u16,
    cigar: Vec<Cigar>,
) -> Record {
    let mut record = Record::new();
    record.set(name.as_bytes(), Some(&CigarString(cigar)), seq, qual);
    record.set_tid(0);
    record.set_pos(pos);
    record.set_mapq(mapq);
    record.set_flags(flags);
    record
}
fn matched(name: &str, pos: i64, len: usize, mapq: u8, flags: u16) -> Record {
    record(
        name,
        pos,
        &vec![b'A'; len],
        &vec![30; len],
        mapq,
        flags,
        vec![Cigar::Match(len as u32)],
    )
}
#[derive(Default)]
struct Capture(Vec<EvidenceRow>);
impl EvidenceAnalyzer for Capture {
    fn on_batch(&mut self, batch: &EvidenceBatch) -> Result<(), EvidenceError> {
        self.0.extend(batch.rows.clone());
        Ok(())
    }
}
fn collect(request: EvidenceRequest) -> Vec<EvidenceRow> {
    let mut engine = EvidenceEngine::open(request).unwrap();
    let mut capture = Capture::default();
    engine.run(&mut capture).unwrap();
    capture.0
}

#[test]
fn exact_depth_exceeds_old_cap_and_tile_width_does_not_change_rows() {
    let f = Fixture::new(160);
    f.write(
        (0..1101)
            .map(|n| matched(&format!("r{n}"), 0, 150, 60, 0))
            .collect(),
        bam::index::Type::Bai,
    );
    let small = collect(f.request(7));
    let large = collect(f.request(16384));
    assert_eq!(small, large);
    assert_eq!(small[0].callable_depth, 1101);
    assert_eq!(small[149].callable_depth, 1101);
    assert_eq!(small[150].callable_depth, 0);
    assert_eq!(small.len(), 160);
    assert_eq!(small[149].base_quality_histogram[30], 1101);
    assert_eq!(small[149].read_position_sum, 149 * 1101);
}

#[test]
fn independent_hand_calculated_cigar_and_filter_oracle() {
    let f = Fixture::new(16);
    let reads = vec![
        record(
            "cigar",
            5,
            b"AACCGGTT",
            &[30; 8],
            60,
            0,
            vec![
                Cigar::SoftClip(2),
                Cigar::Match(3),
                Cigar::Ins(1),
                Cigar::Match(2),
                Cigar::Del(1),
            ],
        ),
        record(
            "base-filters",
            7,
            b"AGNT",
            &[25, 5, 30, 255],
            60,
            16,
            vec![Cigar::Match(4)],
        ),
        matched("low-mapq", 7, 4, 5, 0),
        matched("duplicate", 7, 4, 60, 0x400),
        matched("secondary-first", 7, 1, 60, 0x900),
        matched("qc-fail", 7, 1, 60, 0x200),
        matched("missing-mapq", 7, 1, 255, 0),
    ];
    f.write(reads, bam::index::Type::Bai);
    let rows = collect(f.request(1));
    assert_eq!(rows, collect(f.request(13)));
    let row = &rows[7];
    assert_eq!(row.prefilter_depth, 7);
    assert_eq!(row.aligned_depth, 2);
    assert_eq!(row.callable_depth, 2);
    assert_eq!(row.allele_counts, [1, 0, 1, 0]);
    assert_eq!(row.strand_counts[0], [0, 1]);
    assert_eq!(row.strand_counts[2], [1, 0]);
    assert_eq!(row.filters.secondary, 1);
    assert_eq!(row.filters.supplementary, 0);
    assert_eq!(row.filters.qc_fail, 1);
    assert_eq!(row.filters.duplicate, 1);
    assert_eq!(row.filters.low_mapq, 1);
    assert_eq!(row.filters.unavailable_mapq, 1);
    assert_eq!(row.base_quality_sum, 55);
    assert_eq!(row.mapping_quality_sum, 120);
    assert_eq!(row.read_position_sum, 7);
    assert_eq!(row.read_length_sum, 12);
    assert_eq!(rows[8].filters.low_base_quality, 1);
    assert_eq!(rows[9].filters.ambiguous_base, 1);
    assert_eq!(rows[10].filters.unavailable_base_quality, 1);
    assert_eq!(rows[10].aligned_depth, 1);
    assert_eq!(rows[10].callable_depth, 0);
    assert_eq!(rows[11].prefilter_depth, 0);
}

#[test]
fn snv_sites_validate_reference_and_merge_duplicate_alternates() {
    let f = Fixture::new(20);
    f.write(vec![matched("r", 0, 10, 60, 0)], bam::index::Type::Bai);
    let vcf = f.root.join("sites.vcf");
    std::fs::write(&vcf,"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\nchr1\t2\t.\tA\tT,C\nchr1\t2\t.\tA\tG\nchr1\t18\t.\tA\tC\n").unwrap();
    let mut engine = EvidenceEngine::open(f.request(3)).unwrap();
    let selection = EvidenceSelection::from_vcf(&vcf, engine.contigs()).unwrap();
    engine.set_selection(selection).unwrap();
    assert_eq!(engine.plan().selected_loci, 2);
    let mut capture = Capture::default();
    engine.run(&mut capture).unwrap();
    assert_eq!(capture.0.len(), 2);
    assert_eq!(capture.0[0].requested_alts, b"CGT");
    assert_eq!(capture.0[1].callable_depth, 0);
    let bad = EvidenceSelection::Sites(vec![SnvSite {
        contig: 0,
        position: 1,
        reference: b'C',
        alternates: vec![b'T'],
    }]);
    assert!(engine.set_selection(bad).is_err());
}

#[test]
fn arrow_bytes_are_canonical_and_reader_roundtrips_all_fields() {
    let f = Fixture::new(1100);
    f.write(vec![matched("r", 0, 150, 60, 0)], bam::index::Type::Bai);
    let encode = |width| {
        let mut engine = EvidenceEngine::open(f.request(width)).unwrap();
        let mut writer = EvidenceArrowWriter::new(Vec::new());
        engine.run(&mut writer).unwrap();
        writer.into_inner().unwrap()
    };
    let a = encode(7);
    let b = encode(1100);
    assert_eq!(a, b);
    let engine = EvidenceEngine::open(f.request(100)).unwrap();
    let mut rows = Vec::new();
    read_evidence_batches(a.as_slice(), engine.contigs(), |batch| {
        assert!(batch.rows.len() <= 1024);
        rows.extend(batch.rows.clone());
        Ok(())
    })
    .unwrap();
    assert_eq!(rows, collect(f.request(100)));
}

#[test]
fn overlapping_panel_targets_keep_full_length_denominators_and_fuse() {
    let f = Fixture::new(20);
    f.write(vec![matched("r", 2, 4, 60, 0)], bam::index::Type::Bai);
    let mut panel = PanelQcAnalyzer::new(vec![
        PanelTarget {
            id: "a".into(),
            contig: 0,
            start: 0,
            end: 10,
        },
        PanelTarget {
            id: "b".into(),
            contig: 0,
            start: 4,
            end: 12,
        },
    ])
    .unwrap()
    .with_min_callable_depth(1);
    let mut request = f.request(3);
    request.selection = panel.selection();
    let mut engine = EvidenceEngine::open(request).unwrap();
    let mut writer = EvidenceTsvWriter::new(Vec::new());
    {
        let mut fused = FusedAnalyzers::new(vec![&mut panel, &mut writer]);
        engine.run(&mut fused).unwrap();
    }
    assert_eq!(panel.summaries()[0].callable_depth_sum, 4);
    assert_eq!(panel.summaries()[0].target_length(), 10);
    assert_eq!(panel.summaries()[0].callable_positions, 4);
    assert_eq!(panel.summaries()[1].callable_depth_sum, 2);
    assert_eq!(panel.summaries()[1].target_length(), 8);
    assert_eq!(panel.summaries()[1].base_quality_sum, 60);
    let text = String::from_utf8(writer.into_inner()).unwrap();
    assert_eq!(text.lines().count(), 13);
}

#[test]
fn csi_and_relocated_fai_index_inputs_work_without_adjacent_sidecars() {
    let f = Fixture::new(30);
    f.write(vec![matched("r", 0, 10, 60, 0)], bam::index::Type::Csi(14));
    let moved_index = f.root.join("index-by-hash");
    let moved_fai = f.root.join("fai-by-hash");
    std::fs::rename(f.root.join("reads.bam.csi"), &moved_index).unwrap();
    std::fs::rename(f.root.join("reference.fa.fai"), &moved_fai).unwrap();
    let mut request = f.request(11);
    request.alignment_index = Some(moved_index);
    request.reference_fai = Some(moved_fai);
    let rows = collect(request);
    assert_eq!(rows[0].callable_depth, 1);
    assert_eq!(rows[20].callable_depth, 0);
}

#[test]
fn coverage_without_reference_and_tiny_budget_refusal_are_explicit() {
    let f = Fixture::new(20);
    f.write(vec![matched("r", 0, 10, 60, 0)], bam::index::Type::Bai);
    let request = EvidenceRequest::coverage(&f.bam);
    let rows = collect(request);
    assert_eq!(rows[0].reference, b'N');
    assert_eq!(rows[0].callable_depth, 1);
    let mut request = f.request(4);
    request.execution.memory_budget_bytes = Some(1);
    assert!(matches!(
        EvidenceEngine::open(request),
        Err(EvidenceError::Refused { .. })
    ));
}

#[test]
fn worker_factory_shares_open_reference_without_reopening_its_path() {
    let f = Fixture::new(30);
    f.write(vec![matched("r", 0, 10, 60, 0)], bam::index::Type::Bai);
    let engine = EvidenceEngine::open(f.request(5)).unwrap();
    let factory = engine.worker_factory();
    std::fs::rename(&f.fasta, f.root.join("reference-moved.fa")).unwrap();
    let handles: Vec<_> = (0..2)
        .map(|i| {
            let factory = factory.clone();
            std::thread::spawn(move || {
                let mut engine = factory
                    .open(EvidenceSelection::Intervals(vec![GenomicInterval {
                        contig: 0,
                        start: i * 15,
                        end: (i + 1) * 15,
                    }]))
                    .unwrap();
                let mut capture = Capture::default();
                engine.run(&mut capture).unwrap();
                capture.0
            })
        })
        .collect();
    let rows: Vec<_> = handles
        .into_iter()
        .flat_map(|h| h.join().unwrap())
        .collect();
    assert_eq!(rows.len(), 30);
    assert_eq!(rows[0].callable_depth, 1);
    assert_eq!(rows[29].callable_depth, 0);
}

#[test]
fn cram_uses_explicit_local_fasta_and_matches_bam_evidence() {
    let f = Fixture::new(100);
    f.write(vec![matched("r", 5, 20, 60, 0)], bam::index::Type::Bai);
    let expected = collect(f.request(7));
    let cram = f.root.join("reads.cram");
    let mut writer = bam::Writer::from_path(&cram, &f.header(), bam::Format::Cram).unwrap();
    writer.set_reference(&f.fasta).unwrap();
    let mut source = bam::Reader::from_path(&f.bam).unwrap();
    for record in source.records() {
        writer.write(&record.unwrap()).unwrap();
    }
    drop(writer);
    bam::index::build(&cram, None::<&PathBuf>, bam::index::Type::Bai, 1).unwrap();
    let mut request = EvidenceRequest::new(&cram, &f.fasta);
    request.execution.max_microtile_bases = 7;
    assert_eq!(collect(request), expected);
    let request = EvidenceRequest::coverage(&cram);
    assert!(matches!(
        EvidenceEngine::open(request),
        Err(EvidenceError::InvalidRequest(_))
    ));
    struct RequiresReference(usize);
    impl EvidenceAnalyzer for RequiresReference {
        fn requirements(&self) -> EvidenceRequirements {
            EvidenceRequirements {
                fields: EvidenceFields::DEPTHS,
                requires_reference: true,
                context_bases: 0,
                retained_bytes: Some(0),
            }
        }
        fn on_batch(&mut self, batch: &EvidenceBatch) -> Result<(), EvidenceError> {
            assert!(batch.rows.iter().all(|row| row.reference == b'A'));
            self.0 += batch.rows.len();
            Ok(())
        }
    }
    let mut coverage_bam = EvidenceEngine::open(EvidenceRequest::coverage(&f.bam)).unwrap();
    let mut analyzer = RequiresReference(0);
    assert!(matches!(
        coverage_bam.run(&mut analyzer),
        Err(EvidenceError::InvalidRequest(_))
    ));
    assert_eq!(analyzer.0, 0);
    let mut cram_request = EvidenceRequest::coverage(&cram);
    cram_request.cram_reference = Some(f.fasta.clone());
    let mut coverage_cram = EvidenceEngine::open(cram_request).unwrap();
    coverage_cram.run(&mut analyzer).unwrap();
    assert_eq!(analyzer.0, 100);
}

#[test]
fn malformed_cigar_and_declared_read_envelope_fail_instead_of_truncating() {
    let f = Fixture::new(300);
    f.write(vec![matched("r", 0, 251, 60, 0)], bam::index::Type::Bai);
    let mut engine = EvidenceEngine::open(f.request(8)).unwrap();
    let mut sink = EvidenceTsvWriter::new(std::io::sink());
    assert!(matches!(
        engine.run(&mut sink),
        Err(EvidenceError::RecordLimit(_))
    ));
}

#[test]
fn analyzer_capabilities_and_unknown_budget_bounds_are_validated_before_callbacks() {
    struct NeedsContext;
    impl EvidenceAnalyzer for NeedsContext {
        fn requirements(&self) -> EvidenceRequirements {
            EvidenceRequirements {
                fields: EvidenceFields::DEPTHS,
                requires_reference: false,
                context_bases: 1,
                retained_bytes: Some(0),
            }
        }
        fn on_batch(&mut self, _: &EvidenceBatch) -> Result<(), EvidenceError> {
            panic!("unsupported analyzer must never execute")
        }
    }
    let f = Fixture::new(20);
    f.write(vec![matched("r", 0, 10, 60, 0)], bam::index::Type::Bai);
    let mut engine = EvidenceEngine::open(f.request(4)).unwrap();
    assert!(engine.run(&mut NeedsContext).is_err());
    let mut projected = f.request(4);
    projected.fields = EvidenceFields::DEPTHS;
    assert!(EvidenceEngine::open(projected).is_err());
    let mut budgeted = f.request(4);
    budgeted.execution.memory_budget_bytes = Some(1 << 30);
    let mut engine = EvidenceEngine::open(budgeted).unwrap();
    assert!(engine.run(&mut Capture::default()).is_err());
}

#[test]
fn selection_identity_is_independent_of_input_interval_order_and_duplicates() {
    let f = Fixture::new(20);
    f.write(vec![matched("r", 0, 10, 60, 0)], bam::index::Type::Bai);
    let mut engine = EvidenceEngine::open(f.request(4)).unwrap();
    engine
        .set_selection(EvidenceSelection::Intervals(vec![
            GenomicInterval {
                contig: 0,
                start: 4,
                end: 10,
            },
            GenomicInterval {
                contig: 0,
                start: 0,
                end: 5,
            },
            GenomicInterval {
                contig: 0,
                start: 0,
                end: 5,
            },
        ]))
        .unwrap();
    let first = engine.selection_digest();
    engine
        .set_selection(EvidenceSelection::Intervals(vec![GenomicInterval {
            contig: 0,
            start: 0,
            end: 10,
        }]))
        .unwrap();
    assert_eq!(first, engine.selection_digest());
}

#[test]
fn repeated_planning_preserves_admitted_baseline_and_charges_new_selection() {
    let f = Fixture::new(10);
    f.write(Vec::new(), bam::index::Type::Bai);
    let mut engine = EvidenceEngine::open(f.request(10)).unwrap();
    let baseline = engine.plan().baseline_rss_bytes;
    engine
        .set_selection(EvidenceSelection::Sites(vec![SnvSite {
            contig: 0,
            position: 1,
            reference: b'A',
            alternates: vec![b'C'],
        }]))
        .unwrap();
    assert_eq!(engine.plan().baseline_rss_bytes, baseline);
    assert!(engine.plan().selection_bytes >= 256);
    let analyzer = EvidenceCallback::new(|_: &EvidenceBatch| Ok(()), 1024);
    let admitted = engine
        .plan_for_analyzer(&analyzer)
        .unwrap()
        .predicted_peak_rss_bytes;
    let allocation = vec![1u8; 32 << 20];
    std::hint::black_box(&allocation);
    assert_eq!(
        engine
            .plan_for_analyzer(&analyzer)
            .unwrap()
            .predicted_peak_rss_bytes,
        admitted
    );
    assert_eq!(engine.plan().baseline_rss_bytes, baseline);
}

#[test]
fn arrow_reader_rejects_oversized_frames_before_decoder_allocation() {
    let f = Fixture::new(10);
    f.write(Vec::new(), bam::index::Type::Bai);
    let mut engine = EvidenceEngine::open(f.request(10)).unwrap();
    let mut writer = EvidenceArrowWriter::new(Vec::new());
    engine.run(&mut writer).unwrap();
    let bytes = writer.into_inner().unwrap();
    let mut huge_metadata = bytes.clone();
    huge_metadata[4..8].copy_from_slice(&i32::MAX.to_le_bytes());
    let error =
        read_evidence_batches(huge_metadata.as_slice(), engine.contigs(), |_| Ok(())).unwrap_err();
    assert!(error.to_string().contains("metadata exceeds"), "{error}");

    // FlatBuffers table fields are located via their verified vtable offsets.
    // Mutate only bodyLength in an otherwise valid record batch message.
    let first_metadata = i32::from_le_bytes(bytes[4..8].try_into().unwrap()) as usize;
    let second_start = 8 + first_metadata;
    assert_eq!(&bytes[second_start..second_start + 4], &[255; 4]);
    let metadata_start = second_start + 8;
    let mut huge_body = bytes.clone();
    let root = u32::from_le_bytes(
        huge_body[metadata_start..metadata_start + 4]
            .try_into()
            .unwrap(),
    ) as usize
        + metadata_start;
    let vtable = (root as isize
        - i32::from_le_bytes(huge_body[root..root + 4].try_into().unwrap()) as isize)
        as usize;
    let field_slot = vtable + arrow_ipc::Message::VT_BODYLENGTH as usize;
    let field = root
        + u16::from_le_bytes(huge_body[field_slot..field_slot + 2].try_into().unwrap()) as usize;
    assert!(field > root);
    huge_body[field..field + 8].copy_from_slice(&i64::MAX.to_le_bytes());
    let error =
        read_evidence_batches(huge_body.as_slice(), engine.contigs(), |_| Ok(())).unwrap_err();
    assert!(error.to_string().contains("body exceeds"), "{error}");
}

#[test]
fn equal_position_read_order_does_not_change_exact_integer_summaries() {
    let f = Fixture::new(20);
    let records = [
        record("a", 0, b"AAAA", &[30; 4], 60, 0, vec![Cigar::Match(4)]),
        record("c", 0, b"CCCC", &[40; 4], 30, 16, vec![Cigar::Match(4)]),
        record("g", 0, b"GGGG", &[20; 4], 20, 0, vec![Cigar::Match(4)]),
        record("t", 0, b"TTTT", &[25; 4], 50, 16, vec![Cigar::Match(4)]),
        record(
            "excluded",
            0,
            b"AAAA",
            &[30; 4],
            60,
            1024,
            vec![Cigar::Match(4)],
        ),
    ];
    let mut canonical = None;
    for (order, width) in [
        ([0, 1, 2, 3, 4], 1),
        ([4, 3, 2, 1, 0], 3),
        ([2, 0, 4, 1, 3], 16384),
    ] {
        f.write(
            order
                .into_iter()
                .map(|index| records[index].clone())
                .collect(),
            bam::index::Type::Bai,
        );
        let rows = collect(f.request(width));
        assert_eq!(rows[0].allele_counts, [1; 4]);
        assert_eq!(rows[0].filters.duplicate, 1);
        if let Some(canonical) = &canonical {
            assert_eq!(&rows, canonical);
        } else {
            canonical = Some(rows);
        }
    }
}

#[test]
fn file_backed_pack_and_legacy_reference_windows_match_fasta_across_word_and_contig_boundaries() {
    use rosalind::genomics::{GenomeIndex, IndexWriter, ReferencePackBuilder};
    let f = Fixture::new(150);
    let sequences = vec![
        ("chr1".to_string(), b"ACGTN".repeat(15)),
        ("chr2".to_string(), b"NNTGCA".repeat(13)),
    ];
    let first = String::from_utf8(sequences[0].1.clone()).unwrap();
    let second = String::from_utf8(sequences[1].1.clone()).unwrap();
    std::fs::write(&f.fasta, format!(">chr1\n{first}\n>chr2\n{second}\n")).unwrap();
    std::fs::write(
        f.root.join("reference.fa.fai"),
        "chr1\t75\t6\t75\t76\nchr2\t78\t88\t78\t79\n",
    )
    .unwrap();
    let pack = f.root.join("reference.rref");
    let index = f.root.join("reference.idx");
    ReferencePackBuilder::build(&f.fasta, &pack, false).unwrap();
    let genome = GenomeIndex::from_named_sequences(&sequences).unwrap();
    IndexWriter::create(&index)
        .unwrap()
        .write_genome_index(&genome)
        .unwrap();
    for path in [&f.fasta, &pack, &index] {
        let reference = EvidenceReference::open(path).unwrap();
        for (contig, (_, expected)) in sequences.iter().enumerate() {
            for start in 0..expected.len() {
                for width in [1, 7, 31, 32, 63, 64, 79] {
                    let end = (start + width).min(expected.len());
                    assert_eq!(
                        reference
                            .read_window(contig as u32, start as u32, end as u32)
                            .unwrap(),
                        expected[start..end],
                        "{} contig={contig} start={start} end={end}",
                        path.display()
                    );
                }
            }
        }
    }
}

#[test]
fn sam_reference_equality_bases_resolve_with_reference_and_refuse_without_it() {
    let f = Fixture::new(10);
    std::fs::write(&f.fasta, ">chr1\nACGTNAAAAA\n").unwrap();
    f.write(
        vec![record(
            "equal",
            0,
            b"=====",
            &[30; 5],
            60,
            16,
            vec![Cigar::Match(5)],
        )],
        bam::index::Type::Bai,
    );
    let rows = collect(f.request(1));
    for (position, allele) in [0, 1, 2, 3].into_iter().enumerate() {
        assert_eq!(rows[position].callable_depth, 1);
        assert_eq!(rows[position].allele_counts[allele], 1);
        assert_eq!(rows[position].strand_counts[allele], [0, 1]);
    }
    assert_eq!(rows[4].filters.ambiguous_base, 1);
    assert_eq!(rows[4].callable_depth, 0);
    let mut engine = EvidenceEngine::open(EvidenceRequest::coverage(&f.bam)).unwrap();
    let error = engine.run(&mut Capture::default()).unwrap_err();
    assert!(
        error
            .to_string()
            .contains("SEQ '=' requires a local analysis reference"),
        "{error}"
    );
}
