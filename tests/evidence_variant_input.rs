use rosalind::core::ContigSet;
use rosalind::evidence::{EvidenceError, EvidenceSelection};
use rosalind::variant_io::{
    parse_snv_record, CheckedVariantWriter, VariantFormat, VariantLimits, VariantReader,
};
use rust_htslib::bcf;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

static FIXTURE_ID: AtomicU64 = AtomicU64::new(0);

struct Fixture(PathBuf);
impl Fixture {
    fn new() -> Self {
        let path = std::env::temp_dir().join(format!(
            "rosalind-variant-io-{}-{}",
            std::process::id(),
            FIXTURE_ID.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir_all(&path).unwrap();
        Self(path)
    }
    fn write(&self, name: &str, text: impl AsRef<[u8]>) -> PathBuf {
        let path = self.0.join(name);
        std::fs::write(&path, text).unwrap();
        path
    }
}
impl Drop for Fixture {
    fn drop(&mut self) {
        std::fs::remove_dir_all(&self.0).unwrap();
    }
}

fn header() -> &'static str {
    "##fileformat=VCFv4.2\n\
     ##contig=<ID=chr1,length=100>\n\
     ##INFO=<ID=OLD,Number=1,Type=String,Description=\"Original annotation\">\n\
     ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
     ##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Original depth\">\n\
     #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
}

fn source(f: &Fixture) -> PathBuf {
    f.write(
        "input.vcf",
        format!(
            "{}chr1\t8\tsecond\tA\tT,C\t31\tPASS\tOLD=keep\tGT:DP\t2|1:7\n\
             chr1\t2\tfirst\tA\tG\t.\t.\tOLD=two\tGT:DP\t0/1:3\n\
             chr1\t8\tduplicate\tA\tG,T\t9\tPASS\t.\tGT:DP\t1/2:4\n",
            header()
        ),
    )
}

fn contigs() -> ContigSet {
    let mut contigs = ContigSet::new();
    contigs.push("chr1", 100);
    contigs
}

fn sites(path: &Path) -> Vec<rosalind::evidence::SnvSite> {
    match EvidenceSelection::from_vcf(path, &contigs()).unwrap() {
        EvidenceSelection::Sites(sites) => sites,
        _ => panic!("expected sites"),
    }
}

#[test]
fn vcf_bgzf_and_bcf_have_identical_normalized_unsorted_duplicate_selection() {
    let f = Fixture::new();
    let input = source(&f);
    let expected = sites(&input);
    assert_eq!(expected.len(), 2);
    assert_eq!(expected[0].position, 1);
    assert_eq!(expected[1].position, 7);
    assert_eq!(expected[1].alternates, b"CGT");
    for (name, format) in [
        ("out.vcf", VariantFormat::Vcf),
        ("out.vcf.gz", VariantFormat::VcfGz),
        ("out.bcf", VariantFormat::Bcf),
    ] {
        let mut reader = VariantReader::open(&input, VariantLimits::default()).unwrap();
        let output = f.0.join(name);
        let mut writer = CheckedVariantWriter::create(
            &output,
            &bcf::Header::from_template(reader.header()),
            format,
            VariantLimits::default(),
        )
        .unwrap();
        while let Some(record) = reader.read().unwrap() {
            writer.write(record).unwrap();
        }
        writer.finish().unwrap();
        assert_eq!(sites(&output), expected);
        let mut reader = VariantReader::open(output, VariantLimits::default()).unwrap();
        assert_eq!(reader.header().samples(), vec![b"SAMPLE".as_slice()]);
        for (position, id, alternates, genotype, depth) in [
            (7, "second", "TC", "2|1", 7),
            (1, "first", "G", "0/1", 3),
            (7, "duplicate", "GT", "1/2", 4),
        ] {
            let record = reader.read().unwrap().unwrap();
            let site = parse_snv_record(record, &contigs()).unwrap();
            assert_eq!(site.position, position);
            assert_eq!(site.alternates, alternates.as_bytes());
            assert_eq!(record.id(), id.as_bytes());
            assert_eq!(record.genotypes().unwrap().get(0).to_string(), genotype);
            assert_eq!(record.format(b"DP").integer().unwrap()[0], &[depth]);
        }
        assert!(reader.read().unwrap().is_none());
        assert_eq!(reader.record_number(), 3);
    }
}

#[test]
fn ordinary_gzip_vcf_is_supported_without_a_tabix_index() {
    let f = Fixture::new();
    let input = source(&f);
    let output = f.0.join("ordinary.vcf.gz");
    let mut gzip = flate2::write::GzEncoder::new(
        std::fs::File::create(&output).unwrap(),
        flate2::Compression::default(),
    );
    gzip.write_all(&std::fs::read(&input).unwrap()).unwrap();
    gzip.finish().unwrap();
    assert_eq!(sites(&output), sites(&input));
    let file = std::fs::OpenOptions::new()
        .write(true)
        .open(&output)
        .unwrap();
    file.set_len(file.metadata().unwrap().len() - 8).unwrap();
    drop(file);
    assert!(
        EvidenceSelection::from_vcf(&output, &contigs()).is_err(),
        "ordinary gzip truncation must be rejected"
    );
}

#[test]
fn sample_free_vcf_and_empty_stream_roundtrip_without_null_sample_access() {
    let f = Fixture::new();
    for (name, row) in [("site", "chr1\t1\t.\tA\tC\t.\t.\t.\n"), ("empty", "")] {
        let input = f.write(name, format!("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{row}"));
        let mut reader = VariantReader::open(&input, VariantLimits::default()).unwrap();
        assert_eq!(reader.header().sample_count(), 0);
        let output = f.0.join(format!("{name}.bcf"));
        let mut writer = CheckedVariantWriter::create(
            &output,
            &bcf::Header::from_template(reader.header()),
            VariantFormat::Bcf,
            VariantLimits::default(),
        )
        .unwrap();
        while let Some(record) = reader.read().unwrap() {
            writer.write(record).unwrap();
        }
        writer.finish().unwrap();
        assert_eq!(sites(&input), sites(&output));
    }
}

#[test]
fn annotation_updates_copy_and_preserves_original_header_record_and_allele_order() {
    let f = Fixture::new();
    let input = source(&f);
    let mut reader = VariantReader::open(&input, VariantLimits::default()).unwrap();
    let mut output_header = bcf::Header::from_template(reader.header());
    output_header
        .push_record(b"##INFO=<ID=COUNT,Number=R,Type=Integer,Description=\"Read counts\">");
    output_header
        .push_record(b"##INFO=<ID=SUM,Number=R,Type=String,Description=\"Exact integer sums\">");
    let output = f.0.join("annotated.bcf");
    let mut writer = CheckedVariantWriter::create(
        &output,
        &output_header,
        VariantFormat::Bcf,
        VariantLimits::default(),
    )
    .unwrap();
    let record = reader.read().unwrap().unwrap();
    writer
        .write_with_info(
            record,
            &[(b"COUNT", &[12, 7, 2])],
            &[(b"SUM", &[b"18446744073709551615", b"42", b"0"])],
        )
        .unwrap();
    assert!(record.header().info_type(b"COUNT").is_err());
    assert!(record.info(b"COUNT").integer().is_err());
    assert_eq!(record.alleles(), [b"A", b"T", b"C"]);
    assert_eq!(record.genotypes().unwrap().get(0).to_string(), "2|1");
    assert_eq!(record.info(b"OLD").string().unwrap().unwrap()[0], b"keep");
    writer.finish().unwrap();
    let mut reader = VariantReader::open(output, VariantLimits::default()).unwrap();
    let record = reader.read().unwrap().unwrap();
    assert_eq!(
        &*record.info(b"COUNT").integer().unwrap().unwrap(),
        &[12, 7, 2]
    );
    assert_eq!(
        &*record.info(b"SUM").string().unwrap().unwrap(),
        &[b"18446744073709551615".as_slice(), b"42", b"0"]
    );
    assert_eq!(record.id(), b"second");
    assert_eq!(record.genotypes().unwrap().get(0).to_string(), "2|1");
    assert_eq!(record.info(b"OLD").string().unwrap().unwrap()[0], b"keep");
}

#[test]
fn parser_warning_recovery_is_rejected_and_cannot_be_resumed() {
    let f = Fixture::new();
    for (name, row) in [
        ("info", "chr1\t1\t.\tA\tC\t.\t.\tUNDEFINED=4\tGT\t0/1\n"),
        ("format", "chr1\t1\t.\tA\tC\t.\t.\t.\tUNKNOWN\t4\n"),
        ("contig", "chr2\t1\t.\tA\tC\t.\t.\t.\tGT\t0/1\n"),
        ("columns", "chr1\t1\t.\tA\tC\t.\t.\t.\tGT\n"),
    ] {
        let path = f.write(name, format!("{}{row}", header()));
        let mut reader = VariantReader::open(path, VariantLimits::default()).unwrap();
        assert!(reader.read().is_err(), "{name}");
        assert!(
            reader.read().is_err(),
            "reader must remain poisoned: {name}"
        );
    }
}

#[test]
fn malformed_or_missing_header_fails_without_native_null_header_access() {
    let f = Fixture::new();
    for (name, contents) in [
        ("empty", ""),
        ("missing", "chr1\t1\t.\tA\tC\t.\t.\t.\n"),
        ("partial", "##fileformat=VCFv4.2\n"),
        ("bad", "##fileformat=VCFv4.2\n#CHROM\tPOS\n"),
    ] {
        let path = f.write(name, contents);
        assert!(
            VariantReader::open(path, VariantLimits::default()).is_err(),
            "{name}"
        );
    }
}

#[test]
fn snv_selection_rejects_non_snvs_coordinates_and_conflicting_duplicate_ref() {
    let f = Fixture::new();
    for (name, pos, reference, alternate) in [
        ("deletion", "1", "AT", "A"),
        ("insertion", "1", "A", "AT"),
        ("symbolic", "1", "A", "<DEL>"),
        ("missing", "1", "A", "."),
        ("n", "1", "N", "C"),
        ("same", "1", "A", "A"),
        ("duplicate_alt", "1", "A", "C,C"),
        ("zero", "0", "A", "C"),
        ("past_end", "101", "A", "C"),
    ] {
        let path = f.write(
            name,
            format!(
                "{}chr1\t{pos}\t.\t{reference}\t{alternate}\t.\t.\t.\tGT\t0/1\n",
                header()
            ),
        );
        assert!(
            EvidenceSelection::from_vcf(path, &contigs()).is_err(),
            "{name}"
        );
    }
    let conflict = f.write(
        "conflict",
        format!(
            "{}chr1\t1\t.\tA\tC\t.\t.\t.\tGT\t0/1\nchr1\t1\t.\tG\tT\t.\t.\t.\tGT\t0/1\n",
            header()
        ),
    );
    assert!(EvidenceSelection::from_vcf(conflict, &contigs())
        .unwrap_err()
        .to_string()
        .contains("conflicting REF"));
}

#[test]
fn header_record_and_annotation_envelopes_fail_closed() {
    let f = Fixture::new();
    let input = source(&f);
    let tiny_header = VariantLimits {
        max_header_bytes: 10,
        ..VariantLimits::default()
    };
    assert!(matches!(
        VariantReader::open(&input, tiny_header),
        Err(EvidenceError::RecordLimit(_))
    ));
    let tiny_record = VariantLimits {
        max_record_bytes: 10,
        ..VariantLimits::default()
    };
    let mut reader = VariantReader::open(&input, tiny_record).unwrap();
    assert!(matches!(reader.read(), Err(EvidenceError::RecordLimit(_))));
    let mut reader = VariantReader::open(&input, VariantLimits::default()).unwrap();
    let mut output_header = bcf::Header::from_template(reader.header());
    output_header.push_record(b"##INFO=<ID=LONG,Number=1,Type=String,Description=\"Large value\">");
    assert!(matches!(
        CheckedVariantWriter::create(
            f.0.join("bad-header"),
            &output_header,
            VariantFormat::Bcf,
            tiny_header
        ),
        Err(EvidenceError::RecordLimit(_))
    ));
    assert!(!f.0.join("bad-header").exists());
    let mut writer = CheckedVariantWriter::create(
        f.0.join("large.bcf"),
        &output_header,
        VariantFormat::Bcf,
        VariantLimits {
            max_record_bytes: 4096,
            ..VariantLimits::default()
        },
    )
    .unwrap();
    let record = reader.read().unwrap().unwrap();
    let large = vec![b'x'; 5000];
    assert!(matches!(
        writer.write_with_info(record, &[], &[(b"LONG", &[&large])]),
        Err(EvidenceError::RecordLimit(_))
    ));
    assert!(writer.write(record).is_err());
    assert!(writer.finish().is_err());
}

#[test]
fn output_type_errors_do_not_mutate_source_and_prevent_successful_finish() {
    let f = Fixture::new();
    let input = source(&f);
    let mut reader = VariantReader::open(&input, VariantLimits::default()).unwrap();
    let header = bcf::Header::from_template(reader.header());
    let record = reader.read().unwrap().unwrap();
    for key in [b"UNKNOWN".as_slice(), b"OLD"] {
        let mut writer = CheckedVariantWriter::create(
            f.0.join(String::from_utf8_lossy(key).as_ref()),
            &header,
            VariantFormat::Vcf,
            VariantLimits::default(),
        )
        .unwrap();
        assert!(writer.write_with_info(record, &[(key, &[1])], &[]).is_err());
        assert!(writer.finish().is_err());
        assert_eq!(record.info(b"OLD").string().unwrap().unwrap()[0], b"keep");
    }
}

#[test]
fn translation_cache_is_isolated_between_writers_using_the_same_input_record() {
    let f = Fixture::new();
    let input = source(&f);
    let mut reader = VariantReader::open(&input, VariantLimits::default()).unwrap();
    let original_header = bcf::Header::from_template(reader.header());
    let mut reordered = bcf::Header::new();
    reordered.push_record(b"##contig=<ID=chr1,length=100>");
    reordered
        .push_record(b"##INFO=<ID=BEFORE,Number=1,Type=String,Description=\"Changes native IDs\">");
    reordered
        .push_record(b"##INFO=<ID=OLD,Number=1,Type=String,Description=\"Original annotation\">");
    reordered.push_record(b"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    reordered.push_record(b"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Original depth\">");
    reordered.push_sample(b"SAMPLE");
    let record = reader.read().unwrap().unwrap();
    for (name, header) in [
        ("translated.bcf", reordered),
        ("original.bcf", original_header),
    ] {
        let output = f.0.join(name);
        let mut writer = CheckedVariantWriter::create(
            &output,
            &header,
            VariantFormat::Bcf,
            VariantLimits::default(),
        )
        .unwrap();
        writer.write(record).unwrap();
        writer.finish().unwrap();
        let mut output = VariantReader::open(output, VariantLimits::default()).unwrap();
        let roundtrip = output.read().unwrap().unwrap();
        assert_eq!(
            roundtrip.info(b"OLD").string().unwrap().unwrap()[0],
            b"keep"
        );
        assert_eq!(roundtrip.genotypes().unwrap().get(0).to_string(), "2|1");
        assert_eq!(roundtrip.format(b"DP").integer().unwrap()[0], &[7]);
    }
    assert_eq!(record.info(b"OLD").string().unwrap().unwrap()[0], b"keep");
}

#[test]
fn output_cannot_silently_omit_input_definitions_or_relabel_samples() {
    let f = Fixture::new();
    let input = source(&f);
    let mut reader = VariantReader::open(&input, VariantLimits::default()).unwrap();
    let mut missing = bcf::Header::from_template(reader.header());
    missing.remove_info(b"OLD");
    let mut samples = bcf::Header::from_template(reader.header());
    samples.push_sample(b"EXTRA");
    let mut length = bcf::Header::from_template(reader.header());
    length.remove_contig(b"chr1");
    length.push_record(b"##contig=<ID=chr1,length=200>");
    let record = reader.read().unwrap().unwrap();
    for (name, header) in [
        ("missing.bcf", missing),
        ("sample.bcf", samples),
        ("length.bcf", length),
    ] {
        let mut writer = CheckedVariantWriter::create(
            f.0.join(name),
            &header,
            VariantFormat::Bcf,
            VariantLimits::default(),
        )
        .unwrap();
        assert!(writer.write(record).is_err());
        assert!(writer.finish().is_err());
    }
}

#[test]
fn missing_bgzf_end_marker_cannot_be_accepted_as_successful_eof() {
    let f = Fixture::new();
    let input = source(&f);
    for format in [VariantFormat::VcfGz, VariantFormat::Bcf] {
        let mut reader = VariantReader::open(&input, VariantLimits::default()).unwrap();
        let output = f.0.join(format!("truncated-{format:?}"));
        let mut writer = CheckedVariantWriter::create(
            &output,
            &bcf::Header::from_template(reader.header()),
            format,
            VariantLimits::default(),
        )
        .unwrap();
        while let Some(record) = reader.read().unwrap() {
            writer.write(record).unwrap();
        }
        writer.finish().unwrap();
        let file = std::fs::OpenOptions::new()
            .write(true)
            .open(&output)
            .unwrap();
        let length = file.metadata().unwrap().len();
        file.set_len(length - 28).unwrap();
        drop(file);
        assert!(VariantReader::open(output, VariantLimits::default())
            .unwrap_err()
            .to_string()
            .contains("BGZF end marker"));
    }
}

#[cfg(target_os = "linux")]
#[test]
fn native_header_and_final_flush_errors_are_reported() {
    let mut header = bcf::Header::new();
    let writer = CheckedVariantWriter::create(
        "/dev/full",
        &header,
        VariantFormat::Vcf,
        VariantLimits::default(),
    )
    .unwrap();
    assert!(
        writer.finish().is_err(),
        "buffered header flush must fail explicitly"
    );
    header.push_record(format!("##long={}", "x".repeat(100_000)).as_bytes());
    assert!(
        CheckedVariantWriter::create(
            "/dev/full",
            &header,
            VariantFormat::Vcf,
            VariantLimits::default()
        )
        .is_err(),
        "large header must surface a write failure"
    );
}

#[cfg(target_os = "linux")]
#[test]
fn native_record_write_error_prevents_successful_finish() {
    let f = Fixture::new();
    let input = f.write(
        "large.vcf",
        format!(
            "{}chr1\t1\t.\tA\tC\t.\t.\tOLD={}\tGT\t0/1\n",
            header(),
            "x".repeat(100_000)
        ),
    );
    let mut reader = VariantReader::open(input, VariantLimits::default()).unwrap();
    let mut writer = CheckedVariantWriter::create(
        "/dev/full",
        &bcf::Header::from_template(reader.header()),
        VariantFormat::Vcf,
        VariantLimits::default(),
    )
    .unwrap();
    let record = reader.read().unwrap().unwrap();
    assert!(matches!(writer.write(record), Err(EvidenceError::Io(_))));
    assert!(writer.finish().is_err());
}
