//! Record-preserving SNV annotation from verified canonical evidence partitions.
//!
//! INFO annotations describe alignment read evidence for the explicitly resolved
//! sample scope; they never overwrite VCF genotypes, FILTER, or existing fields.

use std::path::Path;

use rust_htslib::bcf;

use crate::core::governor::checkpoint;
use crate::dataset::{DatasetError, DatasetOutcome, VerifiedEvidenceLookup};
use crate::evidence::{EvidenceEngine, EvidenceError, EvidenceFields, EvidenceProfile};
use crate::variant_io::{
    parse_snv_record, CheckedVariantWriter, VariantFormat, VariantLimits, VariantReader,
};

/// Version of the INFO annotations and record-preservation contract.
pub const ANNOTATION_SEMANTICS: &str = "snv-read-evidence-info-v1";

/// Infer the requested encoding from the final (not staging) filename.
pub fn annotation_format(path: &Path) -> Result<VariantFormat, EvidenceError> {
    let name = path.to_string_lossy();
    if name.ends_with(".vcf.gz") {
        Ok(VariantFormat::VcfGz)
    } else if name.ends_with(".vcf") {
        Ok(VariantFormat::Vcf)
    } else if name.ends_with(".bcf") {
        Ok(VariantFormat::Bcf)
    } else {
        Err(EvidenceError::InvalidRequest(
            "annotated variant output must end in .vcf, .vcf.gz, or .bcf".into(),
        ))
    }
}

/// Additional conservative memory reservation, including verified lookup state,
/// native headers/records, formatted VCF scratch and compression buffers. Native
/// decoding is checked cooperatively after HTSlib allocation.
pub fn annotation_memory_bytes(engine: &EvidenceEngine, limits: VariantLimits) -> u64 {
    VerifiedEvidenceLookup::memory_bytes(engine)
        .saturating_add((limits.max_header_bytes as u64).saturating_mul(8))
        .saturating_add((limits.max_record_bytes as u64).saturating_mul(8))
        .saturating_add(8 << 20)
}

struct Tag {
    id: &'static str,
    number: &'static str,
    kind: &'static str,
    description: &'static str,
}

fn tags(fields: EvidenceFields) -> Vec<Tag> {
    let mut tags = vec![
        Tag {
            id: "RSL_DP",
            number: "1",
            kind: "Integer",
            description: "Rosalind callable A/C/G/T read depth",
        },
        Tag {
            id: "RSL_PF",
            number: "1",
            kind: "Integer",
            description:
                "Rosalind prefilter aligned base observations excluding deletion and reference skip",
        },
        Tag {
            id: "RSL_ED",
            number: "1",
            kind: "Integer",
            description:
                "Rosalind eligible aligned base observations before base-quality filtering",
        },
        Tag {
            id: "RSL_AD",
            number: "R",
            kind: "Integer",
            description: "Rosalind callable read counts in original REF and ALT order",
        },
    ];
    if fields.contains(EvidenceFields::STRANDS) {
        tags.extend([
            Tag {
                id: "RSL_ADF",
                number: "R",
                kind: "Integer",
                description: "Rosalind forward callable read counts in REF and ALT order",
            },
            Tag {
                id: "RSL_ADR",
                number: "R",
                kind: "Integer",
                description: "Rosalind reverse callable read counts in REF and ALT order",
            },
        ]);
    }
    if fields.contains(EvidenceFields::ALLELE_QUALITY) {
        tags.extend([
            Tag {
                id: "RSL_BQS",
                number: "R",
                kind: "String",
                description: "Rosalind exact unsigned base-quality sums in REF and ALT order",
            },
            Tag {
                id: "RSL_MQS",
                number: "R",
                kind: "String",
                description: "Rosalind exact unsigned mapping-quality sums in REF and ALT order",
            },
            Tag {
                id: "RSL_RPS",
                number: "R",
                kind: "String",
                description:
                    "Rosalind exact unsigned zero-based sequencing-cycle sums in REF and ALT order",
            },
            Tag {
                id: "RSL_RLS",
                number: "R",
                kind: "String",
                description: "Rosalind exact unsigned stored-read-length sums in REF and ALT order",
            },
        ]);
    }
    tags
}

/// Check required fields and header collisions before starting extraction.
pub fn preflight_annotation(
    input: &Path,
    fields: EvidenceFields,
    limits: VariantLimits,
) -> Result<(), EvidenceError> {
    if !fields.contains(EvidenceFields::DEPTHS.union(EvidenceFields::ALLELES)) {
        return Err(EvidenceError::InvalidRequest(
            "variant annotation requires --fields to include depths and alleles".into(),
        ));
    }
    let reader = VariantReader::open(input, limits)?;
    for tag in tags(fields) {
        if reader.header().info_type(tag.id.as_bytes()).is_ok() {
            return Err(EvidenceError::InvalidInput(format!(
                "variant INFO {} already exists; annotation never overwrites existing fields",
                tag.id
            )));
        }
    }
    Ok(())
}

fn integer(value: u64) -> Result<i32, EvidenceError> {
    i32::try_from(value).map_err(|_| {
        EvidenceError::InvalidInput(
            "read count exceeds the VCF Integer range; retain exact evidence in Arrow/TSV".into(),
        )
    })
}

fn allele_index(base: u8) -> usize {
    match base {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => unreachable!("validated SNV"),
    }
}

/// Write annotations to a staging path. The caller must publish it atomically
/// with its receipt only after this function succeeds and inputs are rechecked.
/// Record order, repeated sites, original allele order and sample data survive.
#[allow(clippy::too_many_arguments)]
pub fn annotate_variants(
    input: &Path,
    staging_output: &Path,
    format: VariantFormat,
    limits: VariantLimits,
    engine: &EvidenceEngine,
    dataset: &DatasetOutcome,
    science_digest: &str,
) -> Result<u64, DatasetError> {
    if science_digest.len() != 64 || !science_digest.bytes().all(|byte| byte.is_ascii_hexdigit()) {
        return Err(EvidenceError::InvalidRequest(
            "annotation science identity must be a 64-character digest".into(),
        )
        .into());
    }
    let fields = engine.request().fields;
    preflight_annotation(input, fields, limits)?;
    let mut lookup = VerifiedEvidenceLookup::new(engine, dataset)?;
    let mut reader = VariantReader::open(input, limits)?;
    let mut header = bcf::Header::from_template(reader.header());
    for tag in tags(fields) {
        header.push_record(
            format!(
                "##INFO=<ID={},Number={},Type={},Description=\"{}\">",
                tag.id, tag.number, tag.kind, tag.description
            )
            .as_bytes(),
        );
    }
    // JSON quoting escapes sample names without putting execution-dependent
    // paths, timestamps, worker counts or receipt hashes into deterministic VCF.
    let scope = format!(
        "\"{}\"",
        engine
            .sample_scope()
            .canonical_json()
            .replace('\\', "\\\\")
            .replace('\"', "\\\"")
    );
    header.push_record(format!("##rosalind_evidence=<Semantics={ANNOTATION_SEMANTICS},Profile={},Fields={},Science={science_digest},SampleScope={scope}>", EvidenceProfile::ID, fields.bits()).as_bytes());
    let mut writer = CheckedVariantWriter::create(staging_output, &header, format, limits)?;
    // The writer owns a checked header copy; avoid retaining another copy while
    // decoding records and evidence partitions.
    drop(header);
    let mut count = 0u64;
    while let Some(record) = reader.read()? {
        checkpoint().map_err(EvidenceError::from)?;
        let site = parse_snv_record(record, engine.contigs())?;
        let row = lookup.get(site.contig, site.position)?;
        if row.reference != site.reference
            || site
                .alternates
                .iter()
                .any(|alt| !row.requested_alts.contains(alt))
        {
            return Err(EvidenceError::InvalidInput(format!(
                "variant and evidence alleles differ at contig {} position {}",
                site.contig,
                site.position + 1
            ))
            .into());
        }
        let depths = row.depths.ok_or_else(|| {
            EvidenceError::InvalidInput("annotation evidence has no depths".into())
        })?;
        let alleles = row.alleles.ok_or_else(|| {
            EvidenceError::InvalidInput("annotation evidence has no allele counts".into())
        })?;
        let indices: Vec<_> = std::iter::once(site.reference)
            .chain(site.alternates)
            .map(allele_index)
            .collect();
        let mut ints: Vec<(&[u8], Vec<i32>)> = vec![
            (b"RSL_DP", vec![integer(depths.callable_depth)?]),
            (b"RSL_PF", vec![integer(depths.prefilter_depth)?]),
            (b"RSL_ED", vec![integer(depths.aligned_depth)?]),
            (
                b"RSL_AD",
                indices
                    .iter()
                    .map(|&i| integer(alleles.allele_counts[i]))
                    .collect::<Result<_, _>>()?,
            ),
        ];
        if let Some(strands) = row.strands {
            for (tag, strand) in [(b"RSL_ADF".as_slice(), 0), (b"RSL_ADR".as_slice(), 1)] {
                ints.push((
                    tag,
                    indices
                        .iter()
                        .map(|&i| integer(strands.strand_counts[i][strand]))
                        .collect::<Result<_, _>>()?,
                ));
            }
        }
        let mut strings: Vec<(&[u8], Vec<String>)> = Vec::new();
        if let Some(q) = row.allele_quality {
            for (tag, sums) in [
                (b"RSL_BQS".as_slice(), &q.base_quality_sum),
                (b"RSL_MQS".as_slice(), &q.mapping_quality_sum),
                (b"RSL_RPS".as_slice(), &q.read_position_sum),
                (b"RSL_RLS".as_slice(), &q.read_length_sum),
            ] {
                strings.push((tag, indices.iter().map(|&i| sums[i].to_string()).collect()));
            }
        }
        let int_refs: Vec<_> = ints.iter().map(|(tag, v)| (*tag, v.as_slice())).collect();
        let values: Vec<Vec<&[u8]>> = strings
            .iter()
            .map(|(_, v)| v.iter().map(|s| s.as_bytes()).collect())
            .collect();
        let string_refs: Vec<_> = strings
            .iter()
            .zip(&values)
            .map(|((tag, _), v)| (*tag, v.as_slice()))
            .collect();
        writer.write_with_info(record, &int_refs, &string_refs)?;
        count = count.checked_add(1).ok_or(EvidenceError::CounterOverflow)?;
    }
    writer.finish()?;
    lookup.verify_unchanged()?;
    checkpoint().map_err(EvidenceError::from)?;
    Ok(count)
}
