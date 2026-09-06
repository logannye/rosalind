//! Verified random access retaining at most one canonical projected partition.

use super::*;
use crate::core::ContigSet;
use crate::evidence::{EvidenceBatch, EvidenceRowRef};

const PARTITION_RECEIPT_BYTES: u64 = 64 << 10;
const PARENT_BYTES_PER_PARTITION: u64 = 4096;
const LOOKUP_PATH_BYTES: usize = 4096;

#[derive(Debug)]
struct LookupPartition {
    ownership: Partition,
    receipt_hash: String,
    artifact_hash: String,
    observed: Option<InputSnapshot>,
}

/// Verified coordinate lookup into one completed dataset. Queries may be
/// unsorted or repeated. Only the current canonical partition's projected rows
/// are retained; returning to another partition re-verifies and decodes it.
///
/// Dataset files must remain immutable throughout the lookup session. Ordinary
/// replacements and timestamp/size changes are detected; unsigned receipts do
/// not establish authorship, and metadata stamps are not cryptographic identity.
#[derive(Debug)]
pub struct VerifiedEvidenceLookup {
    root: PathBuf,
    contigs: ContigSet,
    fields: EvidenceFields,
    science_digest: String,
    request_digest: String,
    parent_snapshot: InputSnapshot,
    partitions: Vec<LookupPartition>,
    loaded: Option<(usize, EvidenceBatch)>,
}

impl VerifiedEvidenceLookup {
    /// Conservative additional working set, excluding the already-open engine
    /// and caller's output/record buffers. Includes canonical descriptors and
    /// selections, receipt parsing, tracked file stamps, one projected 16,384-row
    /// partition, and the bounded Arrow decoder. No allocations are made here.
    pub fn memory_bytes(engine: &EvidenceEngine) -> u64 {
        let count = partition_count(engine);
        let fields = engine.request().fields;
        let parent_bytes = parent_receipt_limit(count);
        // ContigSet currently shares Arc names, but include copied names as well
        // so the bound also covers the active batch/decoder's owned strings.
        let dictionary = engine.contigs().iter().fold(0u64, |bytes, contig| {
            bytes.saturating_add(128 + (contig.name.len() as u64).saturating_mul(2))
        });
        let metadata = engine
            .plan()
            .selection_bytes
            .saturating_mul(3)
            .saturating_add(count.saturating_mul(12 << 10))
            // Receipt integrity verification clones the parsed claim and
            // serializes it again; reserve those transients as well as parsing.
            .saturating_add(parent_bytes.saturating_mul(8))
            .saturating_add(PARTITION_RECEIPT_BYTES.saturating_mul(8))
            .saturating_add(dictionary);
        metadata
            .saturating_add(
                fields
                    .storage_bytes_per_locus()
                    .saturating_add(32)
                    .saturating_mul(u64::from(CANONICAL_TILE_BASES)),
            )
            .saturating_add(evidence_reader_memory_bytes(fields))
            .saturating_add(64 << 10)
    }

    /// Verify the completed parent receipt against the engine's exact request.
    /// Partition payloads are verified when first accessed and on every reload.
    /// Reserve [`Self::memory_bytes`] before constructing a budgeted lookup.
    pub fn new(engine: &EvidenceEngine, outcome: &DatasetOutcome) -> Result<Self, DatasetError> {
        let fields = engine.request().fields;
        let digest = request_digest(engine);
        let count = partition_count(engine);
        require_path_size(&outcome.dataset_manifest)?;
        let parent_snapshot = InputSnapshot::capture([outcome.dataset_manifest.clone()])?;
        require_receipt_size(&outcome.dataset_manifest, parent_receipt_limit(count))?;
        let manifest = read_receipt(&outcome.dataset_manifest)?;
        let expected = [
            ("evidence.science_blake3", outcome.science_digest.clone()),
            ("dataset.request_blake3", digest.clone()),
            ("dataset.partition_count", count.to_string()),
            ("dataset.partition_bases", CANONICAL_TILE_BASES.to_string()),
            ("outcome.rows", engine.plan().selected_loci.to_string()),
            ("run_status", "completed".to_string()),
            ("pileup.semantics", "exact-or-fail-v1".to_string()),
        ];
        if expected
            .iter()
            .any(|(key, value)| manifest.params.get(*key) != Some(value))
            || manifest.inputs.len() as u64 != count
            || manifest.outputs.len() as u64 != count
            || receipt_fields(&manifest)? != fields
        {
            return Err(DatasetError::Corrupt(
                "lookup parent receipt differs from the completed evidence request".into(),
            ));
        }
        if outcome.science_digest.len() != 64
            || !outcome
                .science_digest
                .bytes()
                .all(|byte| byte.is_ascii_hexdigit())
        {
            return Err(DatasetError::Incompatible(
                "lookup science identity must be a 64-character digest".into(),
            ));
        }
        let ownership = partitions(engine);
        if ownership.len() as u64 != count {
            return Err(DatasetError::Corrupt(
                "lookup partition count is inconsistent".into(),
            ));
        }
        let mut verified = Vec::with_capacity(ownership.len());
        for ((part, receipt), artifact) in ownership
            .into_iter()
            .zip(&manifest.inputs)
            .zip(&manifest.outputs)
        {
            if part.row_count() > u64::from(CANONICAL_TILE_BASES) {
                return Err(DatasetError::Corrupt(
                    "lookup ownership exceeds a canonical partition".into(),
                ));
            }
            verified.push(LookupPartition {
                ownership: part,
                receipt_hash: receipt.blake3.clone(),
                artifact_hash: artifact.blake3.clone(),
                observed: None,
            });
        }
        parent_snapshot.verify()?;
        let root = outcome
            .dataset_manifest
            .parent()
            .ok_or_else(|| {
                DatasetError::Incompatible("dataset manifest has no parent directory".into())
            })?
            .to_path_buf();
        Ok(Self {
            root,
            contigs: engine.contigs().clone(),
            fields,
            science_digest: outcome.science_digest.clone(),
            request_digest: digest,
            parent_snapshot,
            partitions: verified,
            loaded: None,
        })
    }

    /// Borrow the requested selected locus. This never returns a synthetic zero
    /// for a missing position. The borrow prevents a reload until it is released.
    pub fn get(&mut self, contig: u32, position: u32) -> Result<EvidenceRowRef<'_>, DatasetError> {
        crate::core::governor::checkpoint().map_err(EvidenceError::from)?;
        self.parent_snapshot.verify()?;
        let canonical = position / CANONICAL_TILE_BASES * CANONICAL_TILE_BASES;
        let index = self
            .partitions
            .binary_search_by_key(&(contig, canonical), |partition| {
                (partition.ownership.contig, partition.ownership.start)
            })
            .map_err(|_| missing_locus(contig, position))?;
        let part = &self.partitions[index].ownership;
        let interval = part
            .intervals
            .partition_point(|interval| interval.end <= position);
        if part
            .intervals
            .get(interval)
            .is_none_or(|interval| interval.start > position)
        {
            return Err(missing_locus(contig, position));
        }
        if let Some(observed) = &self.partitions[index].observed {
            observed.verify()?;
        }
        if self.loaded.as_ref().map(|(loaded, _)| *loaded) != Some(index) {
            self.load(index)?;
        }
        let batch = &self.loaded.as_ref().expect("successful partition load").1;
        let row = batch
            .loci()
            .binary_search_by_key(&position, |locus| locus.position)
            .map_err(|_| {
                DatasetError::Corrupt("verified partition lost a selected locus".into())
            })?;
        Ok(batch.row(row).expect("verified locus index"))
    }

    /// Recheck ordinary mutation of the parent and every partition observed so
    /// far. Call this immediately before publishing an annotation artifact.
    pub fn verify_unchanged(&self) -> Result<(), DatasetError> {
        crate::core::governor::checkpoint().map_err(EvidenceError::from)?;
        self.parent_snapshot.verify()?;
        for partition in &self.partitions {
            if let Some(observed) = &partition.observed {
                observed.verify()?;
            }
        }
        Ok(())
    }

    fn load(&mut self, index: usize) -> Result<(), DatasetError> {
        // Release the previous canonical partition before retaining or decoding
        // another; a failed replacement cannot leave stale rows addressable.
        self.loaded = None;
        self.parent_snapshot.verify()?;
        let entry = &self.partitions[index];
        let part = &entry.ownership;
        let directory = self.root.join(part.name());
        let receipt_path = directory.join("manifest.json");
        let artifact_path = directory.join("evidence.arrow");
        require_path_size(&receipt_path)?;
        require_path_size(&artifact_path)?;
        let observed = InputSnapshot::capture([receipt_path.clone(), artifact_path.clone()])?;
        require_receipt_size(&receipt_path, PARTITION_RECEIPT_BYTES)?;
        if blake3_file(&receipt_path)? != entry.receipt_hash {
            return Err(DatasetError::Corrupt(
                "lookup partition receipt differs from its completed parent".into(),
            ));
        }
        let receipt = verify_partition(
            &directory,
            part,
            &self.science_digest,
            &self.request_digest,
            self.fields,
        )?;
        if receipt.outputs[0].blake3 != entry.artifact_hash {
            return Err(DatasetError::Corrupt(
                "lookup partition artifact differs from its completed parent".into(),
            ));
        }
        observed.verify()?;
        let name = self
            .contigs
            .by_id(part.contig)
            .ok_or_else(|| DatasetError::Corrupt("unknown lookup ownership contig".into()))?
            .name
            .to_string();
        let mut retained =
            EvidenceBatch::new(part.contig, name, part.start, self.fields, Vec::new());
        retained.reserve_exact(part.row_count() as usize);
        let mut interval_index = 0;
        let mut next_position = part.intervals.first().map_or(0, |interval| interval.start);
        let mut row_index = 0;
        let metadata = read_evidence_batches_expected_fields(
            BufReader::new(File::open(&artifact_path)?),
            &self.contigs,
            self.fields,
            |batch| {
                crate::core::governor::checkpoint()?;
                if batch.fields() != self.fields
                    || batch.contig_id != part.contig
                    || batch.canonical_tile_start != part.start
                {
                    return Err(EvidenceError::InvalidInput(
                        "lookup batch fields or canonical ownership differ".into(),
                    ));
                }
                for row in batch.rows() {
                    let interval = part.intervals.get(interval_index).ok_or_else(|| {
                        EvidenceError::InvalidInput("lookup partition contains extra loci".into())
                    })?;
                    if row.position != next_position {
                        return Err(EvidenceError::InvalidInput(
                            "lookup partition omits or duplicates a selected locus".into(),
                        ));
                    }
                    if let EvidenceSelection::Sites(sites) = &part.selection {
                        let site = sites.get(row_index).ok_or_else(|| {
                            EvidenceError::InvalidInput(
                                "lookup has an unselected variant locus".into(),
                            )
                        })?;
                        if row.position != site.position
                            || row.reference != site.reference
                            || row.requested_alts != site.alternates
                        {
                            return Err(EvidenceError::InvalidInput(
                                "lookup variant REF/ALT annotation differs from selection".into(),
                            ));
                        }
                    } else if !row.requested_alts.is_empty() {
                        return Err(EvidenceError::InvalidInput(
                            "lookup interval evidence contains unexpected ALT annotations".into(),
                        ));
                    }
                    retained.push_row(row)?;
                    row_index += 1;
                    next_position += 1;
                    if next_position == interval.end {
                        interval_index += 1;
                        if let Some(next) = part.intervals.get(interval_index) {
                            next_position = next.start;
                        }
                    }
                }
                Ok(())
            },
        )
        .map_err(|error| match error {
            error @ EvidenceError::Core(_) => DatasetError::Evidence(error),
            error => DatasetError::Corrupt(format!("invalid lookup Arrow evidence: {error}")),
        })?;
        if metadata.fields != self.fields
            || metadata.schema_version != self.fields.schema_version()
            || interval_index != part.intervals.len()
            || retained.len() as u64 != part.row_count()
        {
            return Err(DatasetError::Corrupt(
                "lookup partition lacks its complete selected field/locus set".into(),
            ));
        }
        observed.verify()?;
        self.parent_snapshot.verify()?;
        crate::core::governor::checkpoint().map_err(EvidenceError::from)?;
        self.partitions[index].observed = Some(observed);
        self.loaded = Some((index, retained));
        Ok(())
    }
}

fn partition_count(engine: &EvidenceEngine) -> u64 {
    let mut count = 0u64;
    let mut last = None;
    for interval in engine.intervals() {
        if interval.start == interval.end {
            continue;
        }
        let first = interval.start / CANONICAL_TILE_BASES;
        let final_tile = (interval.end - 1) / CANONICAL_TILE_BASES;
        count +=
            u64::from(final_tile - first + 1) - u64::from(last == Some((interval.contig, first)));
        last = Some((interval.contig, final_tile));
    }
    count
}

fn parent_receipt_limit(partitions: u64) -> u64 {
    (64u64 << 10).saturating_add(partitions.saturating_mul(PARENT_BYTES_PER_PARTITION))
}

fn require_receipt_size(path: &Path, maximum: u64) -> Result<(), DatasetError> {
    if fs::metadata(path)?.len() > maximum {
        return Err(DatasetError::Corrupt(format!(
            "lookup receipt exceeds its {maximum}-byte metadata envelope: {}",
            path.display()
        )));
    }
    Ok(())
}

fn require_path_size(path: &Path) -> Result<(), DatasetError> {
    // The per-partition metadata reservation includes two cloned stamp paths
    // plus ownership/hash metadata. Refuse a larger path before cloning it.
    if path.as_os_str().len() > LOOKUP_PATH_BYTES {
        return Err(DatasetError::Incompatible(format!(
            "lookup path exceeds the {LOOKUP_PATH_BYTES}-byte path envelope"
        )));
    }
    Ok(())
}

fn missing_locus(contig: u32, position: u32) -> DatasetError {
    DatasetError::Incompatible(format!(
        "locus {contig}:{position} is absent from the selected evidence dataset"
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::evidence::{EvidenceCallback, EvidenceLocus, EvidenceRequest, SnvSite};
    use rust_htslib::bam::{self, header::HeaderRecord, record::Cigar, record::CigarString};

    struct Fixture {
        root: PathBuf,
        engine: EvidenceEngine,
        outcome: DatasetOutcome,
    }

    impl Fixture {
        fn new() -> Self {
            let root = std::env::temp_dir().join(format!(
                "rosalind-evidence-lookup-{}-{}",
                std::process::id(),
                TEMP_SEQUENCE.fetch_add(1, Ordering::Relaxed)
            ));
            fs::create_dir(&root).unwrap();
            let length = CANONICAL_TILE_BASES + 4;
            let reference = root.join("reference.fa");
            fs::write(
                &reference,
                format!(">chr1\n{}\n>chr2\nCCCC\n", "A".repeat(length as usize)),
            )
            .unwrap();
            fs::write(
                root.join("reference.fa.fai"),
                format!(
                    "chr1\t{length}\t6\t{length}\t{}\nchr2\t4\t{}\t4\t5\n",
                    length + 1,
                    length + 13
                ),
            )
            .unwrap();
            let bam_path = root.join("reads.bam");
            let mut header = bam::Header::new();
            header.push_record(HeaderRecord::new(b"HD").push_tag(b"SO", "coordinate"));
            for (name, length) in [("chr1", length), ("chr2", 4)] {
                header.push_record(
                    HeaderRecord::new(b"SQ")
                        .push_tag(b"SN", name)
                        .push_tag(b"LN", length),
                );
            }
            let mut writer = bam::Writer::from_path(&bam_path, &header, bam::Format::Bam).unwrap();
            for (tid, position, base, quality, mapq) in [
                (0, 1, b'C', 30, 60),
                (0, CANONICAL_TILE_BASES + 1, b'G', 35, 42),
                (1, 1, b'T', 40, 55),
            ] {
                let mut record = bam::Record::new();
                record.set(
                    format!("r{tid}-{position}").as_bytes(),
                    Some(&CigarString(vec![Cigar::Match(1)])),
                    &[base],
                    &[quality],
                );
                record.set_tid(tid);
                record.set_pos(i64::from(position));
                record.set_flags(0);
                record.set_mapq(mapq);
                writer.write(&record).unwrap();
            }
            drop(writer);
            bam::index::build(&bam_path, None::<&PathBuf>, bam::index::Type::Bai, 1).unwrap();
            let mut request = EvidenceRequest::new(&bam_path, &reference);
            request.fields = EvidenceFields::DEPTHS
                .union(EvidenceFields::ALLELES)
                .union(EvidenceFields::ALLELE_QUALITY);
            request.selection = EvidenceSelection::Sites(vec![
                SnvSite {
                    contig: 0,
                    position: 1,
                    reference: b'A',
                    alternates: vec![b'C'],
                },
                SnvSite {
                    contig: 0,
                    position: 2,
                    reference: b'A',
                    alternates: vec![b'G'],
                },
                SnvSite {
                    contig: 0,
                    position: CANONICAL_TILE_BASES + 1,
                    reference: b'A',
                    alternates: vec![b'G', b'T'],
                },
                SnvSite {
                    contig: 1,
                    position: 1,
                    reference: b'C',
                    alternates: vec![b'T'],
                },
            ]);
            let mut engine = EvidenceEngine::open(request).unwrap();
            let mut consumer = EvidenceCallback::with_fields(
                |_: &EvidenceBatch| Ok(()),
                0,
                engine.request().fields,
            );
            let outcome = run_dataset(
                &mut engine,
                &"a".repeat(64),
                &DatasetOptions {
                    cache_dir: root.join("cache"),
                    resume: false,
                    workers: 1,
                },
                &mut consumer,
            )
            .unwrap();
            Self {
                root,
                engine,
                outcome,
            }
        }

        fn partition_path(&self, index: usize) -> PathBuf {
            self.outcome
                .dataset_manifest
                .parent()
                .unwrap()
                .join(partitions(&self.engine)[index].name())
        }

        fn replace_partition(&self, replacement: &EvidenceBatch) {
            let directory = self.partition_path(0);
            let artifact = directory.join("evidence.arrow");
            let mut writer = EvidenceArrowWriter::with_fields(
                File::create(&artifact).unwrap(),
                replacement.fields(),
            );
            writer.on_batch(replacement).unwrap();
            writer.finish().unwrap();
            drop(writer);
            let receipt_path = directory.join("manifest.json");
            let mut receipt = read_receipt(&receipt_path).unwrap();
            receipt.outputs[0].blake3 = blake3_file(&artifact).unwrap();
            receipt.finalize();
            fs::write(&receipt_path, receipt.to_canonical_json()).unwrap();
            let mut parent = read_receipt(&self.outcome.dataset_manifest).unwrap();
            parent.inputs[0].blake3 = blake3_file(&receipt_path).unwrap();
            parent.outputs[0].blake3 = receipt.outputs[0].blake3.clone();
            parent.finalize();
            fs::write(&self.outcome.dataset_manifest, parent.to_canonical_json()).unwrap();
        }
    }

    impl Drop for Fixture {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.root);
        }
    }

    #[test]
    fn unsorted_repeated_lookups_preserve_projected_allele_quality_and_zero_depth() {
        let fixture = Fixture::new();
        let mut lookup = VerifiedEvidenceLookup::new(&fixture.engine, &fixture.outcome).unwrap();
        for (contig, position, allele, bq, mq) in [
            (0, CANONICAL_TILE_BASES + 1, 2, 35, 42),
            (0, 1, 1, 30, 60),
            (1, 1, 3, 40, 55),
            (0, CANONICAL_TILE_BASES + 1, 2, 35, 42),
            (0, 1, 1, 30, 60),
        ] {
            let row = lookup.get(contig, position).unwrap();
            assert_eq!(row.depths.unwrap().callable_depth, 1);
            assert_eq!(row.alleles.unwrap().allele_counts[allele], 1);
            assert_eq!(row.allele_quality.unwrap().base_quality_sum[allele], bq);
            assert_eq!(row.allele_quality.unwrap().mapping_quality_sum[allele], mq);
            assert!(row.quality_histograms.is_none());
            assert!(lookup.loaded.as_ref().unwrap().1.len() <= CANONICAL_TILE_BASES as usize);
        }
        assert_eq!(lookup.get(0, 2).unwrap().depths.unwrap().callable_depth, 0);
        assert!(matches!(
            lookup.get(0, 0),
            Err(DatasetError::Incompatible(_))
        ));
        assert!(matches!(
            lookup.get(7, 1),
            Err(DatasetError::Incompatible(_))
        ));
        lookup.verify_unchanged().unwrap();
        assert!(
            VerifiedEvidenceLookup::memory_bytes(&fixture.engine)
                > fixture.engine.request().fields.storage_bytes_per_locus()
                    * u64::from(CANONICAL_TILE_BASES)
        );
    }

    #[test]
    fn parent_request_and_partition_hash_anchors_are_required() {
        let fixture = Fixture::new();
        let mut changed = fixture.engine.request().clone();
        changed.profile.min_mapq += 1;
        let changed = EvidenceEngine::open(changed).unwrap();
        assert!(VerifiedEvidenceLookup::new(&changed, &fixture.outcome).is_err());
        let path = fixture.partition_path(0).join("manifest.json");
        let mut receipt = read_receipt(&path).unwrap();
        receipt.params.insert("extra".into(), "resealed".into());
        receipt.finalize();
        fs::write(&path, receipt.to_canonical_json()).unwrap();
        let mut lookup = VerifiedEvidenceLookup::new(&fixture.engine, &fixture.outcome).unwrap();
        assert!(matches!(lookup.get(0, 1), Err(DatasetError::Corrupt(_))));
        assert!(lookup.loaded.is_none());
    }

    #[test]
    fn resealed_rows_cannot_change_selected_order_annotations_or_fields() {
        for mutation in [
            "missing",
            "duplicate",
            "alt",
            "reference",
            "fields",
            "contig",
        ] {
            let fixture = Fixture::new();
            let mut loci = vec![
                EvidenceLocus {
                    position: 1,
                    reference: b'A',
                    requested_alts: vec![b'C'],
                },
                EvidenceLocus {
                    position: 2,
                    reference: b'A',
                    requested_alts: vec![b'G'],
                },
            ];
            let mut fields = fixture.engine.request().fields;
            let mut contig = 0;
            let mut name = "chr1";
            match mutation {
                "missing" => {
                    loci.pop();
                }
                "duplicate" => {
                    loci[1] = loci[0].clone();
                }
                "alt" => {
                    loci[0].requested_alts = vec![b'T'];
                }
                "reference" => {
                    loci[0].reference = b'G';
                }
                "fields" => {
                    fields = EvidenceFields::DEPTHS;
                }
                "contig" => {
                    contig = 1;
                    name = "chr2";
                }
                _ => unreachable!(),
            }
            fixture.replace_partition(&EvidenceBatch::new(contig, name, 0, fields, loci));
            let mut lookup =
                VerifiedEvidenceLookup::new(&fixture.engine, &fixture.outcome).unwrap();
            assert!(
                matches!(lookup.get(0, 1), Err(DatasetError::Corrupt(_))),
                "{mutation}"
            );
            assert!(lookup.loaded.is_none(), "{mutation}");
        }
    }

    #[test]
    fn final_guard_checks_previously_loaded_partitions_and_parent_mutation() {
        let fixture = Fixture::new();
        let mut lookup = VerifiedEvidenceLookup::new(&fixture.engine, &fixture.outcome).unwrap();
        lookup.get(0, 1).unwrap();
        lookup.get(1, 1).unwrap();
        fs::write(fixture.partition_path(0).join("evidence.arrow"), b"changed").unwrap();
        assert!(lookup.verify_unchanged().is_err());
        assert!(lookup.get(0, 1).is_err());
        let fixture = Fixture::new();
        let mut lookup = VerifiedEvidenceLookup::new(&fixture.engine, &fixture.outcome).unwrap();
        fs::write(&fixture.outcome.dataset_manifest, b"changed").unwrap();
        assert!(lookup.get(0, 1).is_err());
    }

    #[test]
    fn oversized_parent_receipt_refuses_before_parsing() {
        let fixture = Fixture::new();
        File::options()
            .write(true)
            .open(&fixture.outcome.dataset_manifest)
            .unwrap()
            .set_len(parent_receipt_limit(partition_count(&fixture.engine)) + 1)
            .unwrap();
        match VerifiedEvidenceLookup::new(&fixture.engine, &fixture.outcome) {
            Err(DatasetError::Corrupt(message)) => assert!(message.contains("metadata envelope")),
            _ => panic!("oversized metadata was not refused"),
        }
    }
}
