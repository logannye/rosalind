use super::descriptor::*;
use super::store::*;
use crate::dataset::{
    publish_evidence_dataset, run_dataset_with_snapshot, DatasetDescriptor, DatasetOptions,
    DatasetReadLimits, DescriptorLimits, VerifiedEvidenceDataset, VerifiedInputSession,
};
use crate::evidence::{
    EvidenceBatch, EvidenceCallback, EvidenceEngine, EvidenceFields, EvidenceProfile,
    EvidenceRequest, EvidenceSelection,
};
use crate::selection::GenomicInterval;
use rust_htslib::bam::{
    self,
    header::HeaderRecord,
    record::{Aux, Cigar, CigarString},
};
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

static FIXTURE_ID: AtomicU64 = AtomicU64::new(0);

#[derive(Debug, Clone)]
pub(crate) struct FixtureOptions {
    pub sample: Option<String>,
    pub fields: EvidenceFields,
    pub profile: EvidenceProfile,
    pub reference_base: u8,
    pub reference_length: u32,
    pub cram: bool,
    pub cram_reference_only: bool,
}
impl Default for FixtureOptions {
    fn default() -> Self {
        Self {
            sample: Some("sample-A".into()),
            fields: EvidenceFields::DEPTHS.union(EvidenceFields::ALLELES),
            profile: EvidenceProfile::default(),
            reference_base: b'A',
            reference_length: 32,
            cram: false,
            cram_reference_only: false,
        }
    }
}

#[derive(Debug)]
pub(crate) struct Fixture {
    pub root: PathBuf,
    pub manifest: PathBuf,
    request: EvidenceRequest,
    inputs: Vec<(String, PathBuf)>,
}
impl Fixture {
    pub fn new(options: FixtureOptions) -> Self {
        let root = temp("fixture");
        fs::create_dir(&root).unwrap();
        let reference = root.join("reference.fa");
        fs::write(
            &reference,
            format!(
                ">chr1\n{}\n",
                (options.reference_base as char)
                    .to_string()
                    .repeat(options.reference_length as usize)
            ),
        )
        .unwrap();
        let fai = root.join("reference.fa.fai");
        fs::write(
            &fai,
            format!(
                "chr1\t{}\t6\t{}\t{}\n",
                options.reference_length,
                options.reference_length,
                options.reference_length + 1
            ),
        )
        .unwrap();
        let alignment = root.join(if options.cram {
            "reads.cram"
        } else {
            "reads.bam"
        });
        let mut header = bam::Header::new();
        header.push_record(HeaderRecord::new(b"HD").push_tag(b"SO", "coordinate"));
        header.push_record(
            HeaderRecord::new(b"SQ")
                .push_tag(b"SN", "chr1")
                .push_tag(b"LN", options.reference_length),
        );
        if let Some(sample) = &options.sample {
            header.push_record(
                HeaderRecord::new(b"RG")
                    .push_tag(b"ID", "rg1")
                    .push_tag(b"SM", sample),
            );
        }
        let mut writer = bam::Writer::from_path(
            &alignment,
            &header,
            if options.cram {
                bam::Format::Cram
            } else {
                bam::Format::Bam
            },
        )
        .unwrap();
        if options.cram {
            writer.set_reference(&reference).unwrap();
        }
        for position in [1, 3, 8] {
            let mut record = bam::Record::new();
            record.set(
                format!("read-{position}").as_bytes(),
                Some(&CigarString(vec![Cigar::Match(1)])),
                b"C",
                &[35],
            );
            record.set_tid(0);
            record.set_pos(position);
            record.set_flags(0);
            record.set_mapq(60);
            if options.sample.is_some() {
                record.push_aux(b"RG", Aux::String("rg1")).unwrap();
            }
            writer.write(&record).unwrap();
        }
        drop(writer);
        let index = root.join(if options.cram {
            "reads.cram.crai"
        } else {
            "reads.bam.bai"
        });
        bam::index::build(&alignment, Some(&index), bam::index::Type::Bai, 1).unwrap();
        let mut request = EvidenceRequest::new(&alignment, &reference);
        request.alignment_index = Some(index.clone());
        request.fields = options.fields;
        request.profile = options.profile;
        request.selection = intervals(0, 4);
        let mut inputs = vec![
            ("alignments".into(), alignment),
            ("alignment-index".into(), index),
        ];
        if options.cram_reference_only {
            assert!(options.cram);
            request.reference = None;
            request.cram_reference = Some(reference.clone());
            request.cram_reference_fai = Some(fai.clone());
            inputs.extend([
                ("cram-reference".into(), reference),
                ("cram-reference-fai".into(), fai),
            ]);
        } else {
            request.reference_fai = Some(fai.clone());
            inputs.extend([
                ("reference".into(), reference),
                ("reference-fai".into(), fai),
            ]);
        }
        let manifest = make_dataset(&root, &request, &inputs);
        Self {
            root,
            manifest,
            request,
            inputs,
        }
    }
    pub fn descriptor(&self) -> DatasetDescriptor {
        VerifiedEvidenceDataset::open(&self.manifest, DatasetReadLimits::default())
            .unwrap()
            .descriptor()
            .clone()
    }
    pub fn extra_leaf(&self, start: u32, end: u32, fields: EvidenceFields) -> PathBuf {
        let mut request = self.request.clone();
        request.selection = intervals(start, end);
        request.fields = fields;
        make_dataset(&self.root, &request, &self.inputs)
    }
    pub fn member(&self, id: &str) -> ImportMember {
        ImportMember {
            metadata: metadata(id),
            manifests: vec![self.manifest.clone()],
        }
    }
}
impl Drop for Fixture {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.root);
    }
}

fn make_dataset(root: &Path, request: &EvidenceRequest, inputs: &[(String, PathBuf)]) -> PathBuf {
    let mut engine = EvidenceEngine::open(request.clone()).unwrap();
    let session = VerifiedInputSession::open(inputs.to_vec()).unwrap();
    let namespace = session.dataset_namespace(&engine).unwrap();
    let mut callback = EvidenceCallback::with_fields(|_: &EvidenceBatch| Ok(()), 0, request.fields);
    let outcome = run_dataset_with_snapshot(
        &mut engine,
        &namespace,
        &DatasetOptions {
            cache_dir: root.join(format!(
                "cache-{}",
                FIXTURE_ID.fetch_add(1, Ordering::Relaxed)
            )),
            resume: false,
            workers: 1,
        },
        &mut callback,
        session.snapshot(),
    )
    .unwrap();
    publish_evidence_dataset(&engine, &outcome, &session, DescriptorLimits::default()).unwrap()
}
fn intervals(start: u32, end: u32) -> EvidenceSelection {
    EvidenceSelection::Intervals(vec![GenomicInterval {
        contig: 0,
        start,
        end,
    }])
}
pub(crate) fn metadata(id: &str) -> MemberMetadata {
    MemberMetadata {
        id: id.into(),
        group: None,
        subject: None,
        timepoint: None,
    }
}
fn temp(label: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "rosalind-cohort-{label}-{}-{}",
        std::process::id(),
        FIXTURE_ID.fetch_add(1, Ordering::Relaxed)
    ))
}
struct Root(PathBuf);
impl Root {
    fn new() -> Self {
        Self(temp("store"))
    }
}
impl Drop for Root {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

#[test]
fn imports_complete_portable_inventory_relocates_and_deduplicates() {
    let fixture = Fixture::new(FixtureOptions::default());
    let root = Root::new();
    let first = create_snapshot(
        &root.0,
        &[fixture.member("A")],
        None,
        CohortLimits::default(),
    )
    .unwrap();
    let second = create_snapshot(
        &root.0,
        &[fixture.member("A")],
        None,
        CohortLimits::default(),
    )
    .unwrap();
    assert_eq!(first.id, second.id);
    assert_eq!(fs::read_dir(root.0.join("objects")).unwrap().count(), 1);
    let object = root
        .0
        .join("objects")
        .join(&first.descriptor.members[0].leaves[0].object_id);
    assert!(!object.join("dataset.manifest.json").exists());
    assert_eq!(
        fs::read(&fixture.manifest).unwrap(),
        fs::read(object.join("evidence-dataset.manifest.json")).unwrap()
    );
    let relocated = Root::new();
    copy_tree(&root.0, &relocated.0);
    drop(fixture); // All original source alignments/reference/cache paths disappear.
    fs::remove_dir_all(&root.0).unwrap();
    verify_snapshot(&relocated.0, &first.id, CohortLimits::default()).unwrap();
}

#[test]
fn named_samples_fields_and_disjoint_same_partition_leaves_are_independent() {
    let a = Fixture::new(FixtureOptions::default());
    let b = Fixture::new(FixtureOptions {
        sample: Some("sample-B".into()),
        fields: EvidenceFields::ALL,
        ..FixtureOptions::default()
    });
    let mut member = a.member("A");
    member
        .manifests
        .push(a.extra_leaf(4, 9, EvidenceFields::ALL));
    let root = Root::new();
    let snapshot = create_snapshot(
        &root.0,
        &[b.member("B"), member],
        None,
        CohortLimits::default(),
    )
    .unwrap();
    assert_eq!(snapshot.descriptor.members[0].metadata.id, "A");
    verify_snapshot(&root.0, &snapshot.id, CohortLimits::default()).unwrap();
}

#[test]
fn rejects_duplicate_members_aliases_and_overlapping_ownership() {
    let fixture = Fixture::new(FixtureOptions::default());
    let root = Root::new();
    assert!(create_snapshot(
        &root.0,
        &[fixture.member("A"), fixture.member("A")],
        None,
        CohortLimits::default()
    )
    .is_err());
    assert!(!root.0.exists());
    assert!(create_snapshot(
        &root.0,
        &[fixture.member("A"), fixture.member("alias")],
        None,
        CohortLimits::default()
    )
    .is_err());
    let mut overlap = fixture.member("A");
    overlap
        .manifests
        .push(fixture.extra_leaf(3, 8, EvidenceFields::ALL));
    assert!(create_snapshot(&root.0, &[overlap], None, CohortLimits::default()).is_err());
    assert_eq!(fs::read_dir(root.0.join("snapshots")).unwrap().count(), 0);
}

#[test]
fn parent_identity_is_verified_and_existing_snapshots_survive_failed_import() {
    let fixture = Fixture::new(FixtureOptions::default());
    let root = Root::new();
    let first = create_snapshot(
        &root.0,
        &[fixture.member("A")],
        None,
        CohortLimits::default(),
    )
    .unwrap();
    let child = create_snapshot(
        &root.0,
        &[fixture.member("A")],
        Some(&first.id),
        CohortLimits::default(),
    )
    .unwrap();
    assert_ne!(first.id, child.id);
    let bad = Fixture::new(FixtureOptions {
        sample: Some("B".into()),
        ..FixtureOptions::default()
    });
    let arrow = bad
        .manifest
        .parent()
        .unwrap()
        .join(&bad.descriptor().partitions[0].arrow.path);
    fs::write(arrow, b"corrupt").unwrap();
    assert!(create_snapshot(
        &root.0,
        &[bad.member("B")],
        Some(&child.id),
        CohortLimits::default()
    )
    .is_err());
    verify_snapshot(&root.0, &first.id, CohortLimits::default()).unwrap();
    verify_snapshot(&root.0, &child.id, CohortLimits::default()).unwrap();
    fs::write(
        root.0.join("snapshots").join(&first.id).join(SNAPSHOT_FILE),
        b"{}",
    )
    .unwrap();
    assert!(open_snapshot(&root.0, &child.id, CohortLimits::default()).is_err());
}

#[test]
fn corrupt_existing_object_is_never_overwritten_or_deduplicated() {
    let fixture = Fixture::new(FixtureOptions::default());
    let root = Root::new();
    let snapshot = create_snapshot(
        &root.0,
        &[fixture.member("A")],
        None,
        CohortLimits::default(),
    )
    .unwrap();
    let object = root
        .0
        .join("objects")
        .join(&snapshot.descriptor.members[0].leaves[0].object_id);
    let arrow = object.join(&fixture.descriptor().partitions[0].arrow.path);
    fs::write(&arrow, b"changed").unwrap();
    assert!(create_snapshot(
        &root.0,
        &[fixture.member("A")],
        None,
        CohortLimits::default()
    )
    .is_err());
    assert_eq!(fs::read(&arrow).unwrap(), b"changed");
    assert!(verify_snapshot(&root.0, &snapshot.id, CohortLimits::default()).is_err());
}

#[test]
fn metadata_and_budget_envelopes_refuse_without_success() {
    let fixture = Fixture::new(FixtureOptions::default());
    let root = Root::new();
    let limits = CohortLimits {
        memory_budget_bytes: Some(1),
        ..CohortLimits::default()
    };
    assert!(create_snapshot(&root.0, &[fixture.member("A")], None, limits).is_err());
    assert!(!root.0.exists());
    let limits = CohortLimits {
        max_snapshot_bytes: 1,
        ..CohortLimits::default()
    };
    assert!(create_snapshot(&root.0, &[fixture.member("A")], None, limits).is_err());
    assert_eq!(fs::read_dir(root.0.join("snapshots")).unwrap().count(), 0);
}

#[test]
fn nested_dataset_budget_never_weakens_storage_admission() {
    for (outer, inner, expected) in [
        (None, Some(1), Some(1)),
        (Some(1), None, Some(1)),
        (Some(100), Some(1), Some(1)),
        (Some(1), Some(100), Some(1)),
        (None, None, None),
    ] {
        let limits = CohortLimits {
            memory_budget_bytes: outer,
            dataset: DatasetReadLimits {
                memory_budget_bytes: inner,
                ..DatasetReadLimits::default()
            },
            ..CohortLimits::default()
        };
        assert_eq!(limits.dataset_limits().memory_budget_bytes, expected);
        if expected.is_some() {
            let fixture = Fixture::new(FixtureOptions::default());
            let root = Root::new();
            assert!(create_snapshot(&root.0, &[fixture.member("A")], None, limits).is_err());
            assert!(!root.0.exists());
        }
    }
}

#[cfg(unix)]
#[test]
fn import_rejects_file_directory_and_manifest_symlinks() {
    use std::os::unix::fs::symlink;
    for name in ["manifest", "payload", "directory"] {
        let fixture = Fixture::new(FixtureOptions::default());
        let descriptor = fixture.descriptor();
        let cache = fixture.manifest.parent().unwrap();
        let path = match name {
            "manifest" => fixture.manifest.clone(),
            "payload" => cache.join(&descriptor.partitions[0].arrow.path),
            _ => cache.join(
                Path::new(&descriptor.partitions[0].arrow.path)
                    .parent()
                    .unwrap(),
            ),
        };
        let original = path.with_extension("original");
        fs::rename(&path, &original).unwrap();
        symlink(&original, &path).unwrap();
        let root = Root::new();
        assert!(
            create_snapshot(
                &root.0,
                &[fixture.member("A")],
                None,
                CohortLimits::default()
            )
            .is_err(),
            "{name}"
        );
        assert_eq!(fs::read_dir(root.0.join("snapshots")).unwrap().count(), 0);
    }
}

#[test]
fn snapshot_schema_and_metadata_are_canonical_bounded_and_versioned() {
    let limits = CohortLimits::default();
    let descriptor = SnapshotDescriptor {
        version: 1,
        comparison_version: 1,
        parent: None,
        members: vec![],
    };
    let bytes = descriptor.to_bytes(&limits).unwrap();
    assert_eq!(
        SnapshotDescriptor::from_bytes(&bytes, &limits).unwrap(),
        descriptor
    );
    let mut whitespace = bytes.clone();
    whitespace.push(b'\n');
    assert!(SnapshotDescriptor::from_bytes(&whitespace, &limits).is_err());
    assert!(SnapshotDescriptor::from_bytes(
        br#"{"version":2,"comparison_version":1,"parent":null,"members":[]}"#,
        &limits
    )
    .is_err());
    assert!(SnapshotDescriptor::from_bytes(
        br#"{"version":1,"comparison_version":1,"parent":null,"members":[],"other":1}"#,
        &limits
    )
    .is_err());
    let mut invalid = metadata("A");
    invalid.subject = Some("private\nnewline".into());
    assert!(invalid.validate().is_err());
}

#[test]
fn concurrent_identical_imports_publish_one_snapshot_without_replacement() {
    let fixture = Fixture::new(FixtureOptions::default());
    let root = Root::new();
    fs::create_dir(&root.0).unwrap();
    std::thread::scope(|scope| {
        let first = scope.spawn(|| {
            create_snapshot(
                &root.0,
                &[fixture.member("A")],
                None,
                CohortLimits::default(),
            )
        });
        let second = scope.spawn(|| {
            create_snapshot(
                &root.0,
                &[fixture.member("A")],
                None,
                CohortLimits::default(),
            )
        });
        assert_eq!(
            first.join().unwrap().unwrap().id,
            second.join().unwrap().unwrap().id
        );
    });
    assert_eq!(fs::read_dir(root.0.join("snapshots")).unwrap().count(), 1);
    assert_eq!(fs::read_dir(root.0.join("objects")).unwrap().count(), 1);
}

#[test]
fn unlisted_empty_directories_are_rejected_in_objects_and_snapshots() {
    let fixture = Fixture::new(FixtureOptions::default());
    let root = Root::new();
    let snapshot = create_snapshot(
        &root.0,
        &[fixture.member("A")],
        None,
        CohortLimits::default(),
    )
    .unwrap();
    let object_extra = root
        .0
        .join("objects")
        .join(&snapshot.descriptor.members[0].leaves[0].object_id)
        .join("undeclared");
    fs::create_dir(&object_extra).unwrap();
    assert!(open_snapshot(&root.0, &snapshot.id, CohortLimits::default()).is_err());
    fs::remove_dir(&object_extra).unwrap();
    let snapshot_extra = root
        .0
        .join("snapshots")
        .join(&snapshot.id)
        .join("undeclared");
    fs::create_dir(&snapshot_extra).unwrap();
    assert!(open_snapshot(&root.0, &snapshot.id, CohortLimits::default()).is_err());
}

#[test]
fn failure_at_each_publication_boundary_never_publishes_a_snapshot() {
    for stage in ["before_object_publication", "before_snapshot_publication"] {
        let fixture = Fixture::new(FixtureOptions::default());
        let root = Root::new();
        let result = with_test_hook(
            stage,
            || Err(std::io::Error::other("injected publication failure").into()),
            || {
                create_snapshot(
                    &root.0,
                    &[fixture.member("A")],
                    None,
                    CohortLimits::default(),
                )
            },
        );
        assert!(result.is_err());
        assert_eq!(fs::read_dir(root.0.join("snapshots")).unwrap().count(), 0);
        assert_eq!(
            fs::read_dir(root.0.join("objects")).unwrap().count(),
            usize::from(stage == "before_snapshot_publication")
        );
        let retry = create_snapshot(
            &root.0,
            &[fixture.member("A")],
            None,
            CohortLimits::default(),
        )
        .unwrap();
        verify_snapshot(&root.0, &retry.id, CohortLimits::default()).unwrap();
    }
}

#[test]
fn mutation_after_verification_is_detected_before_snapshot_publication() {
    for which in ["payload", "descriptor", "staged-snapshot"] {
        let fixture = Fixture::new(FixtureOptions::default());
        let root = Root::new();
        let root_path = root.0.clone();
        let payload = fixture.descriptor().partitions[0].arrow.path.clone();
        let result = with_test_hook(
            "before_snapshot_publication",
            move || {
                let path = if which == "staged-snapshot" {
                    fs::read_dir(root_path.join("snapshots"))?
                        .next()
                        .unwrap()?
                        .path()
                        .join(SNAPSHOT_FILE)
                } else {
                    fs::read_dir(root_path.join("objects"))?
                        .next()
                        .unwrap()?
                        .path()
                        .join(if which == "payload" {
                            payload.as_str()
                        } else {
                            "dataset.descriptor.json"
                        })
                };
                fs::write(path, b"changed after verification")?;
                Ok(())
            },
            || {
                create_snapshot(
                    &root.0,
                    &[fixture.member("A")],
                    None,
                    CohortLimits::default(),
                )
            },
        );
        assert!(result.is_err(), "{which}");
        assert_eq!(fs::read_dir(root.0.join("snapshots")).unwrap().count(), 0);
    }
}

#[test]
fn parent_mutation_during_child_creation_is_detected() {
    let fixture = Fixture::new(FixtureOptions::default());
    let root = Root::new();
    let first = create_snapshot(
        &root.0,
        &[fixture.member("A")],
        None,
        CohortLimits::default(),
    )
    .unwrap();
    let parent_file = root.0.join("snapshots").join(&first.id).join(SNAPSHOT_FILE);
    let parent_bytes = fs::read(&parent_file).unwrap();
    let changed_path = parent_file.clone();
    let result = with_test_hook(
        "before_snapshot_publication",
        move || fs::write(changed_path, b"{}").map_err(Into::into),
        || {
            create_snapshot(
                &root.0,
                &[fixture.member("A")],
                Some(&first.id),
                CohortLimits::default(),
            )
        },
    );
    assert!(result.is_err());
    assert_eq!(fs::read_dir(root.0.join("snapshots")).unwrap().count(), 1);
    fs::write(parent_file, parent_bytes).unwrap();
    verify_snapshot(&root.0, &first.id, CohortLimits::default()).unwrap();
}

#[cfg(unix)]
#[test]
fn copied_payload_has_its_own_inode_and_rejects_late_symlink() {
    use std::os::unix::fs::{symlink, MetadataExt};
    let fixture = Fixture::new(FixtureOptions::default());
    let root = Root::new();
    let first = create_snapshot(
        &root.0,
        &[fixture.member("A")],
        None,
        CohortLimits::default(),
    )
    .unwrap();
    let relative = fixture.descriptor().partitions[0].arrow.path.clone();
    let original = fixture.manifest.parent().unwrap().join(&relative);
    let copied = root
        .0
        .join("objects")
        .join(&first.descriptor.members[0].leaves[0].object_id)
        .join(&relative);
    assert_ne!(
        fs::metadata(&original).unwrap().ino(),
        fs::metadata(&copied).unwrap().ino()
    );
    let target = copied.clone();
    let result = with_test_hook(
        "before_snapshot_publication",
        move || {
            let moved = target.with_extension("moved");
            fs::rename(&target, &moved)?;
            symlink(moved, target)?;
            Ok(())
        },
        || {
            create_snapshot(
                &root.0,
                &[fixture.member("A")],
                Some(&first.id),
                CohortLimits::default(),
            )
        },
    );
    assert!(result.is_err());
    assert_eq!(fs::read_dir(root.0.join("snapshots")).unwrap().count(), 1);
}

fn copy_tree(source: &Path, target: &Path) {
    fs::create_dir(target).unwrap();
    for entry in fs::read_dir(source).unwrap() {
        let entry = entry.unwrap();
        if entry.file_type().unwrap().is_dir() {
            copy_tree(&entry.path(), &target.join(entry.file_name()));
        } else {
            fs::copy(entry.path(), target.join(entry.file_name())).unwrap();
        }
    }
}
