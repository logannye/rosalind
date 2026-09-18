use super::descriptor::*;
use super::{CohortError, Result};
use crate::dataset::{
    DatasetDescriptor, DatasetQuery, InputSnapshot, VerifiedEvidenceDataset,
    DATASET_DESCRIPTOR_NAME, EVIDENCE_DATASET_MANIFEST_NAME,
};
use crate::evidence::{EvidenceBatch, EvidenceCallback, EvidenceExecution};
use crate::provenance::{FileHash, RunManifest};
use std::collections::{BTreeMap, BTreeSet};
use std::fs::{self, File, OpenOptions};
use std::io::{Read, Write};
use std::path::{Component, Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

static TEMP_ID: AtomicU64 = AtomicU64::new(0);

#[derive(Debug, Clone)]
pub(crate) struct ImportMember {
    pub metadata: MemberMetadata,
    pub manifests: Vec<PathBuf>,
}

#[derive(Debug)]
pub(crate) struct SnapshotHandle {
    pub root: PathBuf,
    pub id: String,
    pub descriptor: SnapshotDescriptor,
    guards: Vec<RetainedGuard>,
}

impl SnapshotHandle {
    pub fn verify_unchanged(&self) -> Result<()> {
        verify_guards(&self.guards)
    }

    pub fn open_leaf(
        &self,
        leaf: &LeafRef,
        limits: CohortLimits,
    ) -> Result<VerifiedEvidenceDataset> {
        if !self
            .descriptor
            .members
            .iter()
            .any(|member| member.leaves.contains(leaf))
        {
            return Err(CohortError::Incompatible(
                "leaf is not part of the opened snapshot".into(),
            ));
        }
        let (dataset, _) = verify_object(&self.root, leaf, &limits, false)?;
        Ok(dataset)
    }
}

/// Copy each declared portable inventory, verify all bytes/rows, then publish
/// the snapshot last. Completed objects can remain unreferenced after failure.
pub(crate) fn create_snapshot(
    root: &Path,
    members: &[ImportMember],
    parent: Option<&str>,
    limits: CohortLimits,
) -> Result<SnapshotHandle> {
    limits.validate()?;
    let mut ids = BTreeSet::new();
    let mut leaf_count = 0usize;
    let mut input_bytes = 0usize;
    if members.len() > limits.max_members {
        return Err(CohortError::Limit("too many members".into()));
    }
    limits.admit(
        (members.len() as u64)
            .saturating_mul(8192)
            .saturating_add(128 << 10),
    )?;
    for member in members {
        member.metadata.validate()?;
        if !ids.insert(&member.metadata.id) || member.manifests.is_empty() {
            return Err(CohortError::Incompatible(
                "duplicate member ID or empty member".into(),
            ));
        }
        leaf_count = leaf_count
            .checked_add(member.manifests.len())
            .ok_or_else(|| CohortError::Limit("leaf count overflow".into()))?;
        for path in &member.manifests {
            if path.as_os_str().len() > 4096 {
                return Err(CohortError::Limit("input path exceeds 4096 bytes".into()));
            }
            input_bytes = input_bytes
                .checked_add(path.as_os_str().len() + 256)
                .ok_or_else(|| CohortError::Limit("input metadata overflow".into()))?;
        }
    }
    if leaf_count > limits.max_leaf_references {
        return Err(CohortError::Limit("too many leaf references".into()));
    }
    limits.admit(
        (input_bytes as u64)
            .saturating_mul(16)
            .saturating_add(128 << 10),
    )?;
    if parent.is_some_and(|id| !valid_hash(id)) {
        return Err(CohortError::Corrupt("invalid parent snapshot ID".into()));
    }
    let parent_guard = if let Some(id) = parent {
        let root = checked_root(root)?;
        Some(validate_ancestry(&root, id, &limits)?)
    } else {
        None
    };
    if !root.exists() {
        fs::create_dir(root)?;
    }
    let root = checked_root(root)?;
    for directory in ["objects", "snapshots"] {
        let path = root.join(directory);
        match fs::create_dir(&path) {
            Ok(()) => {}
            Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => {
                checked_directory(&path)?;
            }
            Err(e) => return Err(e.into()),
        }
    }
    let mut snapshot = SnapshotDescriptor {
        version: SNAPSHOT_VERSION,
        comparison_version: COMPARISON_VERSION,
        parent: parent.map(str::to_owned),
        members: Vec::with_capacity(members.len()),
    };
    for member in members {
        let mut leaves = Vec::with_capacity(member.manifests.len());
        for source in &member.manifests {
            leaves.push(import_object(&root, source, &limits)?);
        }
        leaves.sort_by(|a, b| a.object_id.cmp(&b.object_id));
        snapshot.members.push(SnapshotMember {
            metadata: member.metadata.clone(),
            leaves,
        });
    }
    snapshot
        .members
        .sort_by(|a, b| a.metadata.id.cmp(&b.metadata.id));
    snapshot.validate(&limits)?;
    let object_guards = validate_members(&root, &snapshot, &limits, true)?;
    let bytes = snapshot.to_bytes(&limits)?;
    let id = blake3::hash(&bytes).to_hex().to_string();
    let final_path = root.join("snapshots").join(&id);
    let staging = Staging::new(&root.join("snapshots"))?;
    write_new(&staging.path.join(SNAPSHOT_FILE), &bytes)?;
    limits.admit((bytes.len() as u64).saturating_mul(16))?;
    let receipt = snapshot_receipt(&snapshot, &id);
    let receipt_bytes = receipt.to_canonical_json();
    if receipt_bytes.len() > limits.max_snapshot_bytes.saturating_mul(4) {
        return Err(CohortError::Limit(
            "snapshot receipt exceeds envelope".into(),
        ));
    }
    write_new(
        &staging.path.join(SNAPSHOT_RECEIPT),
        receipt_bytes.as_bytes(),
    )?;
    let staged_guard = RetainedGuard::capture(
        vec![
            staging.path.join(SNAPSHOT_FILE),
            staging.path.join(SNAPSHOT_RECEIPT),
        ],
        &limits,
    )?;
    test_checkpoint("before_snapshot_publication")?;
    if let Some(guards) = parent_guard {
        verify_guards(&guards)?;
    }
    // Keep every source/object observation guarded through the final boundary.
    verify_guards(&object_guards)?;
    staged_guard.verify()?;
    limits.admit(0)?;
    match publish_directory(&staging.path, &final_path) {
        Ok(()) => {}
        Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => {
            // Idempotent publication only accepts an intact existing snapshot.
            let existing = open_snapshot(&root, &id, limits)?;
            if existing.descriptor != snapshot {
                return Err(CohortError::Corrupt(
                    "existing snapshot identity differs".into(),
                ));
            }
        }
        Err(e) => return Err(e.into()),
    }
    open_snapshot(&root, &id, limits)
}

/// Verify snapshot/ancestor metadata and each referenced leaf's metadata and
/// ownership. Payload bytes are checked by verify_snapshot, not this opener.
pub(crate) fn open_snapshot(root: &Path, id: &str, limits: CohortLimits) -> Result<SnapshotHandle> {
    limits.validate()?;
    let root = checked_root(root)?;
    let mut ancestry_guards = validate_ancestry(&root, id, &limits)?;
    let descriptor = read_snapshot(&root, id, &limits)?;
    let object_guards = validate_members(&root, &descriptor, &limits, false)?;
    verify_guards(&ancestry_guards)?;
    verify_guards(&object_guards)?;
    ancestry_guards.extend(object_guards);
    Ok(SnapshotHandle {
        root,
        id: id.into(),
        descriptor,
        guards: ancestry_guards,
    })
}

/// Fully stream and verify the current snapshot's unique leaf datasets. Parent
/// snapshot metadata is checked; unrelated ancestor payloads are not consumed.
pub(crate) fn verify_snapshot(
    root: &Path,
    id: &str,
    limits: CohortLimits,
) -> Result<SnapshotHandle> {
    let mut snapshot = open_snapshot(root, id, limits)?;
    let ancestry_guards = validate_ancestry(&snapshot.root, id, &limits)?;
    let object_guards = validate_members(&snapshot.root, &snapshot.descriptor, &limits, true)?;
    verify_guards(&ancestry_guards)?;
    verify_guards(&object_guards)?;
    snapshot.guards.extend(object_guards);
    snapshot.verify_unchanged()?;
    Ok(snapshot)
}

#[derive(Debug)]
pub(crate) struct ExtensionPublication {
    pub snapshot: SnapshotHandle,
    /// Bytes read while hashing each unique parent's declared portable inventory
    /// once. Separate parsing, new-object work and source hashing are excluded.
    pub verified_existing_bytes: u64,
}

/// Preserve every parent member/leaf, append fully verified missing-only leaves,
/// and publish an immutable child last. Existing Arrow bytes are rehashed, never
/// decoded here. The caller checks retained original-source sessions in the final
/// callback; returning an error leaves every existing snapshot untouched.
pub(crate) fn publish_extension(
    parent: &SnapshotHandle,
    additions: &[(usize, PathBuf)],
    limits: CohortLimits,
    before_publish: impl FnOnce() -> Result<()>,
) -> Result<ExtensionPublication> {
    limits.validate()?;
    parent.verify_unchanged()?;
    let root = checked_root(&parent.root)?;
    let ancestry_guards = validate_ancestry(&root, &parent.id, &limits)?;
    if read_snapshot(&root, &parent.id, &limits)? != parent.descriptor {
        return Err(CohortError::Corrupt(
            "parent handle differs from its immutable descriptor".into(),
        ));
    }
    let old_count: usize = parent
        .descriptor
        .members
        .iter()
        .map(|member| member.leaves.len())
        .sum();
    if old_count
        .checked_add(additions.len())
        .is_none_or(|count| count > limits.max_leaf_references)
    {
        return Err(CohortError::Limit(
            "extension exceeds leaf-reference envelope".into(),
        ));
    }
    for (member, source) in additions {
        if *member >= parent.descriptor.members.len() {
            return Err(CohortError::Incompatible(
                "extension names an unknown parent member".into(),
            ));
        }
        if source.as_os_str().len() > 4096 {
            return Err(CohortError::Limit(
                "extension input path exceeds 4096 bytes".into(),
            ));
        }
    }
    let parent_bytes = parent.descriptor.to_bytes(&limits)?;
    limits.admit(
        (parent_bytes.len() as u64)
            .saturating_mul(16)
            .saturating_add((additions.len() as u64).saturating_mul(16_384))
            .saturating_add((old_count as u64).saturating_mul(256)),
    )?;
    let mut guarded_objects = Vec::new();
    let mut verified_existing_bytes = 0u64;
    let mut existing = BTreeSet::new();
    for leaf in parent
        .descriptor
        .members
        .iter()
        .flat_map(|member| &member.leaves)
    {
        if existing.insert(leaf.object_id.clone()) {
            let (bytes, guards) = verify_object_bytes(&root, leaf, &limits)?;
            verified_existing_bytes =
                verified_existing_bytes.checked_add(bytes).ok_or_else(|| {
                    CohortError::Limit("existing verification byte count overflow".into())
                })?;
            guarded_objects.extend(guards);
        }
    }
    if additions.is_empty() {
        before_publish()?;
        parent.verify_unchanged()?;
        verify_guards(&ancestry_guards)?;
        verify_guards(&guarded_objects)?;
        return Ok(ExtensionPublication {
            snapshot: open_snapshot(&root, &parent.id, limits)?,
            verified_existing_bytes,
        });
    }
    let mut child = parent.descriptor.clone();
    child.parent = Some(parent.id.clone());
    let mut imported = Vec::with_capacity(additions.len());
    for (member_index, manifest) in additions {
        let leaf = import_object(&root, manifest, &limits)?;
        child.members[*member_index].leaves.push(leaf.clone());
        imported.push(leaf);
    }
    for member in &mut child.members {
        member.leaves.sort_by(|a, b| a.object_id.cmp(&b.object_id));
    }
    child.validate(&limits)?;
    // This metadata check enforces the original source key and no overlapping
    // ownership, while keeping arbitrary mixtures of adequate physical masks.
    guarded_objects.extend(validate_members(&root, &child, &limits, false)?);
    for leaf in &imported {
        // Close the gap after import's full verification and before retaining
        // this operation's guards, without decoding already validated new rows.
        let (_, guards) = verify_object_bytes(&root, leaf, &limits)?;
        guarded_objects.extend(guards);
    }
    let bytes = child.to_bytes(&limits)?;
    let id = blake3::hash(&bytes).to_hex().to_string();
    let staging = Staging::new(&root.join("snapshots"))?;
    write_new(&staging.path.join(SNAPSHOT_FILE), &bytes)?;
    limits.admit((bytes.len() as u64).saturating_mul(16))?;
    let receipt_bytes = snapshot_receipt(&child, &id).to_canonical_json();
    if receipt_bytes.len() > limits.max_snapshot_bytes.saturating_mul(4) {
        return Err(CohortError::Limit(
            "extension receipt exceeds envelope".into(),
        ));
    }
    write_new(
        &staging.path.join(SNAPSHOT_RECEIPT),
        receipt_bytes.as_bytes(),
    )?;
    let staged_guard = RetainedGuard::capture(
        vec![
            staging.path.join(SNAPSHOT_FILE),
            staging.path.join(SNAPSHOT_RECEIPT),
        ],
        &limits,
    )?;
    test_checkpoint("before_snapshot_publication")?;
    before_publish()?;
    parent.verify_unchanged()?;
    verify_guards(&ancestry_guards)?;
    verify_guards(&guarded_objects)?;
    staged_guard.verify()?;
    limits.admit(0)?;
    match publish_directory(&staging.path, &root.join("snapshots").join(&id)) {
        Ok(()) => {}
        Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => {
            let existing = open_snapshot(&root, &id, limits)?;
            if existing.descriptor != child {
                return Err(CohortError::Corrupt(
                    "existing extension snapshot differs".into(),
                ));
            }
        }
        Err(error) => return Err(error.into()),
    }
    Ok(ExtensionPublication {
        snapshot: open_snapshot(&root, &id, limits)?,
        verified_existing_bytes,
    })
}

fn verify_object_bytes(
    root: &Path,
    leaf: &LeafRef,
    limits: &CohortLimits,
) -> Result<(u64, Vec<RetainedGuard>)> {
    let (dataset, guards) = verify_object(root, leaf, limits, false)?;
    let files = inventory(&dataset, limits)?;
    let object = root.join("objects").join(&leaf.object_id);
    limits.admit(64 << 10)?;
    let mut buffer = [0u8; 64 << 10];
    let mut total = 0u64;
    for artifact in files {
        let mut file = open_read(&checked_file(&object, &artifact.path)?)?;
        let mut hash = blake3::Hasher::new();
        loop {
            limits.admit(0)?;
            let count = file.read(&mut buffer)?;
            if count == 0 {
                break;
            }
            hash.update(&buffer[..count]);
            total = total
                .checked_add(count as u64)
                .ok_or_else(|| CohortError::Limit("verification byte count overflow".into()))?;
        }
        if hash.finalize().to_hex().as_str() != artifact.blake3 {
            return Err(CohortError::Corrupt(
                "existing cohort inventory hash differs".into(),
            ));
        }
    }
    dataset.verify_unchanged()?;
    verify_guards(&guards)?;
    Ok((total, guards))
}

fn import_object(root: &Path, manifest: &Path, limits: &CohortLimits) -> Result<LeafRef> {
    if manifest.file_name().and_then(|name| name.to_str()) != Some(EVIDENCE_DATASET_MANIFEST_NAME) {
        return Err(CohortError::Incompatible(
            "import requires evidence-dataset.manifest.json".into(),
        ));
    }
    let source_root = checked_root(
        manifest
            .parent()
            .ok_or_else(|| CohortError::Corrupt("manifest has no parent".into()))?,
    )?;
    checked_file(&source_root, EVIDENCE_DATASET_MANIFEST_NAME)?;
    checked_file(&source_root, DATASET_DESCRIPTOR_NAME)?;
    let dataset = VerifiedEvidenceDataset::open(
        source_root.join(EVIDENCE_DATASET_MANIFEST_NAME),
        limits.dataset_limits(),
    )?;
    let inventory = inventory(&dataset, limits)?;
    let leaf = LeafRef {
        object_id: inventory[0].blake3.clone(),
        descriptor_blake3: inventory[1].blake3.clone(),
    };
    let destination = root.join("objects").join(&leaf.object_id);
    if fs::symlink_metadata(&destination).is_ok() {
        verify_object(root, &leaf, limits, true)?;
        return Ok(leaf);
    }
    limits.admit(
        (inventory.len() as u64)
            .saturating_mul(4096)
            .saturating_add(128 << 10),
    )?;
    let source_paths = inventory
        .iter()
        .map(|file| checked_file(&source_root, &file.path))
        .collect::<Result<Vec<_>>>()?;
    let guard = InputSnapshot::capture(source_paths)?;
    let staging = Staging::new(&root.join("objects"))?;
    for file in &inventory {
        let source = checked_file(&source_root, &file.path)?;
        let target = staging.path.join(&file.path);
        if let Some(parent) = target.parent() {
            fs::create_dir_all(parent)?;
        }
        copy_verified(&source, &target, &file.blake3, limits)?;
    }
    guard.verify()?;
    dataset.verify_unchanged()?;
    for file in &inventory {
        checked_file(&source_root, &file.path)?;
    }
    let (_, staged_guards) = verify_object_at(&staging.path, &leaf, limits, true)?;
    guard.verify()?;
    limits.admit(0)?;
    test_checkpoint("before_object_publication")?;
    guard.verify()?;
    verify_guards(&staged_guards)?;
    match publish_directory(&staging.path, &destination) {
        Ok(()) => {}
        Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => {
            verify_object(root, &leaf, limits, true)?;
        }
        Err(e) => return Err(e.into()),
    }
    Ok(leaf)
}

fn inventory(dataset: &VerifiedEvidenceDataset, limits: &CohortLimits) -> Result<Vec<FileHash>> {
    let count = dataset
        .descriptor()
        .partitions
        .len()
        .saturating_mul(2)
        .saturating_add(2);
    limits.admit((count as u64).saturating_mul(8192))?;
    let hashes = dataset.source_hashes();
    let mut files = vec![
        FileHash {
            path: EVIDENCE_DATASET_MANIFEST_NAME.into(),
            blake3: hashes[0].blake3.clone(),
        },
        FileHash {
            path: DATASET_DESCRIPTOR_NAME.into(),
            blake3: hashes[1].blake3.clone(),
        },
    ];
    for part in &dataset.descriptor().partitions {
        for file in [&part.receipt, &part.arrow] {
            files.push(FileHash {
                path: file.path.clone(),
                blake3: file.blake3.clone(),
            });
        }
    }
    Ok(files)
}

fn verify_object(
    root: &Path,
    leaf: &LeafRef,
    limits: &CohortLimits,
    payload: bool,
) -> Result<(VerifiedEvidenceDataset, Vec<RetainedGuard>)> {
    checked_directory(&root.join("objects"))?;
    verify_object_at(
        &root.join("objects").join(&leaf.object_id),
        leaf,
        limits,
        payload,
    )
}

fn verify_object_at(
    path: &Path,
    leaf: &LeafRef,
    limits: &CohortLimits,
    payload: bool,
) -> Result<(VerifiedEvidenceDataset, Vec<RetainedGuard>)> {
    checked_directory(path)?;
    checked_file(path, EVIDENCE_DATASET_MANIFEST_NAME)?;
    checked_file(path, DATASET_DESCRIPTOR_NAME)?;
    let metadata_guard = RetainedGuard::capture(
        vec![
            path.join(EVIDENCE_DATASET_MANIFEST_NAME),
            path.join(DATASET_DESCRIPTOR_NAME),
        ],
        limits,
    )?;
    let mut dataset = VerifiedEvidenceDataset::open(
        path.join(EVIDENCE_DATASET_MANIFEST_NAME),
        limits.dataset_limits(),
    )?;
    if dataset.source_hashes()[0].blake3 != leaf.object_id
        || dataset.source_hashes()[1].blake3 != leaf.descriptor_blake3
    {
        return Err(CohortError::Corrupt(
            "object manifest or descriptor hash differs from snapshot".into(),
        ));
    }
    let files = inventory(&dataset, limits)?;
    let mut expected = BTreeSet::new();
    for file in &files {
        checked_file(path, &file.path)?;
        expected.insert(file.path.clone());
    }
    let (actual, directories) = directory_inventory(path, files.len())?;
    let expected_directories = declared_directories(&expected);
    if actual != expected || directories != expected_directories {
        return Err(CohortError::Corrupt(
            "object contains undeclared files".into(),
        ));
    }
    let inventory_guard = RetainedGuard::capture(
        files.iter().map(|file| path.join(&file.path)).collect(),
        limits,
    )?;
    metadata_guard.verify()?;
    if payload {
        test_checkpoint("before_arrow_decode")?;
        let query = DatasetQuery {
            selection: dataset.selection(),
            fields: dataset.fields(),
        };
        let mut callback =
            EvidenceCallback::with_fields(|_: &EvidenceBatch| Ok(()), 0, query.fields);
        let execution = EvidenceExecution {
            memory_budget_bytes: limits.memory_budget_bytes,
            ..EvidenceExecution::default()
        };
        dataset.visit_batches(&query, &mut callback, &execution)?;
        dataset.verify_unchanged()?;
    }
    metadata_guard.verify()?;
    inventory_guard.verify()?;
    Ok((dataset, vec![metadata_guard, inventory_guard]))
}

fn validate_members(
    root: &Path,
    snapshot: &SnapshotDescriptor,
    limits: &CohortLimits,
    payload: bool,
) -> Result<Vec<RetainedGuard>> {
    let leaf_count: usize = snapshot.members.iter().map(|m| m.leaves.len()).sum();
    limits.admit(
        (snapshot.members.len() as u64)
            .saturating_mul(2048)
            .saturating_add((leaf_count as u64).saturating_mul(512)),
    )?;
    let mut retained_guards = Vec::new();
    let mut aliases = BTreeMap::<String, String>::new();
    let mut verified = BTreeSet::new();
    for member in &snapshot.members {
        let mut compatibility = None;
        let mut source_alias = None;
        let mut ownership = Vec::new();
        for leaf in &member.leaves {
            let (dataset, guards) = verify_object(
                root,
                leaf,
                limits,
                payload && verified.insert(&leaf.object_id),
            )?;
            retained_guards.extend(guards);
            let descriptor = dataset.descriptor();
            if compatibility
                .as_ref()
                .is_some_and(|key| key != &descriptor.compatibility_blake3)
            {
                return Err(CohortError::Incompatible(format!(
                    "member {} leaves do not share source-reuse identity",
                    member.metadata.id
                )));
            }
            compatibility = Some(descriptor.compatibility_blake3.clone());
            source_alias = Some(source_scope_identity(descriptor, limits)?);
            let count = ownership
                .len()
                .checked_add(descriptor.selection.intervals.len())
                .ok_or_else(|| CohortError::Limit("ownership interval count overflow".into()))?;
            if count > limits.max_ownership_intervals {
                return Err(CohortError::Limit(
                    "member ownership exceeds interval envelope".into(),
                ));
            }
            limits.admit((count as u64).saturating_mul(64))?;
            ownership.extend(
                descriptor
                    .selection
                    .intervals
                    .iter()
                    .map(|i| (i.contig, i.start, i.end)),
            );
            dataset.verify_unchanged()?;
        }
        ownership.sort_unstable();
        if ownership
            .windows(2)
            .any(|pair| pair[0].0 == pair[1].0 && pair[0].2 > pair[1].1)
        {
            return Err(CohortError::Incompatible(format!(
                "member {} has overlapping leaf ownership",
                member.metadata.id
            )));
        }
        if let Some(previous) = aliases.insert(
            source_alias.expect("validated nonempty member"),
            member.metadata.id.clone(),
        ) {
            return Err(CohortError::Incompatible(format!(
                "members {previous} and {} alias the same alignment and sample scope",
                member.metadata.id
            )));
        }
    }
    verify_guards(&retained_guards)?;
    Ok(retained_guards)
}

fn source_scope_identity(descriptor: &DatasetDescriptor, limits: &CohortLimits) -> Result<String> {
    let alignment = descriptor
        .sources
        .iter()
        .find(|source| source.role == "alignments")
        .ok_or_else(|| CohortError::Corrupt("dataset has no alignment identity".into()))?;
    let scope_bytes = descriptor
        .sample_scope
        .read_groups
        .iter()
        .fold(1024u64, |sum, group| {
            sum.saturating_add(group.id.len() as u64)
                .saturating_add(group.sample.as_ref().map_or(0, |name| name.len()) as u64)
                .saturating_add(128)
        });
    limits.admit(scope_bytes.saturating_mul(16))?;
    let scope = descriptor.sample_scope.canonical_json()?;
    Ok(blake3::hash(
        format!(
            "cohort-source-scope-v1\n{}\n{}\n{scope}",
            alignment.blake3, alignment.bytes
        )
        .as_bytes(),
    )
    .to_hex()
    .to_string())
}

fn validate_ancestry(root: &Path, id: &str, limits: &CohortLimits) -> Result<Vec<RetainedGuard>> {
    let mut current = Some(id.to_owned());
    let mut seen = BTreeSet::new();
    let mut guards = Vec::new();
    while let Some(id) = current {
        if !seen.insert(id.clone()) || seen.len() > limits.max_parent_depth {
            return Err(CohortError::Limit(
                "snapshot parent cycle or depth envelope exceeded".into(),
            ));
        }
        if !valid_hash(&id) {
            return Err(CohortError::Corrupt("invalid ancestor ID".into()));
        }
        limits.admit((seen.len() as u64).saturating_mul(16384))?;
        let directory = root.join("snapshots").join(&id);
        checked_directory(&directory)?;
        let paths = [
            checked_file(&directory, SNAPSHOT_FILE)?,
            checked_file(&directory, SNAPSHOT_RECEIPT)?,
        ];
        let guard = RetainedGuard::capture(paths.to_vec(), limits)?;
        let descriptor = read_snapshot(root, &id, limits)?;
        guard.verify()?;
        guards.push(guard);
        current = descriptor.parent;
    }
    verify_guards(&guards)?;
    Ok(guards)
}

fn read_snapshot(root: &Path, id: &str, limits: &CohortLimits) -> Result<SnapshotDescriptor> {
    if !valid_hash(id) {
        return Err(CohortError::Corrupt("invalid snapshot ID".into()));
    }
    checked_directory(&root.join("snapshots"))?;
    let path = root.join("snapshots").join(id);
    checked_directory(&path)?;
    let bytes = read_bounded(
        &checked_file(&path, SNAPSHOT_FILE)?,
        limits.max_snapshot_bytes,
        limits,
    )?;
    if blake3::hash(&bytes).to_hex().as_str() != id {
        return Err(CohortError::Corrupt(
            "snapshot content hash differs from its ID".into(),
        ));
    }
    let descriptor = SnapshotDescriptor::from_bytes(&bytes, limits)?;
    let receipt_bytes = read_bounded(
        &checked_file(&path, SNAPSHOT_RECEIPT)?,
        limits.max_snapshot_bytes.saturating_mul(4),
        limits,
    )?;
    let text =
        std::str::from_utf8(&receipt_bytes).map_err(|e| CohortError::Corrupt(e.to_string()))?;
    let receipt =
        RunManifest::from_canonical_json(text).map_err(|e| CohortError::Corrupt(e.to_string()))?;
    let expected = snapshot_receipt(&descriptor, id);
    if receipt.self_hash_ok() != Some(true)
        || receipt.measurement_hash_ok() == Some(false)
        || (receipt.claims_measurements() && receipt.measurement_hash_ok() != Some(true))
        || receipt.subcommand != expected.subcommand
        // Producer/toolchain provenance belongs to the recorded producer, not
        // the reader's build. Only the versioned cohort lineage must agree with
        // this immutable descriptor; the receipt self-hash binds the
        // original claim without rewriting it for the current executable.
        || ["cohort.snapshot_version", "cohort.snapshot_blake3", "run_status"]
            .iter()
            .any(|key| receipt.params.get(*key) != expected.params.get(*key))
        || inventory_map(&receipt.inputs)? != inventory_map(&expected.inputs)?
        || inventory_map(&receipt.outputs)? != inventory_map(&expected.outputs)?
    {
        return Err(CohortError::Corrupt(
            "snapshot receipt lineage differs".into(),
        ));
    }
    if directory_inventory(&path, 2)?
        != (
            BTreeSet::from([SNAPSHOT_FILE.into(), SNAPSHOT_RECEIPT.into()]),
            BTreeSet::new(),
        )
    {
        return Err(CohortError::Corrupt(
            "snapshot contains undeclared files".into(),
        ));
    }
    Ok(descriptor)
}

fn snapshot_receipt(snapshot: &SnapshotDescriptor, id: &str) -> RunManifest {
    let mut receipt = RunManifest::new("cohort snapshot");
    receipt.params.extend([
        (
            "cohort.snapshot_version".into(),
            SNAPSHOT_VERSION.to_string(),
        ),
        ("cohort.snapshot_blake3".into(), id.into()),
        ("run_status".into(), "completed".into()),
    ]);
    let mut inputs = BTreeMap::new();
    for leaf in snapshot.members.iter().flat_map(|member| &member.leaves) {
        inputs.insert(
            format!(
                "objects/{}/{}",
                leaf.object_id, EVIDENCE_DATASET_MANIFEST_NAME
            ),
            leaf.object_id.clone(),
        );
        inputs.insert(
            format!("objects/{}/{}", leaf.object_id, DATASET_DESCRIPTOR_NAME),
            leaf.descriptor_blake3.clone(),
        );
    }
    if let Some(parent) = &snapshot.parent {
        inputs.insert(
            format!("snapshots/{parent}/{SNAPSHOT_FILE}"),
            parent.clone(),
        );
    }
    receipt.inputs = inputs
        .into_iter()
        .map(|(path, blake3)| FileHash { path, blake3 })
        .collect();
    receipt.outputs.push(FileHash {
        path: SNAPSHOT_FILE.into(),
        blake3: id.into(),
    });
    receipt.finalize();
    receipt
}

fn inventory_map(files: &[FileHash]) -> Result<BTreeMap<&str, &str>> {
    let mut map = BTreeMap::new();
    for file in files {
        if map
            .insert(file.path.as_str(), file.blake3.as_str())
            .is_some()
        {
            return Err(CohortError::Corrupt(
                "duplicate receipt inventory path".into(),
            ));
        }
    }
    Ok(map)
}

fn checked_root(path: &Path) -> Result<PathBuf> {
    checked_directory(path)?;
    Ok(fs::canonicalize(path)?)
}
fn checked_directory(path: &Path) -> Result<()> {
    let metadata = fs::symlink_metadata(path)?;
    if !metadata.is_dir() || metadata.file_type().is_symlink() {
        return Err(CohortError::Corrupt(format!(
            "expected nonsymlink directory: {}",
            path.display()
        )));
    }
    Ok(())
}
fn checked_file(root: &Path, relative: &str) -> Result<PathBuf> {
    checked_directory(root)?;
    let components: Vec<_> = Path::new(relative).components().collect();
    if components.is_empty()
        || components
            .iter()
            .any(|c| !matches!(c, Component::Normal(_)))
    {
        return Err(CohortError::Corrupt("nonlocal inventory path".into()));
    }
    let mut path = root.to_path_buf();
    for (index, component) in components.iter().enumerate() {
        path.push(component.as_os_str());
        let metadata = fs::symlink_metadata(&path)?;
        if metadata.file_type().is_symlink()
            || (index + 1 == components.len() && !metadata.is_file())
            || (index + 1 != components.len() && !metadata.is_dir())
        {
            return Err(CohortError::Corrupt(
                "inventory contains a symlink or nonregular file".into(),
            ));
        }
    }
    if path.as_os_str().len() > 4096 {
        return Err(CohortError::Limit(
            "inventory path exceeds 4096 bytes".into(),
        ));
    }
    Ok(path)
}

fn directory_inventory(
    root: &Path,
    maximum: usize,
) -> Result<(BTreeSet<String>, BTreeSet<String>)> {
    let mut files = BTreeSet::new();
    let mut found_directories = BTreeSet::new();
    let mut directories = vec![PathBuf::new()];
    let mut entries = 0usize;
    while let Some(relative) = directories.pop() {
        for entry in fs::read_dir(root.join(&relative))? {
            let entry = entry?;
            entries = entries
                .checked_add(1)
                .ok_or_else(|| CohortError::Limit("inventory count overflow".into()))?;
            if entries > maximum.saturating_mul(2).saturating_add(2) {
                return Err(CohortError::Limit(
                    "undeclared object inventory exceeds envelope".into(),
                ));
            }
            let path = relative.join(entry.file_name());
            let kind = entry.file_type()?;
            if kind.is_dir() && relative.as_os_str().is_empty() {
                found_directories.insert(
                    path.to_str()
                        .ok_or_else(|| {
                            CohortError::Corrupt("non-UTF-8 inventory directory".into())
                        })?
                        .to_owned(),
                );
                directories.push(path);
            } else if kind.is_file() {
                files.insert(
                    path.to_str()
                        .ok_or_else(|| CohortError::Corrupt("non-UTF-8 inventory path".into()))?
                        .to_owned(),
                );
            } else {
                return Err(CohortError::Corrupt(
                    "symlink, nested directory or special file in object".into(),
                ));
            }
        }
    }
    Ok((files, found_directories))
}

fn declared_directories(files: &BTreeSet<String>) -> BTreeSet<String> {
    files
        .iter()
        .filter_map(|file| Path::new(file).parent())
        .filter(|path| !path.as_os_str().is_empty())
        .map(|path| path.to_str().expect("validated UTF-8 path").to_owned())
        .collect()
}

fn verify_guards(guards: &[RetainedGuard]) -> Result<()> {
    for guard in guards {
        guard.verify()?;
    }
    Ok(())
}

fn read_bounded(path: &Path, maximum: usize, limits: &CohortLimits) -> Result<Vec<u8>> {
    let mut file = open_read(path)?;
    let size = file.metadata()?.len();
    if size > maximum as u64 {
        return Err(CohortError::Limit("metadata file exceeds envelope".into()));
    }
    limits.admit(size.saturating_mul(16))?;
    let mut bytes = Vec::with_capacity(size as usize);
    Read::by_ref(&mut file)
        .take(size.saturating_add(1))
        .read_to_end(&mut bytes)?;
    if bytes.len() as u64 != size {
        return Err(CohortError::Corrupt(
            "metadata changed while reading".into(),
        ));
    }
    Ok(bytes)
}

fn open_read(path: &Path) -> std::io::Result<File> {
    let mut options = OpenOptions::new();
    options.read(true);
    #[cfg(unix)]
    {
        use std::os::unix::fs::OpenOptionsExt;
        options.custom_flags(libc::O_NOFOLLOW);
    }
    options.open(path)
}

fn copy_verified(
    source: &Path,
    target: &Path,
    expected: &str,
    limits: &CohortLimits,
) -> Result<()> {
    let mut input = open_read(source)?;
    let mut output = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(target)?;
    let mut hash = blake3::Hasher::new();
    let mut buffer = [0u8; 64 << 10];
    loop {
        limits.admit(0)?;
        let count = input.read(&mut buffer)?;
        if count == 0 {
            break;
        }
        hash.update(&buffer[..count]);
        output.write_all(&buffer[..count])?;
    }
    if hash.finalize().to_hex().as_str() != expected {
        return Err(CohortError::Corrupt("copied inventory hash differs".into()));
    }
    output.sync_all()?;
    Ok(())
}
fn write_new(path: &Path, bytes: &[u8]) -> std::io::Result<()> {
    let mut file = OpenOptions::new().write(true).create_new(true).open(path)?;
    file.write_all(bytes)?;
    file.sync_all()
}

#[derive(Debug)]
struct Staging {
    path: PathBuf,
}
impl Staging {
    fn new(parent: &Path) -> Result<Self> {
        checked_directory(parent)?;
        for _ in 0..128 {
            let path = parent.join(format!(
                ".staging-{}-{}",
                std::process::id(),
                TEMP_ID.fetch_add(1, Ordering::Relaxed)
            ));
            match fs::create_dir(&path) {
                Ok(()) => return Ok(Self { path }),
                Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => continue,
                Err(e) => return Err(e.into()),
            }
        }
        Err(CohortError::Incompatible(
            "cannot reserve cohort staging directory".into(),
        ))
    }
}
impl Drop for Staging {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.path);
    }
}

/// Same-filesystem, create-new directory publication. No replacement fallback.
fn publish_directory(source: &Path, destination: &Path) -> std::io::Result<()> {
    #[cfg(any(target_os = "linux", target_os = "macos"))]
    {
        use std::ffi::CString;
        use std::os::unix::ffi::OsStrExt;
        let source = CString::new(source.as_os_str().as_bytes())
            .map_err(|_| std::io::Error::other("NUL directory path"))?;
        let destination = CString::new(destination.as_os_str().as_bytes())
            .map_err(|_| std::io::Error::other("NUL directory path"))?;
        // SAFETY: both NUL-terminated paths live through this synchronous syscall;
        // the platform flag atomically refuses every existing destination.
        let result = unsafe {
            #[cfg(target_os = "linux")]
            {
                // Call the kernel directly: libc's renameat2 wrapper would add
                // a GLIBC_2.28 symbol requirement to older supported bundles.
                libc::syscall(
                    libc::SYS_renameat2,
                    libc::AT_FDCWD,
                    source.as_ptr(),
                    libc::AT_FDCWD,
                    destination.as_ptr(),
                    libc::RENAME_NOREPLACE,
                )
            }
            #[cfg(target_os = "macos")]
            {
                libc::renamex_np(source.as_ptr(), destination.as_ptr(), libc::RENAME_EXCL)
            }
        };
        if result != 0 {
            return Err(std::io::Error::last_os_error());
        }
        Ok(())
    }
    #[cfg(not(any(target_os = "linux", target_os = "macos")))]
    {
        let _ = (source, destination);
        Err(std::io::Error::new(
            std::io::ErrorKind::Unsupported,
            "atomic cohort publication requires Linux or macOS",
        ))
    }
}

#[derive(Debug)]
struct RetainedGuard {
    snapshot: InputSnapshot,
    paths: Vec<PathBuf>,
    directories: Vec<DirectoryStamp>,
}
impl RetainedGuard {
    fn capture(paths: Vec<PathBuf>, limits: &CohortLimits) -> Result<Self> {
        // Inventory membership is immutable too: guard each declared file's
        // containing directory so late extra entries cannot evade file guards.
        limits.admit((paths.len() as u64).saturating_mul(3 * 16_384))?;
        let directory_paths: BTreeSet<_> = paths.iter().filter_map(|path| path.parent()).collect();
        let directories = directory_paths
            .into_iter()
            .map(DirectoryStamp::read)
            .collect::<Result<Vec<_>>>()?;
        let snapshot = InputSnapshot::capture(paths.iter().cloned())?;
        let guard = Self {
            snapshot,
            paths,
            directories,
        };
        guard.verify()?;
        Ok(guard)
    }
    fn verify(&self) -> Result<()> {
        for path in &self.paths {
            let mut current = PathBuf::new();
            for component in path.components() {
                current.push(component.as_os_str());
                if fs::symlink_metadata(&current)?.file_type().is_symlink() {
                    return Err(CohortError::Corrupt(
                        "retained inventory path became a symlink".into(),
                    ));
                }
            }
        }
        for directory in &self.directories {
            if DirectoryStamp::read(&directory.path)? != *directory {
                return Err(CohortError::Corrupt(
                    "declared inventory directory changed during the operation".into(),
                ));
            }
        }
        self.snapshot.verify()?;
        Ok(())
    }
}

#[derive(Debug, PartialEq, Eq)]
struct DirectoryStamp {
    path: PathBuf,
    modified: Option<std::time::SystemTime>,
    bytes: u64,
    #[cfg(unix)]
    identity: (u64, u64, i64, i64, i64, i64),
}
impl DirectoryStamp {
    fn read(path: &Path) -> Result<Self> {
        let metadata = fs::symlink_metadata(path)?;
        if !metadata.is_dir() || metadata.file_type().is_symlink() {
            return Err(CohortError::Corrupt(
                "inventory parent is no longer a regular directory".into(),
            ));
        }
        Ok(Self {
            path: path.to_path_buf(),
            modified: metadata.modified().ok(),
            bytes: metadata.len(),
            #[cfg(unix)]
            identity: {
                use std::os::unix::fs::MetadataExt;
                (
                    metadata.dev(),
                    metadata.ino(),
                    metadata.mtime(),
                    metadata.mtime_nsec(),
                    metadata.ctime(),
                    metadata.ctime_nsec(),
                )
            },
        })
    }
}

#[cfg(not(test))]
fn test_checkpoint(_: &str) -> Result<()> {
    Ok(())
}

#[cfg(test)]
type TestHook = (&'static str, Box<dyn FnOnce() -> Result<()>>);
#[cfg(test)]
thread_local! { static TEST_HOOK: std::cell::RefCell<Option<TestHook>> = std::cell::RefCell::new(None); }
#[cfg(test)]
fn test_checkpoint(stage: &str) -> Result<()> {
    let hook = TEST_HOOK.with(|slot| {
        let mut slot = slot.borrow_mut();
        if slot.as_ref().is_some_and(|(name, _)| *name == stage) {
            slot.take()
        } else {
            None
        }
    });
    if let Some((_, hook)) = hook {
        hook()?;
    }
    Ok(())
}
#[cfg(test)]
pub(crate) fn with_test_hook<T>(
    stage: &'static str,
    hook: impl FnOnce() -> Result<()> + 'static,
    run: impl FnOnce() -> T,
) -> T {
    struct Reset;
    impl Drop for Reset {
        fn drop(&mut self) {
            TEST_HOOK.with(|slot| {
                slot.borrow_mut().take();
            });
        }
    }
    TEST_HOOK.with(|slot| {
        assert!(slot.borrow().is_none());
        *slot.borrow_mut() = Some((stage, Box::new(hook)));
    });
    let _reset = Reset;
    run()
}
