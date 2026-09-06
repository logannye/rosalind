//! Portable dataset metadata and a source-hashed publication boundary.
//!
//! The wire structs describe scientific metadata, not trusted executions. Only
//! VerifiedInputSession can mint the new source-bound cache namespace. Publishing
//! requires that namespace and rechecks every partition and source snapshot.

use super::{DatasetError, DatasetOutcome, InputSnapshot};
use crate::core::ContigSet;
use crate::evidence::{
    self, EvidenceEngine, EvidenceError, EvidenceFields, EvidenceProfile, EvidenceReadGroup,
    EvidenceSampleMode, EvidenceSampleScope, EvidenceSelection, SnvSite, CANONICAL_TILE_BASES,
};
use crate::provenance::{FileHash, RunManifest};
use crate::selection::GenomicInterval;
use crate::util::atomic::{commit_group, AtomicFile};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::fs::{self, File};
use std::io::{BufReader, Read, Write};
use std::path::{Path, PathBuf};

/// Portable descriptor filename, relative to the dataset root.
pub const DATASET_DESCRIPTOR_NAME: &str = "dataset.descriptor.json";
/// Receipt binding the portable descriptor and every partition artifact.
pub const EVIDENCE_DATASET_MANIFEST_NAME: &str = "evidence-dataset.manifest.json";

/// Bounded metadata encoding and parsing controls.
#[derive(Debug, Clone, Copy)]
pub struct DescriptorLimits {
    /// Maximum serialized descriptor or receipt size; defaults to 32 MiB.
    pub max_bytes: usize,
}
impl Default for DescriptorLimits {
    fn default() -> Self {
        Self {
            max_bytes: 32 << 20,
        }
    }
}

/// Source identity verified once by a live session; paths are provenance only.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SourceIdentity {
    /// Stable input role, for example alignments, reference, or alignment-index.
    pub role: String,
    /// Original canonical path; offline readers do not need or rehash this path.
    pub path: String,
    /// BLAKE3 of the original file bytes.
    pub blake3: String,
    /// Original file size in bytes.
    pub bytes: u64,
}

/// A single immutable input session with private verified content identities.
#[derive(Debug, Clone)]
pub struct VerifiedInputSession {
    identities: Vec<SourceIdentity>,
    snapshot: InputSnapshot,
    access_paths: Vec<(PathBuf, PathBuf)>,
}
impl VerifiedInputSession {
    /// Hash each unique local file once, bracketed by metadata/inode guards.
    /// Roles must be unique. All inputs must stay immutable throughout the run.
    pub fn open(inputs: impl IntoIterator<Item = (String, PathBuf)>) -> Result<Self, DatasetError> {
        let mut paths = BTreeMap::new();
        let mut access_paths = Vec::new();
        for (role, path) in inputs {
            if paths.len() >= 256 || !valid_role(&role) || paths.contains_key(&role) {
                return Err(incompatible(
                    "source roles must be unique, bounded ASCII names (at most 256)",
                ));
            }
            let access = if path.is_absolute() {
                path
            } else {
                std::env::current_dir()?.join(path)
            };
            let canonical = fs::canonicalize(&access)?;
            access_paths.push((access, canonical.clone()));
            paths.insert(role, canonical);
        }
        if paths
            .iter()
            .map(|(role, path)| role.len() + path.as_os_str().len())
            .sum::<usize>()
            .saturating_add(
                access_paths
                    .iter()
                    .map(|(access, canonical)| {
                        access
                            .as_os_str()
                            .len()
                            .saturating_add(canonical.as_os_str().len())
                    })
                    .sum::<usize>(),
            )
            > 1 << 20
        {
            return Err(incompatible("verified source metadata exceeds 1 MiB"));
        }
        let snapshot = InputSnapshot::capture(
            access_paths
                .iter()
                .flat_map(|(access, canonical)| [access.clone(), canonical.clone()]),
        )?;
        let mut hashes = BTreeMap::<PathBuf, (String, u64)>::new();
        let mut identities = Vec::with_capacity(paths.len());
        for (role, path) in paths {
            let (blake3, bytes) = match hashes.get(&path) {
                Some(identity) => identity.clone(),
                None => {
                    let identity = (hash_file(&path)?, fs::metadata(&path)?.len());
                    hashes.insert(path.clone(), identity.clone());
                    identity
                }
            };
            let path = path
                .to_str()
                .ok_or_else(|| incompatible("source path is not UTF-8"))?
                .to_owned();
            identities.push(SourceIdentity {
                role,
                path,
                blake3,
                bytes,
            });
        }
        let session = Self {
            identities,
            snapshot,
            access_paths,
        };
        session.verify()?;
        Ok(session)
    }
    /// Read verified source identities without allowing callers to replace them.
    pub fn identities(&self) -> &[SourceIdentity] {
        &self.identities
    }
    /// Copy the role-to-content-hash mapping for receipt construction.
    pub fn hashes(&self) -> BTreeMap<String, String> {
        self.identities
            .iter()
            .map(|s| (s.role.clone(), s.blake3.clone()))
            .collect()
    }
    /// The guard captured before hashing; use it for worker execution.
    pub fn snapshot(&self) -> &InputSnapshot {
        &self.snapshot
    }
    /// Refuse ordinary mutations since this session was opened.
    pub fn verify(&self) -> Result<(), DatasetError> {
        self.snapshot.verify()?;
        for (access, target) in &self.access_paths {
            if fs::canonicalize(access)? != *target {
                return Err(incompatible(
                    "input access path was retargeted during the verified session",
                ));
            }
        }
        Ok(())
    }
    /// Compatibility identity excluding selection, projection, execution and code.
    pub fn compatibility_key(&self, engine: &EvidenceEngine) -> Result<String, DatasetError> {
        self.assert_engine(engine)?;
        engine_compatibility(engine, &self.identities)
    }
    /// New source-bound namespace for run_dataset_with_snapshot. Legacy arbitrary
    /// caller namespaces cannot be promoted by the portable publisher.
    pub fn dataset_namespace(&self, engine: &EvidenceEngine) -> Result<String, DatasetError> {
        self.assert_engine(engine)?;
        let compatibility = engine_compatibility(engine, &self.identities)?;
        Ok(namespace(&compatibility, &super::request_digest(engine)))
    }
    fn assert_engine(&self, engine: &EvidenceEngine) -> Result<(), DatasetError> {
        self.verify()?;
        let request = engine.request();
        for (role, path) in [
            ("alignments", Some(&request.alignments)),
            ("reference", request.reference.as_ref()),
            ("alignment-index", request.alignment_index.as_ref()),
            ("reference-fai", request.reference_fai.as_ref()),
            ("cram-reference", request.cram_reference.as_ref()),
            ("cram-reference-fai", request.cram_reference_fai.as_ref()),
        ] {
            if let Some(path) = path {
                let path = fs::canonicalize(path)?;
                if !self
                    .identities
                    .iter()
                    .any(|source| source.role == role && Path::new(&source.path) == path)
                {
                    return Err(incompatible(&format!(
                        "verified input session did not capture engine {role}"
                    )));
                }
            }
        }
        // Include auto-discovered sidecars as well. Explicit paths are preferred
        // when a caller wants only one index to participate in compatibility.
        let actual = InputSnapshot::for_engine(engine)?;
        for (path, _) in actual.files.iter() {
            let path = fs::canonicalize(path)?;
            let expected_role = if path == fs::canonicalize(&request.alignments)? {
                "alignments"
            } else if request.reference.as_ref().is_some_and(|reference| {
                fs::canonicalize(format!("{}.fai", reference.display()))
                    .is_ok_and(|fai| fai == path)
            }) {
                "reference-fai"
            } else if request.cram_reference.as_ref().is_some_and(|reference| {
                fs::canonicalize(format!("{}.fai", reference.display()))
                    .is_ok_and(|fai| fai == path)
            }) {
                "cram-reference-fai"
            } else {
                // Explicit engine dependencies were checked by exact role above.
                // Remaining discovered sidecars must participate in scientific
                // compatibility; a reserved selection role must never hide one.
                if [
                    request.reference.as_ref(),
                    request.reference_fai.as_ref(),
                    request.cram_reference.as_ref(),
                    request.cram_reference_fai.as_ref(),
                    request.alignment_index.as_ref(),
                ]
                .into_iter()
                .flatten()
                .any(|explicit| fs::canonicalize(explicit).is_ok_and(|p| p == path))
                {
                    continue;
                }
                "alignment-index"
            };
            if !self
                .identities
                .iter()
                .any(|source| source.role == expected_role && Path::new(&source.path) == path)
            {
                return Err(incompatible("verified input session omitted an engine dependency under its canonical source role; choose an explicit index when multiple sidecars exist"));
            }
        }
        Ok(())
    }
}

/// Ordered reference dictionary entry.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DescriptorContig {
    /// Contiguous zero-based dictionary identifier.
    pub id: u32,
    /// Exact reference contig name.
    pub name: String,
    /// Reference length in bases.
    pub length: u32,
}
/// Half-open normalized genomic interval.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DescriptorInterval {
    /// Dictionary contig identifier.
    pub contig: u32,
    /// Inclusive zero-based start.
    pub start: u32,
    /// Exclusive zero-based end.
    pub end: u32,
}
/// Canonical SNV selection annotation.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DescriptorSite {
    /// Dictionary contig identifier.
    pub contig: u32,
    /// Zero-based coordinate.
    pub position: u32,
    /// ASCII A/C/G/T REF byte.
    pub reference: u8,
    /// Sorted, distinct ASCII A/C/G/T ALT bytes.
    pub alternates: Vec<u8>,
}
/// Canonical selection, including full denominators and optional SNV annotations.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DescriptorSelection {
    /// Ordered disjoint intervals; adjacent intervals have already been merged.
    pub intervals: Vec<DescriptorInterval>,
    /// Some for SNV selection; None for BED/whole-reference coverage selection.
    pub sites: Option<Vec<DescriptorSite>>,
}
impl DescriptorSelection {
    /// Convert validated metadata to the engine's scientific selection type.
    pub fn to_evidence_selection(&self) -> EvidenceSelection {
        match &self.sites {
            Some(sites) => EvidenceSelection::Sites(
                sites
                    .iter()
                    .map(|s| SnvSite {
                        contig: s.contig,
                        position: s.position,
                        reference: s.reference,
                        alternates: s.alternates.clone(),
                    })
                    .collect(),
            ),
            None => EvidenceSelection::Intervals(
                self.intervals
                    .iter()
                    .map(|i| GenomicInterval {
                        contig: i.contig,
                        start: i.start,
                        end: i.end,
                    })
                    .collect(),
            ),
        }
    }
}
/// Exact read-filter parameters; execution limits never appear here.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DescriptorProfile {
    /// Minimum available mapping quality.
    pub min_mapq: u8,
    /// Minimum available base quality.
    pub min_base_quality: u8,
    /// Exclude secondary alignments.
    pub exclude_secondary: bool,
    /// Exclude supplementary alignments.
    pub exclude_supplementary: bool,
    /// Exclude QC-failed reads.
    pub exclude_qc_fail: bool,
    /// Exclude duplicate-flagged reads.
    pub exclude_duplicates: bool,
}
impl DescriptorProfile {
    /// Restore the exact scientific filter profile.
    pub fn to_profile(&self) -> EvidenceProfile {
        EvidenceProfile {
            min_mapq: self.min_mapq,
            min_base_quality: self.min_base_quality,
            exclude_secondary: self.exclude_secondary,
            exclude_supplementary: self.exclude_supplementary,
            exclude_qc_fail: self.exclude_qc_fail,
            exclude_duplicates: self.exclude_duplicates,
        }
    }
}
impl From<&EvidenceProfile> for DescriptorProfile {
    fn from(p: &EvidenceProfile) -> Self {
        Self {
            min_mapq: p.min_mapq,
            min_base_quality: p.min_base_quality,
            exclude_secondary: p.exclude_secondary,
            exclude_supplementary: p.exclude_supplementary,
            exclude_qc_fail: p.exclude_qc_fail,
            exclude_duplicates: p.exclude_duplicates,
        }
    }
}
/// One declared alignment read group.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DescriptorReadGroup {
    /// Read-group identifier.
    pub id: String,
    /// Declared sample, if any.
    pub sample: Option<String>,
}
/// Resolved sample scope; Auto and equivalent Named requests serialize alike.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DescriptorSampleScope {
    /// Sample-scope schema version (1).
    pub version: u32,
    /// unknown, named, or pooled.
    pub mode: String,
    /// Selected named sample when applicable.
    pub selected_sample: Option<String>,
    /// Ordered distinct declared samples.
    pub declared_samples: Vec<String>,
    /// Ordered distinct read groups.
    pub read_groups: Vec<DescriptorReadGroup>,
    /// Whether records without a named assignment can contribute.
    pub allows_unassigned: bool,
}
impl DescriptorSampleScope {
    /// Restore a validated resolved scope.
    pub fn to_scope(&self) -> Result<EvidenceSampleScope, DatasetError> {
        let mode = match self.mode.as_str() {
            "unknown" => EvidenceSampleMode::Unknown,
            "named" => EvidenceSampleMode::Named,
            "pooled" => EvidenceSampleMode::Pooled,
            _ => return Err(corrupt("unknown sample scope mode")),
        };
        if self.version != 1
            || !strictly_sorted(&self.declared_samples)
            || !self.read_groups.windows(2).all(|w| w[0].id < w[1].id)
            || self
                .read_groups
                .iter()
                .any(|r| r.id.is_empty() || r.sample.as_ref().is_some_and(String::is_empty))
        {
            return Err(corrupt("invalid canonical sample scope"));
        }
        let declared: Vec<_> = self
            .read_groups
            .iter()
            .filter_map(|g| g.sample.clone())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        if declared != self.declared_samples
            || self.allows_unassigned != (mode != EvidenceSampleMode::Named)
            || (mode == EvidenceSampleMode::Unknown && !declared.is_empty())
            || (mode == EvidenceSampleMode::Named
                && !self
                    .selected_sample
                    .as_ref()
                    .is_some_and(|s| declared.contains(s)))
            || (mode != EvidenceSampleMode::Named && self.selected_sample.is_some())
        {
            return Err(corrupt("inconsistent sample scope"));
        }
        Ok(EvidenceSampleScope {
            mode,
            selected_sample: self.selected_sample.clone(),
            declared_samples: self.declared_samples.clone(),
            read_groups: self
                .read_groups
                .iter()
                .map(|g| EvidenceReadGroup {
                    id: g.id.clone(),
                    sample: g.sample.clone(),
                })
                .collect(),
            allows_unassigned: self.allows_unassigned,
        })
    }
    /// The same stable scope encoding used by extraction and legacy cache keys.
    pub fn canonical_json(&self) -> Result<String, DatasetError> {
        Ok(self.to_scope()?.canonical_json())
    }
}
impl From<&EvidenceSampleScope> for DescriptorSampleScope {
    fn from(s: &EvidenceSampleScope) -> Self {
        Self {
            version: 1,
            mode: match s.mode {
                EvidenceSampleMode::Unknown => "unknown",
                EvidenceSampleMode::Named => "named",
                EvidenceSampleMode::Pooled => "pooled",
            }
            .into(),
            selected_sample: s.selected_sample.clone(),
            declared_samples: s.declared_samples.clone(),
            read_groups: s
                .read_groups
                .iter()
                .map(|g| DescriptorReadGroup {
                    id: g.id.clone(),
                    sample: g.sample.clone(),
                })
                .collect(),
            allows_unassigned: s.allows_unassigned,
        }
    }
}
/// A content-bound artifact under the portable dataset root.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DescriptorArtifact {
    /// Canonical relative path; absolute paths and traversal are rejected.
    pub path: String,
    /// BLAKE3 of the complete artifact bytes.
    pub blake3: String,
}
/// Complete ownership and integrity descriptor for one canonical partition.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DescriptorPartition {
    /// Dictionary contig identifier.
    pub contig: u32,
    /// Canonical tile start, divisible by CANONICAL_TILE_BASES.
    pub start: u32,
    /// Ordered selected intervals within the tile.
    pub intervals: Vec<DescriptorInterval>,
    /// Exact selected SNV REF/ALT annotations, when applicable.
    pub sites: Option<Vec<DescriptorSite>>,
    /// Full expected row denominator, including zero-depth loci.
    pub expected_rows: u64,
    /// Canonical Arrow artifact.
    pub arrow: DescriptorArtifact,
    /// Self-verifying partition receipt.
    pub receipt: DescriptorArtifact,
}
/// Versioned portable dataset description; deserialize with from_bytes to enforce
/// the envelope, ownership and scientific-identity invariants.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DatasetDescriptor {
    /// Descriptor schema version (1).
    pub version: u32,
    /// Source-bound execution namespace used by all partition receipts.
    pub dataset_namespace: String,
    /// Scientific reuse compatibility, excluding fields and selection.
    pub compatibility_blake3: String,
    /// Existing request digest, preserved without changing legacy cache identity.
    pub request_blake3: String,
    /// Extractor package version that generated the legacy request digest.
    pub extractor_version: String,
    /// Exact extractor scientific semantics identifier.
    pub semantics: String,
    /// Read, rather than fragment, counting.
    pub counting_unit: String,
    /// Exact enumeration, always none.
    pub sampling: String,
    /// Whether actual reference sequence was available during extraction.
    pub has_reference: bool,
    /// Verified original sources; source paths remain provenance only.
    pub sources: Vec<SourceIdentity>,
    /// Exact ordered reference dictionary.
    pub contigs: Vec<DescriptorContig>,
    /// Scientific filter profile.
    pub profile: DescriptorProfile,
    /// Resolved sample selection.
    pub sample_scope: DescriptorSampleScope,
    /// Physical field mask.
    pub fields: u32,
    /// Arrow evidence schema version.
    pub schema_version: u32,
    /// Field-mask definition version.
    pub fields_version: u32,
    /// Canonical scientific selection.
    pub selection: DescriptorSelection,
    /// Complete ordered canonical partition inventory.
    pub partitions: Vec<DescriptorPartition>,
}

impl DatasetDescriptor {
    /// Decode bounded JSON, rejecting unknown keys, malformed ownership and
    /// inconsistent scientific identities before exposing a usable descriptor.
    pub fn from_bytes(bytes: &[u8], limits: DescriptorLimits) -> Result<Self, DatasetError> {
        check_size(bytes.len(), limits)?;
        let descriptor: Self = serde_json::from_slice(bytes)
            .map_err(|e| corrupt(&format!("invalid dataset descriptor: {e}")))?;
        descriptor.validate()?;
        Ok(descriptor)
    }
    /// Encode deterministic bounded JSON after checking all structural invariants.
    pub fn to_bytes(&self, limits: DescriptorLimits) -> Result<Vec<u8>, DatasetError> {
        let mut out = LimitedBytes {
            bytes: Vec::new(),
            limit: limits.max_bytes,
        };
        serde_json::to_writer(&mut out, self).map_err(|e| {
            incompatible(&format!(
                "descriptor encoding exceeds envelope or is invalid: {e}"
            ))
        })?;
        self.validate()?;
        Ok(out.bytes)
    }
    /// Reconstruct the ordered dictionary without consulting original inputs.
    pub fn contig_set(&self) -> ContigSet {
        let mut set = ContigSet::new();
        for contig in &self.contigs {
            set.push(contig.name.clone(), contig.length);
        }
        set
    }
    /// Compute the scientific compatibility key independently of stored claims.
    pub fn computed_compatibility_key(&self) -> Result<String, DatasetError> {
        compatibility(
            &self.sources,
            &self.contigs,
            &self.profile,
            &self.sample_scope,
            self.has_reference,
            &self.semantics,
        )
    }
    /// Check canonical structure and all redundant identities. Hashes establish
    /// internal consistency, not third-party authenticity of source claims.
    pub fn validate(&self) -> Result<(), DatasetError> {
        if self.version != 1
            || self.semantics != evidence::EVIDENCE_SEMANTICS_VERSION
            || self.counting_unit != "read"
            || self.sampling != "none"
            || self.extractor_version.is_empty()
            || self.extractor_version.len() > 128
        {
            return Err(corrupt("unsupported dataset descriptor semantics"));
        }
        let fields = EvidenceFields::from_bits(self.fields)?;
        if self.profile.min_base_quality > 93 || self.profile.min_mapq == 255 {
            return Err(corrupt("unsupported evidence quality thresholds"));
        }
        if self.schema_version != fields.schema_version()
            || self.fields_version != fields.mask_version()
        {
            return Err(corrupt("inconsistent descriptor field/schema versions"));
        }
        if self.sources.len() > 256
            || !self.sources.windows(2).all(|w| w[0].role < w[1].role)
            || self.sources.iter().any(|s| {
                !valid_role(&s.role)
                    || !valid_hash(&s.blake3)
                    || s.path.len() > 32768
                    || s.path.contains('\0')
                    || !Path::new(&s.path).is_absolute()
            })
            || !self.sources.iter().any(|s| s.role == "alignments")
        {
            return Err(corrupt("invalid canonical descriptor sources"));
        }
        if self.has_reference
            != self
                .sources
                .iter()
                .any(|s| matches!(s.role.as_str(), "reference" | "cram-reference"))
        {
            return Err(corrupt(
                "descriptor reference capability differs from source roles",
            ));
        }
        let mut names = BTreeSet::new();
        for (index, contig) in self.contigs.iter().enumerate() {
            if contig.id as usize != index
                || contig.name.is_empty()
                || contig.name.len() > 4096
                || contig.name.contains('\0')
                || !names.insert(&contig.name)
            {
                return Err(corrupt("invalid ordered descriptor dictionary"));
            }
        }
        self.sample_scope.to_scope()?;
        validate_selection(&self.selection, &self.contigs)?;
        if self
            .selection
            .sites
            .as_ref()
            .is_some_and(|sites| !sites.is_empty())
            && !self.has_reference
        {
            return Err(corrupt("SNV descriptor has no reference"));
        }
        let expected_count = partition_count(&self.selection.intervals)?;
        if expected_count != self.partitions.len() {
            return Err(corrupt("descriptor partition count differs from selection"));
        }
        let expected = expected_partitions(&self.selection);
        for (part, expected) in self.partitions.iter().zip(expected) {
            if part.contig != expected.contig
                || part.start != expected.start
                || part.intervals != expected.intervals
                || part.sites != expected.sites
                || part.expected_rows != expected.expected_rows
                || !valid_hash(&part.arrow.blake3)
                || !valid_hash(&part.receipt.blake3)
                || part.arrow.path != expected.arrow.path
                || part.receipt.path != expected.receipt.path
            {
                return Err(corrupt("descriptor partition ownership, annotation, or path differs from canonical selection"));
            }
        }
        if self.compatibility_blake3 != self.computed_compatibility_key()?
            || self.request_blake3 != self.computed_request_digest()?
            || self.dataset_namespace != namespace(&self.compatibility_blake3, &self.request_blake3)
        {
            return Err(corrupt("descriptor scientific identity mismatch"));
        }
        Ok(())
    }
    fn computed_request_digest(&self) -> Result<String, DatasetError> {
        let mut hash = blake3::Hasher::new();
        let mut field = |value: &[u8]| {
            hash.update(&(value.len() as u64).to_le_bytes());
            hash.update(value);
        };
        field(b"rosalind-dataset-request-v1");
        field(self.semantics.as_bytes());
        field(self.extractor_version.as_bytes());
        field(&self.schema_version.to_le_bytes());
        field(&self.fields_version.to_le_bytes());
        field(&self.fields.to_le_bytes());
        field(EvidenceProfile::ID.as_bytes());
        field(self.sample_scope.canonical_json()?.as_bytes());
        let p = &self.profile;
        field(&[
            p.min_mapq,
            p.min_base_quality,
            u8::from(p.exclude_secondary),
            u8::from(p.exclude_supplementary),
            u8::from(p.exclude_qc_fail),
            u8::from(p.exclude_duplicates),
        ]);
        field(selection_digest(&self.selection).as_bytes());
        field(&[u8::from(self.has_reference)]);
        for c in &self.contigs {
            field(&c.id.to_le_bytes());
            field(c.name.as_bytes());
            field(&c.length.to_le_bytes());
        }
        Ok(hash.finalize().to_hex().to_string())
    }
}

fn descriptor_base(engine: &EvidenceEngine, session: &VerifiedInputSession) -> DatasetDescriptor {
    let fields = engine.request().fields;
    DatasetDescriptor {
        version: 1,
        dataset_namespace: String::new(),
        compatibility_blake3: String::new(),
        request_blake3: super::request_digest(engine),
        extractor_version: env!("CARGO_PKG_VERSION").into(),
        semantics: evidence::EVIDENCE_SEMANTICS_VERSION.into(),
        counting_unit: "read".into(),
        sampling: "none".into(),
        has_reference: engine.request().reference.is_some()
            || engine.request().cram_reference.is_some(),
        sources: session.identities.clone(),
        contigs: dictionary(engine),
        profile: (&engine.request().profile).into(),
        sample_scope: engine.sample_scope().into(),
        fields: fields.bits(),
        schema_version: fields.schema_version(),
        fields_version: fields.mask_version(),
        selection: DescriptorSelection {
            intervals: engine.intervals().iter().map(interval).collect(),
            sites: match &engine.request().selection {
                EvidenceSelection::Sites(sites) => Some(sites.iter().map(site).collect()),
                _ => None,
            },
        },
        partitions: Vec::new(),
    }
}
fn dictionary(engine: &EvidenceEngine) -> Vec<DescriptorContig> {
    engine
        .contigs()
        .iter()
        .map(|c| DescriptorContig {
            id: c.id,
            name: c.name.to_string(),
            length: c.length,
        })
        .collect()
}
fn interval(i: &GenomicInterval) -> DescriptorInterval {
    DescriptorInterval {
        contig: i.contig,
        start: i.start,
        end: i.end,
    }
}
fn site(s: &SnvSite) -> DescriptorSite {
    DescriptorSite {
        contig: s.contig,
        position: s.position,
        reference: s.reference,
        alternates: s.alternates.clone(),
    }
}
fn engine_compatibility(
    engine: &EvidenceEngine,
    sources: &[SourceIdentity],
) -> Result<String, DatasetError> {
    compatibility(
        sources,
        &dictionary(engine),
        &(&engine.request().profile).into(),
        &engine.sample_scope().into(),
        engine.request().reference.is_some() || engine.request().cram_reference.is_some(),
        evidence::EVIDENCE_SEMANTICS_VERSION,
    )
}
fn compatibility(
    sources: &[SourceIdentity],
    contigs: &[DescriptorContig],
    profile: &DescriptorProfile,
    scope: &DescriptorSampleScope,
    has_reference: bool,
    semantics: &str,
) -> Result<String, DatasetError> {
    let identities: Vec<_> = sources
        .iter()
        .filter(|s| !matches!(s.role.as_str(), "sites" | "regions"))
        .map(|s| (&s.role, &s.blake3, s.bytes))
        .collect();
    let bytes = serde_json::to_vec(&(
        "rosalind-evidence-compatibility-v1",
        identities,
        contigs,
        profile,
        scope,
        has_reference,
        semantics,
        "read",
        "none",
    ))
    .map_err(|e| corrupt(&e.to_string()))?;
    Ok(blake3::hash(&bytes).to_hex().to_string())
}
fn namespace(compatibility: &str, request: &str) -> String {
    blake3::hash(format!("rosalind-verified-dataset-v1\n{compatibility}\n{request}\n").as_bytes())
        .to_hex()
        .to_string()
}
fn selection_digest(selection: &DescriptorSelection) -> String {
    let mut hash = blake3::Hasher::new();
    hash.update(b"rosalind-evidence-selection-v1\0");
    for i in &selection.intervals {
        hash.update(&i.contig.to_le_bytes());
        hash.update(&i.start.to_le_bytes());
        hash.update(&i.end.to_le_bytes());
    }
    hash.update(b"\0sites\0");
    for s in selection.sites.iter().flatten() {
        hash.update(&s.contig.to_le_bytes());
        hash.update(&s.position.to_le_bytes());
        hash.update(&[s.reference, s.alternates.len() as u8]);
        hash.update(&s.alternates);
    }
    hash.finalize().to_hex().to_string()
}
fn validate_selection(
    selection: &DescriptorSelection,
    contigs: &[DescriptorContig],
) -> Result<(), DatasetError> {
    let mut before: Option<&DescriptorInterval> = None;
    let mut rows = 0u64;
    for i in &selection.intervals {
        if i.start >= i.end
            || contigs
                .get(i.contig as usize)
                .is_none_or(|c| i.end > c.length)
            || before
                .is_some_and(|p| p.contig > i.contig || (p.contig == i.contig && p.end >= i.start))
        {
            return Err(corrupt("selection intervals are not normalized"));
        }
        rows = rows
            .checked_add(u64::from(i.end - i.start))
            .ok_or_else(|| corrupt("selection count overflow"))?;
        before = Some(i);
    }
    if let Some(sites) = &selection.sites {
        if rows != sites.len() as u64 {
            return Err(corrupt("site selection denominator differs from intervals"));
        }
        let mut interval_index = 0;
        let mut previous = None;
        for s in sites {
            while selection
                .intervals
                .get(interval_index)
                .is_some_and(|i| (i.contig, i.end) <= (s.contig, s.position))
            {
                interval_index += 1;
            }
            if !selection.intervals.get(interval_index).is_some_and(|i| {
                i.contig == s.contig && i.start <= s.position && s.position < i.end
            }) || previous.is_some_and(|p| p >= (s.contig, s.position))
                || !b"ACGT".contains(&s.reference)
                || s.alternates.is_empty()
                || s.alternates.len() > 3
                || !strictly_sorted(&s.alternates)
                || s.alternates
                    .iter()
                    .any(|b| !b"ACGT".contains(b) || *b == s.reference)
            {
                return Err(corrupt("invalid canonical SNV selection"));
            }
            previous = Some((s.contig, s.position));
        }
    }
    Ok(())
}
fn partition_count(intervals: &[DescriptorInterval]) -> Result<usize, DatasetError> {
    let mut count = 0usize;
    let mut previous = None;
    for i in intervals {
        let first = i.start / CANONICAL_TILE_BASES;
        let last = (i.end - 1) / CANONICAL_TILE_BASES;
        count = count
            .checked_add(
                (last - first + 1) as usize - usize::from(previous == Some((i.contig, first))),
            )
            .ok_or_else(|| corrupt("partition count overflow"))?;
        previous = Some((i.contig, last));
    }
    Ok(count)
}
fn expected_partitions(selection: &DescriptorSelection) -> Vec<DescriptorPartition> {
    let mut parts = BTreeMap::<(u32, u32), DescriptorPartition>::new();
    for i in &selection.intervals {
        let mut start = i.start;
        while start < i.end {
            let canonical = start / CANONICAL_TILE_BASES * CANONICAL_TILE_BASES;
            let end = i.end.min(canonical.saturating_add(CANONICAL_TILE_BASES));
            let p = parts.entry((i.contig, canonical)).or_insert_with(|| {
                let name = format!("c{:08}-p{:010}", i.contig, canonical);
                DescriptorPartition {
                    contig: i.contig,
                    start: canonical,
                    intervals: Vec::new(),
                    sites: selection.sites.as_ref().map(|_| Vec::new()),
                    expected_rows: 0,
                    arrow: DescriptorArtifact {
                        path: format!("{name}/evidence.arrow"),
                        blake3: String::new(),
                    },
                    receipt: DescriptorArtifact {
                        path: format!("{name}/manifest.json"),
                        blake3: String::new(),
                    },
                }
            });
            p.intervals.push(DescriptorInterval {
                contig: i.contig,
                start,
                end,
            });
            p.expected_rows += u64::from(end - start);
            start = end;
        }
    }
    for s in selection.sites.iter().flatten() {
        parts
            .get_mut(&(
                s.contig,
                s.position / CANONICAL_TILE_BASES * CANONICAL_TILE_BASES,
            ))
            .expect("validated selection")
            .sites
            .as_mut()
            .unwrap()
            .push(s.clone());
    }
    parts.into_values().collect()
}

fn valid_role(role: &str) -> bool {
    !role.is_empty()
        && role.len() <= 128
        && role
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || b"-_".contains(&b))
}
fn valid_hash(hash: &str) -> bool {
    hash.len() == 64
        && hash
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn strictly_sorted<T: Ord>(values: &[T]) -> bool {
    values.windows(2).all(|w| w[0] < w[1])
}
fn corrupt(message: &str) -> DatasetError {
    DatasetError::Corrupt(message.into())
}
fn incompatible(message: &str) -> DatasetError {
    DatasetError::Incompatible(message.into())
}
fn check_size(size: usize, limits: DescriptorLimits) -> Result<(), DatasetError> {
    if limits.max_bytes == 0 || size > limits.max_bytes {
        Err(incompatible(
            "dataset metadata exceeds declared byte envelope",
        ))
    } else {
        Ok(())
    }
}
fn hash_file(path: &Path) -> Result<String, DatasetError> {
    let mut file = File::open(path)?;
    let mut hash = blake3::Hasher::new();
    let mut buffer = [0u8; 65536];
    loop {
        crate::core::governor::checkpoint().map_err(EvidenceError::from)?;
        let n = file.read(&mut buffer)?;
        if n == 0 {
            break;
        }
        hash.update(&buffer[..n]);
    }
    Ok(hash.finalize().to_hex().to_string())
}
struct LimitedBytes {
    bytes: Vec<u8>,
    limit: usize,
}
impl Write for LimitedBytes {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        if self
            .bytes
            .len()
            .checked_add(bytes.len())
            .is_none_or(|n| n > self.limit)
        {
            return Err(std::io::Error::other("metadata byte limit exceeded"));
        }
        self.bytes.extend_from_slice(bytes);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

/// Conservative incremental publication reservation derived from actual selected
/// loci, intervals, dictionary and canonical partition counts, not the 32 MiB cap.
pub fn publication_memory_bytes(
    engine: &EvidenceEngine,
    limits: DescriptorLimits,
) -> Result<u64, DatasetError> {
    if limits.max_bytes == 0 {
        return Err(incompatible("descriptor byte limit must be positive"));
    }
    publication_metadata_bytes(engine)?
        .checked_mul(6)
        .and_then(|n| {
            n.checked_add(evidence::evidence_reader_memory_bytes(
                engine.request().fields,
            ))
        })
        .ok_or_else(|| incompatible("publication memory envelope overflow"))
}

fn publication_metadata_bytes(engine: &EvidenceEngine) -> Result<u64, DatasetError> {
    let mut count = 0u64;
    let mut previous = None;
    for i in engine.intervals() {
        let first = i.start / CANONICAL_TILE_BASES;
        let last = (i.end - 1) / CANONICAL_TILE_BASES;
        count = count
            .checked_add(
                u64::from(last - first + 1) - u64::from(previous == Some((i.contig, first))),
            )
            .ok_or_else(|| incompatible("publication metadata count overflow"))?;
        previous = Some((i.contig, last));
    }
    let sites = match &engine.request().selection {
        EvidenceSelection::Sites(sites) => sites.len() as u64,
        _ => 0,
    };
    let dictionary = engine.contigs().iter().try_fold(0u64, |sum, c| {
        sum.checked_add(c.name.len() as u64 + 128)
            .ok_or_else(|| incompatible("publication dictionary overflow"))
    })?;
    count
        .checked_mul(16 << 10)
        .and_then(|n| n.checked_add(sites.checked_mul(256)?))
        .and_then(|n| n.checked_add((engine.intervals().len() as u64).checked_mul(128)?))
        .and_then(|n| n.checked_add(dictionary))
        .and_then(|n| n.checked_add(engine.plan().sample_scope_bytes.checked_mul(4)?))
        .and_then(|n| n.checked_add((1 << 20) + 4096))
        .ok_or_else(|| incompatible("publication metadata envelope overflow"))
}

/// Publish a portable descriptor and new receipt alongside an untouched legacy
/// cache manifest. Only a complete source-bound session outcome is accepted.
/// Every Arrow partition, ownership annotation, receipt and source snapshot is
/// checked before the two new files are committed as one transactional group.
pub fn publish_evidence_dataset(
    engine: &EvidenceEngine,
    outcome: &DatasetOutcome,
    session: &VerifiedInputSession,
    limits: DescriptorLimits,
) -> Result<PathBuf, DatasetError> {
    session.assert_engine(engine)?;
    let expected_namespace = session.dataset_namespace(engine)?;
    if outcome.science_digest != expected_namespace {
        return Err(incompatible(
            "dataset outcome was not produced under this verified input session namespace",
        ));
    }
    let memory = publication_memory_bytes(engine, limits)?;
    if let Some(budget) = engine.request().execution.memory_budget_bytes {
        let needed = crate::util::rss::peak_rss_bytes().saturating_add(memory);
        if needed > budget {
            return Err(EvidenceError::Refused { needed, budget }.into());
        }
    }
    if outcome
        .dataset_manifest
        .file_name()
        .and_then(|n| n.to_str())
        != Some(super::MANIFEST_NAME)
    {
        return Err(incompatible("outcome is not a legacy dataset manifest"));
    }
    let root = fs::canonicalize(
        outcome
            .dataset_manifest
            .parent()
            .ok_or_else(|| incompatible("dataset manifest has no directory"))?,
    )?;
    let metadata_limit = DescriptorLimits {
        max_bytes: limits
            .max_bytes
            .min(usize::try_from(publication_metadata_bytes(engine)?).unwrap_or(usize::MAX)),
    };
    let legacy = bounded_receipt(&root.join(super::MANIFEST_NAME), metadata_limit)?;
    let request = super::request_digest(engine);
    let fields = engine.request().fields;
    if legacy.subcommand != "dataset evidence"
        || legacy.params.get("run_status").map(String::as_str) != Some("completed")
        || legacy.params.get("evidence.science_blake3") != Some(&expected_namespace)
        || legacy.params.get("dataset.request_blake3") != Some(&request)
        || super::receipt_fields(&legacy)? != fields
    {
        return Err(incompatible(
            "complete legacy dataset claim differs from verified session",
        ));
    }
    let mut descriptor = descriptor_base(engine, session);
    descriptor.compatibility_blake3 = descriptor.computed_compatibility_key()?;
    descriptor.dataset_namespace = expected_namespace;
    descriptor.partitions = expected_partitions(&descriptor.selection);
    if legacy.outputs.len() != descriptor.partitions.len()
        || legacy.inputs.len() != descriptor.partitions.len()
        || legacy.params.get("dataset.partition_count")
            != Some(&descriptor.partitions.len().to_string())
    {
        return Err(corrupt("legacy dataset inventory differs from selection"));
    }
    let legacy_outputs = resolved_inventory(&legacy.outputs)?;
    let legacy_inputs = resolved_inventory(&legacy.inputs)?;
    let parts = super::partitions(engine);
    let mut artifact_guards = Vec::with_capacity(parts.len() * 2 + 1);
    artifact_guards.push(root.join(super::MANIFEST_NAME));
    for part in &descriptor.partitions {
        artifact_guards.push(root.join(&part.arrow.path));
        artifact_guards.push(root.join(&part.receipt.path));
    }
    let partition_snapshot = InputSnapshot::capture(artifact_guards)?;
    let mut total_rows = 0u64;
    for (part, published) in parts.iter().zip(&mut descriptor.partitions) {
        crate::core::governor::checkpoint().map_err(EvidenceError::from)?;
        let part_root = root.join(part.name());
        let receipt = bounded_receipt(
            &part_root.join("manifest.json"),
            DescriptorLimits {
                max_bytes: limits.max_bytes.min(64 << 10),
            },
        )?;
        let expected = [
            ("dataset.request_blake3", request.clone()),
            (
                "evidence.science_blake3",
                descriptor.dataset_namespace.clone(),
            ),
            ("evidence.schema", fields.schema_version().to_string()),
            ("pileup.semantics", "exact-or-fail-v1".into()),
            ("dataset.contig", part.contig.to_string()),
            ("dataset.partition_start", part.start.to_string()),
            ("dataset.partition_bases", CANONICAL_TILE_BASES.to_string()),
            ("dataset.selection_blake3", part.selection_hash()),
            ("outcome.rows", part.row_count().to_string()),
            ("run_status", "completed".into()),
        ];
        if receipt.outputs.len() != 1
            || super::receipt_fields(&receipt)? != fields
            || expected
                .iter()
                .any(|(key, value)| receipt.params.get(*key) != Some(value))
        {
            return Err(corrupt("partition claim differs from verified session"));
        }
        let original_arrow = root.join(&published.arrow.path);
        let original_receipt = root.join(&published.receipt.path);
        published.arrow.blake3 = hash_file(&original_arrow)?;
        published.receipt.blake3 = hash_file(&original_receipt)?;
        if receipt.outputs[0].blake3 != published.arrow.blake3
            || legacy_outputs.get(&original_arrow) != Some(&published.arrow.blake3)
            || legacy_inputs.get(&original_receipt) != Some(&published.receipt.blake3)
        {
            return Err(corrupt(
                "legacy dataset inventory does not bind a partition",
            ));
        }
        let mut interval_index = 0usize;
        let mut next = published.intervals.first().map_or(0, |i| i.start);
        let mut rows = 0u64;
        evidence::read_evidence_batches_expected_fields(
            BufReader::new(File::open(&original_arrow)?),
            engine.contigs(),
            fields,
            |batch| {
                if batch.contig_id != published.contig
                    || batch.canonical_tile_start != published.start
                {
                    return Err(EvidenceError::InvalidInput(
                        "partition batch ownership differs".into(),
                    ));
                }
                for row in batch.rows() {
                    let interval = published.intervals.get(interval_index).ok_or_else(|| {
                        EvidenceError::InvalidInput("partition contains extra rows".into())
                    })?;
                    if row.position != next {
                        return Err(EvidenceError::InvalidInput(
                            "partition omits, repeats or reorders selected loci".into(),
                        ));
                    }
                    if let Some(sites) = &published.sites {
                        let site = sites.get(rows as usize).ok_or_else(|| {
                            EvidenceError::InvalidInput("partition has extra SNV rows".into())
                        })?;
                        if row.reference != site.reference || row.requested_alts != site.alternates
                        {
                            return Err(EvidenceError::InvalidInput(
                                "partition SNV annotations differ".into(),
                            ));
                        }
                    } else if !row.requested_alts.is_empty() {
                        return Err(EvidenceError::InvalidInput(
                            "interval partition has unexpected SNV annotations".into(),
                        ));
                    }
                    rows += 1;
                    next += 1;
                    if next == interval.end {
                        interval_index += 1;
                        if let Some(i) = published.intervals.get(interval_index) {
                            next = i.start;
                        }
                    }
                }
                Ok(())
            },
        )?;
        if rows != published.expected_rows || interval_index != published.intervals.len() {
            return Err(corrupt("partition row denominator differs"));
        }
        total_rows = total_rows
            .checked_add(rows)
            .ok_or_else(|| corrupt("dataset row count overflow"))?;
    }
    if legacy.params.get("outcome.rows") != Some(&total_rows.to_string()) {
        return Err(corrupt("legacy dataset row denominator differs"));
    }
    let bytes = descriptor.to_bytes(limits)?;
    let descriptor_hash = blake3::hash(&bytes).to_hex().to_string();
    let mut manifest = RunManifest::new("dataset evidence portable");
    manifest.tool_version = env!("CARGO_PKG_VERSION").into();
    manifest.params.extend([
        ("dataset.descriptor_version".into(), "1".into()),
        (
            "dataset.namespace".into(),
            descriptor.dataset_namespace.clone(),
        ),
        (
            "dataset.compatibility_blake3".into(),
            descriptor.compatibility_blake3.clone(),
        ),
        (
            "dataset.request_blake3".into(),
            descriptor.request_blake3.clone(),
        ),
        (
            "dataset.partition_count".into(),
            descriptor.partitions.len().to_string(),
        ),
        (
            "evidence.schema".into(),
            fields.schema_version().to_string(),
        ),
        ("evidence.fields".into(), fields.bits().to_string()),
        (
            "evidence.fields_version".into(),
            fields.mask_version().to_string(),
        ),
        ("evidence.semantics".into(), descriptor.semantics.clone()),
        ("run_status".into(), "completed".into()),
        ("outcome.rows".into(), total_rows.to_string()),
    ]);
    manifest.inputs = session
        .identities
        .iter()
        .map(|s| FileHash {
            path: s.path.clone(),
            blake3: s.blake3.clone(),
        })
        .collect();
    manifest.outputs.push(FileHash {
        path: DATASET_DESCRIPTOR_NAME.into(),
        blake3: descriptor_hash,
    });
    for part in &descriptor.partitions {
        manifest.inputs.push(FileHash {
            path: part.receipt.path.clone(),
            blake3: part.receipt.blake3.clone(),
        });
        manifest.outputs.push(FileHash {
            path: part.arrow.path.clone(),
            blake3: part.arrow.blake3.clone(),
        });
    }
    manifest.finalize();
    let receipt_bytes = manifest.to_canonical_json();
    check_size(receipt_bytes.len(), limits)?;
    let descriptor_path = root.join(DATASET_DESCRIPTOR_NAME);
    let receipt_path = root.join(EVIDENCE_DATASET_MANIFEST_NAME);
    let mut descriptor_file = AtomicFile::create(&descriptor_path)?;
    descriptor_file.file_mut().write_all(&bytes)?;
    let mut receipt_file = AtomicFile::create(&receipt_path)?;
    receipt_file
        .file_mut()
        .write_all(receipt_bytes.as_bytes())?;
    session.verify()?;
    partition_snapshot.verify()?;
    crate::core::governor::checkpoint().map_err(EvidenceError::from)?;
    commit_group(
        vec![
            (descriptor_file, descriptor_path),
            (receipt_file, receipt_path.clone()),
        ],
        true,
    )?;
    Ok(receipt_path)
}

fn bounded_receipt(path: &Path, limits: DescriptorLimits) -> Result<RunManifest, DatasetError> {
    let mut bytes = Vec::new();
    File::open(path)?
        .take((limits.max_bytes as u64).saturating_add(1))
        .read_to_end(&mut bytes)?;
    check_size(bytes.len(), limits)?;
    let text = std::str::from_utf8(&bytes).map_err(|_| corrupt("receipt is not UTF-8"))?;
    let receipt = RunManifest::from_canonical_json(text).map_err(|e| corrupt(&e.to_string()))?;
    if receipt.self_hash_ok() != Some(true) || receipt.measurement_hash_ok() == Some(false) {
        return Err(corrupt("dataset receipt self-hash failed"));
    }
    Ok(receipt)
}
fn resolved_inventory(entries: &[FileHash]) -> Result<BTreeMap<PathBuf, String>, DatasetError> {
    let mut inventory = BTreeMap::new();
    for entry in entries {
        if entry.path.len() > 4096 {
            return Err(corrupt(
                "legacy artifact path exceeds 4096-byte publication envelope",
            ));
        }
        if inventory
            .insert(fs::canonicalize(&entry.path)?, entry.blake3.clone())
            .is_some()
        {
            return Err(corrupt("legacy dataset inventory repeats an artifact"));
        }
    }
    Ok(inventory)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dataset::{run_dataset_with_snapshot, DatasetOptions};
    use crate::evidence::{EvidenceCallback, EvidenceRequest};
    use rust_htslib::bam::{
        self,
        record::{Cigar, CigarString},
    };
    use std::sync::atomic::{AtomicU64, Ordering};
    static SEQUENCE: AtomicU64 = AtomicU64::new(0);
    struct Fixture {
        root: PathBuf,
        bam: PathBuf,
        reference: PathBuf,
    }
    impl Fixture {
        fn new() -> Self {
            let root = std::env::temp_dir().join(format!(
                "rosalind-descriptor-{}-{}",
                std::process::id(),
                SEQUENCE.fetch_add(1, Ordering::Relaxed)
            ));
            fs::create_dir_all(&root).unwrap();
            let reference = root.join("reference.fa");
            fs::write(&reference, format!(">chr1\n{}\n", "A".repeat(20000))).unwrap();
            fs::write(
                root.join("reference.fa.fai"),
                "chr1\t20000\t6\t20000\t20001\n",
            )
            .unwrap();
            let bam = root.join("reads.bam");
            let fixture = Self {
                root,
                bam,
                reference,
            };
            fixture.write_bam(b"ACAA");
            fixture
        }
        fn write_bam(&self, seq: &[u8]) {
            let mut header = bam::Header::new();
            header.push_record(
                bam::header::HeaderRecord::new(b"SQ")
                    .push_tag(b"SN", "chr1")
                    .push_tag(b"LN", 20000),
            );
            let mut writer = bam::Writer::from_path(&self.bam, &header, bam::Format::Bam).unwrap();
            let mut record = bam::Record::new();
            record.set(
                b"read",
                Some(&CigarString(vec![Cigar::Match(4)])),
                seq,
                &[30; 4],
            );
            record.set_tid(0);
            record.set_pos(0);
            record.set_mapq(60);
            record.set_flags(0);
            writer.write(&record).unwrap();
            drop(writer);
            bam::index::build(&self.bam, None::<&PathBuf>, bam::index::Type::Bai, 1).unwrap();
        }
        fn paths(&self) -> Vec<(String, PathBuf)> {
            vec![
                ("alignments".into(), self.bam.clone()),
                ("alignment-index".into(), self.root.join("reads.bam.bai")),
                ("reference".into(), self.reference.clone()),
                ("reference-fai".into(), self.root.join("reference.fa.fai")),
            ]
        }
        fn session(&self) -> VerifiedInputSession {
            VerifiedInputSession::open(self.paths()).unwrap()
        }
        fn engine(&self) -> EvidenceEngine {
            let mut request = EvidenceRequest::new(&self.bam, &self.reference);
            request.alignment_index = Some(self.root.join("reads.bam.bai"));
            request.reference_fai = Some(self.root.join("reference.fa.fai"));
            request.fields = EvidenceFields::DEPTHS.union(EvidenceFields::ALLELES);
            request.selection = EvidenceSelection::Sites(vec![
                SnvSite {
                    contig: 0,
                    position: 1,
                    reference: b'A',
                    alternates: vec![b'C'],
                },
                SnvSite {
                    contig: 0,
                    position: 18000,
                    reference: b'A',
                    alternates: vec![b'T'],
                },
            ]);
            EvidenceEngine::open(request).unwrap()
        }
        fn run(
            &self,
            engine: &mut EvidenceEngine,
            session: &VerifiedInputSession,
            namespace: &str,
        ) -> DatasetOutcome {
            let fields = engine.request().fields;
            run_dataset_with_snapshot(
                engine,
                namespace,
                &DatasetOptions {
                    cache_dir: self.root.join("cache"),
                    resume: false,
                    workers: 1,
                },
                &mut EvidenceCallback::with_fields(|_: &evidence::EvidenceBatch| Ok(()), 0, fields),
                session.snapshot(),
            )
            .unwrap()
        }
    }
    impl Drop for Fixture {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.root);
        }
    }

    #[test]
    fn verified_publication_preserves_legacy_and_binds_relative_partition_inventory() {
        let f = Fixture::new();
        let session = f.session();
        let mut engine = f.engine();
        let namespace = session.dataset_namespace(&engine).unwrap();
        let outcome = f.run(&mut engine, &session, &namespace);
        let before = fs::read(&outcome.dataset_manifest).unwrap();
        let published =
            publish_evidence_dataset(&engine, &outcome, &session, DescriptorLimits::default())
                .unwrap();
        assert_eq!(fs::read(&outcome.dataset_manifest).unwrap(), before);
        let root = published.parent().unwrap();
        let descriptor = DatasetDescriptor::from_bytes(
            &fs::read(root.join(DATASET_DESCRIPTOR_NAME)).unwrap(),
            DescriptorLimits::default(),
        )
        .unwrap();
        assert_eq!(descriptor.dataset_namespace, namespace);
        assert_eq!(
            descriptor.request_blake3,
            super::super::request_digest(&engine)
        );
        assert_eq!(descriptor.partitions.len(), 2);
        assert_eq!(
            descriptor
                .partitions
                .iter()
                .map(|p| p.expected_rows)
                .sum::<u64>(),
            2
        );
        let receipt = bounded_receipt(&published, DescriptorLimits::default()).unwrap();
        assert_eq!(receipt.outputs.len(), 3);
        for output in &receipt.outputs {
            assert!(!Path::new(&output.path).is_absolute());
            assert_eq!(hash_file(&root.join(&output.path)).unwrap(), output.blake3);
        }
        for p in &descriptor.partitions {
            assert_eq!(
                hash_file(&root.join(&p.receipt.path)).unwrap(),
                p.receipt.blake3
            );
        }
        let mut bad = descriptor.clone();
        bad.partitions[0].arrow.path = "../escape.arrow".into();
        assert!(bad.to_bytes(DescriptorLimits::default()).is_err());
        let mut value = serde_json::to_value(&descriptor).unwrap();
        value["unexpected"] = true.into();
        assert!(DatasetDescriptor::from_bytes(
            &serde_json::to_vec(&value).unwrap(),
            DescriptorLimits::default()
        )
        .is_err());
        assert!(DatasetDescriptor::from_bytes(
            &fs::read(root.join(DATASET_DESCRIPTOR_NAME)).unwrap(),
            DescriptorLimits { max_bytes: 10 }
        )
        .is_err());
    }

    #[test]
    fn arbitrary_old_key_and_changed_input_session_cannot_promote_stale_partitions() {
        let f = Fixture::new();
        let old_session = f.session();
        let mut engine = f.engine();
        let legacy = f.run(&mut engine, &old_session, &"a".repeat(64));
        assert!(publish_evidence_dataset(
            &engine,
            &legacy,
            &old_session,
            DescriptorLimits::default()
        )
        .unwrap_err()
        .to_string()
        .contains("namespace"));
        let namespace = old_session.dataset_namespace(&engine).unwrap();
        let mut outcome = f.run(&mut engine, &old_session, &namespace);
        f.write_bam(b"AGAA");
        assert!(old_session.verify().is_err());
        let new_session = f.session();
        let new_engine = f.engine();
        let new_namespace = new_session.dataset_namespace(&new_engine).unwrap();
        assert_ne!(namespace, new_namespace);
        assert!(publish_evidence_dataset(
            &new_engine,
            &outcome,
            &new_session,
            DescriptorLimits::default()
        )
        .is_err());
        // Public outcome fields cannot substitute for the receipts' old identity.
        outcome.science_digest = new_namespace;
        assert!(publish_evidence_dataset(
            &new_engine,
            &outcome,
            &new_session,
            DescriptorLimits::default()
        )
        .is_err());
        assert!(!outcome
            .dataset_manifest
            .parent()
            .unwrap()
            .join(DATASET_DESCRIPTOR_NAME)
            .exists());
    }

    #[test]
    fn compatibility_ignores_selection_projection_and_locations_but_binds_filters() {
        let f = Fixture::new();
        let source = f.root.join("selection.vcf");
        fs::write(&source, b"first selection source").unwrap();
        let mut paths = f.paths();
        paths.push(("sites".into(), source.clone()));
        let first = VerifiedInputSession::open(paths.clone()).unwrap();
        let engine = f.engine();
        let key = first.compatibility_key(&engine).unwrap();
        let mut request = engine.request().clone();
        request.fields = EvidenceFields::ALLELES;
        request.selection = EvidenceSelection::Intervals(vec![GenomicInterval {
            contig: 0,
            start: 0,
            end: 10,
        }]);
        request.execution.max_microtile_bases = 1;
        let projected = EvidenceEngine::open(request.clone()).unwrap();
        assert_eq!(key, first.compatibility_key(&projected).unwrap());
        assert_ne!(
            first.dataset_namespace(&engine).unwrap(),
            first.dataset_namespace(&projected).unwrap()
        );
        request.profile.min_mapq = 30;
        let changed = EvidenceEngine::open(request).unwrap();
        assert_ne!(key, first.compatibility_key(&changed).unwrap());
        fs::write(source, b"reordered equivalent raw selection source").unwrap();
        let second = VerifiedInputSession::open(paths).unwrap();
        assert_eq!(key, second.compatibility_key(&engine).unwrap());
        assert_eq!(first.hashes()["alignments"], second.hashes()["alignments"]);
        let mut missing = f.paths();
        missing.retain(|(r, _)| r != "alignment-index");
        assert!(VerifiedInputSession::open(missing)
            .unwrap()
            .dataset_namespace(&engine)
            .is_err());
    }
    #[test]
    fn discovered_index_cannot_hide_under_a_selection_role() {
        let f = Fixture::new();
        let mut request = f.engine().request().clone();
        request.alignment_index = None;
        let engine = EvidenceEngine::open(request).unwrap();
        let mut paths = f.paths();
        paths
            .iter_mut()
            .find(|(role, _)| role == "alignment-index")
            .unwrap()
            .0 = "regions".into();
        let session = VerifiedInputSession::open(paths).unwrap();
        assert!(session
            .compatibility_key(&engine)
            .unwrap_err()
            .to_string()
            .contains("canonical source role"));
    }

    #[test]
    fn publication_memory_reserves_read_group_and_sample_metadata() {
        let f = Fixture::new();
        let before = publication_memory_bytes(&f.engine(), DescriptorLimits::default()).unwrap();
        let mut header = bam::Header::new();
        header.push_record(
            bam::header::HeaderRecord::new(b"SQ")
                .push_tag(b"SN", "chr1")
                .push_tag(b"LN", 20000),
        );
        for index in 0..32 {
            header.push_record(
                bam::header::HeaderRecord::new(b"RG")
                    .push_tag(b"ID", format!("{index}-{}", "x".repeat(1024)))
                    .push_tag(b"SM", "sample-name"),
            );
        }
        drop(bam::Writer::from_path(&f.bam, &header, bam::Format::Bam).unwrap());
        bam::index::build(&f.bam, None::<&PathBuf>, bam::index::Type::Bai, 1).unwrap();
        let engine = f.engine();
        let after = publication_memory_bytes(&engine, DescriptorLimits::default()).unwrap();
        assert!(after > before + 32 * 1024 * 6);
    }
    #[cfg(unix)]
    #[test]
    fn verified_session_detects_symlink_retargeting() {
        use std::os::unix::fs::symlink;
        let f = Fixture::new();
        let first = f.root.join("target-a");
        let second = f.root.join("target-b");
        let access = f.root.join("access");
        fs::write(&first, "same length A").unwrap();
        fs::write(&second, "same length B").unwrap();
        symlink(&first, &access).unwrap();
        let session =
            VerifiedInputSession::open(vec![("alignments".into(), access.clone())]).unwrap();
        session.verify().unwrap();
        fs::remove_file(&access).unwrap();
        symlink(&second, &access).unwrap();
        assert!(session.verify().is_err());
        assert_eq!(fs::read_to_string(first).unwrap(), "same length A");
    }
}
