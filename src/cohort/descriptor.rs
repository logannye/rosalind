use super::{CohortError, Result};
use crate::dataset::DatasetReadLimits;
use serde::{Deserialize, Serialize};
use std::io::Write;

pub(crate) const SNAPSHOT_VERSION: u32 = 1;
pub(crate) const COMPARISON_VERSION: u32 = 1;
pub(crate) const SNAPSHOT_FILE: &str = "snapshot.json";
pub(crate) const SNAPSHOT_RECEIPT: &str = "manifest.json";

/// Operational envelopes, not scientific settings or universal allocation caps.
#[derive(Debug, Clone, Copy)]
pub(crate) struct CohortLimits {
    pub max_snapshot_bytes: usize,
    pub max_members: usize,
    pub max_leaf_references: usize,
    pub max_ownership_intervals: usize,
    pub max_parent_depth: usize,
    pub memory_budget_bytes: Option<u64>,
    pub dataset: DatasetReadLimits,
}

impl Default for CohortLimits {
    fn default() -> Self {
        Self {
            max_snapshot_bytes: 8 << 20,
            max_members: 65_536,
            max_leaf_references: 262_144,
            max_ownership_intervals: 262_144,
            max_parent_depth: 1_024,
            memory_budget_bytes: None,
            dataset: DatasetReadLimits::default(),
        }
    }
}

impl CohortLimits {
    pub fn validate(&self) -> Result<()> {
        if self.max_snapshot_bytes == 0
            || self.max_members == 0
            || self.max_leaf_references == 0
            || self.max_ownership_intervals == 0
            || self.max_parent_depth == 0
        {
            return Err(CohortError::Limit("envelopes must be positive".into()));
        }
        Ok(())
    }

    pub fn admit(&self, additional: u64) -> Result<()> {
        crate::core::governor::checkpoint().map_err(crate::evidence::EvidenceError::from)?;
        let needed = crate::util::rss::peak_rss_bytes().saturating_add(additional);
        if self
            .memory_budget_bytes
            .is_some_and(|budget| needed > budget)
        {
            return Err(CohortError::Limit(format!(
                "needs at least {needed} process bytes; budget is {}",
                self.memory_budget_bytes.unwrap()
            )));
        }
        Ok(())
    }

    pub fn dataset_limits(&self) -> DatasetReadLimits {
        DatasetReadLimits {
            memory_budget_bytes: self
                .memory_budget_bytes
                .or(self.dataset.memory_budget_bytes),
            ..self.dataset
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct MemberMetadata {
    pub id: String,
    pub group: Option<String>,
    pub subject: Option<String>,
    pub timepoint: Option<String>,
}

impl MemberMetadata {
    pub fn validate(&self) -> Result<()> {
        if !valid_text(&self.id, 256)
            || [&self.group, &self.subject, &self.timepoint]
                .into_iter()
                .flatten()
                .any(|value| !valid_text(value, 1_024))
        {
            return Err(CohortError::Corrupt(
                "member IDs must be nonempty, at most 256 UTF-8 bytes; optional metadata at most 1024; control characters are forbidden".into(),
            ));
        }
        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct LeafRef {
    pub object_id: String,
    pub descriptor_blake3: String,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct SnapshotMember {
    pub metadata: MemberMetadata,
    pub leaves: Vec<LeafRef>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct SnapshotDescriptor {
    pub version: u32,
    pub comparison_version: u32,
    pub parent: Option<String>,
    pub members: Vec<SnapshotMember>,
}

impl SnapshotDescriptor {
    pub fn validate(&self, limits: &CohortLimits) -> Result<()> {
        limits.validate()?;
        if self.version != SNAPSHOT_VERSION || self.comparison_version != COMPARISON_VERSION {
            return Err(CohortError::Corrupt(
                "unsupported snapshot/comparison version".into(),
            ));
        }
        if self.parent.as_ref().is_some_and(|id| !valid_hash(id)) {
            return Err(CohortError::Corrupt("invalid parent snapshot ID".into()));
        }
        if self.members.len() > limits.max_members {
            return Err(CohortError::Limit("too many members".into()));
        }
        let mut count = 0usize;
        let mut previous: Option<&str> = None;
        for member in &self.members {
            member.metadata.validate()?;
            if previous.is_some_and(|id| id >= member.metadata.id.as_str()) {
                return Err(CohortError::Corrupt(
                    "members must have unique sorted IDs".into(),
                ));
            }
            previous = Some(&member.metadata.id);
            if member.leaves.is_empty() {
                return Err(CohortError::Corrupt(
                    "a member requires at least one leaf".into(),
                ));
            }
            count = count
                .checked_add(member.leaves.len())
                .ok_or_else(|| CohortError::Limit("leaf count overflow".into()))?;
            if count > limits.max_leaf_references {
                return Err(CohortError::Limit("too many leaf references".into()));
            }
            if member
                .leaves
                .iter()
                .any(|leaf| !valid_hash(&leaf.object_id) || !valid_hash(&leaf.descriptor_blake3))
                || !member
                    .leaves
                    .windows(2)
                    .all(|pair| pair[0].object_id < pair[1].object_id)
            {
                return Err(CohortError::Corrupt(
                    "leaf identities must be valid, unique and sorted".into(),
                ));
            }
        }
        Ok(())
    }

    pub fn to_bytes(&self, limits: &CohortLimits) -> Result<Vec<u8>> {
        self.validate(limits)?;
        // Counting first bounds both the serialized allocation and parser copies.
        let mut count = LimitedWriter {
            count: 0,
            max: limits.max_snapshot_bytes,
        };
        serde_json::to_writer(&mut count, self).map_err(|e| CohortError::Limit(e.to_string()))?;
        limits.admit((count.count as u64).saturating_mul(16))?;
        serde_json::to_vec(self).map_err(|e| CohortError::Corrupt(e.to_string()))
    }

    pub fn from_bytes(bytes: &[u8], limits: &CohortLimits) -> Result<Self> {
        if bytes.len() > limits.max_snapshot_bytes {
            return Err(CohortError::Limit("snapshot exceeds byte envelope".into()));
        }
        limits.admit((bytes.len() as u64).saturating_mul(16))?;
        let descriptor: Self = serde_json::from_slice(bytes)
            .map_err(|e| CohortError::Corrupt(format!("invalid snapshot JSON: {e}")))?;
        descriptor.validate(limits)?;
        if descriptor.to_bytes(limits)? != bytes {
            return Err(CohortError::Corrupt(
                "snapshot is not canonical JSON".into(),
            ));
        }
        Ok(descriptor)
    }
}

struct LimitedWriter {
    count: usize,
    max: usize,
}
impl Write for LimitedWriter {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        self.count = self
            .count
            .checked_add(bytes.len())
            .filter(|n| *n <= self.max)
            .ok_or_else(|| std::io::Error::other("snapshot exceeds byte envelope"))?;
        Ok(bytes.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

pub(crate) fn valid_hash(hash: &str) -> bool {
    hash.len() == 64
        && hash
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn valid_text(value: &str, max: usize) -> bool {
    !value.is_empty() && value.len() <= max && !value.chars().any(char::is_control)
}
