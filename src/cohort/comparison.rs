//! Cross-sample comparison is deliberately distinct from same-source reuse.
//!
//! Descriptors still validate their source-bound compatibility identities. A
//! comparison contract excludes alignment bytes, sample labels, paths, physical
//! projections and execution choices. It never relaxes a saved dataset's filters.

use super::descriptor::COMPARISON_VERSION;
use super::{CohortError, Result};
use crate::dataset::{DatasetDescriptor, DescriptorContig, DescriptorProfile};
use crate::evidence::EvidenceFields;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct ComparisonReference {
    pub blake3: String,
    pub bytes: u64,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct ComparisonContract {
    pub version: u32,
    pub semantics: String,
    pub counting_unit: String,
    pub sampling: String,
    pub profile: DescriptorProfile,
    pub contigs: Vec<DescriptorContig>,
    pub reference: ComparisonReference,
    /// FASTA coordinate mapping is an input to the observed reference bases.
    /// Native packed references have no FAI; presence must also agree.
    pub reference_index: Option<ComparisonReference>,
}

/// Stable codes let a planner explain incompatibility without parsing prose.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub(crate) enum ComparisonMismatchCode {
    ContractVersion,
    EvidenceSemantics,
    CountingUnit,
    SamplingPolicy,
    FilterProfile,
    ReferenceDictionary,
    ReferenceContent,
    ReferenceIndex,
    RequiredFields,
    InvalidDescriptor,
    SampleScope,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct ComparisonMismatch {
    pub code: ComparisonMismatchCode,
    pub expected: String,
    pub observed: String,
}

impl std::fmt::Display for ComparisonMismatch {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{:?}: expected {}; observed {}",
            self.code, self.expected, self.observed
        )
    }
}

impl ComparisonContract {
    /// The caller admits the descriptor and this compact comparison copy. No
    /// original reference/alignment path is opened or trusted as an identity.
    pub fn from_descriptor(descriptor: &DatasetDescriptor) -> Result<Self> {
        descriptor.validate()?;
        if !descriptor.has_reference {
            return Err(CohortError::Incompatible(
                "candidate comparison requires stored reference evidence".into(),
            ));
        }
        // Match the engine's explicit analysis-reference precedence. A BAM and
        // CRAM can name the same FASTA under different roles; role and pathname
        // are not biological incompatibilities. Byte identity remains strict.
        let source = descriptor
            .sources
            .iter()
            .find(|source| source.role == "reference")
            .or_else(|| {
                descriptor
                    .sources
                    .iter()
                    .find(|source| source.role == "cram-reference")
            })
            .ok_or_else(|| {
                CohortError::Incompatible("missing analysis reference identity".into())
            })?;
        Ok(Self {
            version: COMPARISON_VERSION,
            semantics: descriptor.semantics.clone(),
            counting_unit: descriptor.counting_unit.clone(),
            sampling: descriptor.sampling.clone(),
            profile: descriptor.profile.clone(),
            contigs: descriptor.contigs.clone(),
            reference: ComparisonReference {
                blake3: source.blake3.clone(),
                bytes: source.bytes,
            },
            reference_index: descriptor
                .sources
                .iter()
                .find(|index| {
                    index.role
                        == if source.role == "reference" {
                            "reference-fai"
                        } else {
                            "cram-reference-fai"
                        }
                })
                .map(|index| ComparisonReference {
                    blake3: index.blake3.clone(),
                    bytes: index.bytes,
                }),
        })
    }

    /// Canonical identity of this versioned comparison contract, not an identity
    /// for reusing one sample's evidence to represent another sample.
    pub fn digest(&self) -> Result<String> {
        let bytes = serde_json::to_vec(self)
            .map_err(|error| CohortError::Corrupt(format!("comparison encoding: {error}")))?;
        let mut hash = blake3::Hasher::new();
        hash.update(b"rosalind-cohort-comparison-v1\0");
        hash.update(&bytes);
        Ok(hash.finalize().to_hex().to_string())
    }

    /// Differences are deterministic and bounded by the contract fields, not by
    /// the number of differing contigs. Report only the first dictionary change.
    pub fn compare(&self, observed: &Self) -> Vec<ComparisonMismatch> {
        let mut failures = Vec::new();
        let mut add = |code, expected: String, observed: String| {
            failures.push(ComparisonMismatch {
                code,
                expected,
                observed,
            });
        };
        if self.version != observed.version {
            add(
                ComparisonMismatchCode::ContractVersion,
                self.version.to_string(),
                observed.version.to_string(),
            );
        }
        for (code, expected, actual) in [
            (
                ComparisonMismatchCode::EvidenceSemantics,
                &self.semantics,
                &observed.semantics,
            ),
            (
                ComparisonMismatchCode::CountingUnit,
                &self.counting_unit,
                &observed.counting_unit,
            ),
            (
                ComparisonMismatchCode::SamplingPolicy,
                &self.sampling,
                &observed.sampling,
            ),
        ] {
            if expected != actual {
                add(code, expected.clone(), actual.clone());
            }
        }
        if self.profile != observed.profile {
            add(
                ComparisonMismatchCode::FilterProfile,
                format!("{:?}", self.profile),
                format!("{:?}", observed.profile),
            );
        }
        if self.contigs.len() != observed.contigs.len() {
            add(
                ComparisonMismatchCode::ReferenceDictionary,
                format!("{} contigs", self.contigs.len()),
                format!("{} contigs", observed.contigs.len()),
            );
        } else if let Some((expected, actual)) = self
            .contigs
            .iter()
            .zip(&observed.contigs)
            .find(|(a, b)| a != b)
        {
            add(
                ComparisonMismatchCode::ReferenceDictionary,
                format!("{expected:?}"),
                format!("{actual:?}"),
            );
        }
        if self.reference != observed.reference {
            add(
                ComparisonMismatchCode::ReferenceContent,
                format!("{} ({} bytes)", self.reference.blake3, self.reference.bytes),
                format!(
                    "{} ({} bytes)",
                    observed.reference.blake3, observed.reference.bytes
                ),
            );
        }
        if self.reference_index != observed.reference_index {
            add(
                ComparisonMismatchCode::ReferenceIndex,
                format!("{:?}", self.reference_index),
                format!("{:?}", observed.reference_index),
            );
        }
        failures
    }
}

/// A different physical mask/schema can be adequate. Omitted fields can never
/// become zero evidence or a partial cell; this is a planning error.
pub(crate) fn check_required_fields(
    descriptor: &DatasetDescriptor,
    requested: EvidenceFields,
) -> std::result::Result<(), ComparisonMismatch> {
    descriptor.validate().map_err(|error| ComparisonMismatch {
        code: ComparisonMismatchCode::InvalidDescriptor,
        expected: "a valid supported portable dataset descriptor".into(),
        observed: error.to_string(),
    })?;
    let available =
        EvidenceFields::from_bits(descriptor.fields).map_err(|error| ComparisonMismatch {
            code: ComparisonMismatchCode::InvalidDescriptor,
            expected: "a supported physical field mask".into(),
            observed: error.to_string(),
        })?;
    if available.contains(requested) {
        Ok(())
    } else {
        Err(ComparisonMismatch {
            code: ComparisonMismatchCode::RequiredFields,
            expected: requested.names().join(","),
            observed: available.names().join(","),
        })
    }
}

/// Metadata is inspectable for unknown/pooled inputs, but a first-party cohort
/// summary must not present either scope as one identified analysis sample.
pub(crate) fn require_named_sample(
    descriptor: &DatasetDescriptor,
) -> std::result::Result<&str, ComparisonMismatch> {
    descriptor.validate().map_err(|error| ComparisonMismatch {
        code: ComparisonMismatchCode::InvalidDescriptor,
        expected: "a valid supported portable dataset descriptor".into(),
        observed: error.to_string(),
    })?;
    match (
        &*descriptor.sample_scope.mode,
        descriptor.sample_scope.selected_sample.as_deref(),
    ) {
        ("named", Some(name)) if !name.is_empty() => Ok(name),
        (mode, _) => Err(ComparisonMismatch {
            code: ComparisonMismatchCode::SampleScope,
            expected: "one explicitly named analysis sample".into(),
            observed: mode.into(),
        }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cohort::tests::{Fixture, FixtureOptions};
    use crate::evidence::EvidenceProfile;

    #[test]
    fn different_sources_samples_and_adequate_schema_projections_compare() {
        let a = Fixture::new(FixtureOptions::default());
        let b = Fixture::new(FixtureOptions {
            sample: Some("sample-B".into()),
            fields: EvidenceFields::FULL_V1,
            ..FixtureOptions::default()
        });
        let left = a.descriptor();
        let right = b.descriptor();
        assert_ne!(left.compatibility_blake3, right.compatibility_blake3);
        assert_eq!(left.schema_version, 2);
        assert_eq!(right.schema_version, 1);
        // Source paths become unavailable before any comparison work.
        drop(a);
        drop(b);
        let first = ComparisonContract::from_descriptor(&left).unwrap();
        let second = ComparisonContract::from_descriptor(&right).unwrap();
        assert!(first.compare(&second).is_empty());
        assert_eq!(first.digest().unwrap(), second.digest().unwrap());
        let requested = EvidenceFields::DEPTHS.union(EvidenceFields::ALLELES);
        check_required_fields(&left, requested).unwrap();
        check_required_fields(&right, requested).unwrap();
        check_required_fields(&right, EvidenceFields::STRANDS).unwrap();
        let missing = check_required_fields(&left, EvidenceFields::STRANDS).unwrap_err();
        assert_eq!(missing.code, ComparisonMismatchCode::RequiredFields);
        assert_eq!(missing.expected, "strands");
        assert_eq!(require_named_sample(&left).unwrap(), "sample-A");
        assert_eq!(require_named_sample(&right).unwrap(), "sample-B");
    }

    #[test]
    fn identical_reference_bytes_compare_across_bam_and_cram_roles() {
        let bam = Fixture::new(FixtureOptions::default());
        let cram = Fixture::new(FixtureOptions {
            sample: Some("CRAM sample".into()),
            cram: true,
            cram_reference_only: true,
            ..FixtureOptions::default()
        });
        let left = bam.descriptor();
        let right = cram.descriptor();
        assert!(left.sources.iter().any(|source| source.role == "reference"));
        assert!(!right
            .sources
            .iter()
            .any(|source| source.role == "reference"));
        assert!(right
            .sources
            .iter()
            .any(|source| source.role == "cram-reference"));
        drop(bam);
        drop(cram);
        let first = ComparisonContract::from_descriptor(&left).unwrap();
        let second = ComparisonContract::from_descriptor(&right).unwrap();
        assert!(first.compare(&second).is_empty());
        assert_eq!(first.digest().unwrap(), second.digest().unwrap());
    }

    #[test]
    fn reference_content_and_filter_changes_have_separate_explanations() {
        let a = Fixture::new(FixtureOptions::default());
        let changed = Fixture::new(FixtureOptions {
            reference_base: b'G',
            profile: EvidenceProfile {
                min_base_quality: 36,
                ..EvidenceProfile::default()
            },
            ..FixtureOptions::default()
        });
        let first = ComparisonContract::from_descriptor(&a.descriptor()).unwrap();
        let second = ComparisonContract::from_descriptor(&changed.descriptor()).unwrap();
        assert_eq!(first.contigs, second.contigs); // An assembly/dictionary match is insufficient.
        let differences = first.compare(&second);
        assert_eq!(
            differences
                .iter()
                .map(|issue| issue.code)
                .collect::<Vec<_>>(),
            [
                ComparisonMismatchCode::FilterProfile,
                ComparisonMismatchCode::ReferenceContent
            ]
        );
        assert_ne!(first.digest().unwrap(), second.digest().unwrap());
        let json = serde_json::to_value(&differences).unwrap();
        assert_eq!(json[0]["code"], "filter-profile");
        assert_eq!(json[1]["code"], "reference-content");
    }

    #[test]
    fn comparison_explains_dictionary_counting_and_semantics_independently() {
        let fixture = Fixture::new(FixtureOptions::default());
        let first = ComparisonContract::from_descriptor(&fixture.descriptor()).unwrap();
        let mut other = first.clone();
        other.contigs[0].length += 1;
        other.counting_unit = "molecule".into();
        other.semantics = "future-profile".into();
        other.sampling = "downsampled".into();
        assert_eq!(
            first
                .compare(&other)
                .iter()
                .map(|issue| issue.code)
                .collect::<Vec<_>>(),
            [
                ComparisonMismatchCode::EvidenceSemantics,
                ComparisonMismatchCode::CountingUnit,
                ComparisonMismatchCode::SamplingPolicy,
                ComparisonMismatchCode::ReferenceDictionary
            ]
        );
    }

    #[test]
    fn unknown_sample_is_inspectable_but_not_an_individual_for_summary() {
        let fixture = Fixture::new(FixtureOptions {
            sample: None,
            ..FixtureOptions::default()
        });
        let descriptor = fixture.descriptor();
        ComparisonContract::from_descriptor(&descriptor).unwrap();
        let refusal = require_named_sample(&descriptor).unwrap_err();
        assert_eq!(refusal.code, ComparisonMismatchCode::SampleScope);
        assert_eq!(refusal.observed, "unknown");
    }

    #[test]
    fn invalid_descriptor_is_not_rescued_by_a_matching_projection() {
        let fixture = Fixture::new(FixtureOptions::default());
        let mut descriptor = fixture.descriptor();
        descriptor.counting_unit = "molecule".into();
        assert!(ComparisonContract::from_descriptor(&descriptor).is_err());
        let issue = check_required_fields(&descriptor, EvidenceFields::DEPTHS).unwrap_err();
        assert_eq!(issue.code, ComparisonMismatchCode::InvalidDescriptor);
    }

    #[test]
    fn same_fasta_and_dictionary_do_not_hide_a_changed_coordinate_index() {
        let fixture = Fixture::new(FixtureOptions::default());
        let descriptor = fixture.descriptor();
        let first = ComparisonContract::from_descriptor(&descriptor).unwrap();
        assert!(first.reference_index.is_some());
        let mut changed = descriptor.clone();
        let index = changed
            .sources
            .iter_mut()
            .find(|source| source.role == "reference-fai")
            .unwrap();
        index.blake3 = blake3::hash(b"different FAI offsets with the same names and lengths")
            .to_hex()
            .to_string();
        changed.compatibility_blake3 = changed.computed_compatibility_key().unwrap();
        changed.dataset_namespace = blake3::hash(
            format!(
                "rosalind-verified-dataset-v1\n{}\n{}\n",
                changed.compatibility_blake3, changed.request_blake3
            )
            .as_bytes(),
        )
        .to_hex()
        .to_string();
        let second = ComparisonContract::from_descriptor(&changed).unwrap();
        assert_eq!(first.reference, second.reference);
        assert_eq!(first.contigs, second.contigs);
        assert_eq!(
            first.compare(&second)[0].code,
            ComparisonMismatchCode::ReferenceIndex
        );
        assert_ne!(first.digest().unwrap(), second.digest().unwrap());
        let mut absent = first.clone();
        absent.reference_index = None;
        assert_eq!(
            first.compare(&absent)[0].code,
            ComparisonMismatchCode::ReferenceIndex
        );
    }
}
