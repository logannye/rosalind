//! Internal implementation of the next-minor cohort CLI/Python preview.
//! Generic cohort consumers are not a public Rust SDK contract.
#![allow(dead_code)]

pub(crate) mod artifact;
pub(crate) mod comparison;
pub(crate) mod descriptor;
pub(crate) mod encoding;
pub(crate) mod extend;
pub(crate) mod query;
pub(crate) mod runtime;
pub(crate) mod store;
pub(crate) mod summary;

#[cfg(test)]
pub(crate) mod tests;

#[derive(Debug, thiserror::Error)]
pub(crate) enum CohortError {
    #[error("cohort I/O: {0}")]
    Io(#[from] std::io::Error),
    #[error("invalid cohort: {0}")]
    Corrupt(String),
    #[error("incompatible cohort: {0}")]
    Incompatible(String),
    #[error("cohort metadata or resource envelope: {0}")]
    Limit(String),
    #[error(transparent)]
    Dataset(#[from] crate::dataset::DatasetError),
    #[error(transparent)]
    Evidence(#[from] crate::evidence::EvidenceError),
    #[error(transparent)]
    Artifact(#[from] crate::evidence::EvidenceArtifactError),
}

pub(crate) type Result<T> = std::result::Result<T, CohortError>;
