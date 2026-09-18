//! Internal next-minor cohort storage preparation. No CLI or public SDK yet.
#![allow(dead_code)]

pub(crate) mod comparison;
pub(crate) mod descriptor;
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
}

pub(crate) type Result<T> = std::result::Result<T, CohortError>;
