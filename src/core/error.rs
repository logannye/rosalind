//! Typed error taxonomy for the Rosalind core/library layer.
//!
//! The CLI boundary maps these into `anyhow`; library code returns `CoreError`.

use thiserror::Error;

/// Errors produced by the Rosalind core/library layer.
#[derive(Debug, Error)]
pub enum CoreError {
    /// The active native analysis was cooperatively cancelled.
    #[error("analysis cancelled")]
    Cancelled,
    /// A streaming stage's working set would exceed the declared budget.
    #[error("memory budget exceeded: working set {needed} bytes > budget {budget} bytes")]
    BudgetExceeded {
        /// Bytes the stage would require.
        needed: u64,
        /// Bytes the caller permitted.
        budget: u64,
    },
    /// Exact pileup would require more simultaneous reads than its declared
    /// capacity. No read is sampled or evicted to make the analysis fit.
    #[error("pileup capacity exceeded at contig {contig}, zero-based position {position}: {required} active reads require capacity greater than --max-depth {capacity}; raise --max-depth and rerun")]
    CapacityExceeded {
        /// Contig identifier of the first overflowing locus.
        contig: u32,
        /// Zero-based position of the first overflowing locus.
        position: u32,
        /// Maximum number of active reads permitted.
        capacity: u32,
        /// Number of active reads required to retain all eligible evidence.
        required: u64,
    },
    /// An input record could not be interpreted.
    #[error("malformed record: {0}")]
    MalformedRecord(String),
    /// A contig id was not present in the active `ContigSet`.
    #[error("invalid contig id {0}")]
    InvalidContig(u32),
    /// A read longer than the declared `--max-read-len` voids the predicted
    /// memory envelope (raised only under `--enforce`).
    #[error("read length {len} exceeds declared --max-read-len {declared}; raise --max-read-len or drop --enforce")]
    ReadExceedsDeclaredLength {
        /// The offending read's sequence length.
        len: u32,
        /// The declared cap.
        declared: u32,
    },
    /// An underlying I/O failure.
    #[error("io error: {0}")]
    Io(#[from] std::io::Error),
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn budget_exceeded_displays_both_numbers() {
        let e = CoreError::BudgetExceeded {
            needed: 4096,
            budget: 1024,
        };
        let msg = e.to_string();
        assert!(
            msg.contains("4096"),
            "message should report needed bytes: {msg}"
        );
        assert!(
            msg.contains("1024"),
            "message should report budget bytes: {msg}"
        );
    }

    #[test]
    fn io_error_converts_via_from() {
        let io = std::io::Error::new(std::io::ErrorKind::NotFound, "missing");
        let e: CoreError = io.into();
        assert!(matches!(e, CoreError::Io(_)));
    }
}
