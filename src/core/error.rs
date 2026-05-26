//! Typed error taxonomy for the Rosalind core/library layer.
//!
//! The CLI boundary maps these into `anyhow`; library code returns `CoreError`.

use thiserror::Error;

/// Errors produced by the Rosalind core/library layer.
#[derive(Debug, Error)]
pub enum CoreError {
    /// A streaming stage's working set would exceed the declared budget.
    #[error("memory budget exceeded: working set {needed} bytes > budget {budget} bytes")]
    BudgetExceeded {
        /// Bytes the stage would require.
        needed: u64,
        /// Bytes the caller permitted.
        budget: u64,
    },
    /// An input record could not be interpreted.
    #[error("malformed record: {0}")]
    MalformedRecord(String),
    /// A contig id was not present in the active `ContigSet`.
    #[error("invalid contig id {0}")]
    InvalidContig(u32),
    /// An underlying I/O failure.
    #[error("io error: {0}")]
    Io(#[from] std::io::Error),
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn budget_exceeded_displays_both_numbers() {
        let e = CoreError::BudgetExceeded { needed: 4096, budget: 1024 };
        let msg = e.to_string();
        assert!(msg.contains("4096"), "message should report needed bytes: {msg}");
        assert!(msg.contains("1024"), "message should report budget bytes: {msg}");
    }
}
