//! Generic compressed evaluation framework built on top of the O(√t) engine.
//!
//! Provides abstractions for block-respecting algorithms that reuse the
//! height-compressed tree, streaming ledger, and space tracker without
//! relying on the Turing machine implementation.

mod compressed_eval;

pub use compressed_eval::{
    BlockContext, BlockProcessor, CompressedEvaluator, EvaluationResult, EvaluatorConfig,
    FrameworkError,
};
