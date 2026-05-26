//! Core types — the lingua franca shared by every Rosalind layer.
//!
//! Note: this module is named `core`; inside the crate always reference it as
//! `crate::core::…`. Reach the std `core` crate (rarely needed) as `::core::…`.

pub mod error;

pub use error::CoreError;
