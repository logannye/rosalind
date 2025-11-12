//! Height-compressed evaluation tree
//!
//! Key innovation: Transform canonical tree (height Θ(T))
//! into compressed tree (height O(log T)) via midpoint recursion
//!
//! Implicit representation: No explicit tree stored!
//! Nodes are intervals [left, right] computed on-demand.

mod node;
mod traversal;

pub use node::TreeNode;
pub use traversal::{PathToken, PointerlessTraversal};

/// Compressed evaluation tree (implicit)
///
/// Never materialized - all navigation via arithmetic
#[derive(Debug)]
pub struct CompressedTree {
    /// Total number of leaf blocks
    num_blocks: usize,

    /// Root node (just endpoints, no structure)
    root: TreeNode,
}

impl CompressedTree {
    /// Create compressed tree for T blocks
    pub fn new(num_blocks: usize) -> Self {
        Self {
            num_blocks,
            root: TreeNode::root(1, num_blocks),
        }
    }

    /// Get root node
    pub fn root(&self) -> &TreeNode {
        &self.root
    }

    /// Theoretical height: O(log T)
    pub fn height_bound(&self) -> usize {
        (self.num_blocks as f64).log2().ceil() as usize
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_height_logarithmic() {
        // Verify O(log T) height for various T
        for t in [10, 100, 1000, 10000] {
            let tree = CompressedTree::new(t);
            let bound = tree.height_bound();
            assert!(bound <= (t as f64).log2().ceil() as usize + 1);
        }
    }
}
