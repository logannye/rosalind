//! Implicit tree node representation
//!
//! Node = interval [left, right] ⊆ [1, T]
//! Children computed via midpoint: m = ⌊(left + right) / 2⌋
//!   Left child: [left, m]
//!   Right child: [m+1, right]

use std::fmt;

/// Tree node (implicit - just an interval)
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct TreeNode {
    /// Left block index (inclusive)
    pub left: usize,

    /// Right block index (inclusive)  
    pub right: usize,
}

impl TreeNode {
    /// Create root spanning [1, T]
    pub fn root(left: usize, right: usize) -> Self {
        Self { left, right }
    }

    /// Check if leaf (unit interval)
    #[inline]
    pub fn is_leaf(&self) -> bool {
        self.left == self.right
    }

    /// Interval length
    #[inline]
    pub fn length(&self) -> usize {
        self.right - self.left + 1
    }

    /// Compute midpoint for split
    ///
    /// Key to height compression: m = ⌊(left + right) / 2⌋
    #[inline]
    pub fn midpoint(&self) -> usize {
        (self.left + self.right) / 2
    }

    /// Get children via midpoint split
    ///
    /// Returns: ([left, mid], [mid+1, right])
    /// Geometric shrinkage: length(child) ≤ ⌈length(parent) / 2⌉
    pub fn children(&self) -> (TreeNode, TreeNode) {
        debug_assert!(!self.is_leaf(), "Leaf has no children");

        let mid = self.midpoint();
        let left_child = TreeNode {
            left: self.left,
            right: mid,
        };
        let right_child = TreeNode {
            left: mid + 1,
            right: self.right,
        };

        (left_child, right_child)
    }

    /// Get leaf block ID
    pub fn leaf_block_id(&self) -> usize {
        debug_assert!(self.is_leaf(), "Only leaves have block IDs");
        self.left
    }

    /// Compute depth from this node to any leaf
    ///
    /// Due to geometric shrinkage: depth = O(log(length))
    pub fn depth_to_leaf(&self) -> usize {
        if self.is_leaf() {
            return 0;
        }

        // Simulate recursive descent, counting splits until length = 1
        let mut node = *self;
        let mut depth = 0;

        while !node.is_leaf() {
            let (left_child, _) = node.children();
            node = left_child; // Follow left path (any path works)
            depth += 1;
        }

        depth
    }
}

impl fmt::Display for TreeNode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_leaf() {
            write!(f, "[{}]", self.left)
        } else {
            write!(f, "[{}, {}]", self.left, self.right)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_midpoint_split() {
        let node = TreeNode::root(1, 100);
        let (left, right) = node.children();

        assert_eq!(left.left, 1);
        assert_eq!(left.right, 50);
        assert_eq!(right.left, 51);
        assert_eq!(right.right, 100);
    }

    #[test]
    fn test_geometric_shrinkage() {
        // Verify geometric shrinkage: length halves at each level
        let root = TreeNode::root(1, 128);
        let mut node = root;
        let mut lengths = vec![node.length()];

        while !node.is_leaf() {
            let (left, _) = node.children();
            node = left;
            lengths.push(node.length());
        }

        // Verify: ℓ(d+1) ≤ ⌈ℓ(d)/2⌉ at each step
        for window in lengths.windows(2) {
            let parent_len = window[0];
            let child_len = window[1];
            let expected_max = (parent_len + 1) / 2; // ⌈parent_len/2⌉
            assert!(
                child_len <= expected_max,
                "Child length {} should be <= ⌈parent {}/2⌉ = {}",
                child_len,
                parent_len,
                expected_max
            );
        }

        // Verify: depth = O(log T)
        let depth = root.depth_to_leaf();
        let log_bound = (128_f64).log2().ceil() as usize;
        assert!(
            depth <= log_bound + 1,
            "Depth {} should be <= log2({}) + 1 = {}",
            depth,
            128,
            log_bound + 1
        );
    }
}
