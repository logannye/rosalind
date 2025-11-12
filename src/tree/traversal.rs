//! Pointerless DFS traversal
//!
//! Key innovation: Store only O(1) bits per level, not O(log b)
//! Path token = (node_type, direction) = 2 bits
//! Endpoints recomputed on-demand from root

use super::TreeNode;

/// Path token stored at each recursion level
///
/// CRITICAL: Only 2 bits per level!
/// This eliminates the O(log b) per-level factor
#[derive(Debug, Clone, Copy)]
pub struct PathToken {
    /// Node type (1 bit)
    pub node_type: NodeType,

    /// Which child we're in (1 bit)
    pub direction: Direction,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum NodeType {
    /// Split node (midpoint recursion)
    Split,

    /// Combiner node (merge operation)
    Combiner,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Direction {
    /// Currently processing left child
    Left,

    /// Currently processing right child
    Right,
}

/// Pointerless traversal state
///
/// Stack depth: O(log T)
/// Per-level: O(1) bits
/// Total stack: O(log T) cells
#[derive(Debug)]
pub struct PointerlessTraversal {
    /// Stack of path tokens
    path_stack: Vec<PathToken>,

    /// Current node (recomputed, not stored long-term)
    #[allow(dead_code)]
    current: Option<TreeNode>,
}

impl PointerlessTraversal {
    /// Create new traversal starting at root
    pub fn new(_root: TreeNode) -> Self {
        Self {
            path_stack: Vec::new(),
            current: None,
        }
    }

    /// Push level onto stack (O(1) bits)
    pub fn push_level(&mut self, node_type: NodeType, direction: Direction) {
        self.path_stack.push(PathToken {
            node_type,
            direction,
        });
    }

    /// Pop level from stack
    pub fn pop_level(&mut self) -> Option<PathToken> {
        self.path_stack.pop()
    }

    /// Recompute current node endpoints from path
    ///
    /// Space: O(log T) scratch (additive, not per-level!)
    /// Time: O(log T) arithmetic operations
    ///
    /// This is the key to constant per-level storage!
    pub fn recompute_endpoints(&self, root: TreeNode) -> TreeNode {
        // Start with root endpoints
        let mut node = root;

        // For each token in path_stack, navigate down the tree
        for token in &self.path_stack {
            if node.is_leaf() {
                break; // Can't go deeper
            }

            let (left_child, right_child) = node.children();

            // Select left or right child based on direction
            match token.direction {
                Direction::Left => {
                    node = left_child;
                }
                Direction::Right => {
                    node = right_child;
                }
            }
        }

        node
    }

    /// Stack depth (number of active levels)
    pub fn depth(&self) -> usize {
        self.path_stack.len()
    }

    /// Space usage: O(1) per level × depth
    pub fn space_usage(&self) -> usize {
        // Each token = 2 bits, but we count in cells
        // Over fixed alphabet: ceiling division
        (self.path_stack.len() * 2 + 7) / 8
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_recompute_endpoints() {
        let root = TreeNode::root(1, 100);
        let mut traversal = PointerlessTraversal::new(root);

        // Push left direction
        traversal.push_level(NodeType::Split, Direction::Left);
        let node = traversal.recompute_endpoints(root);
        assert_eq!(node.left, 1);
        assert_eq!(node.right, 50); // Midpoint split

        // Push left again (going deeper into left subtree)
        traversal.push_level(NodeType::Split, Direction::Left);
        let node = traversal.recompute_endpoints(root);
        let (expected_left, _) = root.children();
        let (expected_left_left, _) = expected_left.children();
        assert_eq!(node.left, expected_left_left.left);
        assert_eq!(node.right, expected_left_left.right);

        // Now test right direction
        let mut traversal2 = PointerlessTraversal::new(root);
        traversal2.push_level(NodeType::Split, Direction::Right);
        let node = traversal2.recompute_endpoints(root);
        let (_, expected_right) = root.children();
        assert_eq!(node.left, expected_right.left);
        assert_eq!(node.right, expected_right.right);
    }
}
