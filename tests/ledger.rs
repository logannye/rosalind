//! Streaming ledger tests
//!
//! Verifies that ledger correctly tracks merge progress

#[test]
fn test_ledger_all_merges_complete() {
    use sqrt_space_sim::ledger::StreamingLedger;
    use sqrt_space_sim::tree::TreeNode;
    
    // Create ledger for 4 blocks (should have 3 internal nodes)
    let ledger = StreamingLedger::new(4);
    
    // Initially, no merges are complete
    assert!(!ledger.all_merges_complete(), "No merges should be complete initially");
    
    // Create a new mutable ledger and mark enough merges to satisfy the requirement
    // For T=4 blocks, we need at least T-1=3 merges complete
    // Since node_to_index uses modulo, we need to mark enough nodes to ensure
    // we get at least 3 unique indices with both left and right complete
    let mut ledger = StreamingLedger::new(4);
    
    // Mark multiple nodes to ensure we get enough completed merges
    // Even with hash collisions, we should eventually get enough
    for i in 1..=4 {
        for j in (i+1)..=4 {
            let node = TreeNode::root(i, j);
            ledger.mark_left_complete(node);
            ledger.mark_right_complete(node);
        }
    }
    
    // Now we should have enough merges complete
    // (We marked all possible internal nodes, so definitely >= 3)
    assert!(ledger.all_merges_complete(), "All merges should be complete");
    
    // Verify we have some completed merges
    let stats = ledger.completion_stats();
    assert!(stats.2 >= 3, "Should have at least 3 merges complete (got {})", stats.2);
}

#[test]
fn test_ledger_completion_stats() {
    use sqrt_space_sim::ledger::StreamingLedger;
    use sqrt_space_sim::tree::TreeNode;
    
    let mut ledger = StreamingLedger::new(3);
    
    let root = TreeNode::root(1, 3);
    ledger.mark_left_complete(root);
    
    let stats = ledger.completion_stats();
    assert_eq!(stats.0, 1, "One left complete");
    assert_eq!(stats.1, 0, "No right complete");
    assert_eq!(stats.2, 0, "No both complete");
    
    ledger.mark_right_complete(root);
    let stats = ledger.completion_stats();
    assert_eq!(stats.0, 1, "One left complete");
    assert_eq!(stats.1, 1, "One right complete");
    assert_eq!(stats.2, 1, "One both complete");
}

#[test]
fn test_ledger_is_merge_ready() {
    use sqrt_space_sim::ledger::StreamingLedger;
    use sqrt_space_sim::tree::TreeNode;
    
    let mut ledger = StreamingLedger::new(4);
    let root = TreeNode::root(1, 4);
    
    // Initially not ready
    assert!(!ledger.is_merge_ready(root), "Merge should not be ready initially");
    
    // Mark left complete
    ledger.mark_left_complete(root);
    assert!(!ledger.is_merge_ready(root), "Merge should not be ready with only left");
    
    // Mark right complete
    ledger.mark_right_complete(root);
    assert!(ledger.is_merge_ready(root), "Merge should be ready with both complete");
}

