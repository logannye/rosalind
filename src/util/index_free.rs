//! Index-free iteration (Lemma 4.2 in paper)
//!
//! Replace O(log b) counters with O(1) marker-based scanning
//! Trade time for space: rescan from marker each iteration

/// Index-free iterator using marker scanning
///
/// Space: O(1) cells (marker position only)
/// Time: O(m^2) for m items (acceptable for space-only bound)
#[allow(dead_code)]
#[derive(Debug)]
pub struct IndexFreeIterator<'a, T> {
    #[allow(dead_code)]
    items: &'a [T],
    #[allow(dead_code)]
    marker_pos: usize,
    #[allow(dead_code)]
    current_item: usize,
}

impl<'a, T> IndexFreeIterator<'a, T> {
    /// Create iterator with marker at start
    pub fn new(items: &'a [T]) -> Self {
        Self {
            items,
            marker_pos: 0,
            current_item: 0,
        }
    }
    
    /// Get next item (rescans from marker)
    ///
    /// Key idea: Don't store counter, reconstruct position by scanning
    pub fn next_item(&mut self) -> Option<&'a T> {
        // TODO:
        // 1. Start at marker_pos
        // 2. Skip current_item delimiters via constant-state recognition
        // 3. Return item
        // 4. Increment current_item (in O(1) state)
        
        // No O(log m) counter stored!
        unimplemented!()
    }
}

/// Delimiter-based scanning (constant-state machine)
#[allow(dead_code)]
fn scan_to_next_delimiter(_pos: usize, _data: &[u8]) -> usize {
    // TODO: Scan using constant-state automaton
    // Recognizes item boundaries without arithmetic
    unimplemented!()
}