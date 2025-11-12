//! Space tracking for verification

/// Tracks space usage during simulation
#[derive(Debug)]
pub struct SpaceTracker {
    /// Current space used
    current: usize,

    /// Maximum seen
    max: usize,

    /// Track frame sizes for correct popping
    frame_sizes: Vec<usize>,

    /// Enable detailed profiling
    #[allow(dead_code)]
    profile_enabled: bool,

    /// Profile data (if enabled)
    profile: Option<super::SpaceProfile>,

    /// Current stack depth (for profiling)
    stack_depth: usize,
}

impl SpaceTracker {
    /// Create new tracker
    pub fn new(profile_enabled: bool) -> Self {
        Self {
            current: 0,
            max: 0,
            frame_sizes: Vec::new(),
            profile_enabled,
            profile: if profile_enabled {
                Some(super::SpaceProfile {
                    max_space: 0,
                    timeline: Vec::new(),
                    leaf_buffer_max: 0,
                    stack_depth_max: 0,
                    ledger_size: 0,
                })
            } else {
                None
            },
            stack_depth: 0,
        }
    }

    /// Allocate space (e.g., leaf buffer)
    pub fn allocate_leaf_buffer(&mut self, size: usize) {
        self.current += size;
        self.update_max();

        if let Some(ref mut p) = self.profile {
            p.leaf_buffer_max = p.leaf_buffer_max.max(size);
        }
    }

    /// Free space
    pub fn free_leaf_buffer(&mut self, size: usize) {
        self.current = self.current.saturating_sub(size);
    }

    /// Push stack frame (O(1) per level)
    pub fn push_stack_frame(&mut self, token_size: usize) {
        self.current += token_size;
        self.frame_sizes.push(token_size);
        self.stack_depth += 1;
        self.update_max();

        if let Some(ref mut p) = self.profile {
            p.stack_depth_max = p.stack_depth_max.max(self.stack_depth);
        }
    }

    /// Pop stack frame
    pub fn pop_stack_frame(&mut self) {
        if let Some(size) = self.frame_sizes.pop() {
            self.current = self.current.saturating_sub(size);
            self.stack_depth = self.stack_depth.saturating_sub(1);
        }
    }

    /// Allocate ledger space
    pub fn allocate_ledger(&mut self, size: usize) {
        self.current += size;
        self.update_max();

        if let Some(ref mut p) = self.profile {
            p.ledger_size = size;
        }
    }

    fn update_max(&mut self) {
        self.max = self.max.max(self.current);

        if let Some(ref mut p) = self.profile {
            p.max_space = self.max;
        }
    }

    /// Get maximum space used
    pub fn max_space_used(&self) -> usize {
        self.max
    }

    /// Take profile (consumes tracker)
    pub fn take_profile(&mut self) -> Option<super::SpaceProfile> {
        self.profile.take()
    }
}
