//! Constant-degree combiner for merge operations
//!
//! F_{parent}(x) = G(F_left(Ax), F_right(Bx), x)
//! where A, B are affine maps and G is constant-degree

use super::{FiniteField, PolynomialEncoding, EvaluationGrid};
use crate::{blocking::{BlockSummary, InterfaceChecker}, space::SpaceTracker, SimulationError};

/// Combiner for merging interval summaries
#[derive(Debug)]
pub struct Combiner {
    field: FiniteField,
}

impl Combiner {
    /// Create combiner for field
    pub fn new(field: &FiniteField) -> Self {
        // Clone field - we need it for operations
        // In practice, we'd store a reference or clone the characteristic
        let characteristic = match field.size() {
            256 => 8,
            size => size.ilog2() as u8,
        };
        Self {
            field: FiniteField::new(characteristic),
        }
    }
    
    /// Merge two adjacent summaries
    ///
    /// Space: O(1) field elements (internal node workspace)
    pub fn merge(
        &self,
        left: &BlockSummary,
        right: &BlockSummary,
        grid: &EvaluationGrid,
        _tracker: &mut SpaceTracker,
    ) -> Result<BlockSummary, SimulationError> {
        // 1. Check interface consistency (exact replay at junction)
        // CRITICAL: This must be done before algebraic merge!
        if !InterfaceChecker::check(left, right)? {
            return Err(SimulationError::InterfaceCheckFailed(left.block_id));
        }
        
        // 2. Extract polynomial encodings from summaries
        // Create encodings from finite-state projections
        let left_encoding = self.extract_encoding(left, grid);
        let right_encoding = self.extract_encoding(right, grid);
        
        // 3. For each grid point x, compute F_parent(x) = G(F_left(Ax), F_right(Bx), x)
        let mut parent_encoding = PolynomialEncoding::new(1); // Output dimension = 1
        
        for x in grid.points() {
            // Apply affine maps
            let ax = self.apply_affine_a(x);
            let bx = self.apply_affine_b(x);
            
            // Evaluate child encodings
            let left_val = left_encoding.eval(&ax);
            let right_val = right_encoding.eval(&bx);
            
            // Compute parent encoding
            let parent_val = self.combine_g(&left_val, &right_val, x);
            parent_encoding.set_eval(x.clone(), parent_val);
        }
        
        // 4. Create merged summary
        // Merge entry/exit states and heads
        let merged = BlockSummary::new(
            left.block_id, // Use left's block ID
            left.entry_state(),
            right.exit_state(),
            left.entry_heads().to_vec(),
            right.exit_heads().to_vec(),
            // Merge movement logs (simplified - in practice would combine)
            left.movement_log().clone(),
            // Merge windows
            self.merge_windows(left.windows(), right.windows()),
        );
        
        Ok(merged)
    }
    
    /// Extract polynomial encoding from block summary
    /// Encodes finite-state projection (control state, flags)
    ///
    /// For each grid point, encodes:
    /// - Entry and exit states (XOR'd together)
    /// - Optionally head positions (truncated to fit in field)
    fn extract_encoding(&self, summary: &BlockSummary, grid: &EvaluationGrid) -> PolynomialEncoding {
        // Create encoding with output dimension 1 (single field element)
        let mut encoding = PolynomialEncoding::new(1);
        
        // Encode summary state information
        let entry_state = summary.entry_state() as u8;
        let exit_state = summary.exit_state() as u8;
        
        // Combine states: entry XOR exit (XOR is addition in GF(2))
        let state_val = self.field.add(entry_state, exit_state);
        
        // Optionally encode head positions (first few positions, truncated)
        // For simplicity, we encode the sum of first head position's low bits
        let head_val = if !summary.entry_heads().is_empty() {
            // Take low 8 bits of first head position
            (summary.entry_heads()[0] & 0xFF) as u8
        } else {
            0
        };
        
        // Combine state and head information
        let combined_val = self.field.add(state_val, head_val);
        
        // For each grid point, set the same encoded value
        // In a more sophisticated implementation, we could use the grid point
        // to create different encodings, but for now we use a constant encoding
        // that captures the essential state information
        for x in grid.points() {
            // Encode: state information + (optionally) grid-dependent component
            // For now, use constant encoding: state_val
            // In full implementation, could interpolate: f(x) = state_val + linear(x)
            let encoded_value = if x.is_empty() {
                combined_val
            } else {
                // Use first component of grid point to add variation
                let grid_component = x[0] % 16; // Keep small to avoid overflow
                self.field.add(combined_val, grid_component)
            };
            
            encoding.set_eval(x.clone(), vec![encoded_value]);
        }
        
        encoding
    }
    
    /// Merge window bounds
    fn merge_windows(
        &self,
        left_windows: &[crate::blocking::WindowBounds],
        right_windows: &[crate::blocking::WindowBounds],
    ) -> Vec<crate::blocking::WindowBounds> {
        left_windows
            .iter()
            .zip(right_windows.iter())
            .map(|(l, r)| crate::blocking::WindowBounds {
                left: l.left.min(r.left),
                right: l.right.max(r.right),
            })
            .collect()
    }
    
    /// Apply affine map A: x → Ax
    /// Simplified: identity map (can be customized)
    fn apply_affine_a(&self, x: &[u8]) -> Vec<u8> {
        x.to_vec() // Identity map for simplicity
    }
    
    /// Apply affine map B: x → Bx
    /// Simplified: identity map (can be customized)
    fn apply_affine_b(&self, x: &[u8]) -> Vec<u8> {
        x.to_vec() // Identity map for simplicity
    }
    
    /// Constant-degree combiner G
    /// G(left, right, x) = left + right + x[0] (simple linear combination)
    fn combine_g(&self, left: &[u8], right: &[u8], x: &[u8]) -> Vec<u8> {
        let left_val = left.get(0).copied().unwrap_or(0);
        let right_val = right.get(0).copied().unwrap_or(0);
        let x_val = x.get(0).copied().unwrap_or(0);
        
        // Simple linear combination: left + right + x[0]
        let result = self.field.add(self.field.add(left_val, right_val), x_val);
        vec![result]
    }
}