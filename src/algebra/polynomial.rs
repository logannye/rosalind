//! Low-degree polynomial extensions
//!
//! For interval summary Σ(I), encode finite-state projection as:
//! F_I: 𝔽^m → 𝔽^q with degree ≤ D in each variable
//!
//! Key: D, m, q all O(1) constants

use super::FiniteField;
use std::collections::HashMap;

/// Polynomial encoding of finite-state projection
///
/// Only encodes: control state + constant-size flags/tags
/// All b-dependent data handled by exact replay!
#[derive(Debug, Clone)]
pub struct PolynomialEncoding {
    /// Evaluation table: grid_point → field_element^q
    /// Grid point is Vec<u8> of length m
    /// Values are Vec<u8> of length q (output dimension)
    evaluations: HashMap<Vec<u8>, Vec<u8>>,

    /// Output dimension q (constant)
    output_dim: usize,
}

impl PolynomialEncoding {
    /// Create new encoding with empty evaluations
    pub fn new(output_dim: usize) -> Self {
        Self {
            evaluations: HashMap::new(),
            output_dim,
        }
    }

    /// Create default encoding (for testing)
    pub fn default() -> Self {
        Self::new(1) // Default: single output dimension
    }

    /// Set evaluation at grid point
    pub fn set_eval(&mut self, point: Vec<u8>, value: Vec<u8>) {
        self.evaluations.insert(point, value);
    }

    /// Evaluate at grid point
    pub fn eval(&self, point: &[u8]) -> Vec<u8> {
        self.evaluations
            .get(point)
            .cloned()
            .unwrap_or_else(|| vec![0; self.output_dim])
    }

    /// Space usage: O(1) cells
    pub fn space_usage(&self) -> usize {
        // |grid| × q × cell_size = O(1) since grid and q are constants
        self.evaluations.len() * self.output_dim
    }
}

/// Constant-size evaluation grid X ⊆ 𝔽^m
///
/// Size: (D+1)^m for degree bound D, dimension m
/// Both D and m are O(1) constants
/// Typical: D=2, m=3 → |X| = 27 points
#[derive(Debug)]
pub struct EvaluationGrid {
    /// Grid points: tensor product {0,1,...,D}^m
    points: Vec<Vec<u8>>,

    /// Degree bound D (constant)
    degree_bound: usize,

    /// Dimension m (constant)
    dimension: usize,
}

impl EvaluationGrid {
    /// Create grid for field
    /// Default: D=2, m=3 → 27 points
    pub fn new(field: &FiniteField) -> Self {
        Self::with_params(field, 2, 3) // D=2, m=3
    }

    /// Create grid with specific parameters
    pub fn with_params(_field: &FiniteField, degree_bound: usize, dimension: usize) -> Self {
        // Generate tensor grid {0,1,...,D}^m
        let mut points = Vec::new();
        Self::generate_grid_points(&mut points, Vec::new(), degree_bound, dimension);

        Self {
            points,
            degree_bound,
            dimension,
        }
    }

    /// Recursively generate grid points
    fn generate_grid_points(
        points: &mut Vec<Vec<u8>>,
        current: Vec<u8>,
        degree_bound: usize,
        dimension: usize,
    ) {
        if current.len() == dimension {
            points.push(current);
            return;
        }

        for i in 0..=degree_bound {
            let mut next = current.clone();
            next.push(i as u8);
            Self::generate_grid_points(points, next, degree_bound, dimension);
        }
    }

    /// Get all grid points
    pub fn points(&self) -> &[Vec<u8>] {
        &self.points
    }

    /// Grid size (constant)
    pub fn size(&self) -> usize {
        (self.degree_bound + 1).pow(self.dimension as u32)
    }
}
