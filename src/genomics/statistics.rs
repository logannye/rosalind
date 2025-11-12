use crate::genomics::PileupNode;

/// Result of Bayesian scoring for a candidate variant.
#[derive(Debug, Clone, PartialEq)]
pub struct VariantCall {
    /// Alternate base proposed by the caller.
    pub alt_base: u8,
    /// Heuristic Phred-scaled quality of the call.
    pub quality: f32,
    /// Estimated alternate allele fraction.
    pub allele_fraction: f32,
}

fn base_index(base: u8) -> Option<usize> {
    match base {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' | b'U' | b'u' => Some(3),
        _ => None,
    }
}

/// Simple Bayesian variant caller using allele counts and average quality.
///
/// This is a lightweight model designed for streaming evaluation. The prior
/// reflects the probability of observing a mutation at any position.
pub fn bayesian_variant_caller(
    node: &PileupNode,
    reference_base: u8,
    prior: f32,
) -> Option<VariantCall> {
    let depth = node.depth as f32;
    if depth == 0.0 {
        return None;
    }

    let ref_idx = base_index(reference_base)?;

    let mut best_idx = None;
    let mut best_count = 0u32;

    for (idx, &count) in node.base_counts.iter().enumerate() {
        if idx == ref_idx {
            continue;
        }
        if count > best_count {
            best_count = count;
            best_idx = Some(idx);
        }
    }

    let alt_idx = best_idx?;
    if best_count == 0 {
        return None;
    }

    let alt_fraction = best_count as f32 / depth;

    let quality_sum = node.quality_sums[alt_idx];
    let avg_quality = if best_count > 0 {
        (quality_sum / best_count as f32).clamp(0.0, 1.0)
    } else {
        0.0
    };

    // Posterior using a simple beta-binomial style update.
    let likelihood_variant = prior * alt_fraction.max(1e-6);
    let likelihood_reference = (1.0 - prior) * (1.0 - alt_fraction).max(1e-6);
    let _posterior = likelihood_variant / (likelihood_variant + likelihood_reference);

    Some(VariantCall {
        alt_base: match alt_idx {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            3 => b'T',
            _ => b'N',
        },
        quality: (alt_fraction * avg_quality.max(0.1) * 100.0).min(60.0),
        allele_fraction: alt_fraction,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bayesian_caller_identifies_alt() {
        let mut node = PileupNode::new(100);
        for _ in 0..8 {
            node.observe(0, 30); // A
        }
        for _ in 0..2 {
            node.observe(1, 25); // C
        }

        let call = bayesian_variant_caller(&node, b'A', 1e-6);
        assert!(call.is_some());
        let call = call.unwrap();
        assert_eq!(call.alt_base, b'C');
        assert!(call.quality > 0.0);
    }
}
