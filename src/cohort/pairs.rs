//! Explicit paired research comparisons with exact right-minus-left fractions.

use super::summary::{CandidateObservation, ExactFraction, ObservationStatus};
use super::{CohortError, Result};
use std::collections::{BTreeMap, BTreeSet};

/// User-supplied direction is scientifically meaningful: right minus left.
/// Names, subjects, groups and timepoints never imply a pairing.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct PairSpec {
    pub id: String,
    pub left: String,
    pub right: String,
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct PairLimits {
    pub max_pairs: usize,
    pub max_members: usize,
    /// Conservative retained pair/member lookup and copied pair-ID reservation.
    /// Managed planning and native traversal include this in admission accounting.
    pub max_metadata_bytes: u64,
}
impl Default for PairLimits {
    fn default() -> Self {
        Self {
            max_pairs: 65_536,
            max_members: 65_536,
            max_metadata_bytes: 64 << 20,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct ResolvedPair {
    pub id: String,
    pub left_member_index: usize,
    pub right_member_index: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct PairPlan {
    /// Preserve explicit table order and left/right direction.
    pub pairs: Vec<ResolvedPair>,
    pub metadata_bytes: u64,
}

/// Resolve explicit unique pair IDs against an already selected member list.
/// A member may participate in several pairs; each output remains a distinct
/// user-requested comparison. Self-pairs refuse as uninformative v1 requests.
pub(crate) fn plan_pairs(
    pairs: &[PairSpec],
    member_ids: &[String],
    limits: PairLimits,
) -> Result<PairPlan> {
    if limits.max_pairs == 0 || limits.max_members == 0 || limits.max_metadata_bytes == 0 {
        return Err(limit("pair envelopes must be positive"));
    }
    if pairs.len() > limits.max_pairs || member_ids.len() > limits.max_members {
        return Err(limit("pair/member cardinality exceeds its envelope"));
    }
    // Check input text and allocations before constructing lookup/output state.
    let mut bytes = 4096u64;
    for id in member_ids {
        validate_id(id)?;
        bytes = bytes.saturating_add(512 + id.len() as u64 * 4);
    }
    for pair in pairs {
        for id in [&pair.id, &pair.left, &pair.right] {
            validate_id(id)?;
            bytes = bytes.saturating_add(id.len() as u64 * 4);
        }
        bytes = bytes.saturating_add(1024);
    }
    if bytes > limits.max_metadata_bytes {
        return Err(limit("pair/member metadata exceeds its byte envelope"));
    }
    let mut members = BTreeMap::new();
    for (index, id) in member_ids.iter().enumerate() {
        if members.insert(id.as_str(), index).is_some() {
            return Err(incompatible("duplicate member ID in pair scope"));
        }
    }
    let mut ids = BTreeSet::new();
    let mut resolved = Vec::with_capacity(pairs.len());
    for pair in pairs {
        if !ids.insert(pair.id.as_str()) {
            return Err(incompatible("pair IDs must be unique"));
        }
        let left = *members
            .get(pair.left.as_str())
            .ok_or_else(|| incompatible("left pair member is absent from the selected cohort"))?;
        let right = *members
            .get(pair.right.as_str())
            .ok_or_else(|| incompatible("right pair member is absent from the selected cohort"))?;
        if left == right {
            return Err(incompatible(
                "self-pairs are unsupported; select two distinct members",
            ));
        }
        resolved.push(ResolvedPair {
            id: pair.id.clone(),
            left_member_index: left,
            right_member_index: right,
        });
    }
    Ok(PairPlan {
        pairs: resolved,
        metadata_bytes: bytes,
    })
}

/// Exact signed rational. Magnitude and denominator can exceed signed 128-bit
/// range; never convert them to i128 or f64. Zero has a nonnegative sign. The
/// denominator is the original depth product, not a hidden normalized fraction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct ExactDifference {
    pub negative: bool,
    pub magnitude: u128,
    pub denominator: u128,
}

/// Both side observations are preserved verbatim after semantic validation.
/// A missing side has null technical eligibility, so the conjunction is also
/// unknown. Low positive depth still permits an exact observed AF difference.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct PairedObservation {
    pub left: CandidateObservation,
    pub right: CandidateObservation,
    pub right_minus_left: Option<ExactDifference>,
    pub both_depth_eligible: Option<bool>,
}

pub(crate) fn compare_observations(
    left: CandidateObservation,
    right: CandidateObservation,
) -> Result<PairedObservation> {
    validate_observation(left)?;
    validate_observation(right)?;
    if left.min_callable_depth != right.min_callable_depth {
        return Err(incompatible(
            "paired observations must use the same technical depth screen",
        ));
    }
    let difference = match (left.observed_alt_fraction, right.observed_alt_fraction) {
        (Some(left), Some(right)) => {
            // u64 x u64 fits u128. Sign+magnitude avoids signed-overflow even
            // when a cross-product is larger than i128::MAX.
            let left_cross = u128::from(left.numerator) * u128::from(right.denominator);
            let right_cross = u128::from(right.numerator) * u128::from(left.denominator);
            Some(ExactDifference {
                negative: right_cross < left_cross,
                magnitude: right_cross.abs_diff(left_cross),
                denominator: u128::from(left.denominator) * u128::from(right.denominator),
            })
        }
        _ => None,
    };
    Ok(PairedObservation {
        left,
        right,
        right_minus_left: difference,
        both_depth_eligible: left
            .depth_eligible
            .zip(right.depth_eligible)
            .map(|(a, b)| a && b),
    })
}

fn validate_observation(observation: CandidateObservation) -> Result<()> {
    if observation.min_callable_depth == 0 {
        return Err(incompatible("technical depth screen must be positive"));
    }
    let invalid = || CohortError::Corrupt("inconsistent paired candidate observation".into());
    match observation.status {
        ObservationStatus::Unmeasured => {
            if observation.callable_depth.is_some()
                || observation.alt_count.is_some()
                || observation.depth_eligible.is_some()
                || observation.alt_supported.is_some()
                || observation.observed_alt_fraction.is_some()
            {
                return Err(invalid());
            }
        }
        ObservationStatus::Observed => {
            let (Some(depth), Some(alt)) = (observation.callable_depth, observation.alt_count)
            else {
                return Err(invalid());
            };
            let eligible = depth >= observation.min_callable_depth;
            let fraction = (depth != 0).then_some(ExactFraction {
                numerator: alt,
                denominator: depth,
            });
            if alt > depth
                || observation.depth_eligible != Some(eligible)
                || observation.alt_supported != Some(eligible && alt > 0)
                || observation.observed_alt_fraction != fraction
            {
                return Err(invalid());
            }
        }
    }
    Ok(())
}

fn validate_id(id: &str) -> Result<()> {
    if id.is_empty() || id.len() > 256 || id.chars().any(char::is_control) {
        return Err(incompatible("pair and member IDs must be nonempty, at most 256 UTF-8 bytes, and contain no controls"));
    }
    Ok(())
}
fn incompatible(message: &str) -> CohortError {
    CohortError::Incompatible(message.into())
}
fn limit(message: &str) -> CohortError {
    CohortError::Limit(message.into())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn observed(depth: u64, alt: u64, threshold: u64) -> CandidateObservation {
        CandidateObservation {
            status: ObservationStatus::Observed,
            callable_depth: Some(depth),
            alt_count: Some(alt),
            depth_eligible: Some(depth >= threshold),
            alt_supported: Some(depth >= threshold && alt > 0),
            observed_alt_fraction: (depth != 0).then_some(ExactFraction {
                numerator: alt,
                denominator: depth,
            }),
            min_callable_depth: threshold,
        }
    }
    fn pair(id: &str, left: &str, right: &str) -> PairSpec {
        PairSpec {
            id: id.into(),
            left: left.into(),
            right: right.into(),
        }
    }

    #[test]
    fn explicit_order_direction_and_repeated_member_use_are_preserved() {
        let plan = plan_pairs(
            &[
                pair("later-first", "A", "B"),
                pair("earlier-second", "B", "A"),
                pair("shared-control", "A", "C"),
            ],
            &["B".into(), "C".into(), "A".into()],
            PairLimits::default(),
        )
        .unwrap();
        assert_eq!(
            plan.pairs,
            vec![
                ResolvedPair {
                    id: "later-first".into(),
                    left_member_index: 2,
                    right_member_index: 0
                },
                ResolvedPair {
                    id: "earlier-second".into(),
                    left_member_index: 0,
                    right_member_index: 2
                },
                ResolvedPair {
                    id: "shared-control".into(),
                    left_member_index: 2,
                    right_member_index: 1
                },
            ]
        );
        assert!(plan.metadata_bytes > 4096);
        assert!(plan_pairs(&[], &[], PairLimits::default())
            .unwrap()
            .pairs
            .is_empty());
    }

    #[test]
    fn ambiguous_invalid_and_unbounded_pairs_refuse() {
        let members = ["A".into(), "B".into()];
        for pairs in [
            vec![pair("same", "A", "B"), pair("same", "B", "A")],
            vec![pair("self", "A", "A")],
            vec![pair("missing", "A", "C")],
            vec![pair("", "A", "B")],
            vec![pair("bad\tname", "A", "B")],
        ] {
            assert!(plan_pairs(&pairs, &members, PairLimits::default()).is_err());
        }
        assert!(plan_pairs(&[], &["A".into(), "A".into()], PairLimits::default()).is_err());
        assert!(plan_pairs(
            &[pair(&"x".repeat(257), "A", "B")],
            &members,
            PairLimits::default()
        )
        .is_err());
        let pairs = [pair("one", "A", "B"), pair("two", "B", "A")];
        assert!(plan_pairs(
            &pairs,
            &members,
            PairLimits {
                max_pairs: 1,
                ..PairLimits::default()
            }
        )
        .is_err());
        assert!(plan_pairs(
            &pairs,
            &members,
            PairLimits {
                max_members: 1,
                ..PairLimits::default()
            }
        )
        .is_err());
        assert!(plan_pairs(
            &pairs,
            &members,
            PairLimits {
                max_metadata_bytes: 4096,
                ..PairLimits::default()
            }
        )
        .is_err());
        assert!(plan_pairs(
            &[],
            &[],
            PairLimits {
                max_pairs: 0,
                ..PairLimits::default()
            }
        )
        .is_err());
    }

    #[test]
    fn independent_fraction_oracles_preserve_low_depth_and_direction() {
        // 3/4 - 1/2 = 2/8: low coverage still has a defined observed difference.
        let low = compare_observations(observed(2, 1, 10), observed(4, 3, 10)).unwrap();
        assert_eq!(
            low.right_minus_left,
            Some(ExactDifference {
                negative: false,
                magnitude: 2,
                denominator: 8
            })
        );
        assert_eq!(low.both_depth_eligible, Some(false));
        let reverse = compare_observations(low.right, low.left).unwrap();
        assert_eq!(
            reverse.right_minus_left,
            Some(ExactDifference {
                negative: true,
                magnitude: 2,
                denominator: 8
            })
        );
        // 5/10 - 4/10 = 10/100; both sides meet the default technical screen.
        let eligible = compare_observations(observed(10, 4, 10), observed(10, 5, 10)).unwrap();
        assert_eq!(
            eligible.right_minus_left,
            Some(ExactDifference {
                negative: false,
                magnitude: 10,
                denominator: 100
            })
        );
        assert_eq!(eligible.both_depth_eligible, Some(true));
        assert_eq!(eligible.left.alt_count, Some(4));
        assert_eq!(eligible.right.alt_count, Some(5));
    }

    #[test]
    fn missing_zero_and_equal_fractions_have_distinct_representations() {
        let missing = CandidateObservation::unmeasured(10).unwrap();
        let positive = observed(12, 4, 10);
        for (left, right) in [(missing, positive), (positive, missing), (missing, missing)] {
            let value = compare_observations(left, right).unwrap();
            assert_eq!(value.right_minus_left, None);
            assert_eq!(value.both_depth_eligible, None);
        }
        for (left, right) in [
            (observed(0, 0, 10), positive),
            (positive, observed(0, 0, 10)),
        ] {
            let value = compare_observations(left, right).unwrap();
            assert_eq!(value.right_minus_left, None);
            assert_eq!(value.both_depth_eligible, Some(false));
        }
        let equal = compare_observations(observed(2, 1, 10), observed(4, 2, 10)).unwrap();
        assert_eq!(
            equal.right_minus_left,
            Some(ExactDifference {
                negative: false,
                magnitude: 0,
                denominator: 8
            })
        );
        let zero_alt = compare_observations(observed(12, 0, 10), observed(10, 0, 10)).unwrap();
        assert_eq!(
            zero_alt.right_minus_left,
            Some(ExactDifference {
                negative: false,
                magnitude: 0,
                denominator: 120
            })
        );
    }

    #[test]
    fn uint64_max_cross_products_are_exact_without_signed_overflow() {
        let max = u64::MAX;
        // (M-2)/(M-1) - (M-1)/M = -1/[M(M-1)]. Both cross products
        // exceed i128::MAX; subtraction must preserve the unit difference.
        let left = observed(max, max - 1, 10);
        let right = observed(max - 1, max - 2, 10);
        let denominator = 340282366920938463408034375210639556610u128;
        let result = compare_observations(left, right).unwrap();
        assert_eq!(
            result.right_minus_left,
            Some(ExactDifference {
                negative: true,
                magnitude: 1,
                denominator
            })
        );
        assert_eq!(
            compare_observations(right, left).unwrap().right_minus_left,
            Some(ExactDifference {
                negative: false,
                magnitude: 1,
                denominator
            })
        );
        let one = compare_observations(observed(max, 0, 10), observed(max, max, 10)).unwrap();
        let square = 340282366920938463426481119284349108225u128;
        assert_eq!(
            one.right_minus_left,
            Some(ExactDifference {
                negative: false,
                magnitude: square,
                denominator: square
            })
        );
        assert_eq!(one.right.callable_depth, Some(max));
    }

    #[test]
    fn forged_observations_and_mixed_thresholds_refuse() {
        let valid = observed(12, 4, 10);
        let mut examples = Vec::new();
        let mut value = valid;
        value.callable_depth = None;
        examples.push(value);
        let mut value = valid;
        value.alt_count = Some(13);
        examples.push(value);
        let mut value = valid;
        value.depth_eligible = Some(false);
        examples.push(value);
        let mut value = valid;
        value.alt_supported = None;
        examples.push(value);
        let mut value = valid;
        value.observed_alt_fraction = Some(ExactFraction {
            numerator: 1,
            denominator: 3,
        });
        examples.push(value);
        let mut value = valid;
        value.status = ObservationStatus::Unmeasured;
        examples.push(value);
        let mut value = valid;
        value.min_callable_depth = 0;
        examples.push(value);
        let mut value = CandidateObservation::unmeasured(10).unwrap();
        value.callable_depth = Some(0);
        examples.push(value);
        for invalid in examples {
            assert!(compare_observations(invalid, valid).is_err());
            assert!(compare_observations(valid, invalid).is_err());
        }
        assert!(compare_observations(valid, observed(12, 4, 11)).is_err());
    }
}
