//! Exact first-party candidate screens. These are read observations, never
//! genotype calls, population frequencies or calibrated confidence estimates.

use super::{CohortError, Result};
use crate::evidence::{EvidenceFields, EvidenceRowRef};
use serde::{Deserialize, Serialize};

pub(crate) const DEFAULT_MIN_CALLABLE_DEPTH: u64 = 10;

pub(crate) fn required_fields() -> EvidenceFields {
    EvidenceFields::DEPTHS.union(EvidenceFields::ALLELES)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub(crate) enum ObservationStatus {
    Observed,
    Unmeasured,
}

/// Exact numerator/denominator, with no loss of uint64 precision from conversion
/// to floating point. A zero denominator is represented by None, not 0/0 or 0.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub(crate) struct ExactFraction {
    pub numerator: u64,
    pub denominator: u64,
}

fn fraction(numerator: u64, denominator: u64) -> Option<ExactFraction> {
    (denominator != 0).then_some(ExactFraction {
        numerator,
        denominator,
    })
}

/// One member/locus/ALT observation. The constructors keep absent loci distinct
/// from observed zero depth and reject absent fields before counting anything.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub(crate) struct CandidateObservation {
    pub status: ObservationStatus,
    pub callable_depth: Option<u64>,
    pub alt_count: Option<u64>,
    pub depth_eligible: Option<bool>,
    pub alt_supported: Option<bool>,
    pub observed_alt_fraction: Option<ExactFraction>,
    pub min_callable_depth: u64,
}

impl CandidateObservation {
    pub fn unmeasured(min_callable_depth: u64) -> Result<Self> {
        threshold(min_callable_depth)?;
        Ok(Self {
            status: ObservationStatus::Unmeasured,
            callable_depth: None,
            alt_count: None,
            depth_eligible: None,
            alt_supported: None,
            observed_alt_fraction: None,
            min_callable_depth,
        })
    }

    pub fn from_row(row: EvidenceRowRef<'_>, alt: u8, min_callable_depth: u64) -> Result<Self> {
        let depths = row.depths.ok_or_else(|| {
            CohortError::Incompatible("candidate summary requires stored depths".into())
        })?;
        let alleles = row.alleles.ok_or_else(|| {
            CohortError::Incompatible("candidate summary requires stored allele counts".into())
        })?;
        let index = match alt {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => {
                return Err(CohortError::Incompatible(
                    "candidate ALT must be A/C/G/T".into(),
                ))
            }
        };
        if !row.requested_alts.contains(&alt) {
            return Err(CohortError::Incompatible(
                "ALT is absent from normalized candidate query".into(),
            ));
        }
        let total = alleles
            .allele_counts
            .iter()
            .try_fold(0u64, |total, &count| {
                checked_add(total, count, "callable allele observations")
            })?;
        if total != depths.callable_depth {
            return Err(CohortError::Corrupt(
                "callable depth differs from A/C/G/T counts".into(),
            ));
        }
        Self::observed(
            depths.callable_depth,
            alleles.allele_counts[index],
            min_callable_depth,
        )
    }

    fn observed(callable_depth: u64, alt_count: u64, min_callable_depth: u64) -> Result<Self> {
        threshold(min_callable_depth)?;
        if alt_count > callable_depth {
            return Err(CohortError::Corrupt(
                "ALT observations exceed callable depth".into(),
            ));
        }
        let eligible = callable_depth >= min_callable_depth;
        Ok(Self {
            status: ObservationStatus::Observed,
            callable_depth: Some(callable_depth),
            alt_count: Some(alt_count),
            depth_eligible: Some(eligible),
            alt_supported: Some(eligible && alt_count != 0),
            observed_alt_fraction: fraction(alt_count, callable_depth),
            min_callable_depth,
        })
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub(crate) struct CandidateSummary {
    pub min_callable_depth: u64,
    pub n_requested: u64,
    pub n_observed: u64,
    pub n_depth_eligible: u64,
    pub n_alt_supported: u64,
    pub callable_total: u64,
    pub alt_total: u64,
    pub eligible_callable_total: u64,
    pub eligible_alt_total: u64,
    /// ALT-supported samples / technically depth-eligible samples.
    pub depth_eligible_support_fraction: Option<ExactFraction>,
    /// Read observations pooled across observed samples, with exact denominator.
    pub observed_alt_fraction: Option<ExactFraction>,
    pub eligible_alt_fraction: Option<ExactFraction>,
}

/// Constant state for one candidate; a runtime owns at most one bounded genomic
/// window of these states and supplies exactly one cell per requested member.
#[derive(Debug, Clone, Copy)]
pub(crate) struct CandidateReducer {
    summary: CandidateSummary,
    received: u64,
}

impl CandidateReducer {
    pub fn new(n_requested: u64, min_callable_depth: u64) -> Result<Self> {
        threshold(min_callable_depth)?;
        Ok(Self {
            summary: CandidateSummary {
                min_callable_depth,
                n_requested,
                n_observed: 0,
                n_depth_eligible: 0,
                n_alt_supported: 0,
                callable_total: 0,
                alt_total: 0,
                eligible_callable_total: 0,
                eligible_alt_total: 0,
                depth_eligible_support_fraction: None,
                observed_alt_fraction: None,
                eligible_alt_fraction: None,
            },
            received: 0,
        })
    }

    /// A rejected cell never partially updates the state. The runtime controls
    /// canonical member order and disjoint leaf ownership; no member matrix is retained.
    pub fn push(&mut self, observation: &CandidateObservation) -> Result<()> {
        if self.received >= self.summary.n_requested {
            return Err(CohortError::Corrupt(
                "more member cells than requested".into(),
            ));
        }
        if observation.min_callable_depth != self.summary.min_callable_depth {
            return Err(CohortError::Incompatible(
                "candidate screen thresholds differ".into(),
            ));
        }
        let mut next = *self;
        next.received = checked_add(next.received, 1, "member cells")?;
        match observation.status {
            ObservationStatus::Unmeasured => {
                if *observation
                    != CandidateObservation::unmeasured(self.summary.min_callable_depth)?
                {
                    return Err(CohortError::Corrupt(
                        "unmeasured cell contains evidence".into(),
                    ));
                }
            }
            ObservationStatus::Observed => {
                let (Some(depth), Some(alt)) = (observation.callable_depth, observation.alt_count)
                else {
                    return Err(CohortError::Corrupt(
                        "observed cell lacks required evidence".into(),
                    ));
                };
                if *observation
                    != CandidateObservation::observed(depth, alt, self.summary.min_callable_depth)?
                {
                    return Err(CohortError::Corrupt(
                        "observed cell has inconsistent eligibility or fraction".into(),
                    ));
                }
                let summary = &mut next.summary;
                summary.n_observed = checked_add(summary.n_observed, 1, "observed samples")?;
                summary.callable_total =
                    checked_add(summary.callable_total, depth, "callable observations")?;
                summary.alt_total = checked_add(summary.alt_total, alt, "ALT observations")?;
                if depth >= summary.min_callable_depth {
                    summary.n_depth_eligible =
                        checked_add(summary.n_depth_eligible, 1, "eligible samples")?;
                    summary.eligible_callable_total = checked_add(
                        summary.eligible_callable_total,
                        depth,
                        "eligible callable observations",
                    )?;
                    summary.eligible_alt_total =
                        checked_add(summary.eligible_alt_total, alt, "eligible ALT observations")?;
                    if alt != 0 {
                        summary.n_alt_supported =
                            checked_add(summary.n_alt_supported, 1, "ALT-supported samples")?;
                    }
                }
            }
        }
        *self = next;
        Ok(())
    }

    pub fn finish(self) -> Result<CandidateSummary> {
        if self.received != self.summary.n_requested {
            return Err(CohortError::Corrupt(format!(
                "received {} of {} requested member cells",
                self.received, self.summary.n_requested
            )));
        }
        let mut summary = self.summary;
        summary.depth_eligible_support_fraction =
            fraction(summary.n_alt_supported, summary.n_depth_eligible);
        summary.observed_alt_fraction = fraction(summary.alt_total, summary.callable_total);
        summary.eligible_alt_fraction =
            fraction(summary.eligible_alt_total, summary.eligible_callable_total);
        Ok(summary)
    }
}

fn threshold(value: u64) -> Result<()> {
    if value == 0 {
        Err(CohortError::Incompatible(
            "min_callable_depth must be positive; it is a technical screen, not confidence".into(),
        ))
    } else {
        Ok(())
    }
}

fn checked_add(left: u64, right: u64, what: &str) -> Result<u64> {
    left.checked_add(right)
        .ok_or_else(|| CohortError::Limit(format!("uint64 overflow in {what}")))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::evidence::{EvidenceAlleles, EvidenceDepths};
    use serde_json::Value;

    #[test]
    fn absent_zero_low_and_eligible_cells_remain_distinct() {
        let missing = CandidateObservation::unmeasured(10).unwrap();
        let zero = CandidateObservation::observed(0, 0, 10).unwrap();
        let low = CandidateObservation::observed(9, 1, 10).unwrap();
        let eligible = CandidateObservation::observed(10, 1, 10).unwrap();
        assert_eq!(missing.status, ObservationStatus::Unmeasured);
        assert_eq!(missing.callable_depth, None);
        assert_eq!(missing.depth_eligible, None);
        assert_eq!(zero.status, ObservationStatus::Observed);
        assert_eq!(zero.callable_depth, Some(0));
        assert_eq!(zero.depth_eligible, Some(false));
        assert_eq!(zero.observed_alt_fraction, None);
        assert_eq!(low.alt_count, Some(1));
        assert_eq!(low.alt_supported, Some(false));
        assert_eq!(low.observed_alt_fraction, fraction(1, 9));
        assert_eq!(eligible.alt_supported, Some(true));
        let mut reducer = CandidateReducer::new(4, 10).unwrap();
        for cell in [missing, zero, low, eligible] {
            reducer.push(&cell).unwrap();
        }
        let result = reducer.finish().unwrap();
        assert_eq!(
            (
                result.n_requested,
                result.n_observed,
                result.n_depth_eligible,
                result.n_alt_supported
            ),
            (4, 3, 1, 1)
        );
        assert_eq!(result.observed_alt_fraction, fraction(2, 19));
        assert_eq!(result.depth_eligible_support_fraction, fraction(1, 1));
        assert_eq!(result.eligible_alt_fraction, fraction(1, 10));
        assert!(serde_json::to_value(missing).unwrap()["callable_depth"].is_null());
        assert_eq!(serde_json::to_value(zero).unwrap()["callable_depth"], 0);
    }

    #[test]
    fn reducer_matches_hand_authored_multi_sample_oracle_before_and_after_fill() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../examples/cohort-reanalysis/expected.json"
        ))
        .unwrap();
        for (complete, key) in [(false, "partial_summary"), (true, "after_fill_summary")] {
            for row in expected[key].as_array().unwrap() {
                let position = row[0].as_u64().unwrap();
                let alt = row[1].as_str().unwrap().as_bytes()[0];
                let index = match alt {
                    b'A' => 4,
                    b'C' => 5,
                    b'G' => 6,
                    b'T' => 7,
                    _ => panic!(),
                };
                let members = expected["members"].as_object().unwrap();
                let mut reducer = CandidateReducer::new(members.len() as u64, 10).unwrap();
                for member in members.values() {
                    let represented = complete
                        || member["stored_positions"]
                            .as_array()
                            .unwrap()
                            .iter()
                            .any(|p| p.as_u64() == Some(position));
                    let cell = if represented {
                        let evidence = member["full_rows"]
                            .as_array()
                            .unwrap()
                            .iter()
                            .find(|r| r[0].as_u64() == Some(position))
                            .unwrap();
                        CandidateObservation::observed(
                            evidence[3].as_u64().unwrap(),
                            evidence[index].as_u64().unwrap(),
                            10,
                        )
                        .unwrap()
                    } else {
                        CandidateObservation::unmeasured(10).unwrap()
                    };
                    reducer.push(&cell).unwrap();
                }
                let actual = reducer.finish().unwrap();
                let counts = [
                    actual.n_requested,
                    actual.n_observed,
                    actual.n_depth_eligible,
                    actual.n_alt_supported,
                    actual.callable_total,
                    actual.alt_total,
                    actual.eligible_callable_total,
                    actual.eligible_alt_total,
                ];
                assert_eq!(
                    counts.to_vec(),
                    row.as_array().unwrap()[2..]
                        .iter()
                        .map(|n| n.as_u64().unwrap())
                        .collect::<Vec<_>>(),
                    "{key} {position}"
                );
            }
        }
    }

    #[test]
    fn uint64_counts_are_exact_and_overflow_is_transactional() {
        let mut reducer = CandidateReducer::new(2, 10).unwrap();
        reducer
            .push(&CandidateObservation::observed(u64::MAX, u64::MAX, 10).unwrap())
            .unwrap();
        assert!(reducer
            .push(&CandidateObservation::observed(1, 1, 10).unwrap())
            .is_err());
        reducer
            .push(&CandidateObservation::observed(0, 0, 10).unwrap())
            .unwrap();
        let result = reducer.finish().unwrap();
        assert_eq!(result.callable_total, u64::MAX);
        assert_eq!(result.alt_total, u64::MAX);
        assert_eq!(result.observed_alt_fraction, fraction(u64::MAX, u64::MAX));
        assert_eq!(
            serde_json::to_value(result).unwrap()["alt_total"].as_u64(),
            Some(u64::MAX)
        );
    }

    #[test]
    fn zero_denominators_missing_members_and_invalid_screens_are_explicit() {
        let empty = CandidateReducer::new(0, 10).unwrap().finish().unwrap();
        assert_eq!(empty.n_requested, 0);
        assert_eq!(empty.depth_eligible_support_fraction, None);
        assert_eq!(empty.observed_alt_fraction, None);
        assert!(CandidateReducer::new(1, 10).unwrap().finish().is_err());
        assert!(CandidateReducer::new(1, 0).is_err());
        assert!(CandidateObservation::unmeasured(0).is_err());
        let mut reducer = CandidateReducer::new(1, 10).unwrap();
        assert!(reducer
            .push(&CandidateObservation::observed(9, 1, 9).unwrap())
            .is_err());
        reducer
            .push(&CandidateObservation::unmeasured(10).unwrap())
            .unwrap();
        assert!(reducer
            .push(&CandidateObservation::unmeasured(10).unwrap())
            .is_err());
        assert_eq!(reducer.finish().unwrap().n_observed, 0);
    }

    #[test]
    fn multiallelic_cells_repeat_depth_without_summing_denominators() {
        let depths = EvidenceDepths {
            callable_depth: 12,
            ..EvidenceDepths::default()
        };
        let alleles = EvidenceAlleles {
            allele_counts: [7, 2, 3, 0],
        };
        let row = EvidenceRowRef {
            position: 9,
            reference: b'A',
            requested_alts: b"CG",
            depths: Some(&depths),
            alleles: Some(&alleles),
            strands: None,
            quality_sums: None,
            quality_histograms: None,
            read_position: None,
            allele_quality: None,
        };
        let c = CandidateObservation::from_row(row, b'C', 10).unwrap();
        let g = CandidateObservation::from_row(row, b'G', 10).unwrap();
        assert_eq!(c.observed_alt_fraction, fraction(2, 12));
        assert_eq!(g.observed_alt_fraction, fraction(3, 12));
        assert!(CandidateObservation::from_row(row, b'T', 10).is_err());
        assert!(CandidateObservation::from_row(
            EvidenceRowRef {
                alleles: None,
                ..row
            },
            b'C',
            10
        )
        .is_err());
        let bad_depths = EvidenceDepths {
            callable_depth: 11,
            ..depths
        };
        assert!(CandidateObservation::from_row(
            EvidenceRowRef {
                depths: Some(&bad_depths),
                ..row
            },
            b'C',
            10
        )
        .is_err());
    }
}
