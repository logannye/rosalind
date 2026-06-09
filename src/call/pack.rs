//! Bin-packing over PREDICTED memory peaks — the contract's prediction turned
//! into a *placement* decision.
//!
//! Each job's peak is a conservative upper bound computed from the index header
//! alone (no run, no read I/O — milliseconds), and peaks are additive. So a
//! scheduler can SUM predicted peaks across co-located jobs and show a node stays
//! within budget before launching a byte. Emergent-peak callers (GATK/DeepVariant)
//! only learn their peak after a possible OOM-kill, so every co-location is a
//! gamble. Here, `Packed(...)` reports a placement whose every node's summed
//! predicted peak is `<=` its capacity, established up front.

/// A job to place: a label and its predicted peak RSS, in bytes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PackJob {
    /// Human-readable label (e.g. the index path or sample id).
    pub label: String,
    /// Predicted peak RSS (a conservative upper bound), in bytes.
    pub predicted_peak_bytes: u64,
}

/// One node's assignment: which jobs land on it and their summed predicted peak
/// (which is `<=` the node capacity — the fit, by predicted peak).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NodeAssignment {
    /// Node index (0-based).
    pub node: usize,
    /// Labels of the jobs placed on this node.
    pub job_labels: Vec<String>,
    /// Summed predicted peak of those jobs, in bytes.
    pub used_bytes: u64,
}

/// The outcome of a packing attempt.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PackOutcome {
    /// Every job placed. Each `NodeAssignment.used_bytes <= node_capacity_bytes`,
    /// so the node is within capacity by predicted peak before any job runs.
    Packed(Vec<NodeAssignment>),
    /// No safe packing exists under the constraints.
    NoFit {
        /// Why no packing was found (an actionable message).
        reason: String,
    },
}

/// Pack `jobs` onto nodes of `node_capacity_bytes` via first-fit-decreasing
/// (largest job first → tighter packing, fewer nodes). With `max_nodes = Some(m)`
/// the pack refuses when the jobs cannot fit within `m` nodes; with `None` it
/// opens as many nodes as needed and reports the count. A single job whose
/// predicted peak exceeds a whole node is an immediate `NoFit` — it can run
/// nowhere. Deterministic: ties break by label, so the same jobs always pack the
/// same way (the contract's reproducibility extends to the schedule).
pub fn first_fit_decreasing(
    jobs: &[PackJob],
    node_capacity_bytes: u64,
    max_nodes: Option<usize>,
) -> PackOutcome {
    // A job larger than a whole node can never be placed — surface it explicitly
    // rather than spinning up unbounded nodes.
    if let Some(j) = jobs
        .iter()
        .find(|j| j.predicted_peak_bytes > node_capacity_bytes)
    {
        return PackOutcome::NoFit {
            reason: format!(
                "job '{}' predicted peak {} B exceeds the node capacity {} B — it cannot run on any node; \
                 use a larger --node-mb or lower --max-depth",
                j.label, j.predicted_peak_bytes, node_capacity_bytes
            ),
        };
    }

    // FFD: descending peak, ties by label (deterministic).
    let mut order: Vec<&PackJob> = jobs.iter().collect();
    order.sort_by(|a, b| {
        b.predicted_peak_bytes
            .cmp(&a.predicted_peak_bytes)
            .then_with(|| a.label.cmp(&b.label))
    });

    let mut nodes: Vec<NodeAssignment> = Vec::new();
    for job in order {
        // Place in the first node that still has room (first-fit).
        let slot = nodes
            .iter_mut()
            .find(|n| n.used_bytes + job.predicted_peak_bytes <= node_capacity_bytes);
        match slot {
            Some(n) => {
                n.used_bytes += job.predicted_peak_bytes;
                n.job_labels.push(job.label.clone());
            }
            None => {
                if let Some(max) = max_nodes {
                    if nodes.len() >= max {
                        return PackOutcome::NoFit {
                            reason: format!(
                                "jobs do not fit in {max} node(s) of {} B each; need more nodes or a larger --node-mb",
                                node_capacity_bytes
                            ),
                        };
                    }
                }
                nodes.push(NodeAssignment {
                    node: nodes.len(),
                    job_labels: vec![job.label.clone()],
                    used_bytes: job.predicted_peak_bytes,
                });
            }
        }
    }
    PackOutcome::Packed(nodes)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn job(label: &str, peak: u64) -> PackJob {
        PackJob {
            label: label.to_string(),
            predicted_peak_bytes: peak,
        }
    }

    #[test]
    fn packs_all_jobs_and_every_node_is_within_capacity() {
        let jobs = vec![job("a", 30), job("b", 40), job("c", 50), job("d", 20)];
        let cap = 100;
        let PackOutcome::Packed(nodes) = first_fit_decreasing(&jobs, cap, None) else {
            panic!("should pack");
        };
        // THE invariant: every node's summed predicted peak is within capacity.
        for n in &nodes {
            assert!(
                n.used_bytes <= cap,
                "node {} sums {} > capacity {cap}",
                n.node,
                n.used_bytes
            );
        }
        // Every job placed exactly once.
        let mut placed: Vec<&String> = nodes.iter().flat_map(|n| &n.job_labels).collect();
        placed.sort();
        assert_eq!(placed, vec!["a", "b", "c", "d"]);
        // FFD on {50,40,30,20} into 100 → 2 nodes (50+40, 30+20).
        assert_eq!(nodes.len(), 2, "FFD should use 2 nodes: {nodes:?}");
    }

    #[test]
    fn refuses_a_job_larger_than_a_node() {
        let jobs = vec![job("ok", 50), job("toobig", 150)];
        match first_fit_decreasing(&jobs, 100, None) {
            PackOutcome::NoFit { reason } => assert!(
                reason.contains("toobig") && reason.contains("cannot run on any node"),
                "unexpected reason: {reason}"
            ),
            other => panic!("a job bigger than a node must NoFit: {other:?}"),
        }
    }

    #[test]
    fn refuses_when_jobs_exceed_the_node_budget() {
        // Three 60 B jobs into 100 B nodes need 3 nodes; cap at 2 → NoFit.
        let jobs = vec![job("a", 60), job("b", 60), job("c", 60)];
        match first_fit_decreasing(&jobs, 100, Some(2)) {
            PackOutcome::NoFit { reason } => assert!(reason.contains("2 node"), "reason: {reason}"),
            other => panic!("should not fit in 2 nodes: {other:?}"),
        }
        // With 3 nodes allowed, it fits.
        assert!(matches!(
            first_fit_decreasing(&jobs, 100, Some(3)),
            PackOutcome::Packed(_)
        ));
    }

    #[test]
    fn is_deterministic_regardless_of_input_order() {
        let a = vec![job("x", 40), job("y", 40), job("z", 40)];
        let b = vec![job("z", 40), job("x", 40), job("y", 40)];
        assert_eq!(
            first_fit_decreasing(&a, 100, None),
            first_fit_decreasing(&b, 100, None),
            "packing must be order-independent (deterministic schedule)"
        );
    }

    #[test]
    fn empty_jobs_pack_to_zero_nodes() {
        assert_eq!(
            first_fit_decreasing(&[], 100, None),
            PackOutcome::Packed(vec![])
        );
    }
}
