//! Pure, `std`-only provenance-DAG walk over a set of run receipts. No filesystem, no
//! htslib — so it is unit-testable and wasm-portable (the future `chain confirm` board /
//! browser verifier reuse it). A node is a receipt (id = its claim `content_hash`); an
//! edge is one input operand, classified by whether its content hash is produced by
//! another node in the set.

use std::collections::HashMap;

use crate::{FileHash, RunManifest};

/// Input operand flags that MUST resolve to a producing node. An unresolved one is a
/// broken chain (a missing/mismatched upstream receipt). Everything else (reads,
/// alignments, reference FASTA) is an external source: integrity-verified by its
/// recorded hash, never a chain failure.
const EXPECTED_INTERNAL: &[&str] = &["--index"];

/// How one input operand resolved against the receipt set.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum EdgeStatus {
    /// The input hash equals some node's output hash — an internal provenance edge.
    Resolved { parent_id: String },
    /// Unresolved, but an expected-external operand (reference/alignments/reads).
    External,
    /// Unresolved AND an expected-internal operand (`--index`) — breaks the chain.
    Broken,
}

/// One classified input operand of one node.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ChainEdge {
    pub child_id: String,
    pub child_subcommand: String,
    pub flag: String,
    pub input_blake3: String,
    pub status: EdgeStatus,
}

/// A receipt as a DAG node.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ChainNode {
    pub id: String,
    pub subcommand: String,
    /// `Some(true/false)` if a claim self-hash is recorded and matches/mismatches;
    /// `None` for a pre-self-hash receipt.
    pub self_hash: Option<bool>,
}

/// The result of walking the set: nodes, classified edges, and the overall verdict.
#[derive(Debug, Clone)]
pub struct ChainReport {
    pub nodes: Vec<ChainNode>,
    pub edges: Vec<ChainEdge>,
    /// `true` iff no node fails its self-hash AND no edge is `Broken`.
    pub intact: bool,
}

impl ChainReport {
    /// A compact, dependency-free JSON summary for `--json` / a scheduler.
    pub fn to_json(&self) -> String {
        let resolved = self
            .edges
            .iter()
            .filter(|e| matches!(e.status, EdgeStatus::Resolved { .. }))
            .count();
        let external = self
            .edges
            .iter()
            .filter(|e| matches!(e.status, EdgeStatus::External))
            .count();
        let broken = self
            .edges
            .iter()
            .filter(|e| matches!(e.status, EdgeStatus::Broken))
            .count();
        format!(
            "{{\"intact\":{},\"nodes\":{},\"edges_resolved\":{},\"edges_external\":{},\"edges_broken\":{}}}",
            self.intact,
            self.nodes.len(),
            resolved,
            external,
            broken
        )
    }
}

/// Recover `(flag, input_blake3)` pairs from a recorded `command` string: each
/// `@in:<hash>` token is preceded by its operand flag.
fn input_operands(command: &str) -> Vec<(String, String)> {
    let toks: Vec<&str> = command.split(' ').collect();
    let mut out = Vec::new();
    for (i, t) in toks.iter().enumerate() {
        if let Some(h) = t.strip_prefix("@in:") {
            let flag = if i > 0 { toks[i - 1] } else { "?" };
            out.push((flag.to_string(), h.to_string()));
        }
    }
    out
}

/// Operands for a node: from its recorded `command` when present (carries the flags),
/// else a flag-less fallback over `inputs[]` (every input treated as external).
fn operands_for(m: &RunManifest) -> Vec<(String, String)> {
    match m.params.get("command") {
        Some(c) if !c.is_empty() => input_operands(c),
        _ => m
            .inputs
            .iter()
            .map(|f: &FileHash| ("?".to_string(), f.blake3.clone()))
            .collect(),
    }
}

/// Walk the receipt set as a provenance DAG. Pure: no I/O.
pub fn walk_chain(receipts: &[RunManifest]) -> ChainReport {
    let ids: Vec<String> = receipts.iter().map(|m| m.content_hash()).collect();

    // output content hash -> producing node id (first producer wins).
    let mut producer: HashMap<String, String> = HashMap::new();
    for (m, id) in receipts.iter().zip(&ids) {
        for o in &m.outputs {
            producer.entry(o.blake3.clone()).or_insert_with(|| id.clone());
        }
    }

    let nodes: Vec<ChainNode> = receipts
        .iter()
        .zip(&ids)
        .map(|(m, id)| ChainNode {
            id: id.clone(),
            subcommand: m.subcommand.clone(),
            self_hash: m.self_hash_ok(),
        })
        .collect();

    let mut edges = Vec::new();
    for (m, id) in receipts.iter().zip(&ids) {
        for (flag, hash) in operands_for(m) {
            let status = if let Some(parent) = producer.get(&hash) {
                EdgeStatus::Resolved {
                    parent_id: parent.clone(),
                }
            } else if EXPECTED_INTERNAL.contains(&flag.as_str()) {
                EdgeStatus::Broken
            } else {
                EdgeStatus::External
            };
            edges.push(ChainEdge {
                child_id: id.clone(),
                child_subcommand: m.subcommand.clone(),
                flag,
                input_blake3: hash,
                status,
            });
        }
    }

    let tampered = nodes.iter().any(|n| n.self_hash == Some(false));
    let broken = edges.iter().any(|e| matches!(e.status, EdgeStatus::Broken));
    ChainReport {
        nodes,
        edges,
        intact: !tampered && !broken,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{FileHash, RunManifest};

    /// Build a finalized manifest with a recorded `command`, inputs, and outputs.
    fn mk(
        sub: &str,
        command: &str,
        inputs: &[(&str, &str)],
        outputs: &[(&str, &str)],
    ) -> RunManifest {
        let mut m = RunManifest::new(sub);
        m.params.insert("command".to_string(), command.to_string());
        m.inputs = inputs
            .iter()
            .map(|(p, h)| FileHash {
                path: p.to_string(),
                blake3: h.to_string(),
            })
            .collect();
        m.outputs = outputs
            .iter()
            .map(|(p, h)| FileHash {
                path: p.to_string(),
                blake3: h.to_string(),
            })
            .collect();
        m.finalize();
        m
    }

    fn index_node() -> RunManifest {
        mk(
            "index",
            "index --reference @in:rh --output @out:ih",
            &[("ref.fa", "rh")],
            &[("ref.idx", "ih")],
        )
    }

    fn variants_node(index_hash: &str) -> RunManifest {
        mk(
            "variants",
            &format!("variants --index @in:{index_hash} --alignments @in:bh -o @out:vh"),
            &[("ref.idx", index_hash), ("s.bam", "bh")],
            &[("calls.vcf", "vh")],
        )
    }

    #[test]
    fn resolves_the_index_edge_and_marks_external_inputs() {
        let report = walk_chain(&[index_node(), variants_node("ih")]);
        assert!(report.intact, "a complete chain is intact");
        // The --index edge resolves; --alignments (bh) and --reference (rh) are external.
        let resolved = report
            .edges
            .iter()
            .filter(|e| matches!(e.status, EdgeStatus::Resolved { .. }))
            .count();
        let external = report
            .edges
            .iter()
            .filter(|e| matches!(e.status, EdgeStatus::External))
            .count();
        assert_eq!(resolved, 1, "exactly the --index edge resolves");
        assert_eq!(external, 2, "--alignments and --reference are external sources");
        assert_eq!(
            report.to_json(),
            "{\"intact\":true,\"nodes\":2,\"edges_resolved\":1,\"edges_external\":2,\"edges_broken\":0}"
        );
    }

    #[test]
    fn an_unresolved_index_input_is_broken() {
        // variants references an --index hash no node produces.
        let report = walk_chain(&[index_node(), variants_node("MISSING")]);
        assert!(!report.intact, "a dangling --index edge breaks the chain");
        assert!(report
            .edges
            .iter()
            .any(|e| e.flag == "--index" && e.status == EdgeStatus::Broken));
    }

    #[test]
    fn an_unresolved_external_input_does_not_break_the_chain() {
        // The index node alone: its --reference (rh) resolves to nothing, but --reference
        // is an external source, so the chain is still intact.
        let report = walk_chain(&[index_node()]);
        assert!(
            report.intact,
            "an unresolved external source is integrity-only, not broken"
        );
        assert!(report
            .edges
            .iter()
            .any(|e| e.flag == "--reference" && e.status == EdgeStatus::External));
    }

    #[test]
    fn a_tampered_node_breaks_the_chain() {
        let mut tampered = index_node();
        // Edit a claim field AFTER finalize → the recorded manifest_blake3 no longer matches.
        tampered
            .params
            .insert("total_bp".to_string(), "999999".to_string());
        let report = walk_chain(&[tampered, variants_node("ih")]);
        assert!(!report.intact, "a self-hash mismatch breaks the chain");
        assert!(report.nodes.iter().any(|n| n.self_hash == Some(false)));
    }
}
