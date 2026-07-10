//! Self-hosted status badges — a shields.io *endpoint*-format JSON and a static SVG —
//! reporting receipt integrity separately from actual reproduction evidence. Self-hosted:
//! no shields.io runtime dependency, so the badge works offline / air-gapped.

use super::{json_escape, TrustReport, TrustState};

/// Evidence level rendered by a badge.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BadgeStatus {
    /// The run receipt is intact, but no linked reproduction certificate was supplied.
    ReceiptIntact,
    /// A valid linked certificate records byte-identical output reproduction.
    Reproduced,
    /// The receipt predates checkable integrity or otherwise lacks necessary evidence.
    EvidenceMissing,
    /// Supplied receipt or reproduction evidence is invalid.
    Invalid,
}

impl BadgeStatus {
    /// Derive badge semantics without equating deterministic claims with reproduction.
    pub fn from_trust(report: &TrustReport) -> Self {
        if report.receipt_integrity.state == TrustState::Failed
            || report.reproduction_evidence.state == TrustState::Failed
        {
            Self::Invalid
        } else if report.reproduction_evidence.state == TrustState::Satisfied {
            Self::Reproduced
        } else if report.receipt_integrity.state == TrustState::Satisfied {
            Self::ReceiptIntact
        } else {
            Self::EvidenceMissing
        }
    }
}

/// The badge message + color for evidence and optional budget fit.
fn badge_message(status: BadgeStatus, fits_mb: Option<u64>) -> (String, &'static str) {
    match (status, fits_mb) {
        (BadgeStatus::ReceiptIntact, Some(mb)) => {
            (format!("receipt intact · fits {mb} MiB"), "blue")
        }
        (BadgeStatus::ReceiptIntact, None) => ("receipt intact".to_string(), "blue"),
        (BadgeStatus::Reproduced, Some(mb)) => {
            (format!("reproduced · fits {mb} MiB"), "brightgreen")
        }
        (BadgeStatus::Reproduced, None) => ("reproduced".to_string(), "brightgreen"),
        (BadgeStatus::EvidenceMissing, _) => ("evidence not supplied".to_string(), "lightgrey"),
        (BadgeStatus::Invalid, _) => ("invalid evidence".to_string(), "red"),
    }
}

/// A shields.io endpoint-format JSON object using explicit evidence semantics.
pub fn badge_json_for(status: BadgeStatus, fits_mb: Option<u64>) -> String {
    let (msg, color) = badge_message(status, fits_mb);
    format!(
        "{{\"schemaVersion\":1,\"label\":\"rosalind\",\"message\":\"{}\",\"color\":\"{}\"}}",
        json_escape(&msg),
        color
    )
}

/// A shields.io endpoint-format JSON object (consume via shields.io's `endpoint` badge,
/// or render the SVG below directly). Kept for source compatibility; new callers should
/// use [`badge_json_for`] so `true` is not confused with mere receipt integrity.
pub fn badge_json(reproducible: bool, fits_mb: Option<u64>) -> String {
    badge_json_for(
        if reproducible {
            BadgeStatus::Reproduced
        } else {
            BadgeStatus::Invalid
        },
        fits_mb,
    )
}

/// A minimal, static, self-contained SVG badge (no external fetch).
pub fn badge_svg_for(status: BadgeStatus, fits_mb: Option<u64>) -> String {
    let (msg, color) = badge_message(status, fits_mb);
    let fill = match color {
        "brightgreen" => "#4c1",
        "blue" => "#007ec6",
        "red" => "#e05d44",
        _ => "#9f9f9f",
    };
    // Width scales loosely with the message length so the text fits.
    let msg_w = 12 + msg.len() as u32 * 7;
    let total = 70 + msg_w;
    let msg_mid = 70 + msg_w / 2;
    format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{total}\" height=\"20\" \
role=\"img\" aria-label=\"rosalind: {msg}\">\
<rect width=\"70\" height=\"20\" fill=\"#555\"/>\
<rect x=\"70\" width=\"{msg_w}\" height=\"20\" fill=\"{fill}\"/>\
<g fill=\"#fff\" font-family=\"Verdana,Geneva,sans-serif\" font-size=\"11\">\
<text x=\"8\" y=\"14\">rosalind</text>\
<text x=\"{msg_mid}\" y=\"14\" text-anchor=\"middle\">{msg}</text>\
</g></svg>"
    )
}

/// Source-compatible wrapper for callers that already possess reproduction evidence.
pub fn badge_svg(reproducible: bool, fits_mb: Option<u64>) -> String {
    badge_svg_for(
        if reproducible {
            BadgeStatus::Reproduced
        } else {
            BadgeStatus::Invalid
        },
        fits_mb,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn badge_json_is_shields_endpoint_shaped() {
        let j = badge_json(true, Some(256));
        assert!(j.contains("\"schemaVersion\":1"), "{j}");
        assert!(j.contains("\"label\":\"rosalind\""), "{j}");
        assert!(j.contains("reproduced"), "{j}");
        assert!(j.contains("256 MiB"), "{j}");
        assert!(j.contains("\"color\":\"brightgreen\""), "{j}");
    }

    #[test]
    fn badge_json_red_when_not_reproducible() {
        let j = badge_json(false, None);
        assert!(j.contains("\"color\":\"red\""), "{j}");
        assert!(j.contains("invalid evidence"), "{j}");
    }

    #[test]
    fn intact_without_certificate_is_blue_not_reproduced() {
        let j = badge_json_for(BadgeStatus::ReceiptIntact, Some(128));
        assert!(j.contains("receipt intact"), "{j}");
        assert!(!j.contains("reproduced"), "{j}");
        assert!(j.contains("\"color\":\"blue\""), "{j}");
    }

    #[test]
    fn badge_svg_is_well_formed() {
        let s = badge_svg(true, Some(256));
        assert!(s.starts_with("<svg"), "{s}");
        assert!(s.trim_end().ends_with("</svg>"), "{s}");
        assert!(s.contains("reproduced"), "{s}");
    }
}
