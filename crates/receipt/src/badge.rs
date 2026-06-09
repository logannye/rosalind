//! Self-hosted status badges — a shields.io *endpoint*-format JSON and a static SVG —
//! asserting "reproducible · fits N MiB" for a run. Self-hosted: no shields.io runtime
//! dependency (we emit the SVG directly), so the badge works offline / air-gapped.

use super::json_escape;

/// The badge message + color for a (reproducible, fits-budget?) status.
fn badge_message(reproducible: bool, fits_mb: Option<u64>) -> (String, &'static str) {
    if !reproducible {
        return ("not reproducible".to_string(), "red");
    }
    match fits_mb {
        Some(mb) => (format!("reproducible, fits {mb} MiB"), "brightgreen"),
        None => ("reproducible".to_string(), "brightgreen"),
    }
}

/// A shields.io endpoint-format JSON object (consume via shields.io's `endpoint` badge,
/// or render the SVG below directly).
pub fn badge_json(reproducible: bool, fits_mb: Option<u64>) -> String {
    let (msg, color) = badge_message(reproducible, fits_mb);
    format!(
        "{{\"schemaVersion\":1,\"label\":\"rosalind\",\"message\":\"{}\",\"color\":\"{}\"}}",
        json_escape(&msg),
        color
    )
}

/// A minimal, static, self-contained SVG badge (no external fetch).
pub fn badge_svg(reproducible: bool, fits_mb: Option<u64>) -> String {
    let (msg, color) = badge_message(reproducible, fits_mb);
    let fill = match color {
        "brightgreen" => "#4c1",
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn badge_json_is_shields_endpoint_shaped() {
        let j = badge_json(true, Some(256));
        assert!(j.contains("\"schemaVersion\":1"), "{j}");
        assert!(j.contains("\"label\":\"rosalind\""), "{j}");
        assert!(j.contains("reproducible"), "{j}");
        assert!(j.contains("256 MiB"), "{j}");
        assert!(j.contains("\"color\":\"brightgreen\""), "{j}");
    }

    #[test]
    fn badge_json_red_when_not_reproducible() {
        let j = badge_json(false, None);
        assert!(j.contains("\"color\":\"red\""), "{j}");
        assert!(j.contains("not reproducible"), "{j}");
    }

    #[test]
    fn badge_svg_is_well_formed() {
        let s = badge_svg(true, Some(256));
        assert!(s.starts_with("<svg"), "{s}");
        assert!(s.trim_end().ends_with("</svg>"), "{s}");
        assert!(s.contains("reproducible"), "{s}");
    }
}
