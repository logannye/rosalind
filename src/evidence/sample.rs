//! Sample assignment is scientific selection, before any locus depth/filter count.

use super::EvidenceError;
use rust_htslib::bam::{record::Aux, HeaderView, Record};
use std::collections::{BTreeMap, BTreeSet};

/// Requested sample scope. The resolved scope, rather than Auto versus Named,
/// defines scientific identity; callers may separately record the original argv.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub enum EvidenceSampleSelection {
    /// Infer one unambiguous named sample, or explicitly unknown identity.
    #[default]
    Auto,
    /// Include only read groups belonging to this declared SM identifier.
    Named(String),
    /// Deliberately aggregate all records, including unassigned read groups.
    Pool,
}

/// Scientific interpretation of all evidence rows in a run.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EvidenceSampleMode {
    /// No sample names are declared; records have unknown sample identity.
    Unknown,
    /// Every counted record is assigned to the selected named sample.
    Named,
    /// All eligible records are deliberately pooled, regardless of assignment.
    Pooled,
}

/// One header-declared read group. Groups are sorted by ID in a resolved scope.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EvidenceReadGroup {
    /// Unique alignment read-group ID.
    pub id: String,
    /// Declared SM identifier, or None when the header does not assign a sample.
    pub sample: Option<String>,
}

/// Canonical, content-serializable sample scope resolved from the alignment header.
/// Auto and Named produce identical scopes when they resolve to the same sample.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EvidenceSampleScope {
    /// Resolved unknown, named, or explicitly pooled interpretation.
    pub mode: EvidenceSampleMode,
    /// Selected SM for named mode; None for unknown and pooled modes.
    pub selected_sample: Option<String>,
    /// All distinct declared SM identifiers, sorted lexically.
    pub declared_samples: Vec<String>,
    /// All declared read groups, sorted by ID, including unassigned groups.
    pub read_groups: Vec<EvidenceReadGroup>,
    /// Whether missing, undeclared, or unnamed read groups can contribute.
    /// True records a possibility, not a claim that such records were observed.
    pub allows_unassigned: bool,
}

impl EvidenceSampleScope {
    /// Stable JSON for recipes, receipts, and cache identities. This versioned
    /// representation includes resolved selection and all header assignments.
    pub fn canonical_json(&self) -> String {
        let mode = match self.mode {
            EvidenceSampleMode::Unknown => "unknown",
            EvidenceSampleMode::Named => "named",
            EvidenceSampleMode::Pooled => "pooled",
        };
        let mut output = format!(
            "{{\"version\":1,\"mode\":\"{mode}\",\"selected_sample\":{},\"declared_samples\":[",
            optional_json(self.selected_sample.as_deref())
        );
        for (i, sample) in self.declared_samples.iter().enumerate() {
            if i != 0 {
                output.push(',');
            }
            output.push_str(&string_json(sample));
        }
        output.push_str("],\"read_groups\":[");
        for (i, group) in self.read_groups.iter().enumerate() {
            if i != 0 {
                output.push(',');
            }
            output.push_str(&format!(
                "{{\"id\":{},\"sample\":{}}}",
                string_json(&group.id),
                optional_json(group.sample.as_deref())
            ));
        }
        output.push_str(&format!(
            "],\"allows_unassigned\":{}}}",
            self.allows_unassigned
        ));
        output
    }

    pub(crate) fn resolve(
        header: &HeaderView,
        selection: &EvidenceSampleSelection,
    ) -> Result<Self, EvidenceError> {
        let mut groups = BTreeMap::new();
        for line in header.as_bytes().split(|byte| *byte == b'\n') {
            if !line.starts_with(b"@RG\t") && line != b"@RG" {
                continue;
            }
            let line = std::str::from_utf8(line).map_err(|_| {
                EvidenceError::InvalidInput("alignment @RG header is not valid UTF-8".into())
            })?;
            let mut id = None;
            let mut sample = None;
            let mut saw_sample = false;
            for field in line.trim_end_matches('\r').split('\t').skip(1) {
                if let Some(value) = field.strip_prefix("ID:") {
                    if id.replace(value.to_owned()).is_some() || value.is_empty() {
                        return Err(EvidenceError::InvalidInput(
                            "alignment @RG requires one nonempty ID".into(),
                        ));
                    }
                } else if let Some(value) = field.strip_prefix("SM:") {
                    if saw_sample {
                        return Err(EvidenceError::InvalidInput(
                            "alignment @RG contains duplicate SM fields".into(),
                        ));
                    }
                    saw_sample = true;
                    sample = (!value.is_empty()).then(|| value.to_owned());
                }
            }
            let id = id.ok_or_else(|| {
                EvidenceError::InvalidInput("alignment @RG is missing its ID".into())
            })?;
            if groups.insert(id.clone(), sample).is_some() {
                return Err(EvidenceError::InvalidInput(format!(
                    "alignment @RG ID {id:?} is declared more than once"
                )));
            }
        }
        let declared_samples: Vec<_> = groups
            .values()
            .filter_map(Option::as_ref)
            .cloned()
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        let has_unnamed = groups.values().any(Option::is_none);
        let (mode, selected_sample) = match selection {
            EvidenceSampleSelection::Auto if declared_samples.is_empty() => {
                (EvidenceSampleMode::Unknown, None)
            }
            EvidenceSampleSelection::Auto if declared_samples.len() == 1 && !has_unnamed => {
                (EvidenceSampleMode::Named, Some(declared_samples[0].clone()))
            }
            EvidenceSampleSelection::Auto => {
                return Err(EvidenceError::InvalidRequest(
                    "alignment sample scope is ambiguous (multiple SM samples or mixed named/unnamed read groups); choose --sample NAME or --pool-samples explicitly".into(),
                ));
            }
            EvidenceSampleSelection::Named(name) => {
                if !declared_samples.contains(name) {
                    return Err(EvidenceError::InvalidRequest(format!(
                        "sample {name:?} is not declared by alignment @RG SM metadata; select a declared --sample or use --pool-samples explicitly"
                    )));
                }
                (EvidenceSampleMode::Named, Some(name.clone()))
            }
            EvidenceSampleSelection::Pool => (EvidenceSampleMode::Pooled, None),
        };
        Ok(Self {
            mode,
            selected_sample,
            declared_samples,
            read_groups: groups
                .into_iter()
                .map(|(id, sample)| EvidenceReadGroup { id, sample })
                .collect(),
            allows_unassigned: mode != EvidenceSampleMode::Named,
        })
    }

    pub(crate) fn includes_record(&self, record: &Record) -> Result<bool, EvidenceError> {
        if self.allows_unassigned {
            return Ok(true);
        }
        let unassigned = |reason: &str| {
            EvidenceError::InvalidInput(format!(
                "read {:?} cannot be assigned to a named sample: {reason}; fix RG/SM metadata or use --pool-samples for deliberate pooling",
                String::from_utf8_lossy(record.qname())
            ))
        };
        let id = match record.aux(b"RG") {
            Ok(Aux::String(id)) if !id.is_empty() => id,
            Ok(_) => return Err(unassigned("RG must be a nonempty string tag")),
            Err(rust_htslib::errors::Error::BamAuxTagNotFound) => {
                return Err(unassigned("missing RG tag"));
            }
            Err(_) => return Err(unassigned("unreadable RG tag")),
        };
        let index = self
            .read_groups
            .binary_search_by(|group| group.id.as_str().cmp(id))
            .map_err(|_| unassigned("RG is absent from the alignment header"))?;
        let sample = self.read_groups[index]
            .sample
            .as_deref()
            .ok_or_else(|| unassigned("RG has no declared SM sample"))?;
        Ok(self.selected_sample.as_deref() == Some(sample))
    }

    /// Conservative additional space for retained scope/request copies and
    /// temporary header resolution during worker creation. Reader header bytes
    /// and decoder allocation remain subject to the existing cooperative model.
    pub(crate) fn memory_bytes(&self) -> u64 {
        let strings = self.selected_sample.as_ref().map_or(0, String::capacity)
            + self
                .declared_samples
                .iter()
                .map(String::capacity)
                .sum::<usize>()
            + self
                .read_groups
                .iter()
                .map(|group| {
                    group.id.capacity() + group.sample.as_ref().map_or(0, String::capacity)
                })
                .sum::<usize>();
        let structures = self.declared_samples.capacity() * std::mem::size_of::<String>()
            + self.read_groups.capacity() * std::mem::size_of::<EvidenceReadGroup>()
            + std::mem::size_of::<Self>();
        (strings as u64)
            .saturating_add(structures as u64)
            .saturating_add((self.read_groups.len() as u64).saturating_mul(256))
            .saturating_mul(4)
    }
}

fn optional_json(value: Option<&str>) -> String {
    value.map_or_else(|| "null".into(), string_json)
}

fn string_json(value: &str) -> String {
    let mut output = String::from("\"");
    for character in value.chars() {
        match character {
            '"' => output.push_str("\\\""),
            '\\' => output.push_str("\\\\"),
            character if character < '\u{20}' => {
                use std::fmt::Write;
                write!(output, "\\u{:04x}", character as u32).unwrap();
            }
            character => output.push(character),
        }
    }
    output.push('"');
    output
}
