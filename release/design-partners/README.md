# Anonymized advisory release feedback

This directory accepts only completed `feedback.json` records produced by
`cargo xtask partners init` and accepted by `cargo xtask partners validate`.

Do not commit names, contact details, organizations, credentials, genomic data, or
raw interview notes. Keep private material outside the repository or under the
ignored `release/private-design-partners/` directory. One accepted record for each
requested persona remains an adoption objective, not a prerequisite for stable
promotion. `partners validate` and `partners report` still reject invalid or
incomplete feedback. Release planning and RC status report those findings and
accepted-record counts as advisory; absence must never be presented as independent
validation. Feedback collection can continue after publication.

A parsed schema-1 record explicitly marked `blocker_severity: release-blocking`
still reports a technical defect that blocks release until resolved. An incomplete
participant session does not erase that known defect. Resolving the defect does
not turn an unfinished session into successful adoption evidence.
