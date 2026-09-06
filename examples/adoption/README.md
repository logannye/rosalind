# Blank adoption session record

Use [the validation kit](../../docs/adoption-validation.md) for the three task
protocols, measurement definitions, privacy rules, and existing release-gate mapping.

Copy [session.template.json](session.template.json) to a private location for each
participant/task attempt. It conforms to [session.schema.json](session.schema.json).
Use [report.template.md](report.template.md) when summarizing actual observations.
Use `record_kind: "session"` for an actual attempt; all unknown observations remain
null. The committed template describes no participant, run, success, or return.
`attempts` starts empty because zero attempts have been collected here.

This is supplemental research information. It is not a `feedback.json` accepted by
the existing release partner gate, and it adds no release requirement. The schema
allows anonymized summaries, not identity or genomic-data fields; schema validation
cannot establish that free text is anonymous. Review consent and content before
sharing. Raw logs/receipts and the participant-to-pseudonym mapping stay private.

For each attempt, record its role (`initial`, `repeat`, `reuse`, `baseline`, or
`verification`), command/tool label, wall time, exit status, cache state, resource
limits/measurement source and available visit/reuse counters. Use the receipt's
claim hash to reference private evidence without adding its local path to the form.
Elapsed installation time includes waiting; active effort and assistance are
separate person-minute measurements. A negative savings result is a valid result.

`followup_30_days.returned_to_real_work` stays null until an observation is recorded.
When the observation is unreachable, do not turn missing evidence into `false`.
Successful return requires an actual useful task at least 30 days after the first
completed task, with a distinct anonymous team ID for the roadmap's team count.
