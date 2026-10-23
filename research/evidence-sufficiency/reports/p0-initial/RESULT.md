# Initial HG002 P0: completed, insufficient assessed calibration data

Hosted run [39878127360](https://github.com/logannye/rosalind/actions/runs/39878127360)
completed successfully at source `c6440019b36c24bc9e1e5f4750c796d79b8aaa09`.
All native suites, 29 existing Python tests, 38 independent oracle tests and
50 new pilot tests passed. Real-data independent counts, native/Parquet feature
equality and saved-only replay passed. No local native build or scientific run
was used.

The scientific gate is **insufficient data**, not evidence of model benefit or
failure. Training contains 492 negative and 2,382 positive assessed candidates.
The calibration and threshold windows each have zero GIAB benchmark-confident
bases: all 5,758 and 9,594 candidates, respectively, have unknown labels.
No models or calibrators were fitted, no test predictions were evaluated, and
model prediction equivalence is therefore unassessed. P1 is disabled.

The analysis took 204.63 seconds, including 88.65 seconds of preparation,
56.54 seconds of candidate scanning and 26.32 seconds of the caller baseline.
Native build/test time is separate. Sampled process-tree peak RSS was
1,565,376,512 bytes and analysis scratch peaked at 5,797,413,211 bytes.
The two native and two Parquet evidence folders total 44,050,502 bytes.
These are observed workflow costs, not a matched storage-engine speed comparison.

All 1,016 files in the retained artifact inventory verified by SHA256 and size,
totaling 215,231,043 bytes including diagnostic features and reports. The full
artifact is retained in Actions for 30 days and in the local ignored
`results/evidence-sufficiency/39878127360` directory. This folder commits the
compact result, exact identities and complete retained-file inventory.

The sole count-driven expansion permitted by the plan is now frozen in
`protocol-p0-expanded.json`. It changes the development blocks and increases the
scope to 15.99 Mb; it is a replacement split, not a strictly nested enlargement.
The old 40–41 Mb test becomes development within the new 39–43 Mb training block.
A fresh 55–59 Mb test stays isolated. Public benchmark-mask coverage is checked
before expensive extraction, without using truth variant labels for candidate
selection. Candidate rules, feature definitions, models and numerical gates
remain unchanged. This expansion was selected from development label counts,
without using test outcomes, and no further count-driven expansion is allowed.
