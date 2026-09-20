# Synthetic paired candidate comparisons

Direction: **right minus left** observed ALT/read-depth fraction. Pairs are explicitly supplied; metadata never infers pairing.

The depth screen is 10 callable reads. It is a technical screen, not confidence or a biological response classification.

| Pair (left → right) | Position / ALT | Left ALT/depth | Right ALT/depth | Exact difference | Both pass depth screen |
|---|---|---|---|---|---|
| a-to-b (specimen-a → specimen-b) | 10 / C | 4/10 | 1/9 (below depth screen) | -13/45 | no |
| a-to-b (specimen-a → specimen-b) | 20 / G | 0/0 (observed; fraction undefined) | 4/12 | undefined | no |
| a-to-b (specimen-a → specimen-b) | 30 / T | 1/10 | 5/10 | 2/5 | yes |
| a-to-b (specimen-a → specimen-b) | 40 / C | 2/12 | unmeasured | undefined | unknown |
| b-to-a (specimen-b → specimen-a) | 10 / C | 1/9 (below depth screen) | 4/10 | 13/45 | no |
| b-to-a (specimen-b → specimen-a) | 20 / G | 4/12 | 0/0 (observed; fraction undefined) | undefined | no |
| b-to-a (specimen-b → specimen-a) | 30 / T | 5/10 | 1/10 | -2/5 | yes |
| b-to-a (specimen-b → specimen-a) | 40 / C | unmeasured | 2/12 | undefined | unknown |
| a-to-c (specimen-a → specimen-c) | 10 / C | 4/10 | 0/10 | -2/5 | yes |
| a-to-c (specimen-a → specimen-c) | 20 / G | 0/0 (observed; fraction undefined) | 0/3 (below depth screen) | undefined | no |
| a-to-c (specimen-a → specimen-c) | 30 / T | 1/10 | unmeasured | undefined | unknown |
| a-to-c (specimen-a → specimen-c) | 40 / C | 2/12 | 0/0 (observed; fraction undefined) | undefined | no |

Unmeasured means no saved observation. Observed zero depth is measured, but has no defined ALT fraction. Positive low-depth fractions remain mathematically defined.

The native TSV and its receipt preserve original integer counts and unreduced difference components. Fractions above are simplified exactly, without floating-point rounding.

[Native paired rows](pairs.tsv) · [Verification and source identity](report.json)

These are authored synthetic inputs, not independent research use or clinical validation.
