# Claim ledger and application wording

The producer, clean replay and independent artifact/cost audits passed at
`401b06b12b875e13e4c101305f41d4bb0ed252fc`,
[run 40151881458](https://github.com/logannye/rosalind/actions/runs/40151881458).
The [retained findings](https://github.com/logannye/rosalind/tree/codex/reproducibility-case-study/docs/findings/reproducibility-case-study-2026-10-24)
contain exact reports, all 27 observations, test logs, archive hashes, and audits.
The source SHA identifies executed code; this later prose is a record of it.

| Claim | Exact metric / denominator | Evidence | Limit |
|---|---|---|---|
| Independent full-field agreement | 1,000 probes × 34 legacy columns; zero differences | Producer full TSV/Arrow and clean replay oracle comparisons | Probe ALTs are not true variants; one specimen/window; shared htslib/specification limits |
| Changed downstream answers | All 27 timed answers, 15,750 total rows, agree | `trials.csv`, retained answer hashes, independent audit | Q1/Q2/Q3 were frozen; no model accuracy claim |
| Clean saved-only replay | Q1/Q2/Q3 on B and D; native Q1 TSV `REPRODUCED`, matching compiled source; 39 original paths absent | `replay-report.json`, original Q1 receipt | Separate clean hosted job, not outside-user adoption; no missing counters interpreted as zero |
| Setup-inclusive three-query cost | A 2.880463 s [2.855964, 2.890056]; B 1.435907 s [1.431209, 1.471755]; D 1.895494 s [1.893129, 1.906313] | Parent-clock prefixes and `cost-audit.json` | Three query rounds, one measured setup; excludes shared acquisition/preparation 33.937795 s and separate evaluation costs |
| Representation tradeoff | Native minimum 3,732,895 bytes; checked Parquet 197,471 bytes; complete bundle 5,840,123 bytes | Inventories and independently summed logical sizes | D shares native extraction and pays conversion; D beat B on Q2/Q3; external binary is another 13,917,792 bytes |
| Observed resources | Preparation peak 297,918,464 bytes; path setup/query A 30,724,096, B 81,440,768, D 117,485,568 bytes | 10 ms simultaneous process-tree RSS samples | Sampled peaks can miss transients; 256 MiB native budget is cooperative, not a whole-workflow cap |
| Refusal and missingness | Consumed-payload corruption: verify/query exit5; absent requests: exit3; REF mismatch: exit2; no successful publication; synthetic zero preserved | Both jobs' controls and independent artifact audit | Corruption tested on copied data; real probes had no zero-callable loci, so zero/missing uses a labeled fixture |
| Historical cohort decision | K=3 analyzer W256 1.251 s; native saved 0.755 s; checked Parquet 0.470 s | [Retained result](https://github.com/logannye/rosalind/blob/70af512ecb803baa61b79abcf68ccb6d42cb990f/docs/findings/bounded-cohort-reanalysis-2026-10-22/README.md), tested merge `a7f3b40daff9aa319a0e089c5af1ef7f37fd5ec4`, run `39816575369` | Synthetic diagnostic; its Parquet cache was independently populated; stopped optimization |
| Historical frozen threshold-transfer decision | False retained 69→53; all retained 5,697→5,672; true retained 4,983→4,981 | [Retained result](https://github.com/logannye/rosalind/blob/70af512ecb803baa61b79abcf68ccb6d42cb990f/research/threshold-transfer/CONTINUATION_RESULTS.md), source `10096419d00c18a03f64098c8c5f35701560bb59`, run `39997188782` | Primary practical gates failed; unknown labels stay distinct; simple baseline retained |

Three-sentence Neuralink application bullet:

> Built a Rust/Python genomic evidence system and independently verified all 34 fields at 1,000 public HG002 loci, then reproduced three changed queries in a clean environment without the original alignments.
> Across three balanced rounds of this fixed three-query workload, native reuse took a median 1.436 s including materialization and verification, versus 2.880 s for fresh extraction and 1.895 s for checked Parquet, excluding shared data acquisition/preparation.
> Verified corrupted-data refusal and preserved earlier no-go decisions, including a richer policy that removed two true candidates and failed practical-benefit gates.

The historical studies are separate runs. No clinical review-time reduction,
BCI benefit, public release, general hard memory limit, cold-cache/p99 result,
raw-to-Parquet speed superiority, or universal speedup is established.
The second sentence concerns the complete fixed sequence, not every individual
query; the measured setup was reused arithmetically, not repeated three times.
