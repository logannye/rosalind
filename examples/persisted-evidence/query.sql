-- DuckDB: replace evidence_export with the committed export directory.
-- UBIGINT columns remain unsigned. SUM uses a wider exact integer accumulator.
SELECT contig,
       count(*) AS requested_positions,
       count(*) FILTER (WHERE callable_depth = 0) AS zero_depth_positions,
       sum(callable_depth) AS callable_read_observations,
       count(*) FILTER (WHERE callable_depth >= 10) AS positions_at_least_10x
FROM read_parquet('evidence_export/part-*.parquet')
GROUP BY contig
ORDER BY contig;

-- Arrays follow A,C,G,T order; DuckDB list indexing is one-based.
SELECT contig, pos,
       c AS c_reads,
       allele_base_quality_sum[2] AS c_quality_sum
FROM read_parquet('evidence_export/part-*.parquet')
WHERE c > 0
ORDER BY contig, pos;
