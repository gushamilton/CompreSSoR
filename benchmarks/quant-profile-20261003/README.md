# Semantic quantisation profile grid (2026-10-03)

Decision record for the default native profile `z10/eaf8/se8+xse`
(format 0.4.6). Runs on BluePebble (`short`/`test`, Cascade Lake for the
grid; 8 threads), CompreSSoR at the PR branch, fastMR 5a189f7.

- `scripts/grid.R`: one dataset x profile cell. Three timed
  `compress_sumstats()` writes from an in-memory table, three timed full
  reads, then error and p-threshold membership against the source.
  Datasets: FinnGen 10M (`finngen_10m_snps.tsv.gz`), simulated `exp01`,
  `out01` (showcase traits, prepared as `convert_tsv_to_cpr()`).
- `scripts/mr.R`: 10 x 10 `fast_mr_compressed()` (IVW, Egger, weighted
  median, simple and weighted mode; nboot 1000, seed 1) with the instruments
  fixed to the showcase 10x10 C8 rep1 lists. The `exact` arm runs the same
  code with the store reader replaced by the TSV values.
- `scripts/compat.R`: same-machine reads of Z9/SE6 stores with a27b32d and
  the PR build, compared by SHA-256 of `serialize()`.

Results: `results/grid_size_time.csv`, `grid_errors.csv`, `grid_flips.csv`
(reference `exact_z_p` = p from source beta/SE; `source_p` = supplied p),
`mr_deviation.csv` (|b_cpr - b_exact| / se_exact per pair and method),
`mr_store_sizes.csv`. Timings are medians of 3; write times vary by node
(the cells did not all land on one node), so only sizes and errors are
compared across profiles.
