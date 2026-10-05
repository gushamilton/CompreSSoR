# Batch strategy benchmark (read_sumstats_batch: one pass per store vs panel sharing)

`bench_share.R` times one `read_sumstats_batch()` call in a fresh R process with
`options(CompreSSoR.batch_share_panels)` set to `FALSE` (one), `TRUE` (share) or unset (auto, the default
rule). It also records the default rule's choice and a digest of the result; the digest is identical across
modes in every cell. `bench_share.sbatch` runs the grids on BluePebble (Gold 6226R, 8 Slurm CPUs, 3 replicates,
modes in random order). `mk_ukb_keys.R` draws the UKB-PPP key lists.

Store sets:
- showcase: 20 or 50 simulated GWAS on one 1000G panel.
- ukbppp: a seeded sample of the 2,940 UKB-PPP protein stores. Their panels repeat in groups; the repeat
  fraction is 0.15, 0.35, 0.56, 0.66 and 0.83 for 20, 60, 160, 300 and 2,940 stores.
- mix: the 20 showcase stores plus 5, 10 or 20 UKB-PPP stores.

Results:
- `results_e551e43_rule075.csv`: full grid on CompreSSoR e551e43 plus the rule with a 0.75 threshold. The one
  and share columns do not depend on the threshold.
- `results_c4f4eb2_confirm.csv`: confirmation of the final rule (0.9).
