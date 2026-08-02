# CompreSSoR real non-EBI end-to-end test

Source: randomly selected FinnGen R2 endpoint `ANTIDEPRESSANTS` from the
[FinnGen public summary-statistics bucket](https://storage.googleapis.com/finngen-public-data-r2/summary_stats/finngen_r2_ANTIDEPRESSANTS.gz).
This is not an EBI/GWAS Catalog file. The source is GRCh38 and uses the
BP/GWASLab-style fields `#chrom`, `pos`, `ref`, `alt`, `rsids`, `pval`,
`beta`, and `sebeta`.

| Measurement | Result |
|---|---:|
| Input gzip | 391,708,712 bytes |
| Input rows | 16,111,549 |
| BP reference rows | 6,271,258 |
| Output rows | 3,294,606 |
| Output Parquet store | 132,687,271 bytes |
| Output/input gzip ratio | 0.338575 |
| End-to-end elapsed time | 68.679 seconds |
| Regional decoded rows | 384 |
| `validate_compressor()` | TRUE |

Rows filtered by the BP fixed-spine process were: 33,372 duplicate rows,
12,302,220 absent from the shared dictionary, 20,074 incompatible allele
rows, and 461,277 unresolved ambiguous rows. The output is the standard
q9/SE10/EAF12 Parquet store with an extra-column sidecar.
