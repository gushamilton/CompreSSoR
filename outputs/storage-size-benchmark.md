# CompreSSoR storage-size benchmark

Measured on the Mac mini using the real 16,111,549-row FinnGen R2 ANTIDEPRESSANTS GWAS.
Sizes are exact byte sums across each representation's files; size_runs records five repeated size reads.
The q8 cache is an optional serving layer and is not a standalone replacement for the q9 variant spine.

| Representation | Bytes | Bytes/row | Size vs source gzip | Size vs raw TSV |
|---|---:|---:|---:|---:|
| Source gzip | 391,708,712 | 24.31 | 1.000x | 0.308x |
| Raw TSV | 1,273,085,087 | 79.02 | 3.250x | 1.000x |
| q9 Parquet (convert-only) | 892,353,845 | 55.39 | 2.278x | 0.701x |
| q9 Parquet (GRCh38 QC) | 899,137,039 | 55.81 | 2.295x | 0.706x |
| Exact Parquet (convert-only) | 916,599,396 | 56.89 | 2.340x | 0.720x |
| q8 cache only (optional serving layer) | 128,881,188 | 8.00 | 0.329x | 0.101x |
| q9 Parquet + q8 cache (convert-only) | 1,021,235,033 | 63.39 | 2.607x | 0.802x |

The temporary exact and q8 outputs were deleted after measurement.
