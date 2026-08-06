# Full FinnGen p-value versus flag-read benchmark

Dataset: FinnGen_chr1_snvs.
Full valid-SNP FinnGen store: 1,124,344 rows; store bytes: 4,375,127; native threads: 4.
Ordinary z/SE read before flag streams: 0.0280 s; after: 0.0290 s; delta: 0.0010 s.
P-filter reads reconstructed p from the full store and filters in R.
The two MR paths fetch the same columns: chromosome, base_pair_location, effect_allele, other_allele, beta, standard_error, effect_allele_frequency.
Flag inspection decodes the full aligned uint8 flag stream; flag + MR fetch then reads selected rows by native row ID.
P filtering + MR fetch reads reconstructed p for the whole store, filters in R, then reads the same selected MR columns.

| Threshold | Source hits | P-filter hits | Flag bytes | P-filter + MR fetch s | Flag inspection + MR fetch s |
|---:|---:|---:|---:|---:|---:|
| p <= 1e-04 | 165 | 165 | 3,795 | 0.1510 | 0.1120 |
| p <= 1e-05 | 10 | 10 | 3,368 | 0.0810 | 0.0380 |
| p <= 1e-06 | 0 | 0 | 3,251 | 0.0340 | 0.0040 |
| p <= 1e-07 | 0 | 0 | 3,251 | 0.0210 | 0.0090 |
| p <= 1e-08 | 0 | 0 | 3,234 | 0.0200 | 0.0050 |

The flag bytes include the Pcodec uint8 payload plus its small metadata/index.
The p-filter path is the fastest current whole-file API path for this question; p is reconstructed because it is not stored in the native payload.
