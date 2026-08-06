# Full FinnGen p-value versus flag-read benchmark

Dataset: FinnGen_full_ANTIDEPRESSANTS.
Full valid-SNP FinnGen store: 14,923,434 rows; store bytes: 57,632,964; native threads: 4.
Ordinary z/SE read before flag streams: 0.4380 s; after: 0.4260 s; delta: -0.0120 s.
P-filter reads reconstructed p from the full store and filters in R.
The two MR paths fetch the same columns: chromosome, base_pair_location, effect_allele, other_allele, beta, standard_error, effect_allele_frequency.
Flag inspection decodes the full aligned uint8 flag stream; flag + MR fetch then reads selected rows by native row ID.
P filtering + MR fetch reads reconstructed p for the whole store, filters in R, then reads the same selected MR columns.

| Threshold | Source hits | P-filter hits | Flag bytes | P-filter + MR fetch s | Flag inspection + MR fetch s |
|---:|---:|---:|---:|---:|---:|
| p <= 1e-04 | 2,728 | 2,729 | 46,247 | 1.7700 | 1.5910 |
| p <= 1e-05 | 432 | 435 | 40,562 | 0.9390 | 0.6630 |
| p <= 1e-06 | 17 | 17 | 38,718 | 0.4240 | 0.1570 |
| p <= 1e-07 | 1 | 1 | 38,466 | 0.3130 | 0.0520 |
| p <= 1e-08 | 0 | 0 | 38,427 | 0.2840 | 0.0690 |

The flag bytes include the Pcodec uint8 payload plus its small metadata/index.
The p-filter path is the fastest current whole-file API path for this question; p is reconstructed because it is not stored in the native payload.
