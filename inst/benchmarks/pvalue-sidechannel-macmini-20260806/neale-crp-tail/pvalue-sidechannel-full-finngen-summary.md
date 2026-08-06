# Full-file p-value versus flag-read benchmark

Dataset: Neale_UKB_CRP_30710_tail.
Full valid-SNP store: 12,524,077 rows; store bytes: 50,706,092; native threads: 4.
Ordinary z/SE read before flag streams: 0.3770 s; after: 0.3640 s; delta: -0.0130 s.
P-filter reads reconstructed p from the full store and filters in R.
The two MR paths fetch the same columns: chromosome, base_pair_location, effect_allele, other_allele, beta, standard_error, effect_allele_frequency.
Flag inspection decodes the full aligned uint8 flag stream; flag + MR fetch then reads selected rows by native row ID.
P filtering + MR fetch reads reconstructed p for the whole store, filters in R, then reads the same selected MR columns.

| Threshold | Source hits | P-filter hits | Flag bytes | P-filter + MR fetch s | Flag inspection + MR fetch s |
|---:|---:|---:|---:|---:|---:|
| p <= 1e-20 | 7,615 | 7,617 | 40,683 | 0.6280 | 0.3830 |
| p <= 1e-30 | 4,160 | 4,163 | 37,166 | 0.4660 | 0.2330 |
| p <= 1e-50 | 2,501 | 2,501 | 35,284 | 0.3770 | 0.1690 |
| p <= 1e-100 | 1,411 | 1,411 | 34,016 | 0.3390 | 0.0890 |
| p <= 1e-150 | 856 | 856 | 33,498 | 0.3370 | 0.0970 |
| p <= 1e-200 | 682 | 682 | 33,297 | 0.3140 | 0.1120 |

The flag bytes include the Pcodec uint8 payload plus its small metadata/index.
The p-filter path is the fastest current whole-file API path for this question; p is reconstructed because it is not stored in the native payload.
