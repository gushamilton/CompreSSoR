# Full-file p-value versus flag-read benchmark

Dataset: Neale_UKB_CRP_30710.
Full valid-SNP store: 12,524,077 rows; store bytes: 50,706,092; native threads: 4.
The main pass covers p <= 1e-4 through 1e-15; a separate five-repeat tail pass adds p <= 1e-20, 1e-30, 1e-50, 1e-100, 1e-150, and 1e-200.
P-filter reads reconstructed p from the full store and filters in R.
The two MR paths fetch the same columns: chromosome, base_pair_location, effect_allele, other_allele, beta, standard_error, effect_allele_frequency.
Flag inspection decodes the full aligned uint8 flag stream; flag + MR fetch then reads selected rows by native row ID.
P filtering + MR fetch reads reconstructed p for the whole store, filters in R, then reads the same selected MR columns.

| Threshold | Source hits | P-filter hits | Flag bytes | P-filter + MR fetch s | Flag inspection + MR fetch s |
|---:|---:|---:|---:|---:|---:|
| p <= 1e-04 | 113,378 | 113,387 | 134,011 | 2.4850 | 2.3880 |
| p <= 1e-05 | 74,928 | 74,930 | 101,134 | 2.0100 | 1.9190 |
| p <= 1e-06 | 54,247 | 54,250 | 83,103 | 1.8630 | 1.6760 |
| p <= 1e-07 | 40,822 | 40,825 | 71,125 | 1.4400 | 1.3040 |
| p <= 1e-08 | 33,095 | 33,095 | 64,218 | 1.2850 | 1.1000 |
| p <= 1e-09 | 27,511 | 27,511 | 59,129 | 1.2590 | 1.0240 |
| p <= 1e-10 | 23,281 | 23,287 | 55,378 | 1.1380 | 0.9340 |
| p <= 1e-11 | 19,599 | 19,599 | 52,066 | 1.0310 | 0.7790 |
| p <= 1e-12 | 16,689 | 16,691 | 49,482 | 0.9140 | 0.7050 |
| p <= 1e-13 | 14,797 | 14,797 | 47,670 | 0.8500 | 0.6410 |
| p <= 1e-14 | 13,170 | 13,171 | 46,108 | 0.8300 | 0.5700 |
| p <= 1e-15 | 11,621 | 11,623 | 44,646 | 0.7840 | 0.5520 |
| p <= 1e-20 | 7,615 | 7,617 | 40,683 | 0.6280 | 0.3830 |
| p <= 1e-30 | 4,160 | 4,163 | 37,166 | 0.4660 | 0.2330 |
| p <= 1e-50 | 2,501 | 2,501 | 35,284 | 0.3770 | 0.1690 |
| p <= 1e-100 | 1,411 | 1,411 | 34,016 | 0.3390 | 0.0890 |
| p <= 1e-150 | 856 | 856 | 33,498 | 0.3370 | 0.0970 |
| p <= 1e-200 | 682 | 682 | 33,297 | 0.3140 | 0.1120 |

The flag bytes include the Pcodec uint8 payload plus its small metadata/index.
The p-filter path is the fastest current whole-file API path for this question; p is reconstructed because it is not stored in the native payload.
