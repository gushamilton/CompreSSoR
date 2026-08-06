# Full FinnGen p-value versus flag-read benchmark

Full valid-SNP FinnGen store: 14,923,434 rows; store bytes: 57,632,964; native threads: 4.
Ordinary z/SE read before flag streams: 0.4800 s; after: 0.4140 s; delta: -0.0660 s.
P-filter reads reconstructed p from the full store and filters in R.
Flag-only reads decode the full aligned uint8 flag stream; flag-fetch additionally reads selected z/SE rows by native row ID.

| Threshold | Source hits | P-filter hits | Flag bytes | P-filter s | P-filter + fetch s | Flag-only s | Flag + fetch s |
|---:|---:|---:|---:|---:|---:|---:|---:|
| p <= 1e-05 | 432 | 435 | 40,562 | 0.3350 | 0.7960 | 0.0590 | 0.5070 |
| p <= 1e-06 | 17 | 17 | 38,718 | 0.3340 | 0.4160 | 0.0450 | 0.1110 |
| p <= 1e-07 | 1 | 1 | 38,466 | 0.3500 | 0.3510 | 0.0670 | 0.0510 |

The flag bytes include the Pcodec uint8 payload plus its small metadata/index.
The p-filter path is the fastest current whole-file API path for this question; p is reconstructed because it is not stored in the native payload.
