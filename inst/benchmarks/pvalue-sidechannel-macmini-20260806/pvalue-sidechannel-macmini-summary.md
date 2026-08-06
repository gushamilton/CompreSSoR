# FinnGen p-value side-channel benchmark

Mac mini; 100,000 valid prepared GRCh38 FinnGen rows.
The fixture includes all source valid-SNP rows through p <= 1e-5.
The p-value filter path reads reconstructed p and filters at p <= 1e-05.

| Scenario | Selected rows | Fraction | Baseline B | Duplicate extra B | Flag extra B | P-filter read+filter s | Duplicate extra write s | Flag extra write s |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| p-1e-07 | 1 | 0.0010% | 416,258 | 32,821 | 666 | 0.0040 | 0.1820 | 0.0270 |
| p-1e-06 | 17 | 0.0170% | 416,258 | 33,274 | 687 | 0.0040 | 0.1770 | 0.0140 |
| p-1e-05 | 432 | 0.4320% | 416,258 | 38,233 | 670 | 0.0040 | 0.1890 | 0.0150 |
| designated-10000 | 10,000 | 10.0000% | 416,258 | 75,471 | 6,531 | 0.0040 | 0.2430 | 0.0150 |

Duplicate extra bytes are a second current native Pcodec store containing the selected rows.
Flag extra bytes are a Pcodec uint8 stream plus its small metadata/index, modelled as an additional main-store stream.
The designated-10000 row scenario is a cardinality stress test, not source-level significance.
Write times are incremental side-channel encoding times; all rows are already prepared and the benchmark uses qc=none.
