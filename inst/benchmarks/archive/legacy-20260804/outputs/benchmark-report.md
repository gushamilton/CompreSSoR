# CompreSSoR benchmark record

This output records the first measured Pareto benchmark used to choose the
standard storage profile. It is a five-run full-read benchmark on a real
10-million-row GWAS.

- Speed: median wall-clock seconds to read the whole dataset; lower is faster.
- Compression: ratio relative to the source gzip file; higher is smaller.
- Pareto: formats retained on the measured speed/compression frontier.

| Format | Median full read (s) | Compression vs source gzip | Pareto |
|---|---:|---:|:---:|
| Source gzip | 0.922 | 1.00× | No |
| Raw float32 binary | 0.0214 | 2.04× | No |
| Raw float32 binary + zstd | 0.1704 | 2.67× | No |
| Current f16/u8 binary | 0.00616 | 5.44× | Yes |
| Current f16/u8 zstd | 0.0744 | 7.71× | No |
| Float32 Parquet | 0.0278 | 2.57× | No |
| Precise Parquet | 0.0308 | 2.37× | No |
| q9 framed binary | 0.0696 | 8.27× | No |
| q9 Parquet | 0.0159 | 7.89× | Yes |
| q8 Parquet | 0.0171 | 8.97× | Yes |
| q8 framed binary | 0.0483 | 9.26× | Yes |

The benchmark is a design record, not a hardware-independent performance
promise. It supports q9/SE10 Parquet as the standard durable store and q8 or
framed representations as optional serving caches.

## Conversion/QC mode benchmark

The mode implementation was checked separately on a balanced synthetic
550,000-row input spanning 22 chromosomes, with three runs per scenario:

| Scenario | Median seconds | Rows retained |
|---|---:|---:|
| Serial QC | 0.374 | 550,000 |
| Four-worker chromosome-parallel QC | 0.956 | 550,000 |
| Four-worker chromosome-parallel core panel | 1.130 | 275,000 |

The fixture is deliberately an edge/mode test rather than a hardware claim.
It verifies that QC and chromosome-parallel output agree and that panel
filtering is accounted for. On this small input, worker startup costs more
than it saves. A full real preserve-all run was also memory/imbalance-bound
because chromosome 1 dominates, so the current chromosome-parallel path is
not yet presented as a large-GWAS acceleration. The next performance step is
a streaming chromosome-batch scheduler.
