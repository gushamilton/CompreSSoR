# CompreSSoR benchmark records

## Current release benchmark

`pcodec-full-api-benchmark.json`, `pcodec-full-api-runs.csv`, and
`pcodec-full-api-roundtrip.json` are the authoritative 0.2.0 release records.
They measure the exported R API on the Mac mini against a real 14,923,434-row
FinnGen SNP GWAS. The source and all temporary data were held on the same
external SSD. Every timing workload has five complete repetitions; byte sizes
are direct measurements and the round-trip file records one deterministic
10,000-row audit.

The source gzip is 198,128,448 bytes and the self-contained Pcodec store is
57,591,123 bytes (3.44x compression). Median end-to-end compression is 62.157
seconds. Median full-core read time is 1.300 seconds versus 1.584 seconds for
TSV.gz; including reconstructed beta and p takes 1.370 versus 1.623 seconds.
A 1 Mb region takes 0.125 versus 1.557 seconds, and 1,000 sparse rows take
0.374 versus 1.498 seconds.

The round-trip record checks 10,000 rows spread across the genome. Full
position/REF/ALT identity is exact. Its observed numeric errors fall within
the declared Z9/EAF8/SE6 profile, and p is derived from decoded Z.

## Pareto design record

`first-pareto.csv` is the data behind `inst/figures/compressor-pareto.svg`.
It records the original five-repeat representation sweep that selected the
independent quantised streams. It is a historical design benchmark rather than
the current exported-package timing. Only its median summary survives here: raw
run logs, source byte accounting, exact projection schema, and complete runtime
metadata were not retained, so it must not be used as an independently
reproducible comparison with the release benchmark.

## Earlier engineering records

The remaining CSV files document earlier Parquet, q8, VCF/Tabix, storage,
conversion-mode, and release-gate experiments. They are retained for
reproducibility and are not the headline benchmark for the current Pcodec
format. Results depend on hardware, filesystem, software versions, dataset,
and the exact access path.
