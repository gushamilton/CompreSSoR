# CompreSSoR benchmark records

## Current v0.3 release evidence

The authoritative format/access evidence is split by the question it answers:

- `cold-mr-final-summary.csv`, `cold-mr-final-runs.csv`, and
  `cold-mr-final-metadata.json`: five randomized, fresh-R-process repetitions
  of one-exposure/one-outcome MR, 5x5 MR, 25x25 MR, and complete loads from a
  real 14,923,434-row FinnGen GWAS. They compare direct CompreSSoR/FastMR,
  explicit CompreSSoR reads, TSV.gz, and indexed VCF.gz. Same-request
  coalescing and the software page cache are disabled. Every logical study in
  every format is a distinct fresh `.noindex` scratch copy written through
  `F_NOCACHE`; Pcodec read descriptors also use `F_NOCACHE`. This is described
  as a symmetric cache-controlled cold approximation, not as proof that macOS
  has no filesystem state. The user-owned `mediaanalysisd` process was paused
  during timed trials and resumed afterwards; no root service or Spotlight
  configuration was changed.
- `reader-profile-r.json`, `reader-profile-python.json`, and
  `reader-runtime-profile.csv`: five-run operating-system-warm component
  profiles. These distinguish analysis-ready R/NumPy data from Python's much
  smaller but still encoded bridge representation.
- `pcodec-canonical-access.json`, `pcodec-canonical-access-runs.csv`, and
  `pcodec-v03-stress.json`: earlier v0.3 warm sparse, regional, and full-reader
  latency checks.
- `pcodec-full-api-roundtrip.json`: the deterministic real-data numeric audit.

The measured release fixture occupies 58,033,297 bytes as a self-contained
Pcodec store, versus 201,658,018 bytes for the eight-column TSV.gz and
228,634,485 bytes for indexed VCF.gz plus `.tbi`. Exact REF/ALT identity is
included in the Pcodec size; there is no uncounted variant spine.

The current native implementation check is recorded in
`native-pcodec-se8-frame.csv` and `native-pcodec-se8-frame-runs.csv`. On a
deterministic 1,000,000-row Mac mini fixture, the 0.4.4 native store is
3,067,151 bytes versus 53,154,119 bytes for the source TSV.gz (17.33x), with
five-run warm medians of 58 ms for a full read, 8 ms for a 10 kb region, and
204 ms for 25 canonical-key reads. It uses 8K identity frames and 65K value
frames.
This is an implementation benchmark, not a replacement for the final
full-FinnGen comparison.

## Pareto design record

`first-pareto.csv` is the data behind `inst/figures/compressor-pareto.svg`.
It records the original five-repeat representation sweep that selected the
independent quantised streams. It is a historical design benchmark rather than
the current exported-package timing. Only its median summary survives here: raw
run logs, source byte accounting, exact projection schema, and complete runtime
metadata were not retained, so it must not be used as an independently
reproducible comparison with the release benchmark.

The plot retains an older diagnostic overlay labelled `Native Pcodec 0.4*`.
That point is not the current implementation result. The current native
evidence is in `native-pcodec-se8-frame.csv` and its companion run table: a
1,000,000-row Mac mini smoke fixture, with the same source written as TSV.gz
and read through `fread` versus the 0.4.4 native `.cpr` reader. It is an
implementation check, not a claim that unlike datasets share one Pareto
frontier.

## Earlier engineering records

`pcodec-full-api-benchmark.json` and `pcodec-full-api-runs.csv` are the old
0.2.0 release records. The remaining CSV files document earlier Parquet, q8,
VCF/Tabix, storage, conversion-mode, and release-gate experiments. They are
retained for reproducibility and are not the headline benchmark for the
current Pcodec format. Results depend on hardware, filesystem, software
versions, dataset, and the exact access path.
