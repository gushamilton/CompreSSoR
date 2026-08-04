# CompreSSoR benchmark records

## Final BP chr1 keyed-format screen

`final-keyed-20260804/` is the authoritative final format comparison. It uses
the real 1,124,344-row FinnGen chr1 core-SNV file on BP, five independent
Slurm tasks, and 34 concrete formats. Every measured format contains the exact
position plus directed REF→ALT identity; no shared spine or external reference
is included in the storage accounting. All 170 observations passed identity and
numeric round-trip checks.

- `final-keyed-runs.csv`: all 170 measured observations;
- `final-keyed-summary.csv`: five-run medians and validation bounds;
- `final-keyed-frontier.csv`: the mathematically calculated Pareto frontier;
- `final-keyed-registry.csv`: measured formats and unavailable bindings;
- `compressor-pareto-final-keyed.png`: the final plot, labelled only at the
  frontier plus the TSV.gz baseline.

The source TSV.gz is 15,186,281 bytes. The default CompreSSoR native Pcodec
store is 4,298,204 bytes (3.53× smaller), with a five-run full-file read median
of 0.321 seconds versus 0.607 seconds for TSV.gz. Its conversion median is
12.178 seconds. The standard store is therefore selected for compression-first
use; the full table shows the speed/storage alternatives.

Zarr, Blosc2, Vortex, and Python-only Pcodec bindings were not installed in the
BP environment. They are recorded as `UNAVAILABLE`, not represented by guessed
numbers.

## Earlier BP chr1 optimisation record

`bp-finngen-chr1-optimization.csv` is the current five-run Linux evidence
from BluePebble compute nodes, using the real 1,124,344-row FinnGen chr1 core
SNV file. It records the identity-frame and Pcodec-page sweep, the selected
131,072-row identity/page geometry. The final key-bearing store is 4,298,204 bytes
versus 15,186,281 bytes for the source TSV.gz (3.53x smaller).

The final same-data Pareto record is split into:

- `pareto-chr1-summary.csv`: all five measured variant-key formats and five access runs;
- `pareto-chr1-access-runs.csv`: individual access timings;
- `pareto-chr1-write.csv`: conversion/write timings;
- `pareto-chr1-validation.csv`: package round-trip and quantisation-bound checks;
- `pareto-chr1-frontier.csv`: the plotted, contract-specific frontier.

The superseded five-format plot is `inst/figures/compressor-pareto.png`. The
authoritative current plot is `inst/figures/compressor-pareto-final-keyed.png`
and the current five-run results are in `final-keyed-20260804/`. Reference-
anchored numeric projections are retained only in older engineering records and
are not used for the headline comparison.

## Earlier release evidence

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

The earlier like-for-like chr1 FinnGen check is retained in
`finngen-chr1-native.csv` and `finngen-chr1-native-runs.csv`. It uses
1,124,344 biallelic SNVs, an eight-column core TSV.gz baseline, and five runs
per workload. Its 4,371,370-byte self-contained store is retained as the
superseded baseline; the current BP result is in
`bp-finngen-chr1-optimization.csv` above.

The clean Linux portability/threading check is recorded in
`pcodec-threads-bp-summary.csv`. It ran as separate one-hour Slurm jobs on
BluePebble compute nodes with R 4.5.1 and the unchanged native Pcodec 1.0.3
backend built using an official Rust 1.87 toolchain. The real SBP GWAS source
yielded 579,586 valid chr1 biallelic SNVs after restricting coordinates to the
GRCh38 chr1 primary range; the 6,181,455-byte core TSV.gz became a
2,584,730-byte self-contained Pcodec store. Five-run medians were:

| Workload | 1 worker | 2 workers | 4 workers |
|---|---:|---:|---:|
| Full core read | 0.156 s | 0.232 s | 0.238 s |
| 1 Mb region | 0.023 s | 0.029 s | 0.029 s |
| 25 canonical keys | 0.832 s | 0.449 s | 0.395 s |
| 1,000 canonical keys | 2.217 s | 1.127 s | 0.844 s |
| Five-store batch, 25 keys | 4.177 s | 2.776 s | 2.441 s |

Serial and parallel reads passed equality checks, with maximum observed Z
reconstruction error 0.00686 and EAF quantisation error 0.00308. Full and
regional reads remain serial because they are already column-streamed; four
workers are useful for sparse and batch access.

## Pareto design record

`first-pareto.csv` is the historical data behind the old Pareto design plot;
the current key-bearing plot is `inst/figures/compressor-pareto-final-keyed.png`.
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
