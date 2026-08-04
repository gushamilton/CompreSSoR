# CompreSSoR

<p align="center">
  <img src="man/figures/logo.png" alt="CompreSSoR logo" width="170">
</p>

<p align="center"><strong>Compact, indexed GWAS summary statistics for R.</strong></p>

[![R-CMD-check](https://github.com/gushamilton/CompreSSoR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gushamilton/CompreSSoR/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

CompreSSoR converts a GWAS summary-statistics table into a self-contained,
indexed `.cpr` store. The standard store is designed for compact storage and
fast whole-genome, regional, and sparse-variant access from R and `fastMR`.

The normal path is:

```text
sumstats → import → optional GRCh38 liftover/QC → Pcodec .cpr store → read/MR
```

![CompreSSoR workflow](docs/figures/compressor-workflow.svg)

## Why this format?

The store keeps the variant identity in every file, without repeating rsIDs or
requiring a study-specific external spine. It stores the core numerical streams
in Pcodec and derives beta and p-value when requested.

On a real 14,923,434-variant FinnGen GWAS:

| Measure | Result |
|---|---:|
| Source TSV.gz | 201,658,018 bytes |
| Self-contained `.cpr` | 58,033,297 bytes |
| Smaller than TSV.gz | **3.47×** |
| Direct 25×25 MR | **1.274 s** |
| Same 25×25 workflow from TSV.gz | 165.266 s |

The headline Pareto result below preserves the measured historical frontier
that motivated the Pcodec representation. The legacy native overlay is kept as
a diagnostic record; the current native implementation benchmark is recorded
separately while the final full-FinnGen comparison is regenerated.

![Compression/access Pareto frontier](inst/figures/compressor-pareto.svg)

The current native check uses a deterministic 1-million-row fixture on the Mac
mini. TSV.gz is 53,154,119 bytes and the self-contained native `.cpr` is
3,067,151 bytes (17.33× smaller), with five-run warm medians of 0.387 s and
0.058 s for full reads. A 10 kb region read is 0.008 s and 25 canonical-key
access is 0.204 s. See the [native SE8 benchmark](inst/benchmarks/native-pcodec-se8-frame.csv).

The native-only smoke benchmark on the Mac mini used one million coherent
rows (`Z ~ N(0, 1)`, `beta = Z × SE`) and five warm full reads:

| Native 0.4.4 measure | Result |
|---|---:|
| Store size | **3,067,151 bytes** |
| Write time | **4.347 s** |
| Full read, all columns | **0.058 s** |
| 10 kb region read | **0.008 s** |
| 25 canonical-key read | **0.204 s** |

Each access number is the median of five runs; the canonical-key test requests
100 keys. The store passed full validation and all checksums. This is a
reproducible engineering check for the installed native backend;
the historical real-GWAS table above predates the native-only format.

The fresh chr1 FinnGen core benchmark is a like-for-like native check on
1,124,344 biallelic SNVs and five repetitions per workload. The baseline is an
eight-column TSV.gz (`chrom`, `pos`, `REF`, `ALT`, `beta`, `SE`, `EAF`, `p`),
not the wider FinnGen annotation table:

| Workload | Pcodec | TSV.gz | Result |
|---|---:|---:|---:|
| Compress end-to-end | 7.038 s | — | — |
| Full core read | 0.073 s | 0.121 s | 1.66× faster |
| Full read with beta/p | 0.075 s | 0.121 s | 1.61× faster |
| 1 Mb region | 0.008 s | 0.132 s | 16.5× faster |
| 25 canonical keys | 0.212 s | 0.123 s | 1.72× slower |
| 1,000 canonical keys | 1.458 s | 0.124 s | 11.8× slower |

The core TSV.gz is 15,186,281 bytes; the self-contained Pcodec store is
4,371,370 bytes (**3.47× smaller**). This exposes the current trade-off
honestly: Pcodec wins whole-file and regional access, while the present
row-by-row sparse-key path still needs optimization. See the [chr1 summary
record](inst/benchmarks/finngen-chr1-native.csv), [five-run record](inst/benchmarks/finngen-chr1-native-runs.csv),
and [reproducible script](scripts/benchmark-finngen-chr1.R).

See [the detailed benchmark record](inst/benchmarks/cold-mr-final-summary.csv)
for all runs and [the technical guide](docs/README-technical.md) for protocol,
limitations, and historical comparisons.

## Installation

The package is native-only: it builds Pcodec through a small Rust-to-R C ABI,
so ordinary reads and writes do not start another process. Install Rust/Cargo
once, then install CompreSSoR:

```bash
# macOS with Homebrew
brew install rust
```

```r
install.packages("remotes")
remotes::install_github("gushamilton/CompreSSoR")
```

If Cargo is missing or the native library cannot be built, installation stops
with an actionable error. The historical Python-backed implementation is
archived in the repository and is not part of the installed package. See the
[technical guide](docs/README-technical.md#native-pcodec-backend) for the
format and build details.

## Quick start

For an already verified GRCh38 REF/ALT table:

```r
library(CompreSSoR)

store <- compress_sumstats(
  "gwas.tsv.gz",
  "gwas.cpr",
  mode = "convert",
  reference = NULL,
  assume_grch38_ref_alt = TRUE,
  overwrite = TRUE
)

validate_compressor(store, full = TRUE)
read_sumstats(
  store,
  columns = c("chromosome", "base_pair_location", "beta", "standard_error")
)
```

For an ordinary third-party GWAS, use the QC path and a configured GRCh38
reference. A GRCh37 input also needs a liftover chain:

```r
store <- compress_sumstats(
  "gwas.tsv.gz",
  "gwas.cpr",
  input_build = "GRCh37",
  chain = "/data/reference/hg19ToHg38.over.chain.gz",
  reference = "GRCh38",
  mode = "qc",
  chrom_threads = 4,
  overwrite = TRUE
)
```

The reference is used during ingestion to establish GRCh38 identity and
allele orientation. The completed standard `.cpr` store carries its own
identity and does not need the reference beside it for reading.

### Optional BP-style variant panels

The default is to keep every valid input variant. To apply a BP common- or
tag-variant panel, pass `variant_set = "common"` or `variant_set = "tag"` and
point the corresponding environment variable at a small compressed panel:

```r
Sys.setenv(COMPRESSOR_COMMON_VARIANTS = "/data/panels/common.parquet")
store <- compress_sumstats(
  "gwas.tsv.gz", "gwas-common.cpr", mode = "convert", reference = NULL,
  assume_grch38_ref_alt = TRUE, variant_set = "common", overwrite = TRUE
)
```

Panels may be Parquet, delimited text, PLINK `.bim`, or a CompreSSoR Pcodec
store. They need either the canonical `variant_id` (`chrom:pos:REF:ALT`) or
`chromosome`, `base_pair_location`, `other_allele` (REF), and `effect_allele`
(ALT). The panel is only a membership filter against the identity key already
stored in each GWAS; it is not a shared spine and is not included in the GWAS
store size. `variant_set = NULL` leaves the main pathway unchanged. The same
flag works with `harmonise_sumstats()` and with `mode = "convert"`, where no
reference harmonisation is performed.

## Reading and FastMR

```r
region <- read_sumstats(
  store,
  region = "chr1:100000000-101000000",
  columns = c("chromosome", "base_pair_location", "beta", "standard_error")
)

instruments <- c("1:109274968:A:G", "6:160589086:C:T")
sparse <- read_sumstats(
  store,
  variants = instruments,
  columns = c("beta", "standard_error")
)
```

`fastMR` can read the requested canonical keys directly from `.cpr` stores and
send the matched beta/SE matrices into its compiled estimator:

```r
remotes::install_github("gushamilton/fastMR")

result <- fastMR::fast_mr_compressed(
  exposure_files = c(exposure = "exposure.cpr"),
  outcome_files = c(outcome = "outcome.cpr"),
  instruments = instruments,
  methods = "ivw"
)
```

The direct path does not load a whole GWAS or reconstruct p-values. It reads
each exposure for its own instruments, reads each outcome for the union, and
uses the shared matrix-grid estimator where possible.

## What is stored?

![CompreSSoR store layout](docs/figures/compressor-store.svg)

Every native 0.4 `.cpr` is a directory with a manifest, checksums,
independently framed Pcodec streams, and a compact exception stream:

```text
gwas.cpr/
├── manifest.json          format, identity, codec, source, and QC metadata
├── manifest.sha256        detached manifest checksum
├── position.pco           global GRCh38 positions
├── native.index.json      block offsets and row counts
├── substitution.pco       four-bit REF→ALT substitution code
├── z.pco                   quantised Z stream
├── eaf.pco                 arcsine-quantised effect-allele frequency
├── se.pco                  block-centred log2-quantised SE stream
└── exceptions.bin          sparse higher-precision numeric exceptions
```

The native 0.4.4 store uses 8,192-row identity frames and 65,536-row numeric
frames by default; `block_rows` changes the numeric frame size when a smaller
random-access frame is preferable. Older 0.2/0.3 stores belong to
the archived pre-native implementation and are not read by this release.

| Logical field | Standard storage | Read-time result |
|---|---|---|
| Chromosome + position | GRCh38 global position, Pcodec `uint32` | Exact chromosome and position |
| REF + ALT | Complete four-bit substitution code | Exact alleles |
| Z | 9-bit semantic code, Pcodec stream | Quantised Z, with exceptions |
| EAF | 8-bit arcsine code, Pcodec stream | Quantised EAF |
| SE | 8-bit block-centred log2 residual | Quantised SE, with exceptions |
| Beta | Not stored per row | Derived as `Z × SE` |
| p-value | Not stored | Derived from Z |
| rsID/text variant ID | Not stored | Use `chrom:pos:REF:ALT` |

Identity is lossless. Standard numeric values are deliberately bounded-lossy;
the manifest records the profile and tolerances. Use the exact Parquet backend
when bitwise-exact doubles or arbitrary non-core fields are required.

### Optional extra columns

The standard Pcodec payload intentionally contains only the core fields above.
This keeps routine MR stores small and predictable. The current released
fallback for retaining arbitrary row-wise columns is:

```r
compress_sumstats(
  input, output,
  backend = "parquet",
  keep_extras = TRUE,
  profile = "exact"
)
```

Those extras are stored as an interoperable Parquet sidecar keyed by the
canonical row number. Pcodec itself can losslessly compress separate numeric
streams, so a typed Pcodec extras layer is feasible; it is not yet part of the
standard `.cpr` contract. We keep that distinction explicit rather than
silently putting arbitrary strings or mixed R columns into a numerical codec.
The design and next implementation step are documented in
[the technical guide](docs/README-technical.md#optional-extra-columns).

## Benchmarks at a glance

Five randomized cold-cache repetitions on the Mac mini, with fresh R
processes and distinct fresh copies for each logical study:

| Input path | 1×1 MR | 5×5 MR | 25×25 MR |
|---|---:|---:|---:|
| **CompreSSoR + FastMR direct** | **0.264 s** | **0.479 s** | **1.274 s** |
| CompreSSoR, explicit reads | 0.387 s | 1.581 s | 7.624 s |
| VCF.gz + Tabix | 0.085 s | 0.330 s | 1.537 s |
| TSV.gz full scans | 6.747 s | 32.830 s | 165.266 s |

Tabix is best for one or a few sparse queries. Direct Pcodec crosses over at
the 25×25 workload. Full-load medians were 3.182 s for Pcodec, 3.954 s for
TSV.gz, and 10.163 s for VCF.gz.

## Further documentation

- [Technical guide, design questions, and full benchmark protocol](docs/README-technical.md)
- [Machine-readable benchmark summary](inst/benchmarks/cold-mr-final-summary.csv)
- [All 75 benchmark repetitions](inst/benchmarks/cold-mr-final-runs.csv)
- [Benchmark metadata and parity gates](inst/benchmarks/cold-mr-final-metadata.json)
- [FastMR direct compressed-MR package](https://github.com/gushamilton/fastMR)

## Verification

```bash
Rscript -e 'testthat::test_local(".")'
R CMD check . --no-manual --as-cran
```

Large GWAS files, reference tables, benchmark stores, and temporary bridges
remain outside the repository. The committed benchmark CSV/JSON/Markdown
files are the durable evidence.

## License

MIT
