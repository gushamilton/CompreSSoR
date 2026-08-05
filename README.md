# CompreSSoR

<p align="center">
  <img src="man/figures/logo.png" alt="CompreSSoR logo" width="170">
</p>

<p align="center"><strong>Compact, indexed GWAS summary statistics for R.</strong></p>

[![R-CMD-check](https://github.com/gushamilton/CompreSSoR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gushamilton/CompreSSoR/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

CompreSSoR converts GWAS summary statistics into a self-contained, indexed
`.cpr` store for compact storage and fast whole-file, regional, and
variant-level access from R. The normal path is:

```text
sumstats → import/QC → optional GRCh38 liftover → native Pcodec .cpr store
```

![CompreSSoR workflow](docs/figures/compressor-workflow.svg)

## Main result

The current headline comparison is the latest five-run BluePebble benchmark
on 10 million real FinnGen SNPs. Every candidate stores its own exact variant
identity key; the shared reference is excluded from storage accounting.

The locked storage contract is documented in
[`docs/pcodec-final-spec.md`](docs/pcodec-final-spec.md).

| Format | Size | Whole-file read |
|---|---:|---:|
| **CompreSSoR native Pcodec, 4 threads** | **38.06 MB** | **0.451 s** |
| Parquet q9 | 69.48 MB | 1.410 s |
| TSV.gz | 135.39 MB | 4.955 s |

The Pcodec store is 3.56× smaller than TSV.gz. The authoritative records and
the archive policy are in [`benchmarks/`](benchmarks/README.md); older plots
and benchmark families are under
[`inst/benchmarks/archive/legacy-20260804/`](inst/benchmarks/archive/legacy-20260804/).

![10m FinnGen keyed-format Pareto frontier](inst/benchmarks/finngen-10m-bp-20260804/format-screen/pareto-10m-bp.png)

## What is stored?

![CompreSSoR store layout](docs/figures/compressor-store.svg)

The standard `.cpr` directory contains a manifest, checksums, an exact
variant-identity stream, independently framed numerical streams, and sparse
exceptions:

```text
gwas.cpr/
├── manifest.json          format, identity, codec, source, and QC metadata
├── manifest.sha256        detached manifest checksum
├── position.pco           compact global GRCh38 position stream
├── native.index.json      block offsets and row counts
├── substitution.pco      four-bit REF→ALT identity code
├── z.pco                  quantised Z stream
├── eaf.pco                arcsine-quantised EAF stream
├── se.pco                 block-centred log2-quantised SE stream
└── exceptions.bin         sparse higher-precision numeric exceptions
```

The identity key is stored in every GWAS, so an external variant spine is not
required for access or comparison. The standard numerical representation is:

| Logical field | Stored representation | Read-time result |
|---|---|---|
| Chromosome + position | GRCh38 global position | Exact chromosome and position |
| REF + ALT | Four-bit directed substitution code | Exact alleles |
| Z | 9-bit semantic code plus exceptions | Quantised Z |
| EAF | 8-bit arcsine code | Quantised EAF |
| SE | 8-bit block-centred log2 residual | Quantised SE |
| Beta | Not stored per row | Derived as `Z × SE` |
| p-value | Not stored | Derived from Z |
| rsID/text ID | Not stored | Use `chrom:pos:REF:ALT` |

This is intentionally a core numerical format. Exact or arbitrary extra
columns can use the Parquet backend; see
[the technical guide](docs/README-technical.md).

## Install

The native backend builds Pcodec through a small Rust-to-R C ABI. Install Rust
once, then install the package:

```bash
# macOS with Homebrew
brew install rust
```

```r
install.packages("remotes")
remotes::install_github("gushamilton/CompreSSoR")
```

The historical Python backend is archived and is not part of the installed
package.

## Quick start

For a GRCh38-aligned, verified table:

```r
library(CompreSSoR)

store <- compress_sumstats(
  "gwas.tsv.gz", "gwas.cpr",
  mode = "convert",
  assume_grch38_ref_alt = TRUE,
  overwrite = TRUE
)

validate_compressor(store, full = TRUE)
read_sumstats(
  store,
  columns = c("chromosome", "base_pair_location", "beta", "standard_error")
)
```

For an ordinary third-party GWAS, use the harmonisation/QC path and provide a
GRCh38 reference or a GRCh37-to-GRCh38 chain:

```r
store <- compress_sumstats(
  "gwas.tsv.gz", "gwas.cpr",
  input_build = "GRCh37",
  chain = "/data/reference/hg19ToHg38.over.chain.gz",
  reference = "GRCh38",
  mode = "qc",
  chrom_threads = 4,
  overwrite = TRUE
)
```

The reference is an ingestion dependency. The completed `.cpr` store carries
its own identity and does not need the reference beside it for reading.

For a core-plus store, harmonise first and then select the canonical core panel
plus a 50 kb window around variants with derived `p < 1e-5`:

```r
harmonised <- harmonise_sumstats(
  "gwas.tsv.gz", reference = "GRCh38", mode = "qc", chrom_threads = 4
)
store <- compress_sumstats(
  harmonised, "gwas-core-plus.cpr", reference = NULL,
  mode = "core_plus", variant_set = "/data/panels/core_by_chrom",
  overwrite = TRUE
)
```

The selection is deterministic, recorded in the manifest, and accompanied by a
small region sidecar. See [the panel preparation guide](docs/variant-panels.md)
for reproducing the canonical core/HM3 inputs and chromosome shards.

## Reading and FastMR

```r
region <- read_sumstats(
  "gwas.cpr",
  region = "chr1:100000000-101000000",
  columns = c("chromosome", "base_pair_location", "beta", "standard_error")
)

instruments <- c("1:109274968:A:G", "6:160589086:C:T")
sparse <- read_sumstats(
  "gwas.cpr", variants = instruments,
  columns = c("beta", "standard_error"), threads = 4
)
```

The companion [fastMR](https://github.com/gushamilton/fastMR) package can use
the same canonical keys to read only the required exposure/outcome values.

## Further documentation

- [Benchmark write-up and reproducible records](benchmarks/README.md)
- [Technical guide and format details](docs/README-technical.md)
- [R package reference](https://gushamilton.github.io/CompreSSoR/)

## Verification

```bash
Rscript -e 'testthat::test_local(".")'
R CMD check . --no-manual --as-cran
```

Large GWAS files, reference tables, benchmark stores, and temporary bridges
remain outside the repository. The committed CSV/JSON/Markdown records are
the durable evidence.

## License

MIT
