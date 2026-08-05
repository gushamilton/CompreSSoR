# CompreSSoR

<p align="center">
  <img src="man/figures/logo.png" alt="CompreSSoR logo" width="170">
</p>

<p align="center"><strong>Compact, indexed GWAS summary statistics for R.</strong></p>

[![R-CMD-check](https://github.com/gushamilton/CompreSSoR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gushamilton/CompreSSoR/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

CompreSSoR converts already-prepared GWAS summary statistics into a
self-contained, indexed `.cpr` store for compact storage and fast whole-file,
regional, and variant-level access from R. The core path is:

```text
prepared sumstats → strict schema/orientation QC → native Pcodec .cpr store
```

![CompreSSoR workflow](docs/figures/compressor-workflow.svg)

## Main result

The current headline comparison is the latest five-run BluePebble benchmark
on 10 million real FinnGen SNPs. Every candidate stores its own exact variant
identity key; no external reference is required by the resulting store.

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
├── position.pco           compact build-specific global position stream
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
| Chromosome + position | GRCh37 or GRCh38 global position | Exact chromosome and position |
| REF + ALT | Four-bit directed substitution code | Exact alleles |
| Z | 9-bit semantic code plus exceptions | Quantised Z |
| EAF | 8-bit arcsine code | Quantised EAF |
| SE | 6-bit semantic block-centred log2 residual in a physical `uint8` stream | Quantised SE |
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

The installed compressor accepts a prepared table only. It must contain:
`chromosome`, `base_pair_location`, explicit `reference_allele`/`REF`,
explicit `alternate_allele`/`ALT`, `effect_allele`, `other_allele`, `beta`,
and `standard_error`/`SE`. The prepared orientation is `other_allele = REF`
and `effect_allele = ALT`; the compressor never flips alleles or resolves an
ID against a reference. Input and store builds must match, and both GRCh37
(`hg19`) and GRCh38 (`hg38`) are supported:

```r
library(CompreSSoR)

store <- compress_sumstats(
  "gwas.tsv.gz", "gwas.cpr",
  input_build = "GRCh38", store_build = "GRCh38",
  threads = 4,
  overwrite = TRUE
)

validate_compressor(store, full = TRUE)
read_sumstats(
  store,
  columns = c("chromosome", "base_pair_location", "beta", "standard_error")
)
```

The boundary importer recognises common aliases through a deterministic
resolution matrix: explicit beta plus `SE` is the primary route; an optional
`Z` is checked against `beta / SE` or derived from it; a positive OR may be
converted to log-OR only at the import boundary when beta is absent. P-values
are not silently converted to SE in the strict writer. The standalone
`import_sumstats(..., allow_p_to_se = TRUE)` opt-in is conflict-checked and
records its bounded conversion rows in resolution provenance; it is not used
by `compress_sumstats()`.

For delimited-file input, the native Pcodec path discovers the header first and
projects only columns that can affect the strict identity/statistical core
before canonicalisation and QC. The complete source header and read-column
counts remain in provenance. Parquet with `keep_extras = TRUE` intentionally
reads and retains the full input table; use that mode when arbitrary source
columns are required.

Native and default compression manifests retain aggregate structural-QC counts,
bounded row examples, source/projection provenance, and phase timings for read,
normalisation, QC, identity/sort, encoding, and the final atomic commit. Full
row-level QC status and canonical-key audits remain available through the
explicit `preflight_sumstats()` path rather than being retained through a
large compression job.

For a prepared GRCh37 table, write a native GRCh37 store directly:

```r
store <- compress_sumstats(
  "gwas.tsv.gz", "gwas.cpr",
  input_build = "hg19", store_build = "GRCh37",
  threads = 4,
  overwrite = TRUE
)
```

`input_build = "GRCh37", store_build = "GRCh38"` is rejected. Perform any
liftover, reference lookup, allele alignment, rsID resolution, and source-file
specific field mapping in an upstream preparation workflow, then pass its
explicit result to CompreSSoR.

### External preparation and archived workflows

The former reference-backed harmonisation and liftover implementation is kept
under [`archive/harmonisation/`](archive/harmonisation/). It is not sourced by
the installed package and is retained only as a migration reference. An
external adapter—such as a GWASLab/MR-Atlas preparation step—should record its
own reference, chain, QC, and timing provenance before calling the strict core.
Panel selection (`core`, `hm3`, and `core_plus`) remains available after the
strict contract has been satisfied; it is filtering, not harmonisation.

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
