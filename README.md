# CompreSSoR <img src="man/figures/logo.png" align="right" height="150" alt="CompreSSoR logo" />

CompreSSoR imports GWAS summary statistics, optionally lifts them to GRCh38,
harmonises their alleles, and writes a compact store with fast whole-genome,
regional, and sparse-row reads.

The default format is a self-contained, block-framed Pcodec store. It does not
write rsIDs, textual variant IDs, chromosome strings, beta, or p-values for
every row. Instead it stores a lossless GRCh38 SNV identity key and three
independent quantised streams: Z9, EAF8, and SE6. Beta and p are reconstructed
when read.

> The identity key is lossless. The standard numeric profile is intentionally
> approximate. Its precision and exception rules are explicit below and in
> every store manifest.

## Installation

```r
# install.packages("remotes")
remotes::install_github("gushamilton/CompreSSoR")
```

The default backend calls a small Python environment containing `numpy`,
`pcodec`, and `zstandard`:

```bash
python3 -m pip install numpy pcodec zstandard
```

Point R at that interpreter if it is not the first `python3` on `PATH`:

```r
options(CompreSSoR.python = "/path/to/python")
# or set COMPRESSOR_PYTHON before starting R
```

The optional Parquet backend and the one-time EBI canonical-reference builder
require the R package `arrow`. Liftover requires the Bioconductor packages
`rtracklayer`, `GenomicRanges`, and `IRanges`. Normal QC mode also requires a
prebuilt, checksum-verified canonical dictionary configured as described below;
the package does not silently manufacture one from a study.
Full Pcodec scans transfer compact integer codes and identity bytes into the
package's compiled reader, which constructs the requested R columns directly;
they do not serialise a decoded GWAS through temporary text.

## Five-function workflow

```r
library(CompreSSoR)

# 1. Import a data.frame, TSV/CSV, .gz/.bgz text file, or single-ALT VCF.
gwas <- import_sumstats("gwas.tsv.gz")

# 2. Optional: lift GRCh37 coordinates and alleles to GRCh38.
gwas38 <- liftover_sumstats(
  gwas,
  input_build = "GRCh37",
  chain = "hg19ToHg38.over.chain.gz"
)

# 3. Match the fixed GRCh38 dictionary and put effects on REF/ALT orientation.
aligned <- harmonise_sumstats(gwas38, reference = "/refs/ebi-grch38")

# 4. Compress. This one function also performs steps 1-3 when requested.
store <- compress_sumstats(
  input = "gwas.tsv.gz",
  output = "gwas.cpr",
  input_build = "GRCh37",
  chain = "hg19ToHg38.over.chain.gz",
  reference = "/refs/ebi-grch38"
)

# 5. Read all rows, a genomic interval, or zero-based stored row IDs.
all_rows <- decompress_sumstats(store)
region <- read_sumstats(store, region = "chr1:1000000-2000000")
sparse <- read_sumstats(store, variants = c(0L, 1000L, 100000L))
validate_compressor(store)
```

For data already harmonised to GRCh38, the minimal conversion path is:

```r
example <- system.file("extdata", "example-grch38.tsv",
                       package = "CompreSSoR")
store <- compress_sumstats(
  example,
  file.path(tempdir(), "example.cpr"),
  mode = "convert",
  reference = NULL,
  assume_grch38_ref_alt = TRUE,
  overwrite = TRUE
)
read_sumstats(store, columns = c("chromosome", "base_pair_location",
                                 "beta", "standard_error", "p_value"))
```

`compress_sumstats()` writes to a staging directory and only replaces the
destination after the complete store succeeds, so a failed overwrite leaves
the previous store intact.

## Inputs and harmonisation

Common aliases are recognised for chromosome/position, effect and other
alleles, beta, odds ratio, SE, Z, EAF, p, -log10(p), and rsID. Odds ratios are
converted to log-odds beta. Missing Z is derived from beta/SE; missing beta is
derived from Z/SE; SE can be inferred from beta and a two-sided Wald p-value.
Unknown columns remain available in memory.

The normal `mode = "qc"` path follows the BP pipeline's conservative rules:

- liftover to GRCh38 when `input_build` is not GRCh38;
- resolve rsID aliases against the fixed reference when coordinates are absent;
- match the canonical `chrom:pos:REF:ALT` allele identity;
- flip beta, Z, and EAF when the effect allele is reversed;
- allow unambiguous strand complements;
- drop unmatched, incompatible, ambiguous, duplicate, zero-mapped, and
  multi-mapped rows by default; and
- record every count and the reference identity/hash in `manifest.json`.

Use `strict = TRUE` to stop instead of dropping unresolved rows. The Parquet
backend can use `drop_unresolved = FALSE` for an explicit audit; canonical
Pcodec stores reject unresolved identities. Pcodec currently stores biallelic
A/C/G/T SNVs on GRCh38 primary chromosomes 1-22, X, and Y. Indels, MNVs,
alternate-contig and mitochondrial rows are dropped and counted by default, or
stop conversion under `strict = TRUE`; multiallelic VCF records must first be
split. The Parquet backend remains available for wider representations.

### Canonical ingestion reference

The harmonisation dictionary is a fixed input, never synthesized from a study.
Build the EBI/Ensembl release-95 GRCh38 dictionary once from downloaded
chromosome VCFs:

```r
build_ebi_reference(
  output = "/refs/ebi-grch38",
  source = "/refs/ebi-harmonisation-vcfs"
)
Sys.setenv(COMPRESSOR_CANONICAL_REFERENCE = "/refs/ebi-grch38")
```

Thereafter `reference = "GRCh38"` resolves that immutable dictionary. The
dictionary is used during ingest but is not needed to decode a Pcodec store:
the per-study key is self-contained.

## Exactly what is stored

Rows are sorted by the identity key before encoding. `A=0`, `C=1`, `G=2`, and
`T=3`.

| Logical field | Physical representation | Read-time result |
|---|---|---|
| chromosome, position, REF, ALT | `key = (global_GRCh38_position << 4) \| (REF << 2) \| ALT`; recursive escaped position gaps plus a 4-bit substitution stream, all block-compressed with Pcodec | Decoded exactly; no external reference |
| Z | 9-bit semantic codes over `[-3.5, 3.5)`; one missing code and one float32 exception code | Central maximum absolute quantisation error about 0.00687; tails retained as float32 |
| EAF | 8-bit `asin(sqrt(EAF))` quantisation; invalid/missing values use float32 exceptions | Absolute error bounded by 0.004 in the standard profile |
| SE | 6-bit block-centred residual of `log2(SE)` after conditioning on decoded EAF; missing and float32 exception codes | Positive SE reconstructed on the original scale; relative precision is data-dependent and recorded by profile |
| beta | Not stored | `Z * SE` |
| p-value | Not stored | Two-sided normal p-value from Z |
| rsID / textual variant ID | Not stored | Not reconstructed; identity is chromosome, position, REF, ALT |
| arbitrary extra columns | Not in the default Pcodec payload | Use `backend = "parquet", keep_extras = TRUE` when required |

Key frames contain 131,072 rows and value frames contain 65,536 rows. Z, EAF,
and SE are independent streams. Sparse exceptional values are held in
value-block-aligned Zstandard-compressed float32 frames. Frame offsets and genomic bounds support
regional and sparse-row reads without scanning the whole file. The complete
layout is specified in [inst/doc/pcodec-format.md](inst/doc/pcodec-format.md).

For a numerically lossless Z/SE/EAF representation, use
`profile = "exact", backend = "parquet"`. This does not preserve an input
p-value that disagrees with Z because p is a derived field in CompreSSoR.

## Why Pcodec?

The historical representation sweep below is the archived median summary from
five full reads of a real 10-million-row GWAS. It was a design experiment—not
an exported-package benchmark—and its raw run logs and full environment record
are not part of this repository. Lower access time and higher compression are
better; the dashed line is its reported Pareto frontier.

![Measured speed versus compression Pareto frontier](inst/figures/compressor-pareto.svg)

That sweep motivated independent quantised numeric streams. Later codec and
identity experiments selected Pcodec and the self-contained key used by the
package.

### Selected-format measurements

The exported R API was measured on the Mac mini against a real 14,923,434-row
FinnGen SNP GWAS. The source was an eight-column gzip containing
`chrom/pos/a1/a2/beta/se/eaf/p`. The source, store, and temporary bridges all
remained on the same external SSD. Timing values below are medians of five
complete runs. Byte sizes are direct file measurements; the compression ratio
is calculated from those sizes, and the 10,000-row round-trip audit was run
once over a deterministic genome-wide sample.

| Storage/write measure | Result |
|---|---:|
| Source TSV.gz | 198,128,448 B |
| Self-contained CompreSSoR store | 57,582,268 B |
| Compression ratio | 3.44x |
| `compress_sumstats()` end-to-end | 63.210 s |

| Access workload | CompreSSoR | TSV.gz | CompreSSoR speedup |
|---|---:|---:|---:|
| Full identity + Z/SE/EAF | 1.487 s | 1.730 s | 1.16x |
| Full identity + Z/SE/EAF + reconstructed beta/p | 1.476 s | 1.790 s | 1.21x |
| chr1 1 Mb region, 6,444 rows | 0.114 s | 1.694 s | 14.86x |
| 1,000 sparse genome-wide rows | 0.387 s | 1.617 s | 4.18x |
| Verify all file checksums | 0.139 s | — | — |
| Decode and validate every frame | 0.430 s | — | — |

The gzip region and sparse workloads necessarily scan the non-indexed source.
The full-core gzip benchmark loads its seven relevant stored columns and
calculates Z; the derived benchmark also loads p. CompreSSoR reconstructs beta
and the two-sided p-value from decoded Z/SE. The Pcodec timings include Python
startup, Pcodec decompression, the compact binary bridge, compiled native
reconstruction, and creation of the final R data frame. These are
hardware-specific engineering measurements, not universal guarantees.

A 10,000-row audit sampled evenly across the same real store. Identity was
exact; maximum observed errors were 0.0069 for Z, 0.0031 for EAF, 1.13% relative
for SE, and 0.0197 absolute for derived beta. P-values exactly followed the
declared function of decoded Z.

Machine-readable records:

- [five-run full API benchmark](inst/benchmarks/pcodec-full-api-benchmark.json)
- [individual benchmark runs](inst/benchmarks/pcodec-full-api-runs.csv)
- [real-data round-trip audit](inst/benchmarks/pcodec-full-api-roundtrip.json)
- `benchmark_table()` for the shipped historical benchmark tables

## Optional Parquet backend

```r
portable <- compress_sumstats(
  aligned, # already harmonised to GRCh38 REF/ALT orientation
  "gwas-parquet.cpr",
  backend = "parquet",
  profile = "exact",
  keep_extras = TRUE,
  reference = NULL,
  mode = "convert",
  assume_grch38_ref_alt = TRUE
)
```

Parquet is useful for Arrow, DuckDB, and Polars interoperability and for exact
double-precision Z/SE/EAF or non-SNV identity. It is not the default compact
format. The legacy q8 serving cache only applies to Parquet stores.

## Development

Reusable code lives in `R/`, `src/`, and `inst/python/`. Large GWAS files,
temporary stores, and timing logs are kept outside the repository; only compact
CSV/JSON summaries and selected figures are versioned.

```bash
R CMD build .
R CMD check CompreSSoR_*.tar.gz --no-manual --as-cran
python3 scripts/test_pcodec_backend.py
```

The package is under active development. Format manifests are versioned; do
not label a standard Pcodec store as a bitwise-exact numerical round trip.
