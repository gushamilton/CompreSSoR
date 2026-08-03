# CompreSSoR <img src="man/figures/logo.png" align="right" height="150" alt="CompreSSoR logo" />

CompreSSoR imports GWAS summary statistics, optionally lifts them to GRCh38,
harmonises their alleles, and writes a compact store with fast whole-genome,
regional, and sparse-row reads.

The default format is a self-contained, page-indexed Pcodec store. It does not
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
`pcodec`, and `zstandard`. Python 3.12 is recommended; the release runtime is
pinned in `inst/python/requirements.txt` in the source repository:

```bash
python3.12 -m venv ~/.virtualenvs/compressor
~/.virtualenvs/compressor/bin/python -m pip install \
  numpy==1.26.4 pcodec==1.0.3 zstandard==0.25.0
```

From a source checkout, `pip install -r inst/python/requirements.txt` installs
the same pins.

Point R at that interpreter if it is not the first `python3` on `PATH`:

```r
options(CompreSSoR.python = "~/.virtualenvs/compressor/bin/python")
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

# Pay codec startup once when extracting instruments from several GWAS files.
panels <- read_sumstats_batch(
  c(exposure = "exposure.cpr", outcome = "outcome.cpr"),
  c("1:100000:G:A", "2:200000:C:T")
)
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

- liftover only from an explicitly declared GRCh37/hg19 build to GRCh38;
- resolve rsID aliases against the fixed reference when coordinates are absent;
- match the canonical `chrom:pos:REF:ALT` allele identity;
- flip beta, Z, and EAF when the effect allele is reversed;
- allow unambiguous strand complements, while always dropping palindromic
  A/T and C/G matches rather than treating population EAF as strand proof;
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
| chromosome, position, REF, ALT | `key = (global_GRCh38_position << 4) \| (REF << 2) \| ALT`; global position is stored as a Pcodec `uint32` stream and the complete substitution as a Pcodec `uint8` stream constrained to four bits | Decoded exactly; no external reference |
| Z | 9-bit semantic codes over `[-3.5, 3.5)`; one missing code and one float32 exception code | Central maximum absolute quantisation error about 0.00687; tails retained as float32 |
| EAF | 8-bit `asin(sqrt(EAF))` quantisation; invalid/missing values use float32 exceptions | Absolute error bounded by 0.004 in the standard profile |
| SE | 6-bit block-centred residual of `log2(SE)` after conditioning on decoded EAF; missing and float32 exception codes | Positive SE reconstructed on the original scale; relative precision is data-dependent and recorded by profile |
| beta | Not stored | `Z * SE` |
| p-value | Not stored | Two-sided normal p-value from Z |
| rsID / textual variant ID | Not stored | Not reconstructed; identity is chromosome, position, REF, ALT |
| arbitrary extra columns | Not in the default Pcodec payload | Use `backend = "parquet", keep_extras = TRUE` when required |

New v0.3 stores use 262,144-row Pcodec chunks divided into independently
readable 4,096-row pages for every stream. Older v0.2 frame stores remain
readable. SE centre blocks remain 65,536 rows independently of physical page
geometry. Z, EAF, and SE are independent streams. Sparse exceptional values
are held in one Zstandard-compressed float32 sidecar. Page offsets and genomic
bounds support
regional, sparse-row, and canonical `chromosome:position:REF:ALT` lookups
without scanning the whole file. Each header, chunk-metadata record and page
has an index-protected CRC32; the manifest and every durable file also have
SHA-256 records. The complete
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

The current exported R API was measured on the Mac mini against a real
14,923,434-row FinnGen SNP GWAS. Its eight-column gzip contains
`chrom/pos/ALT/REF/beta/se/eaf/p`; the canonical Pcodec store was rebuilt with
the same REF/ALT identity as the indexed VCF and passed full validation across
all 18,220 pages. Durable inputs and stores remained on the external SSD;
temporary decoded bridges used the machine's ordinary temporary directory and
were removed after every call. Access timings are medians of five complete
exported-API runs after one warm-up.

| Storage/write measure | Result |
|---|---:|
| Source TSV.gz | 201,658,018 B |
| Self-contained CompreSSoR store | 58,033,297 B |
| Compression ratio | 3.47x |
| Observed canonical rebuild | 62.982 s |

| Access workload | Median |
|---|---:|
| 25 canonical keys, identity + beta/SE | 0.003 s |
| chr1 1 Mb region, 4,325 rows | 0.003 s |
| Full identity + Z/SE/EAF | 0.729 s |
| Full identity + Z/SE/EAF + reconstructed beta/p | 0.728 s |

The persistent reader starts Python once per R session, caches validated store
metadata and a bounded set of decoded pages, and uses a compact binary bridge.
Across 100 random 25-key sets and 50 random 1-Mb regions, both medians were
0.003 s and both p95 values were 0.004 s. Analyses spanning several GWAS should
use `read_sumstats_batch()`, which also coalesces requests that are literally
identical within a batch. Optimized same-file deduplication and independent-read
comparisons are reported separately so codec I/O is not confused with a faster
MR estimator path. In the controlled 5-exposure x 5-outcome IVW benchmark,
ten explicit Pcodec reads plus the same native grid estimator used by every
format took 0.027 s, versus 0.175 s for ten Tabix queries and 18.651 s for ten
TSV.gz scans. The separately labelled optimized same-store batch took 0.006 s
because its ten identical requests were decoded once.

The timings include Pcodec decompression, compact binary bridges, compiled
reconstruction, R data frames, and—for the grid benchmark—the FastMR estimate.
They are hardware-specific engineering measurements, not universal guarantees.

Machine-readable records:

- [five-run canonical access benchmark](inst/benchmarks/pcodec-canonical-access.json)
- [individual canonical access runs](inst/benchmarks/pcodec-canonical-access-runs.csv)
- [random sparse/region stress record](inst/benchmarks/pcodec-v03-stress.json)
- [historical v0.2 full API benchmark](inst/benchmarks/pcodec-full-api-benchmark.json)
- [earlier real-data numerical audit](inst/benchmarks/pcodec-full-api-roundtrip.json)
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
