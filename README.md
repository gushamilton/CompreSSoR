# CompreSSoR <img src="man/figures/logo.png" align="right" height="150" alt="CompreSSoR logo" />

[![R-CMD-check](https://github.com/gushamilton/CompreSSoR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gushamilton/CompreSSoR/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

**Compact, indexed GWAS summary statistics for R.**

CompreSSoR takes heterogeneous GWAS summary-statistics files, optionally lifts
and harmonises them to GRCh38, and writes a compact store designed for fast
whole-genome, regional, and sparse-variant access.

The default `.cpr` format is self-contained. It stores exact GRCh38
`chromosome:position:REF:ALT` identity plus quantised Z, EAF, and SE streams.
It does not repeat rsIDs, textual variant IDs, beta, or p-values in every
study. Beta and p are reconstructed when requested.

> Variant identity is lossless. The standard numeric profile is deliberately
> lossy and has explicit, tested error bounds. Use the Parquet exact profile
> when double-precision numerical round trips are required.

## At a glance

Measured on a real 14,923,434-variant FinnGen GWAS on the Mac mini:

| Result | CompreSSoR v0.3 |
|---|---:|
| Eight-column source TSV.gz | 201,658,018 bytes |
| Self-contained `.cpr` store | 58,033,297 bytes |
| Compression versus TSV.gz | **3.47×** |
| Compression versus indexed VCF.gz + `.tbi` | **3.94×** |
| Compression time | 62.982 s |
| Warm R load, Z/beta/SE/EAF | 0.162 s |
| Warm R load, identity + Z/beta/SE/EAF | 0.463 s |
| Equivalent fully decoded Python/NumPy load | 1.539 s |
| Cache-controlled full load + traversal | **3.182 s** |
| Cache-controlled direct 25x25 MR | **1.274 s** |
| Same 25x25 workflow from TSV.gz | 165.266 s |

The reader figures are medians of five exported-API runs with the operating-
system cache warm; they answer where R/Python materialisation time is spent.
The cache-controlled 1x1, 5x5, 25x25 MR and complete-load benchmark below is
the end-to-end comparison. These are measured engineering results, not
universal guarantees. Full provenance and individual runs are committed under
[`inst/benchmarks`](inst/benchmarks).

## Contents

- [Installation](#installation)
- [Quick start](#quick-start)
- [Choose the right workflow](#choose-the-right-workflow)
- [Input formats and recognised columns](#input-formats-and-recognised-columns)
- [GRCh38 liftover and harmonisation](#grch38-liftover-and-harmonisation)
- [The canonical ingestion reference](#the-canonical-ingestion-reference)
- [Reading compressed GWAS files](#reading-compressed-gwas-files)
- [Exactly what is stored](#exactly-what-is-stored)
- [Accuracy and round-trip guarantees](#accuracy-and-round-trip-guarantees)
- [Integrity and validation](#integrity-and-validation)
- [Benchmarks](#benchmarks)
- [FastMR integration](#fastmr-integration)
- [Optional exact Parquet stores](#optional-exact-parquet-stores)
- [Troubleshooting and FAQ](#troubleshooting-and-faq)
- [Development and reproducibility](#development-and-reproducibility)

## Installation

### 1. Install the R package

The current v0.3 development release is on the pull-request branch:

```r
install.packages("remotes")
remotes::install_github(
  "gushamilton/CompreSSoR@agent/pcodec-r-package",
  upgrade = "never"
)
```

After that branch is merged, install the main branch with:

```r
remotes::install_github("gushamilton/CompreSSoR")
```

The package requires R, a C++17 compiler, and the R packages listed in
[`DESCRIPTION`](DESCRIPTION). The optional capabilities have separate
dependencies:

| Capability | Additional dependency |
|---|---|
| Default Pcodec storage | Python ≥3.10 with `numpy`, `pcodec`, `zstandard` |
| Parquet exact/interoperable storage | R packages `arrow` and `dplyr` |
| GRCh37 → GRCh38 liftover | `rtracklayer`, `GenomicRanges`, `IRanges` |
| Building the EBI reference | R package `arrow` plus downloaded EBI VCFs |

### 2. Install the pinned Pcodec runtime

Python 3.12 is recommended. The exact release pins are in
[`inst/python/requirements.txt`](inst/python/requirements.txt).

```bash
python3.12 -m venv ~/.virtualenvs/compressor
~/.virtualenvs/compressor/bin/python -m pip install \
  numpy==1.26.4 pcodec==1.0.3 zstandard==0.25.0
```

Tell R which interpreter to use:

```r
options(CompreSSoR.python = "~/.virtualenvs/compressor/bin/python")
```

Alternatively, set `COMPRESSOR_PYTHON` before starting R:

```bash
export COMPRESSOR_PYTHON="$PWD/.venv/bin/python"
```

CompreSSoR probes the selected Python version and wrapped Pcodec API before it
writes or reads a Pcodec store. A missing or incompatible runtime fails with an
explicit setup error rather than silently switching formats.

### 3. Smoke-test the installation

This example is intentionally reference-free because the bundled fixture is
already labelled as GRCh38 REF/ALT data:

```r
library(CompreSSoR)

input <- system.file(
  "extdata", "example-grch38.tsv",
  package = "CompreSSoR"
)
output <- file.path(tempdir(), "example.cpr")

store <- compress_sumstats(
  input,
  output,
  mode = "convert",
  reference = NULL,
  assume_grch38_ref_alt = TRUE,
  overwrite = TRUE
)

print(store)
validate_compressor(store, full = TRUE)
read_sumstats(store, columns = c(
  "chromosome", "base_pair_location", "effect_allele", "other_allele",
  "beta", "standard_error", "p_value"
))
```

## Quick start

For an ordinary GWAS that needs proper GRCh38 QC, the complete workflow is:

```r
library(CompreSSoR)

Sys.setenv(
  COMPRESSOR_CANONICAL_REFERENCE = "/data/reference/ebi-ensembl95-grch38"
)

store <- compress_sumstats(
  input = "gwas.tsv.gz",
  output = "gwas.cpr",
  input_build = "GRCh37",
  chain = "/data/reference/hg19ToHg38.over.chain.gz",
  reference = "GRCh38",
  mode = "qc",
  chrom_threads = 4,
  overwrite = TRUE
)

validate_compressor(store)

# Read only what the analysis needs.
region <- read_sumstats(
  store,
  region = "chr1:100000000-101000000",
  columns = c("chromosome", "base_pair_location", "beta", "standard_error")
)

instruments <- read_sumstats(
  store,
  variants = c("1:109274968:A:G", "6:160589086:C:T"),
  columns = c("beta", "standard_error", "effect_allele_frequency")
)
```

`compress_sumstats()` performs import, optional liftover, reference matching,
allele harmonisation, QC, sorting, compression, manifest creation, integrity
sealing, and atomic installation of the completed store.

## Choose the right workflow

Most users should choose `mode = "qc"`.

| Workflow | Use when | Reference matching | Variant filtering |
|---|---|---|---|
| `mode = "qc"` | Ordinary third-party GWAS; recommended default | Required | Drops unresolved and unsupported rows, recording every count |
| `mode = "convert"` | Data are already verified GRCh38 REF/ALT and only need encoding | Skipped | Backend restrictions still apply |
| `mode = "all"` | Backward-compatible spelling of QC mode | Required | Same as `qc` |
| `mode = "core"` | QC plus a caller-supplied core panel | Required | Keeps panel members |
| `mode = "hm3"` | QC plus a caller-supplied HapMap3 panel | Required | Keeps panel members |

### Recommended QC conversion

```r
store <- compress_sumstats(
  "gwas.tsv.gz",
  "gwas.cpr",
  reference = "GRCh38",
  mode = "qc",
  strict = FALSE
)
```

With `strict = FALSE`, unresolved rows are dropped and counted in
`manifest.json`. With `strict = TRUE`, the first unsupported, ambiguous,
duplicated, incompatible, or unmatched condition stops conversion.

### Trusted minimal conversion

```r
store <- compress_sumstats(
  "already-harmonised-grch38.tsv.gz",
  "gwas.cpr",
  mode = "convert",
  reference = NULL,
  assume_grch38_ref_alt = TRUE
)
```

Set `assume_grch38_ref_alt = TRUE` only when all of the following are known:

- coordinates are on GRCh38 primary chromosomes;
- `other_allele` is the GRCh38 REF allele;
- `effect_allele` is ALT; and
- beta, Z, and EAF refer to ALT.

Ordinary effect/other allele columns do not prove REF/ALT orientation. If the
input has explicit `REF` and `ALT` column names, CompreSSoR records that fact
and does not require the override.

### Panel-restricted conversion

`mode = "core"` and `mode = "hm3"` do not use hidden package panels. Supply
the exact panel intended for the project:

```r
store <- compress_sumstats(
  "gwas.tsv.gz",
  "gwas-hm3.cpr",
  reference = "GRCh38",
  mode = "hm3",
  variant_set = "/data/reference/hm3.bim"
)
```

Delimited, Parquet, and PLINK `.bim` panel files are supported. Panel identity
can be expressed as canonical keys, rsIDs, or chromosome/position/alleles.

## Input formats and recognised columns

`import_sumstats()` accepts:

- an R `data.frame`;
- TSV, CSV, and other delimited text accepted by `data.table::fread()`;
- `.gz` or `.bgz` compressed delimited text;
- VCF, VCF.gz, or VCF.bgz with one ALT allele and at most one trait/sample
  column; and
- common GWAS-VCF fields in INFO or FORMAT.

Multiallelic VCF records must be split before import. Malformed delimited rows,
conflicting duplicate aliases, nonnumeric values in numeric columns, invalid
odds ratios, and inconsistent beta/Z/SE values fail explicitly.

Representative recognised aliases are:

| Canonical field | Common input names |
|---|---|
| chromosome | `CHR`, `CHROM`, `#CHROM`, `chrom` |
| base-pair position | `BP`, `POS`, `position`, `base_pair_location` |
| effect allele | `EA`, `A1`, `ALT`, `effect_allele` |
| other allele | `NEA`, `OA`, `A2`, `REF`, `other_allele` |
| beta | `BETA`, `ES`, `LOGOR`, `effect_size`, `estimate` |
| odds ratio | `OR`, `ODDSRATIO` |
| standard error | `SE`, `SEBETA`, `STDERR`, `std_err` |
| Z statistic | `Z`, `ZSCORE`, `ZSTAT`, `z_score` |
| effect-allele frequency | `EAF`, `AF`, `A1FREQ`, `ALT_AF` |
| p-value | `P`, `PVAL`, `PVALUE` |
| −log10(p) | `LP`, `-LOG10P`, `neglog10p` |
| rsID | `rsid`, `RSID`, `dbsnp` |
| variant ID | `SNP`, `SNPID`, `ID`, `variant` |
| sample size | `N`, `SS`, `N_TOTAL`, `TOTALSAMPLESIZE` |
| imputation information | `INFO`, `INFO_SCORE`, `IMPINFO`, `R2`, `RSQ` |

Only effect and other alleles are universally required at import. Coordinates
may be resolved from rsID through the configured canonical reference during QC.
At least one usable association shape must ultimately be available:

- beta and SE;
- Z and SE; or
- beta and a two-sided Wald p-value, from which SE can be inferred.

Odds ratios are converted to log-odds beta. Missing beta or Z is reconstructed
when SE permits it. Input p-values are an ingestion aid and are not preserved
as an independent physical column.

Unrecognised columns survive `import_sumstats()` in memory. The compact Pcodec
backend intentionally stores only the core representation; use the Parquet
backend with `keep_extras = TRUE` when row-level extra columns must remain in
the durable output.

## GRCh38 liftover and harmonisation

### Liftover

CompreSSoR accepts GRCh37/hg19 or GRCh38/hg38 labels. Other source-build labels
are rejected rather than guessed.

```r
imported <- import_sumstats("gwas.tsv.gz")

lifted <- liftover_sumstats(
  imported,
  input_build = "GRCh37",
  chain = "/data/reference/hg19ToHg38.over.chain.gz",
  drop_unmapped = TRUE
)

attr(lifted, "liftover_stats")
```

Liftover uses a caller-supplied chain through `rtracklayer`. It:

- keeps only uniquely mapped primary GRCh38 targets by default;
- records zero maps, multiple maps, and non-primary targets;
- reverse-complements both alleles for reverse-strand mappings; and
- preserves row-level `liftover_status` when called directly.

No jobs or remote services are involved in liftover.

### Allele harmonisation

```r
aligned <- harmonise_sumstats(
  lifted,
  reference = "GRCh38",
  mode = "qc",
  strict = FALSE,
  chrom_threads = 4
)

attr(aligned, "alignment_stats")
attr(aligned, "reference_hash")
```

The QC path follows conservative BP-pipeline behavior:

- coordinates or rsID aliases are resolved against one immutable reference;
- variants are matched as `chromosome:position:REF:ALT`;
- beta, Z, and EAF are flipped when the effect allele is reversed;
- unambiguous strand complements are accepted;
- palindromic A/T and C/G matches fail closed rather than treating population
  frequency as strand proof;
- unmatched, incompatible, ambiguous, duplicated, zero-mapped, and
  multi-mapped rows are dropped by default; and
- the reference identity, checksum, and every QC count are written to the
  store manifest.

The resulting Pcodec contract is always:

```text
identity: chromosome:position:REF:ALT
effect allele: ALT
other allele: REF
beta/Z/EAF allele: ALT
```

This common orientation is what makes sparse keys directly comparable across
independently created stores.

### Current Pcodec variant scope

The default Pcodec backend currently stores biallelic A/C/G/T SNVs on GRCh38
primary chromosomes 1–22, X, and Y. Indels, MNVs, multiallelic records,
mitochondrial variants, alternate contigs, and unresolved identities are
dropped and counted, or cause failure with `strict = TRUE`.

This is an intentional v0.3 scope boundary, not an implication that those
variants are biologically unimportant. Use an exact Parquet store when they
must be retained.

## The canonical ingestion reference

The ingestion reference and the compressed study format solve different
problems:

- the reference is used once during QC to establish GRCh38 identity and
  orientation;
- each `.cpr` store then carries its own exact SNV identity and does not need
  that reference to decode, query, compare, or run MR.

CompreSSoR never manufactures a canonical reference from the study being
converted.

### Build the EBI/Ensembl reference once

The current standard reference is EBI/Ensembl release 95 GRCh38. Download its
per-chromosome VCFs separately, then build the bounded-memory dictionary:

```r
library(CompreSSoR)

build_ebi_reference(
  output = "/data/reference/ebi-ensembl95-grch38",
  source = "/data/reference/ebi-release95-vcfs",
  chromosomes = c(as.character(1:22), "X", "Y", "MT"),
  overwrite = FALSE
)
```

Configure it for subsequent sessions:

```bash
export COMPRESSOR_CANONICAL_REFERENCE=/data/reference/ebi-ensembl95-grch38
```

or inside R:

```r
Sys.setenv(
  COMPRESSOR_CANONICAL_REFERENCE = "/data/reference/ebi-ensembl95-grch38"
)
resolve_reference("GRCh38")
```

The dictionary is checksum-verified before use. It contains one stable global
variant index and is shared by every conversion on that installation. The
default descriptor has no implicit download URL: reference distribution is
kept explicit so a missing or changed upstream asset cannot silently alter a
study conversion.

For another authoritative source, `build_canonical_reference()` accepts VCF,
PLINK BIM/PVAR, delimited, Parquet, or an R data frame and writes an immutable
GRCh38 Parquet reference.

## Reading compressed GWAS files

### Open and inspect

```r
store <- open_compressor("gwas.cpr")
print(store)
store$manifest$n_rows
store$manifest$harmonisation
```

An open `compressor_store` object reuses parsed metadata in hot loops.

### Whole-genome read

```r
all_core <- read_sumstats(
  store,
  columns = c(
    "chromosome", "base_pair_location", "effect_allele", "other_allele",
    "z", "standard_error", "effect_allele_frequency"
  )
)
```

`decompress_sumstats()` is an alias for the same decoded read interface:

```r
all_columns <- decompress_sumstats(store)
```

### Column projection

Independent physical streams mean unrelated values do not need to be opened:

```r
z_only <- read_sumstats(store, columns = "z")
beta_se <- read_sumstats(store, columns = c("beta", "standard_error"))
```

Beta is derived only when requested. The same is true of p-value.

### Regional read

```r
region <- read_sumstats(
  store,
  region = "chr6:25000000-35000000",
  columns = c("chromosome", "base_pair_location", "beta", "standard_error")
)
```

Page-level genomic bounds restrict decoding to candidate pages; a region does
not scan the full study.

### Sparse canonical-key read

```r
keys <- c(
  "1:109274968:A:G",
  "6:160589086:C:T",
  "X:155699751:C:T"
)

hits <- read_sumstats(
  store,
  variants = keys,
  columns = c(
    "chromosome", "base_pair_location", "effect_allele", "other_allele",
    "beta", "standard_error", "effect_allele_frequency"
  )
)
```

Missing keys return no row. Duplicate requested keys return one stored row.
Canonical keys can be constructed safely with:

```r
compressor_variant_key(
  chromosome = "1",
  position = 109274968,
  reference_allele = "A",
  alternate_allele = "G"
)
```

Numeric `variants` are zero-based stored row IDs. Canonical keys are preferred
for portable analysis code because row numbers differ between studies.

### Batch reads across GWAS files

```r
panels <- read_sumstats_batch(
  stores = c(
    exposure = "exposure.cpr",
    outcome = "outcome.cpr"
  ),
  variants = keys,
  columns = c(
    "chromosome", "base_pair_location", "effect_allele", "other_allele",
    "beta", "standard_error"
  ),
  threads = 2
)
```

The batch API starts the codec runtime once, reuses validated indexes and page
caches, and can coalesce requests with exactly the same normalized path, key
vector, and projection. Different stores decode concurrently. For a deliberately
strict reload benchmark, `options(CompreSSoR.coalesce_batch_reads = FALSE,
CompreSSoR.reload_batch_contexts = TRUE)` disables that reuse; ordinary users
normally want the defaults.

## Exactly what is stored

Rows are strictly sorted by a lossless identity key. Bases use `A=0`, `C=1`,
`G=2`, and `T=3`.

```text
key = (global_GRCh38_primary_position << 4) | (REF << 2) | ALT
```

The logical and physical representations are:

| Logical field | Stored representation | Read-time behavior |
|---|---|---|
| chromosome + position | GRCh38 global position, Pcodec `uint32` | Exact |
| REF + ALT | Four-bit complete substitution in Pcodec `uint8` | Exact |
| Z | Z9 semantic code in `uint16` | Central bins plus float32 tail exceptions |
| EAF | EAF8 arcsine-square-root code in `uint8` | Bounded frequency error |
| SE | SE6 block-centred `log2(SE)` residual in `uint8` | Positive reconstructed SE plus float32 exceptions |
| beta | Not stored | `Z × SE` |
| p-value | Not stored | `erfc(abs(Z) / sqrt(2))` |
| rsID / textual variant ID | Not stored | Use the exact REF/ALT key |
| arbitrary extra columns | Not in Pcodec | Use Parquet with `keep_extras = TRUE` |

A `.cpr` store is a directory, not one opaque file:

```text
gwas.cpr/
├── manifest.json
├── manifest.sha256
├── position.pco
├── position.index
├── substitution.pco
├── substitution.index
├── z.pco
├── z.index
├── eaf.pco
├── eaf.index
├── se.pco
├── se.index
└── exceptions.zst
```

Every stream uses 262,144-row Pcodec chunks split into independently decodable
pages of at most 4,096 rows. SE centres use separate 65,536-row blocks. Z, EAF,
SE, position, and substitution are independent streams. Sparse float32
exceptions are held in one Zstandard sidecar.

The full binary contract is documented in
[`inst/doc/pcodec-format.md`](inst/doc/pcodec-format.md).

## Accuracy and round-trip guarantees

### Exact components

- chromosome, position, REF, and ALT;
- row ordering under the canonical identity key;
- missing-key behavior;
- page and file integrity metadata; and
- v0.2 store readability under the v0.3 reader.

### Standard-profile numeric precision

| Stream | Encoding | Declared precision |
|---|---|---|
| Z9 | 510 equal central bins over `[-3.5, 3.5)` | Maximum central absolute error ≈ 0.006863 |
| EAF8 | `asin(sqrt(EAF))` quantisation over 256 codes | Absolute error bound 0.004 |
| SE6 | Block-centred residual of `log2(SE)` conditioned on decoded EAF | Central multiplicative half-step ≈ 1.124% |
| Exceptions | float32 | Finite tail/missing exceptions are not bitwise R doubles |

In the checked 10,000-row real-data audit:

| Audit result | Observed maximum |
|---|---:|
| Identity mismatch | 0 |
| Absolute Z error | 0.0069 |
| Absolute EAF error | 0.0031 |
| Relative SE error | 0.0113 |
| Absolute derived-beta error | 0.0197 |

P-value is exactly derived from decoded Z, not copied from the source. An
input p-value that disagrees with beta/SE/Z is therefore not independently
preserved.

### If exact numbers are required

```r
exact <- compress_sumstats(
  aligned,
  "gwas-exact.cpr",
  backend = "parquet",
  profile = "exact",
  mode = "convert",
  reference = NULL,
  assume_grch38_ref_alt = TRUE
)
```

Do not describe a standard Pcodec store as a bitwise-exact numerical round
trip. It is an explicitly bounded analysis representation.

## Integrity and validation

Stores are written to a staging directory and atomically moved into place only
after the complete write succeeds. A failed `overwrite = TRUE` conversion
therefore leaves the previous completed store intact.

The durable format protects:

- `manifest.json` with a detached SHA-256;
- every durable file with size and SHA-256 in the manifest;
- every Pcodec header, chunk-metadata record, and page with CRC32 in a
  SHA-verified binary index;
- exception flags, ordering, row bounds, and sentinel relationships; and
- position ordering, chromosome bounds, and the complete four-bit
  substitution contract.

Fast validation checks file integrity and index structure:

```r
validate_compressor("gwas.cpr")
```

Full validation additionally decodes and semantically checks every page:

```r
validate_compressor("gwas.cpr", full = TRUE)
```

The 14.9-million-row release store passed full validation across 18,220 pages.

Ordinary reads use a cheaper fail-closed path: the manifest and indexes are
SHA-verified, the exception sidecar is SHA-verified when needed, and each
accessed Pcodec payload is checked by CRC32. They do not hash unrelated stream
payloads on every sparse request. Run `validate_compressor()` when complete
store-level SHA verification is required.

## Benchmarks

### Why Pcodec was selected

The historical representation sweep below reports medians of five full reads
of a real 10-million-row GWAS. Lower read time and higher compression are
better. The dashed line is the measured Pareto frontier.

![Measured speed versus compression Pareto frontier](inst/figures/compressor-pareto.svg)

That experiment motivated independently quantised numeric streams. Subsequent
identity, page-geometry, Pcodec, Zarr, Parquet, Blosc2, byte-shuffle, delta,
and exception-layout experiments selected the current wrapped Pcodec store.
The historical plot is retained as design evidence; the current exported-API
measurements below are the release evidence.

### Release benchmark setup

- machine: Apple Silicon Mac mini;
- data: real FinnGen SNP GWAS;
- rows: 14,923,434;
- source fields: chromosome, position, ALT, REF, beta, SE, EAF, p;
- store: self-contained `0.3.0-pcodec`, 4,096-row pages;
- reader profile: five complete exported-R-API runs with the operating-system
  cache warm, used to isolate decoding and R materialisation;
- end-to-end suite: five randomized fresh-R-process repetitions of 1x1, 5x5,
  25x25 MR and complete loads, with same-request coalescing and the Pcodec page
  cache disabled;
- cold-state control: every logical study in every format is a distinct fresh
  `.noindex` scratch copy written through `F_NOCACHE`; Pcodec read descriptors
  also use `F_NOCACHE` and fresh reader contexts;
- MR panel: 25 real FinnGen index variants selected at p < 1e-5 and clumped
  against GRCh38 1000 Genomes EUR at r2 < 0.001 within 10 Mb;
- Pcodec direct I/O: 20 workers, selected by a fresh-copy 12/16/18/20-worker
  five-run plateau sweep;
- system isolation: the user-owned macOS media-analysis daemon was paused for
  timed trials and resumed afterwards; no root service or Spotlight setting
  was changed;
- durable data: external SSD;
- temporary binary bridges: local temporary storage; and
- validation: all 18,220 pages plus deterministic cross-version parity.

### Storage and write time

| Measure | Result |
|---|---:|
| Source TSV.gz | 201,658,018 B |
| Self-contained CompreSSoR | 58,033,297 B |
| Compression ratio | 3.47× |
| Compression time | 62.982 s |
| Indexed VCF.gz + `.tbi` comparison | 228,634,485 B |

The compressed identity is included in the CompreSSoR size. No external
variant spine is being hidden from this comparison.

### Warm reader profile

| Exported API workload | Median | Range |
|---|---:|---:|
| R: full Z/beta/SE/EAF | 0.162 s | 0.161–0.330 s |
| R: full identity + Z/beta/SE/EAF | 0.463 s | 0.462–0.673 s |
| Python/NumPy: full identity + Z/beta/SE/EAF | 1.539 s | 1.535–1.624 s |
| Python: compact identity/numeric bridge only | 0.141 s | 0.140–0.142 s |

The Python bridge row is deliberately separated: it stops at packed identity
and semantic codes and is not an analysis-ready table. For comparable fully
decoded data, the public R API is 3.3x faster than the NumPy path because the
native C++ bridge reconstructs directly into final R vectors. The analysis-only
projection avoids opening identity streams and materialises a 478 MB R object;
the full public projection materialises 895 MB.

Earlier sparse-access stress testing used 100 random 25-key requests and 50
random 1-Mb regions. Both warm medians were 0.003 s and both p95 values were
0.004 s; these are latency microbenchmarks, not substitutes for the
cache-controlled MR suite.

### End-to-end sparse MR access

The controlled benchmark reuses the association values from one FinnGen GWAS,
but stages each logical exposure/outcome as a distinct file or store. This
models different studies without allowing same-inode cache warmth or identical
request coalescing to flatter a format.

| Workload | Exposure reads | Outcome reads | IVW estimates | Variants per read |
|---|---:|---:|---:|---:|
| 1x1 MR | 1 | 1 | 1 | 25 |
| 5x5 MR | 5 | 5 | 25 | 25 |
| 25x25 MR | 25 | 25 | 625 | 25 |
| Complete load | 0 | 0 | 0 | all 14,923,434 |

Staging is outside the timer; reading, key matching, R materialisation, and the
common FastMR IVW estimator are inside it. The protocol is a symmetric
cache-controlled cold approximation: it is deliberately stronger than an
ordinary warm production loop, but does not claim the operating system has no
other state.

| Input path | 1x1 MR | 5x5 MR | 25x25 MR |
|---|---:|---:|---:|
| **CompreSSoR + FastMR direct** | **0.264 s** (0.255–0.309) | **0.479 s** (0.461–0.626) | **1.274 s** (1.249–1.378) |
| CompreSSoR, explicit sequential reads | 0.387 s (0.379–0.389) | 1.581 s (1.535–1.612) | 7.624 s (7.535–7.653) |
| VCF.gz + Tabix | 0.085 s (0.084–0.097) | 0.330 s (0.326–0.332) | 1.537 s (1.518–1.544) |
| TSV.gz full scans | 6.747 s (6.692–6.777) | 32.830 s (32.412–33.883) | 165.266 s (163.466–173.130) |

Values are five-run median seconds with the complete range in parentheses.
Direct Pcodec is 25.6x, 68.5x, and 129.7x faster than TSV.gz as the grid grows.
Tabix is still the right answer for one or a handful of sparse queries: it is
3.1x faster at 1x1 and 1.45x faster at 5x5. At 25x25, direct Pcodec crosses
over and is 1.21x faster than Tabix. The direct path is 6.0x faster than
explicit sequential Pcodec at 25x25.

### Cache-controlled complete loads

Every returned vector was traversed to compute an identity/numeric checksum;
the total therefore includes more than merely opening a handle.

| Format | Median total | Range | Median read/materialise | Median peak RSS |
|---|---:|---:|---:|---:|
| **CompreSSoR Pcodec** | **3.182 s** | 3.139–3.200 | **0.845 s** | 2.17 GB |
| TSV.gz | 3.954 s | 3.938–4.036 | 1.690 s | 2.09 GB |
| VCF.gz | 10.163 s | 9.834–10.356 | 7.623 s | 2.01 GB |

Pcodec is 1.24x faster than TSV.gz and 3.19x faster than VCF.gz for the full
load-and-traverse endpoint. Its full-load peak memory is about 4% above TSV
because the compact bridge and final R vectors briefly coexist.

The 12/16/18/20-worker 25x25 plateau medians were 1.265, 1.270, 1.266, and
1.258 s. Twenty workers was nominally fastest, but only 0.55% ahead of 12;
12 workers reduced median peak RSS from 270 MB to 219 MB and is the sensible
lower-memory setting.

All 75 trials passed exact cross-format identity/full-load missingness checks,
exact TSV/VCF numeric checks, exact direct/explicit Pcodec MR checks, and the
declared lossy-profile comparison against raw beta/SE/Z. The observed panel
maxima were 0.00827 for absolute beta error, 0.00186 for absolute SE error, and
2.24e-7 for absolute Z error.

### Machine-readable evidence

- [authoritative cold-suite summary](inst/benchmarks/cold-mr-final-summary.csv)
- [all 75 randomized repetitions](inst/benchmarks/cold-mr-final-runs.csv)
- [protocol, revisions, parity gates, and source sizes](inst/benchmarks/cold-mr-final-metadata.json)
- [current Pcodec I/O-worker plateau](inst/benchmarks/pcodec-io-thread-sweep.csv)
- [R/Python reader profile summary](inst/benchmarks/reader-runtime-profile.csv)
- [R reader component runs](inst/benchmarks/reader-profile-r.json)
- [Python reader component runs](inst/benchmarks/reader-profile-python.json)
- [100-request sparse and 50-region warm stress test](inst/benchmarks/pcodec-v03-stress.json)
- [real-data numeric audit](inst/benchmarks/pcodec-full-api-roundtrip.json)
- [historical v0.2 full API benchmark](inst/benchmarks/pcodec-full-api-benchmark.json)
- `benchmark_table("pcodec_access")` for the shipped current access table
- `benchmark_table()` for deliberately retained historical design tables

## FastMR integration

[`fastMR`](https://github.com/gushamilton/fastMR) can run MR directly from
canonical-key CompreSSoR stores without loading complete studies:

```r
library(fastMR)

instruments <- list(
  bmi = c("1:109274968:A:G", "6:160589086:C:T"),
  crp = c("1:154426970:C:T", "2:203867921:G:A")
)

result <- fast_mr_compressed(
  exposure_files = c(
    bmi = "bmi.cpr",
    crp = "crp.cpr"
  ),
  outcome_files = c(
    cad = "cad.cpr",
    stroke = "stroke.cpr"
  ),
  instruments = instruments,
  methods = "ivw",
  threads = 4,
  io_threads = 2
)
```

Each exposure is read for its own instruments; each outcome is read for their
union. The common REF/ALT key and ALT-effect contract remove the need for rsID
lookup or another allele-orientation pass. Instrument discovery, p-value
selection, LD clumping, and biological harmonisation decisions remain explicit
upstream analysis steps.

See the [fastMR compressed-input documentation](https://github.com/gushamilton/fastMR/tree/agent/compressed-gwas-mr#direct-mr-from-compressed-gwas-files)
for the full interface and benchmark.

## Optional exact Parquet stores

Pcodec is the default because it occupies the preferred compression/access
region for the core GWAS representation. Parquet remains useful when a project
needs exact doubles, arbitrary extra columns, non-SNV identity, or direct
Arrow/DuckDB/Polars interoperability.

```r
portable <- compress_sumstats(
  aligned,
  "gwas-parquet.cpr",
  backend = "parquet",
  profile = "exact",
  keep_extras = TRUE,
  reference = NULL,
  mode = "convert",
  assume_grch38_ref_alt = TRUE
)
```

The legacy q8 serving cache applies only to Parquet stores:

```r
build_cache(portable)
```

Pcodec stores are already independently paged; `build_cache()` rejects them.

## Public API

| Function | Purpose |
|---|---|
| `import_sumstats()` | Parse and normalise heterogeneous GWAS input |
| `liftover_sumstats()` | Lift GRCh37 coordinates and alleles to GRCh38 |
| `harmonise_sumstats()` | Match the canonical reference and align effects |
| `compress_sumstats()` | Run the pipeline and write a `.cpr` store |
| `open_compressor()` | Open and inspect a store manifest |
| `read_sumstats()` | Whole, projected, regional, or sparse read |
| `read_sumstats_batch()` | Reuse one codec process across several GWAS files |
| `decompress_sumstats()` | Decode through the same read interface |
| `compressor_variant_key()` | Construct validated canonical REF/ALT keys |
| `validate_compressor()` | Verify file/index integrity and semantic correctness |
| `build_ebi_reference()` | Build the standard partitioned EBI dictionary |
| `build_canonical_reference()` | Build a custom immutable GRCh38 dictionary |
| `resolve_reference()` | Resolve and checksum a configured reference |
| `benchmark_table()` | Read shipped measured benchmark tables |

Use `help(package = "CompreSSoR")` or `?compress_sumstats` for complete
arguments and return values.

## Troubleshooting and FAQ

### Does every `.cpr` need the large canonical reference beside it?

No. The canonical reference is an ingestion dependency, analogous to the
reference used to establish a CRAM-compatible coordinate/allele contract. A
Pcodec `.cpr` carries exact position, REF, and ALT identity and is independently
decodable after conversion.

### Why is there no rsID column?

rsIDs are aliases rather than unique variant identity. They can be absent,
merged, or shared by multiple alleles. `chromosome:position:REF:ALT` is the
stable cross-study key. rsID can be joined from a reference when needed for
presentation.

### Why are beta and p absent from the stored payload?

They are deterministic from the stored representation:

```text
beta = Z × SE
p = erfc(abs(Z) / sqrt(2))
```

Omitting them removes repeated information and guarantees that decoded beta
and p agree with decoded Z and SE.

### Can I preserve N, INFO, phenotype-specific fields, or arbitrary columns?

Yes, with `backend = "parquet", keep_extras = TRUE`. The v0.3 Pcodec core
deliberately stores only exact identity plus Z/SE/EAF. Additional compact
sidecars are a future extension and are not silently created today.

### Can I preserve indels?

Use the Parquet backend. Pcodec v0.3 is deliberately restricted to biallelic
primary-chromosome SNVs.

### Can I preserve the original p-value exactly?

Not in the standard Pcodec representation. P is derived from decoded Z. Keep
the original field in an exact Parquet extras store if it is analytically
important as an independent value.

### Why was a palindromic variant dropped despite having EAF?

Population EAF is not reliable strand proof across ancestries and cohorts.
CompreSSoR therefore drops palindromic A/T and C/G matches rather than risking
a silent sign inversion.

### How do I move temporary bridge files?

Use an R option or environment variable pointing to writable local storage:

```r
options(CompreSSoR.tempdir = "/fast/local/tmp")
```

```bash
export COMPRESSOR_TMPDIR=/fast/local/tmp
```

The bridge is removed after every call and is not part of the `.cpr` format.

### Should I increase `chrom_threads` or batch `threads`?

`chrom_threads` parallelises chromosome-level harmonisation. Batch `threads`
can decode different stores concurrently. Same-store requests are serialized
for decoder safety, and very small sparse reads are often latency-bound, so
more threads are not automatically faster.

### How do I diagnose a store that will not open?

Run:

```r
validate_compressor("gwas.cpr", full = TRUE)
```

Do not bypass checksum or version failures. Readers fail closed on unknown
format versions, altered headers/pages, malformed indexes, invalid semantic
codes, and impossible identity values.

## Development and reproducibility

Repository layout:

```text
R/                  exported R pipeline and reader APIs
src/                compiled semantic reconstruction and bridge reader
inst/python/        Pcodec writer/reader and pinned Python runtime
inst/doc/           durable binary format specification
inst/extdata/       tiny installed examples
inst/benchmarks/    compact, versioned benchmark records
scripts/            reproducible development and benchmark programs
tests/testthat/     R unit, integration, and adversarial tests
```

Large source GWAS files, completed benchmark stores, and temporary decoded
bridges are intentionally kept outside Git. Only compact CSV/JSON evidence and
selected figures are versioned.

Release checks:

```bash
python3 scripts/test_pcodec_backend.py
R CMD build .
R CMD check CompreSSoR_*.tar.gz --as-cran --no-manual
```

The v0.3 release gate includes:

- 20 Python backend corruption and edge cases;
- complete R package tests;
- clean source-tarball installation and checking;
- full validation of the 14.9-million-row FinnGen store;
- exact old-v0.2/new-v0.3 identity and decoded-value parity on a deterministic
  cross-genome sample;
- 100 random sparse requests and 50 random regions; and
- five independent adversarial reviews covering format integrity, ingestion,
  performance, release engineering, and FastMR integration.

Current writers emit `0.3.0-pcodec`. Current readers retain support for
`0.2.0-pcodec`; they fail closed on unknown versions.

## Project status

CompreSSoR is an active pre-1.0 research software project. The core v0.3 format,
reader, harmonisation pipeline, benchmark evidence, and FastMR integration are
implemented and tested. Before treating it as long-term archival infrastructure,
pin the package version, Python requirements, reference manifest, and store
checksums in the analysis environment.

Issues and reproducible examples are welcome at
[github.com/gushamilton/CompreSSoR/issues](https://github.com/gushamilton/CompreSSoR/issues).

## License

CompreSSoR is released under the [MIT License](LICENSE).
