# CompreSSoR agent instructions

## Read first

Before changing the package, read:

- `README.md`
- `docs/pcodec-final-spec.md`
- `benchmarks/README.md`

The native Pcodec specification is locked. Treat the documents above as the
source of truth rather than resurrecting archived experiments.

## Current implementation

- The current writer is native Pcodec format `0.4.5-pcodec-native`.
- The store is self-contained: its identity is the GRCh38 global position plus
  directed REF→ALT substitution code. It does not depend on a shared spine or
  external reference when being read.
- The standard streams are position, substitution, semantic Z9, arcsine EAF8,
  semantic SE6, and Zstandard-compressed exceptions. SE6 is carried in a
  physical `uint8` stream; the physical byte container is not an SE8 semantic
  profile. `beta` and `p` are
  reconstructed on demand.
- The normal ingestion path is import/QC/harmonisation/liftover to GRCh38,
  followed by compression. The current identity path supports biallelic A/C/G/T
  SNVs; indels, ambiguous rows, duplicates, and unsupported alleles are not
  silently treated as valid variants.
- No common-variant or tag-variant filter is applied by default.
- Parquet remains an interoperability backend. The historical Python backend
  is archived and must not become the default again.

## Installation and verification

Rust/Cargo is required to build the native backend. On macOS, install it with
`brew install rust` (or an equivalent Rust toolchain). Arrow is needed for the
Parquet backend and the complete test suite:

```r
install.packages("arrow")
```

The configure script detects the architecture of the running R process. This
matters on Apple Silicon when R is running under Rosetta: do not manually
force Cargo to the host architecture.

Run the focused tests during development:

```sh
Rscript --vanilla -e 'testthat::test_local(".", reporter="summary")'
```

Run the package check before publishing:

```sh
R CMD build --no-manual .
R CMD check --no-manual --as-cran CompreSSoR_*.tar.gz
```

Generated `src/Makevars`, object files, shared libraries, check directories,
and package tarballs are temporary build artifacts and must not be committed.

## Benchmark policy

There is one authoritative benchmark family: the five-run BluePebble
comparison on 10,000,000 real FinnGen SNP rows in
`inst/benchmarks/finngen-10m-bp-20260804/`. Its current headline results are
38.06 MB / 0.451 s for native Pcodec, 69.48 MB / 1.410 s for Parquet q9, and
135.39 MB / 4.955 s for TSV.gz.

Do not rerun, replace, or cite older benchmark families without an explicit
request. Before adding a new benchmark, update `benchmarks/README.md` with its
question, input, access contract, repeats, and storage policy. Keep raw GWAS,
reference files, temporary stores, logs, and scratch data outside the synced
repository. Store only compact summaries, plots, and reproducible scripts.

Run heavy compression, stress, and timing jobs on BluePebble compute nodes or
the Mac mini as appropriate; never run jobs on an HPC login node. Use the
BluePebble skill for BP jobs and keep jobs bounded and reproducible.

## Repository hygiene

- Put reusable R code in `R/`, native code and Rust in `src/`, and maintained
  scripts in `scripts/`.
- Keep current benchmark records under `inst/benchmarks/` and superseded work
  under `inst/benchmarks/archive/legacy-20260804/` or `archive/`.
- Do not move archived results back into the live benchmark paths or use them
  to make current performance claims.
- Preserve the existing clean working tree, use a feature branch and PR for
  changes, and never force-push or reset shared history.
