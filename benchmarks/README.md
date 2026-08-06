# Benchmark policy and current evidence

The locked format contract is [`docs/pcodec-final-spec.md`](../docs/pcodec-final-spec.md).

This repository now has one authoritative benchmark family. Older benchmark
tables are preserved under [`inst/benchmarks/archive/legacy-20260804/`](../inst/benchmarks/archive/legacy-20260804/)
and must not be used for new claims.

## Authoritative current sets

The current evidence is the five-run BluePebble benchmark on the normalized
10-million-row FinnGen SNP file:

1. `inst/benchmarks/finngen-10m-bp-20260804/format-screen/` — the full
   same-data keyed format comparison;
2. `inst/benchmarks/finngen-10m-bp-20260804/pcodec-thread-compare/` — the
   native Pcodec reader at one and four threads.

The source and filtering contract is recorded in
`finngen_10m_metadata.json`. The source has eight columns:
`chrom`, `pos`, `ref`, `alt`, `beta`, `se`, `eaf`, and `p`. Every candidate in
the format screen stores a self-contained variant identity key. No shared
spine or external reference is counted in the candidate size.

The latest headline values are:

| Format | Size | Whole-file read |
|---|---:|---:|
| CompreSSoR native Pcodec, 4 threads | 38.06 MB | 0.451 s |
| Parquet q9 | 69.48 MB | 1.410 s |
| TSV.gz | 135.39 MB | 4.955 s |

The Pcodec store is 3.56× smaller than TSV.gz and the four-thread reader is
18.3% faster than the same store forced to one thread (0.451 s versus
0.552 s). These are the numbers to use in the README and package discussion.

## The SE6/SE8 confusion

The 5.7× screenshot is not the same stored-file contract as the current
10-million-row keyed benchmark. The archived frontier-projection records are
numeric reconstruction experiments: their inputs are `z.u16`, `eaf.u8`, and
`se.u8` plus exception arrays. They do not contain the position and
REF/ALT identity streams required by a self-contained GWAS store.

The current package uses the canonical `z9/eaf8/se6` semantic profile. SE codes
are stored in a physical `uint8` stream, but only 62 central semantic bins plus
missing and exception sentinels are part of the public contract. The historical
full-byte semantic `se8` experiment is separate from the current keyed store
and must not be presented as the 5.7× result.

Until the exact 5.7× file layout—including its identity bytes—is recovered,
the 5.7× screenshot is retained as historical evidence only. It is not a
current compression claim.

## Reproduction rule

Do not start another benchmark without first updating this document and
creating a new dated directory. Any new result must use the same 10m FinnGen
source, self-contained identity key, current package commit, and five-run
protocol. The raw source GWAS, `.cpr` stores, and Slurm logs remain external
to the repository.

## Core-plus acceptance benchmark

The core-plus workstream uses a separate, non-authoritative benchmark question:
how many rows and bytes does the current native store retain when a prepared
UKB-PPP NTRK3 protein file is reduced to the bundled core panel plus associated
regions? The input contract is the prepared GRCh38 file
`UKB-PPP/.../NTRK3_explicit_core_columns.tsv` on BluePebble, with canonical
identity, beta, standard error, and EAF columns. The selection contract is
explicit and locked for this first test: `p <= 1e-5`, inclusive 10,000-bp
padding on each side, core membership unioned row-wise with the regions, and
native Pcodec `Z9/EAF8/SE6`. The p-value is derived from the exact prepared
pre-encoding Z; the input file is not copied into the repository.

Run the maintained `scripts/benchmark_core_plus_ntrk3.R` harness on a Slurm
compute node with a four-thread native install. Keep the raw input, output
store, and Slurm logs under the external BluePebble benchmark workspace. The
harness writes one compact JSON summary containing source/build/commit
provenance, selected rows, region/seed counts, output bytes, manifest phase
timings, and elapsed time; obtain peak RSS from the Slurm `/usr/bin/time -v`
record. This acceptance result is evidence for the core-plus workstream, not a
replacement for the locked five-run FinnGen headline benchmark above.
