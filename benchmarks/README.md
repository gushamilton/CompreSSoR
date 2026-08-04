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

The current package uses `z9/eaf8/se8` because commit `4de9978` deliberately
changed the native SE code from a six-bit semantic domain to the full byte
domain. The code comment gives the reason: six semantic bits caused ordinary
SE values to become exceptions. That change is in the current locked package;
it was not an accidental benchmark switch. However, it is separate from the
old projection/native-C++ storage experiment and should not be presented as
the 5.7× result.

Until the exact 5.7× file layout—including its identity bytes—is recovered,
the 5.7× screenshot is retained as historical evidence only. It is not a
current compression claim.

## Reproduction rule

Do not start another benchmark without first updating this document and
creating a new dated directory. Any new result must use the same 10m FinnGen
source, self-contained identity key, current package commit, and five-run
protocol. The raw source GWAS, `.cpr` stores, and Slurm logs remain external
to the repository.
