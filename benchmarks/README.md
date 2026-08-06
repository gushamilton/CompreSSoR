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
native Pcodec `Z9/EAF8/SE6`. A finite supplied p-value is authoritative; this
NTRK3 input has no p-value column, so its selection falls back to exact
prepared pre-encoding Z. The input file is not copied into the repository.

Run the maintained `scripts/benchmark_core_plus_ntrk3.R` harness on a Slurm
compute node with a four-thread native install. Keep the raw input, output
store, and Slurm logs under the external BluePebble benchmark workspace. The
harness writes one compact JSON summary containing source/build/commit
provenance, selected rows, region/seed counts, output bytes, manifest phase
timings, and elapsed time; obtain peak RSS from the Slurm `/usr/bin/time -v`
record. This acceptance result is evidence for the core-plus workstream, not a
replacement for the locked five-run FinnGen headline benchmark above.

## Issue #27 bounded-memory evidence

This is a separate, non-authoritative acceptance question for the large-file
memory issue: does the prepared UKB-PPP NTRK3 protein complete in the native
Pcodec path without requiring a 32-GiB process? It must not replace the locked
FinnGen benchmark or be quoted as a general scaling law.

The access contract is the read-only BluePebble archive
`/bp1/mrcieu1/data/UKB-PPP/UKB-PPP_pGWAS_summary_statistics/Combined/NTRK3_Q16288_OID21057_v1_Neurology.tar`,
with the externally prepared GRCh38 file
`NTRK3_explicit_core_columns.tsv`. The archive is about 777 MiB and contains
23 compressed chromosome members; the prepared input has 23,745,741 data rows,
1,213,767,088 bytes, and nine canonical columns. Raw input, stores, Slurm
logs, and `/usr/bin/time -v` records stay outside this repository.

The maintained harness is
`scripts/benchmark_issue27_ntrk3.R` plus
`scripts/benchmark_issue27_ntrk3.sbatch`. A benchmark run is one bounded
16-GiB, four-CPU Slurm job per selection/QC mode; the JSON summary records the
source hash and columns, build, commit, requested/effective workers, retained
rows, exception rows, EAF observability, manifest phase timings, output bytes,
and validation status. The Slurm record supplies peak RSS. A native install
requires Cargo/Rust >= 1.87.0; the job fails before reading data if that
toolchain is unavailable rather than silently using a non-native backend.

Commit provenance is resolved before `sbatch`. Use the login-side helper
`scripts/submit_issue27_ntrk3.sh` with `COMPRESSOR_ISSUE27_REPO` and
`COMPRESSOR_ISSUE27_COMMIT` set. It resolves `git -C "$repo" rev-parse HEAD`
on the login node and exports `COMPRESSOR_ISSUE27_ACTUAL_COMMIT` through
Slurm. The compute-side harness uses Git directly when available; when Git is
absent, it requires that pre-resolved value, validates its 40-hex form, and
compares it exactly to the requested commit before creating benchmark
directories or installing the package. Missing, malformed, or mismatched
provenance is a setup failure and is never a benchmark result.

Prior external runs provide the following evidence (the commit and job IDs are
retained here so they are not confused with a fresh exact-base rerun):

| Selection | QC | Source commit | Slurm job | Wall time | MaxRSS | Stored rows | Output |
|---|---|---|---:|---:|---:|---:|---:|
| full | none | `abc7b88` | 18277164 | 11m49s | 11,790,820 KB (~11.2 GiB) | 21,737,357 | 80,339,943 B |
| full | compact | `abc7b88` | 18277165 | 10m53s | 11,430,056 KB (~10.9 GiB) | 21,737,357 | 80,352,487 B |
| core | none | `927535a` | 18277437 | 4m46s | 7,677,716 KB (~7.3 GiB) | 5,578,843 | 22,321,136 B |
| core | compact | `927535a` | 18277438 | 6m25s | 10,261,076 KB (~9.8 GiB) | 5,578,843 | 22,333,731 B |

All four prior runs completed under a 16-GiB allocation, retained complete
EAF coverage, and passed shallow validation. The full `qc="none"` result was
not rerun after the final numeric-QC bypass commit `927535a`; it is therefore
supporting evidence, not final exact-base acceptance evidence. On 2026-08-06,
an exact-`4a0ef5c` rerun was blocked before input read because the current BP
module catalogue exposes Rust 1.78 only, while the native package requires
Rust 1.87.0. No timing or RSS result is claimed for that blocked run.
