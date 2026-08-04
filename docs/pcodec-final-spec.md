# Locked native Pcodec specification

Status: current write format for CompreSSoR, locked after the five-run
BluePebble comparison on 10 million FinnGen SNP rows.

The writer emits `0.4.4-pcodec-native`. Older `0.4.x-pcodec-native` stores are
read-only compatibility inputs; they are not current benchmark targets and
must not be used to make new performance or storage claims.

## Contract

`compress_sumstats()` writes a block-indexed, self-contained Pcodec store with
the standard semantic profile. The normal path is:

```text
sumstats → import/QC/harmonise/liftover → keyed Pcodec store
```

The storage reference is GRCh38. A GRCh37 input therefore needs a
GRCh37-to-GRCh38 chain during ingestion. Once written, the store does not
need a FASTA, a shared spine, or any other external variant reference to be
read or matched.

The default Pcodec store contains only the core GWAS fields and its identity
key. It does not store `rsid`, `p`, or `beta` per row. Extra columns belong in
the Parquet backend and are outside this locked Pcodec format.

## Variant identity

Every row has an unambiguous identity key:

```text
global_position + directed REF→ALT substitution code
```

For GRCh38 primary chromosomes, the zero-based global position is:

```text
chromosome_offset[chromosome] + position - 1
```

The manifest stores the chromosome lengths and offsets. Bases use
`A=0, C=1, G=2, T=3`; the directed substitution code is
`4 * REF_code + ALT_code`. It is stored as a Pcodec `uint8` stream. Positions
are stored as within-block delta-coded `uint32` values, with each block's
first position recorded in the native index.

The identity streams are therefore self-contained. `rsid` is not part of the
key and is not stored.

## Numerical streams

The physical streams are independently framed and compressed with standalone
Pcodec 1.0.3, level 8:

| Logical field | Stored stream | Representation |
|---|---|---|
| Position | `position.pco` | within-key-block delta `uint32` |
| REF→ALT | `substitution.pco` | directed four-bit code in `uint8` |
| Z | `z.pco` | 9-bit semantic code, physically `uint16` |
| EAF | `eaf.pco` | 8-bit arcsine code, physically `uint8` |
| SE | `se.pco` | 8-bit block-centred log2 residual, physically `uint8` |
| Exceptions | `exceptions.bin` | block-local Zstandard level 19 sidecar |

Key frames contain 131,072 rows. Numerical value frames contain 65,536
rows. Pcodec pages use 131,072 rows.

Z uses 510 central bins over `[-3.5, 3.5)`, with reserved missing and exact
exception codes. EAF uses the arcsine transform
`(2/pi) * asin(sqrt(EAF))` over 255 levels. SE uses the residual

```text
log2(SE) + 0.5 * log2(2 * EAF * (1 - EAF))
```

centred by the median of each 65,536-row block. Its central residual range is
`[-4, 4)` with 254 central bins plus missing and exception codes. The current
SE8 choice is intentional: the former six-bit semantic domain produced too
many exceptions.

Exceptions are stored per value block as 17-byte records before Zstandard:
`row:uint32`, `z:float32`, `log2se:float32`, `eaf:float32`, and `flags:uint8`.
Consequently the standard profile is a compact semantic representation, not
a lossless copy of arbitrary input floating-point values.

At read time:

```text
beta = Z * SE
p    = 2 * pnorm(-abs(Z))
```

## Scope and filtering

The Pcodec identity currently supports biallelic A/C/G/T SNVs on chromosomes
1–22, X, and Y. Indels, unsupported alleles, duplicates, unresolved rows,
and ambiguous rows are dropped by the normal ingestion path; `strict=TRUE`
turns these into errors. No common-variant or tag-variant filter is applied by
default. Optional variant panels are an ingestion choice and do not change
the stored format.

## Locked benchmark reference

The authoritative benchmark is the five-run BluePebble comparison in
[`inst/benchmarks/finngen-10m-bp-20260804/`](../inst/benchmarks/finngen-10m-bp-20260804/).
It uses 10,000,000 real, filtered FinnGen SNP rows and requires every format
to store its own identity key. The shared reference is excluded equally.

| Format | Storage | Whole-file read |
|---|---:|---:|
| CompreSSoR native Pcodec, 4 threads | 38.06 MB | 0.451 s |
| Parquet q9 | 69.48 MB | 1.410 s |
| TSV.gz | 135.39 MB | 4.955 s |

The Pcodec store is 3.56× smaller than TSV.gz. The same store read with one
thread takes 0.552 s; four threads take 0.451 s.

The historical 5.7× SE6/native-C++ screenshot is archived separately. Its
hot projection path contained numeric streams but no position/substitution
identity streams, so it is not comparable to this self-contained contract.
