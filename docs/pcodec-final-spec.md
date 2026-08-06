# Locked native Pcodec specification

Status: current write format for CompreSSoR, locked after the five-run
BluePebble comparison on 10 million FinnGen SNP rows.

The writer emits `0.4.5-pcodec-native`. Older `0.4.x-pcodec-native` stores are
read-only compatibility inputs; they are not current benchmark targets and
must not be used to make new performance or storage claims.

## Contract

`compress_sumstats()` writes a block-indexed, self-contained Pcodec store with
the standard semantic profile. The normal prepared-input path is:

```text
prepared sumstats → compact structural QC → keyed Pcodec store
```

Two ingestion modes are part of the current API:

- `qc = "compact"` (default) applies structural QC, duplicate handling, and
  aggregate reporting before encoding.
- `qc = "none"` is an explicit fast path for already-canonical prepared
  input. It skips structural QC, duplicate scans, row-level reports, alias
  resolution, and temporary variant-ID construction. A cheap
  pre-canonicalization identity guard drops rows that cannot be represented by
  the selected native identity: non-primary chromosomes, non-finite or
  non-integer/out-of-range coordinates, and invalid or non-distinct A/C/G/T
  REF/ALT pairs. Malformed canonical statistics are not scanned in this fast
  path; native duplicate checks still apply at write time. The exact canonical
  columns
  remain mandatory. Dropped identity counts and the pre-filter input count are
  recorded in preparation metadata, while selection counts describe the rows
  remaining after that guard. The manifest records the bypass under `qc.mode`,
  `qc.structural_qc`, `qc.statistic_validation`, and `preparation.qc`.

Neither mode performs harmonisation or liftover. Those transformations and
their reference/chain provenance belong to an upstream preparation workflow.

The installed package does not include a position-selective dual-build
reference backend. It does not download or silently substitute dbSNP155 (or
another external table), resolve rsIDs, align alleles, or perform liftover.
An upstream tidyGWAS/dbSNP-style or GWASLab preparation workflow may perform
those operations, but must hand CompreSSoR an explicitly prepared table and
retain its own reference version, build mapping, QC, chain, and timing
provenance. Once written, the native store is self-contained and does not
depend on that upstream reference.

For `selection = "core_plus"`, selection is evaluated before lossy encoding
and retains the core-panel union with inclusive windows around rows at
`p <= pvalue_threshold` (default `1e-5`, with `region_padding = 10000`). A
finite supplied p-value is authoritative, including values from score tests,
mixed models, meta-analysis, or corrected analyses that need not equal a
Wald p-value from beta/SE. When p is absent, the selector derives a two-sided
p-value from the exact prepared Z. A supplied p that is missing or invalid is
never accepted as valid: in compact QC it is rejected before selection and
the malformed/non-finite/out-of-range counts are retained in structural and
selection provenance; in `qc = "none"`, the trusted fast path may fall back
to exact prepared Z and records the fallback and unresolved counts. This is
intentional: compact mode is the safe validation contract, while `qc = "none"`
is an explicitly trusted prepared-input contract. P-values are not stored in
the native payload.

The storage build is explicit. A prepared GRCh37 input may be stored as
GRCh37, and a prepared GRCh38 input may be stored as GRCh38; `input_build` and
`store_build` must match. GRCh37/GRCh38 conversion, liftover, reference lookup,
and allele alignment are external preparation steps, not ingestion features of
the installed compressor. Once written, the store does not need a FASTA, a
shared spine, chain file, or any other external variant reference to be read
or matched.

The default Pcodec store contains only the core GWAS fields and its identity
key. It does not store `rsid`, `p`, or `beta` per row. Extra columns belong in
the Parquet backend and are outside this locked Pcodec format.

The canonical public numerical profile is `Z9/EAF8/SE6`: Z has 9 semantic bits,
EAF has 8 semantic bits, and SE has 6 semantic bits. The SE semantic codes are
carried in a physical `uint8` stream; this byte container must not be confused
with the historical SE8 semantic experiment. Current manifests identify this
distinction with `se_bits = 6`, `se_count = 62`, `se_missing = 62`,
`se_exception = 63`, and `se_physical_dtype = "uint8"`.

The input statistics remain independent: a supplied finite standard error is
the authoritative SE value and is never replaced by a value derived from EAF.
EAF is used both for its own stream and as the predictor in the SE residual
transform. If EAF is missing, the writer uses the quantised decoded fallback
for an internal SE predictor only; it does not write an imputed biological EAF.
The missing EAF is restored on read through an EAF exception and the manifest
records the predictor fallback and affected-row count. Exception flags are
bitwise: `1 = Z`, `2 = SE`, and `4 = EAF`; a missing-EAF row with valid SE
therefore needs only the EAF flag when its SE residual is in range.
For the current EAF8 domain, the `0.5` seed is code `128`, whose decoded
predictor is approximately `0.50307997331906917`; the manifest records this
code and decoded value rather than treating the seed as the stored EAF.

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
| SE | `se.pco` | 6-bit semantic block-centred log2 residual (62 bins plus two sentinels), physically `uint8` |
| Exceptions | `exceptions.bin` | block-local Zstandard level 19 sidecar |

Key frames contain 131,072 rows. Numerical value frames contain 65,536
rows. Pcodec pages use 131,072 rows.

Z uses 510 central bins over `[-3.5, 3.5)`, with reserved missing and exact
exception codes. EAF uses the arcsine transform
`(2/pi) * asin(sqrt(EAF))` over 255 levels. SE uses the residual

```text
log2(SE) + 0.5 * log2(2 * EAF * (1 - EAF))
```

centred by the median of each 65,536-row block. Its standard public residual
range is `[-1, 1)` with 62 central bins plus missing and exception codes. The
former SE8 experiment is historical and is not the public profile.

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

The historical SE8/native-C++ screenshot is archived separately. Its
hot projection path contained numeric streams but no position/substitution
identity streams, so it is not comparable to this self-contained contract.
