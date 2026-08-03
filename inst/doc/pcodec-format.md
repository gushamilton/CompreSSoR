# CompreSSoR Pcodec format 0.2.0-pcodec

A store is a directory containing `manifest.json`, its detached
`manifest.sha256`, two JSON frame indexes, independent Pcodec streams, and
value-block-aligned Zstandard exception frames. Readers fail closed on unknown
format versions and altered decode-critical constants. Every manifest, index,
payload file, and indexed frame has a SHA-256 checksum.

## Identity stream

Rows are sorted by this unsigned integer key:

```text
key = (global_GRCh38_primary_position << 4) | (REF << 2) | ALT
```

`global_GRCh38_primary_position` is zero-based across GRCh38 primary
chromosomes 1-22, X, and Y. Chromosome lengths and offsets are embedded in the
manifest. Bases use `A=0`, `C=1`, `G=2`, and `T=3`. REF and ALT must differ.

This is a lossless biallelic-SNV identity. It needs no rsID, textual variant ID,
FASTA, or shared study reference to decode. The standard encoder rejects
indels, MNVs, multiallelic records, alternate contigs, and mitochondrial rows.

Key frames contain 131,072 rows. Position gaps `0..510` are uint16 values;
`511` escapes into a second uint16 Pcodec stream. Values above `65,534` use
`65,535` there and a uint32 Pcodec stream. The low four substitution bits are
an independent uint8 Pcodec stream.

## Numeric streams

Value frames contain 65,536 rows. Columns are independent so projection can
avoid unrelated streams.

### Z9

- codes `0..509`: equal-width bins over `[-3.5, 3.5)`;
- code `510`: missing;
- code `511`: float32 exception outside the central interval.

Central values decode at bin midpoints. The maximum central absolute error is
`7 / (2 * 510)`, approximately `0.006863`.

### EAF8

For valid `0 <= EAF <= 1`:

```text
code = round(255 * (2 / pi) * asin(sqrt(EAF)))
EAF  = sin(pi * code / (2 * 255))^2
```

Missing EAF is represented in the float32 exception frames. The declared
standard-profile absolute error bound is 0.004.

### SE6

The encoder forms:

```text
residual = log2(SE) + 0.5 * log2(2 * EAF * (1 - EAF))
```

It subtracts the median residual in each value frame and quantises the central
`[-1, 1)` interval into codes `0..61`. Code `62` means missing and code `63`
means a float32 `log2(SE)` exception. A valid SE with missing EAF is always an
SE exception; missing Z never erases SE. Half a residual bin is `1/62`, giving
a central multiplicative half-step of `2^(1/62)-1`, approximately 1.124%.
Because SE is encoded and decoded against the same decoded EAF value, EAF
quantisation does not add a second error term to central SE reconstruction.

### Derived fields

Beta is `Z * SE`. The two-sided p-value is `erfc(abs(Z) / sqrt(2))`. Neither is
stored. Supplied beta, Z, and SE must agree within the importer tolerance before
encoding. Float32 exceptions are full-range exceptions, not bitwise-exact
copies of input R doubles. Use `profile="exact", backend="parquet"` for
lossless double-precision Z/SE/EAF.

## Indexes and integrity

`manifest.sha256` is verified before any decoding metadata is trusted. The
format version fixes chromosome lengths/offsets, stream names, frame sizes,
code widths, sentinels, and Z range; changing them is an incompatible format,
not a manifest option.

`key.index.json` and `value.index.json` store row ranges, byte offsets, byte
lengths, genomic bounds where relevant, and per-frame SHA-256 values. The
manifest stores the byte size and SHA-256 of every index and payload file.

`validate_compressor()` verifies all file checksums and contiguous indexes.
`validate_compressor(full=TRUE)` also decodes every frame, verifies key order,
code widths, recomputed genomic frame bounds, per-frame exception ownership,
exception rows/flags, and sentinel correspondence.

## R reader path

The Python wrapper decompresses only requested Pcodec streams and selected
frames into a temporary,
columnar binary bridge. For a full scan it keeps Z, SE, and EAF as their compact
uint16/uint8 semantic codes. The package's compiled reader loads those files,
applies block centres and exceptions, reconstructs requested identity and
numeric columns directly as R vectors, and derives beta/p in the same pass.
Regional and sparse reads decode only intersecting frames and use a projected
float bridge because their output is small. Canonical
`chromosome:position:REF:ALT` lookups use the key-frame genomic bounds and
binary search within only candidate frames, then decode the corresponding
value rows. Numeric-only whole scans do not
open identity payloads, and Z-only scans do not open EAF or SE payloads. The
bridge is temporary and is not part of the durable file format; its multibyte
fields are little-endian and the compiled reader converts them on big-endian
hosts.
