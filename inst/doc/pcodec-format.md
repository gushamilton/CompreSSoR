# CompreSSoR wrapped Pcodec format 0.3.0-pcodec

A v0.3 store is a directory containing `manifest.json`, its detached
`manifest.sha256`, five independently paged Pcodec streams, five binary stream
indexes, and one Zstandard exception sidecar. Decode-critical constants are
fixed by the format version. Readers fail closed on unknown versions, malformed
indexes, altered headers/pages, and invalid semantic codes.

The Python implementation continues to read `0.2.0-pcodec` frame stores. This
document specifies the format written by current versions.

## Identity

Rows are strictly sorted by this unsigned integer:

```text
key = (global_GRCh38_primary_position << 4) | (REF << 2) | ALT
```

`global_GRCh38_primary_position` is zero-based across GRCh38 primary
chromosomes 1-22, X, and Y. The exact chromosome lengths and offsets are fixed
in the manifest contract. Bases use `A=0`, `C=1`, `G=2`, and `T=3`; REF and ALT
must differ.

The physical identity columns are:

- `position.pco`: global positions as `uint32`;
- `substitution.pco`: `(REF << 2) | ALT` as `uint8`, strictly constrained to
  `0..15`.

This is a lossless biallelic-SNV identity. Decoding needs no rsID, textual
variant ID, FASTA, variant spine, or other external reference. Current ingest
drops indels, MNVs, multiallelic records, alternate contigs, and mitochondrial
rows, or stops under strict mode.

## Numeric semantics

Columns are independent, so projection does not open unrelated streams.

### Z9 (`z.pco`, `uint16`)

- codes `0..509`: equal-width bins over `[-3.5, 3.5)`;
- code `510`: missing;
- code `511`: float32 exception outside the central interval.

Central values decode at bin midpoints. Maximum central absolute error is
`7 / (2 * 510)`, approximately `0.006863`.

### EAF8 (`eaf.pco`, `uint8`)

For valid `0 <= EAF <= 1`:

```text
code = round(255 * (2 / pi) * asin(sqrt(EAF)))
EAF  = sin(pi * code / (2 * 255))^2
```

Missing EAF uses a float32 exception. The declared standard-profile absolute
error bound is 0.004.

### SE6 (`se.pco`, `uint8`)

The encoder forms:

```text
residual = log2(SE) + 0.5 * log2(2 * EAF * (1 - EAF))
```

It subtracts the median residual in each 65,536-row SE-centre block and
quantises the central `[-1, 1)` interval into codes `0..61`. Code `62` is
missing and code `63` is a float32 `log2(SE)` exception. A valid SE with missing
EAF is an SE exception. Half a residual bin is `1/62`, corresponding to a
central multiplicative half-step of `2^(1/62)-1`, approximately 1.124%.

### Exceptions and derived values

`exceptions.zst` is one Zstandard-compressed sorted array with fields:

```text
row:uint32, z:float32, log2se:float32, eaf:float32, flags:uint8
```

Flag bits 0, 1, and 2 select Z, SE, and EAF respectively; all other bits must
be zero. Exception rows are unique, strictly increasing, and below the declared
row count. Stores are capped at `2^31-1` rows by the current bridge contract.

Beta is reconstructed as `Z * SE`; p-value is
`erfc(abs(Z) / sqrt(2))`. Neither is stored. Float32 exceptions are not
bitwise-exact copies of input R doubles. Use `profile="exact",
backend="parquet"` when double-precision losslessness is required.

## Wrapped Pcodec geometry

Every `.pco` stream uses fixed 262,144-row Pcodec chunks. Each chunk is divided
into independently decodable pages of at most 4,096 rows. Constant pages may
have a zero-byte page payload because their value is represented by chunk
metadata.

Each companion `.index` is little-endian and contains:

1. a fixed header: magic, index version, dtype, flags, row count, chunk/page
   geometry, Pcodec header length and CRC32, and chunk/page counts;
2. one record per chunk: metadata offset, length, CRC32, first page, page count;
3. one record per page: payload offset, length, CRC32, row start/count, owner
   chunk;
4. for position pages only, first and last complete identity keys.

Indexes must describe one contiguous physical payload from the Pcodec header to
EOF. Chunk/page counts, ownership, row coverage, offsets, lengths, and identity
bounds are validated. The position bounds allow canonical-key and region
lookups to identify candidate pages without a whole-file scan.

## Integrity

- `manifest.sha256` authenticates `manifest.json` before metadata are trusted.
- The manifest records byte length and SHA-256 for every durable file.
- Ordinary reads SHA-verify the binary indexes before trusting page offsets,
  and SHA-verify the exception sidecar before applying any exception value.
- Every Pcodec header, chunk-metadata record, and page has an index-protected
  CRC32 checked when that payload is read. Ordinary sparse and projected reads
  therefore authenticate metadata and every byte they actually consume, but
  deliberately do not hash unrelated multi-megabyte stream payloads.

`validate_compressor()` is the complete durable-file integrity operation: it
SHA-verifies every file and checks every index and exception record.
`validate_compressor(full=TRUE)` additionally decodes every page, verifies
strict identity order and GRCh38 bounds, enforces the four-bit substitution
contract and numeric widths, and checks exception flags against Z/SE sentinel
codes.

## Reader path

One persistent Python worker caches manifests, binary indexes, Pcodec chunk
decompressors, and a bounded LRU of decoded pages. Sparse canonical-key and
regional reads decode only intersecting pages. Same-store requests are
serialized because cached Pcodec decompressor objects are not assumed
thread-safe; different stores may decode concurrently.

Python writes requested columns to a temporary little-endian binary bridge.
For full scans, Z/SE/EAF remain compact semantic codes and the package's C++
reader allocates the final R vectors, applies centres/exceptions, and derives
beta/p. Full streams decode into preallocated arrays; sparse rows are grouped by
page so each page is decoded once. The bridge is deleted after each call and is
not part of the durable format.

`read_sumstats_batch()` submits several canonical-key requests through the same
worker. By default, requests with exactly the same normalized store path, key
vector, and projection may share one decoded bridge within that batch. The
release benchmark disables coalescing and page caching and requests a fresh
reader context for every logical study, so its 1x1, 5x5, and 25x25 results do
not rely on same-store request deduplication.

## v0.2 compatibility

`0.2.0-pcodec` stores used JSON key/value frame indexes, recursively escaped
position-gap streams, separate substitution frames, and value-frame-aligned
exception frames. Their physical frame sizes are manifest-declared and remain
validated/readable. Current writers never emit that layout.
