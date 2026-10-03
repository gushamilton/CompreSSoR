# CompreSSoR native Pcodec format

CompreSSoR writes `0.4.6-pcodec-native` stores. The package requires
Rust/Cargo at installation time and links a small static Rust Pcodec library
through a narrow C ABI. No Python runtime or external process is used.

Each store is a directory containing:

```text
gwas.cpr/
├── manifest.json + manifest.sha256
├── native.index.json
├── position.pco
├── substitution.pco
├── z.pco
├── eaf.pco
├── se.pco
└── exceptions.bin
```

The streams are independent Pcodec blocks. The standard geometry uses
131,072-row identity frames, 131,072-row Pcodec pages, and 65,536-row value
frames. Position is a GRCh38 global `uint32`; substitution is a four-bit
REF/ALT code; Z is a 10-bit semantic code; EAF is an 8-bit arcsine code; and SE
is an 8-bit semantic block-centred log2 residual carried in a physical `uint8`
stream. The SE domain has 254 central bins plus missing and exception sentinels.
Rows with |Z| >= 3.5 store Z and SE exactly (float32) in their exception
record. Stores written before 0.4.6 use Z9 (510 bins) and SE6 (62 bins) and
keep SE quantised on Z-exception rows; the bit widths and counts are read from
each store's manifest.
Rare out-of-range or missing values
are stored as float32 exceptions. Beta and p-values are derived when read.

The standard profile is bounded-lossy by design: central Z has a 10-bit bin
bound, EAF has a documented absolute bound of 0.004, and SE is encoded as a
semantic SE8 log2 residual with a documented relative bound of `2^(4/254)-1`.
The
manifest records these tolerances so validation can check the representation
against its actual quantisation contract. Use the exact Parquet path when
bit-for-bit preservation of every numeric input is required.

The identity key is self-contained: readers do not need an rsID table, shared
variant spine, or reference file. The installed package accepts prepared input
with an explicit matching build and orientation; any reference lookup,
harmonisation, or liftover needed to prepare a source GWAS belongs to an
external workflow and is not performed by the compressor.

The old Python-backed 0.2/0.3 implementation is retained under
`archive/python-backend/` for historical reproducibility but is not installed,
called, or supported by the current package.
