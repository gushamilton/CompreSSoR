# CompreSSoR native Pcodec format

CompreSSoR writes `0.4.5-pcodec-native` stores. The package requires
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
REF/ALT code; Z is a 9-bit semantic code; EAF is an 8-bit arcsine code; and SE
is a 6-bit semantic block-centred log2 residual carried in a physical `uint8`
stream. The SE domain has 62 central bins plus missing and exception sentinels.
Rare out-of-range or missing values
are stored as float32 exceptions. Beta and p-values are derived when read.

The standard profile is bounded-lossy by design: central Z has a 9-bit bin
bound, EAF has a documented absolute bound of 0.004, and SE is encoded as a
semantic SE6 log2 residual with a documented relative bound of `2^(4/62)-1`.
The
manifest records these tolerances so validation can check the representation
against its actual quantisation contract. Use the exact Parquet path when
bit-for-bit preservation of every numeric input is required.

The identity key is self-contained: readers do not need an rsID table, shared
variant spine, or GRCh38 reference file. The reference is still used by the
ingestion and harmonisation pipeline when a source GWAS needs to be aligned.

The old Python-backed 0.2/0.3 implementation is retained under
`archive/python-backend/` for historical reproducibility but is not installed,
called, or supported by the current package.
