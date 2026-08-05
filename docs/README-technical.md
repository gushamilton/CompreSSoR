# CompreSSoR technical guide

This is the longer companion to the main README. It records the design
decisions, benchmark protocol, limitations, and answers to the questions that
are useful when building or reviewing CompreSSoR.

## Design in one page

CompreSSoR separates three concerns:

1. **Preparation:** an external workflow parses heterogeneous sumstats, applies
   any harmonisation/liftover/reference work, and emits the explicit canonical
   columns required by the compressor.
2. **Storage:** encode one self-contained identity key plus independent core
   numerical streams in a block-framed `.cpr` directory.
3. **Serving:** decode only requested regions, rows, columns, or canonical keys;
   `fastMR` consumes the selected beta/SE rows directly.

The canonical identity is a GRCh38 primary-chromosome global position plus the
complete REF→ALT substitution code. It is not an rsID, and it does not depend
on a study-specific shared spine at read time.

## Why the main store does not contain beta and p

The core numerical contract is:

```text
beta = Z × SE
p    = erfc(abs(Z) / sqrt(2))
```

The standard profile stores Z, EAF, and SE in separate numerical streams:

- Z: 9-bit central quantisation, with sparse float32 exceptions;
- EAF: 8-bit arcsine/square-root quantisation;
- SE: 6-bit semantic block-centred log2 residual in a physical `uint8` stream,
  with sparse exceptions.

This is a bounded-lossy representation. The variant key remains exact, and the
manifest records the numeric profile and tolerances. A Parquet exact store is
available when an analysis needs exact doubles or exact original p-values.
For the standard profile, the manifest records the EAF absolute bound, the
central Z bin bound, the SE relative quantisation bound, and the derived-beta
error formula; the exact Parquet backend is the lossless route.

## Ingestion and GRCh38

The reference is important during ingestion, not because every reader needs to
carry a FASTA or a large variant spine. A normal QC conversion:

```r
store <- compress_sumstats(
  "gwas.tsv.gz",
  "gwas.cpr",
  qc = "compact",
  input_build = "GRCh38",
  store_build = "GRCh38",
  overwrite = TRUE
)
```

The reference establishes the GRCh38 position, REF, ALT, orientation, and
ambiguity decisions. The completed `.cpr` then stores its own exact key and
manifest identity. Reading does not silently consult a different reference.

`mode = "convert"` is intentionally narrower: it is for data already verified
to be GRCh38 with `other_allele = REF`, `effect_allele = ALT`, and matching
beta/Z/EAF orientation. It requires explicit REF/ALT columns or
`assume_grch38_ref_alt = TRUE`.

GWAS Catalog harmonised files commonly contain `hm_*` columns alongside their
ordinary columns. CompreSSoR reports those fields in a structured orientation
diagnostic, but does not infer REF/ALT identity from the prefix alone. Callers
must verify `hm_code` (codes 1--13 are oriented; 14--18 are not) and explicitly
copy `hm_other_allele` to `other_allele`, `hm_effect_allele` to
`effect_allele`, and the matching harmonised coordinates/statistics before
using `assume_grch38_ref_alt = TRUE`. This keeps a malformed or non-oriented
GWAS Catalog row from bypassing the identity safety guard.

Native stores expose two reproducibility hashes in
`manifest$integrity`. `payload_sha256` is deterministic for identical indexed
payload bytes and excludes the manifest and its provenance fields.
`canonical_sha256` is deterministic for the canonical manifest with the
volatile `created_utc` field removed. The detached `manifest.sha256` remains a
byte-integrity check for the actual manifest and may therefore differ between
runs; a whole-directory byte hash is not promised.

The current Pcodec identity scope is biallelic A/C/G/T SNVs on chromosomes
1–22, X, and Y. Indels, unsupported alleles, unresolved reference matches, and
ambiguous rows are handled by the ingestion/QC contract and are not silently
encoded as SNVs.

## On-disk layout

The standard store is a directory rather than an opaque monolithic file so
that the reader can seek to independent streams and pages:

For historical reference, the archived legacy 0.3 layout was:

```text
gwas.cpr/
├── manifest.json
├── manifest.sha256
├── position.pco + position.index
├── substitution.pco + substitution.index
├── z.pco + z.index
├── eaf.pco + eaf.index
├── se.pco + se.index
└── exceptions.zst
```

New native 0.4 stores instead use:

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

`native.index.json` is the only index file: it contains the byte offset and
length of each independently decodable block in every stream.

Each Pcodec stream is column-specific. Pcodec's own documentation warns that
semantically different sequences should not be concatenated into one stream;
CompreSSoR follows that rule. Pages are independently decodable and carry
checksums; the manifest records file checksums, codec constants, identity
constants, source columns, reference metadata, and harmonisation counts.

The current measured geometry is 131,072-row identity frames, 131,072-row
Pcodec pages, and 65,536-row value frames. The identity-frame change was tested
on real FinnGen chr1 rather than inferred from a synthetic estimate: larger
frames reduced the self-contained store until the gain fell below one percent,
while 1,048,576-row identity frames made a 1 Mb regional read more than twice
as slow as the selected 131,072-row frame.

## Benchmark reconciliation

The current headline comparison uses one contract only: every point carries
variant identity in the file. The `.cpr` store intentionally includes its exact
global-position plus REF→ALT key, so it can be moved and read without an
external spine. Earlier reference-anchored experiments remain archived for
historical context but are not presented as competing formats.

On the BP FinnGen chr1 fixture (1,124,344 SNVs; source TSV.gz 15,186,281
bytes), the measured self-contained `.cpr` store is 4,298,204 bytes, or 3.53×
smaller. The complete same-data comparison is in
`inst/benchmarks/archive/legacy-20260804/pareto-chr1-summary.csv`; its five-run access table, write
timings, validation bounds, and plot frontier are stored alongside it. Every
row in that headline record carries variant identity, so the storage numbers
are directly comparable without a separately counted spine.

The native R-facing decoder precomputes the block-centre factors once per
frame, rather than evaluating `exp2()` for every row. A separate direct-file
C++ prototype was benchmarked on BP but was slower than the existing Pcodec
stream path, so it is not enabled in the package; issue #7 remains open for a
proper end-to-end native reader.

## Optional extra columns

The common use case is to keep the store minimal: identity, Z, EAF, and SE are
enough for the main MR and association-serving workflows. N, INFO, case/control
counts, study labels, QC flags, and other fields are useful in some projects,
but they should not enlarge every routine store by default.

### Optional variant-set membership

`variant_set` is an independent, opt-in membership filter. A panel can be a
Parquet or delimited table, a PLINK `.bim`, or a CompreSSoR Pcodec store. Its
canonical identity is the same self-contained `chrom:pos:REF:ALT` key as the
GWAS store. Named panels resolve from `COMPRESSOR_COMMON_VARIANTS`,
`COMPRESSOR_TAG_VARIANTS`, `COMPRESSOR_CORE_VARIANTS`, or
`COMPRESSOR_HM3_VARIANTS` (or from `COMPRESSOR_VARIANT_SET_DIR` using the
corresponding conventional filename).

The default `variant_set = NULL` retains all variants. With
`mode = "convert"`, the panel is matched to the input's existing REF/ALT
identity and no reference table, harmonisation, or shared spine is introduced.
With the QC path, harmonisation happens first and the same identity filter is
then applied. The selected panel is recorded in alignment statistics; its
bytes are external to the per-GWAS store.

### What works today

For exact arbitrary extras, use the Parquet backend:

```r
store <- compress_sumstats(
  "gwas.tsv.gz",
  "gwas.parquet-store",
  backend = "parquet",
  profile = "exact",
  keep_extras = TRUE
)
```

The extras sidecar is keyed by the canonical row number, so it remains aligned
after the same ingestion and sorting steps. It is interoperable with Arrow,
DuckDB, Python, and other Parquet readers.

### What Pcodec can support

Pcodec is a numerical sequence codec. It supports separate integer and floating
point streams, which is a good fit for numeric extras such as N, INFO, and
sample counts. A robust Pcodec extras layer would add, per extra column:

```text
extras/<name>.pco + extras/<name>.index
extras/metadata.json       type, missingness, dictionary, and profile
```

Numeric columns would use a lossless typed stream plus a missingness bitmap.
Character/factor columns would use a dictionary plus a Pcodec integer-code
stream, with a missing-value code. Each extra would remain a separate stream;
we would not interleave columns or claim that text itself is being numerically
compressed.

That layer is feasible, but it changes the `.cpr` manifest contract and needs
its own round-trip, sparse-read, corruption, and size benchmarks. The released
standard profile therefore rejects `keep_extras = TRUE` for Pcodec rather than
silently falling back to a less obvious layout. This is the next sensible
extension if routine numeric extras become important.

## Native Pcodec backend

New stores use a native backend whenever Rust/Cargo is available during
installation. Pcodec itself is Rust; CompreSSoR builds a small static library
with a narrow C ABI and calls it from the package's C++/R bridge. The native
ABI uses caller-allocated buffers and standalone Pcodec streams, which keeps
the package independent of the incomplete upstream wrapped C API while
retaining Pcodec's numerical codec.

The native format is `0.4.5-pcodec-native`. Identity streams use 131,072-row
frames, Pcodec pages use 131,072 rows, and numeric streams use 65,536-row
frames by default (`block_rows` changes the numeric frame size). The five streams are `uint32` global position,
`uint8` substitution code, `uint16` Z code, `uint8` EAF code, and a physical
`uint8` SE code. The public profile is semantic `Z9/EAF8/SE6`: SE centres are
shared across 65,536 rows, with 62 central bins plus missing and exact-exception
sentinels. Exceptions are a small Zstandard-compressed float32 sidecar with
row, Z, log2(SE), EAF, and flags. The index records the byte offset and row
count of every block, so regional and sparse reads do not decode the complete
file.

The native path is the only write and read backend in the current package.
Cargo is required during installation; if the native library cannot be built,
installation stops with the build error. The historical Python-backed 0.2/0.3
implementation remains under `archive/python-backend/` for reference but is
not installed or called.

The upstream Pcodec project documents its standalone C bindings as incomplete;
the CompreSSoR layer therefore owns the ABI, buffer handling, format version,
and round-trip tests rather than exposing that upstream API directly.

On the Mac mini, the historical native-only smoke benchmark uses a deterministic
one-million-row SNV fixture and the self-contained identity key. Its 0.4.4
store is 3,067,151 bytes versus 53,154,119 bytes for source TSV.gz and takes
4.347 s to write. Five-run warm medians are 0.058 s for all-column full load,
0.008 s for a 10 kb region, and 0.204 s for 25 canonical-key reads. Full
validation passed. This historical smoke store used 8K identity frames and 65K
value frames; the current defaults are documented above. This is an engineering
benchmark, not a claim about
every GWAS or sparse workload; the real-GWAS suite should be regenerated for
the native format before making a new production headline.

## Benchmark interpretation

The final cold suite used:

- a real 14,923,434-variant FinnGen GWAS;
- 1×1, 5×5, and 25×25 exposure/outcome grids plus a full-load traversal;
- five randomized repetitions;
- fresh R processes;
- distinct fresh `.noindex` copies for each logical study and format;
- cache-controlled reads, fresh Pcodec reader contexts, and disabled Pcodec
  request coalescing/page cache;
- staging outside the timed interval; and
- exact identity/missingness/round-trip/MR parity gates.

The medians are:

| Input path | 1×1 MR | 5×5 MR | 25×25 MR | Full load |
|---|---:|---:|---:|---:|
| CompreSSoR + FastMR direct | 0.264 s | 0.479 s | **1.274 s** | **3.182 s** |
| CompreSSoR explicit reads | 0.387 s | 1.581 s | 7.624 s | — |
| VCF.gz + Tabix | **0.085 s** | **0.330 s** | 1.537 s | 10.163 s |
| TSV.gz | 6.747 s | 32.830 s | 165.266 s | 3.954 s |

The result is workload-dependent. Tabix wins tiny sparse access; Pcodec wins
the larger sparse grid and is slightly faster for the full load. Pcodec is not
being advertised as universally fastest for every possible query shape.

## FastMR comparison

The compressed-access suite compares data access and MR execution with TSV.gz
and VCF+Tabix. A separate suite compares the compiled FastMR estimator with
native `TwoSampleMR` 0.7.9 on the same in-memory data. Those are different
questions:

- compressed-access benchmark: how quickly can a real GWAS be served into MR?
- TSMR benchmark: how quickly can the estimator run after the data are already
  in R?

FastMR is tens to hundreds of times faster than native TwoSampleMR on the
repeated/grid workloads tested, while direct Pcodec removes the dominant
whole-file scan cost for larger grids.

## Why the reference is not counted in the store size

The default identity key is self-contained. The large canonical GRCh38
reference used to make ingestion decisions is a separate reusable asset and is
not copied into every study. This is analogous to CRAM: the decoder can use a
known reference when required, but the compressed study payload should not
repeat a giant shared asset. In the current standard Pcodec identity design,
the exact REF and ALT identity needed to distinguish SNVs is carried by the
store itself, so reading does not require the external reference.

## Development checks

On the Mac mini:

```bash
Rscript -e 'testthat::test_local(".")'
R CMD check . --no-manual --as-cran
```

The authoritative large files and temporary benchmark stores are on the
external SSD under `/Volumes/crucial_x9/CompreSSoR-benchmarks/`; only compact
benchmark evidence belongs in the repository.
