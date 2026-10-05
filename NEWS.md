# CompreSSoR (development)

P-value flag reads:

- `read_pvalue_flag(as = "row_ids")` and the `pvalue_flag` candidate
  selection decode each flag block natively to row IDs (one C++ pass:
  decode, reject values above 1, emit `row_start + i`), instead of building
  the full row-aligned integer vector in R, checking it with a whole-store
  `%in%` and scanning it again with `which()`. Same rows and validation; on
  UKB-PPP stores (~10M rows) the row-ID read drops from ~0.27 s to ~0.04 s,
  and fastMR's per-store candidate read (which calls it before
  `read_candidates_batch()`) from 150 to 54 CPU-s per 300 stores.
  `as = "row_ids"` also checks the row count against the manifest's
  `hit_rows`. The full-vector path (`as = "logical"`,
  `validate_compressor(full = TRUE)`) checks binary values per block.

Request mapping for dependent packages:

- `read_sumstats_batch(..., request_index = TRUE)` adds an integer
  `request_index` column: for each row, the one-based position in that
  store's `variants` element of the request it answers (first position for a
  repeated key or row ID; `NA` for region-only and whole-store reads). Callers
  no longer rebuild and string-match canonical keys to pair rows with their
  requests. The other columns are unchanged.
- `compressor_identity_code(keys, build)` is exported: the exact numeric
  identity code (`global_position * 16 + substitution`, schema
  `compressor_variant_identity_v1`) of each canonical key, `NA` for keys that
  are not biallelic A/C/G/T SNVs within their primary chromosome.
- `compressor_capabilities()` reports `"key_reads_identity_columns"`,
  `"request_index"` and `"identity_code"`.

Faster, bit-identical reads (adversarial review findings 1, 6-9, 11, 15):

- `read_candidates()` and `read_candidates_batch()` fetch candidate rows
  through the native selective reader: only the key and value blocks that
  hold candidate rows are read and decoded, in native threads, instead of
  decoding whole blocks in R.
- One value decoder for every reader. Row-ID, key, region, region+key,
  candidate and batched reads decode with the same compiled decoder as the
  full read (it now takes a per-row SE-centre index), so every reader returns
  values bit-identical to the full read; previously selective reads decoded
  in R and differed in the last bits (up to ~4.5e-14 relative). The package
  is compiled with `-ffp-contract=off` where the compiler supports it, so the
  decoder's arithmetic does not depend on FMA contraction. Full-read values
  on platforms whose compiler contracted by default (clang on arm64) can
  change in the last bit; stores are unchanged.
- The exception sidecar is read natively in one call with one file handle.
- Opening a store reads `manifest.json` once (hashed in memory and parsed
  from the same bytes). Readers use compact cached block matrices instead of
  the parsed index; parsed indexes are kept for the 16 most recently used
  stores only (`options(CompreSSoR.index_cache_size = )`). Batched readers
  open stores, parse indexes and compute identity signatures in parallel
  workers, and decode the shared keys of separate panel groups in parallel.
- Missing-EAF stores: an exception record that exists only because EAF is
  missing (flags = 4) now stores zero in its unread Z and log2(SE) fields,
  so these records compress (about 4 instead of 9 bytes per missing-EAF row
  in tests). Decoded values are unchanged; stores with complete EAF are
  byte-identical to before.
- Writer: payload files are hashed natively in parallel; input already in
  key order skips the radix order and duplicate scan; only the columns used
  afterwards are reordered. Payloads are unchanged.


Integrity hardening (adversarial review findings 2-5, 10, 12, 13, 16):

- `compress_sumstats(qc = "none")` no longer writes rows whose effect allele
  is REF with their beta attached to ALT (a silent sign error). Such rows,
  and any row whose `effect_allele`/`other_allele` are not ALT/REF, are now
  dropped and counted (`identity_safety$counts$orientation_mismatch`),
  matching compact QC in report mode. Neither mode flips beta or EAF.
- Rows past the end of their chromosome are dropped with a warning in every
  QC mode, counted in a new top-level `manifest$dropped_rows` summary, and
  stop the write above 1% of the input (usually a build mismatch;
  `options(CompreSSoR.max_out_of_range_fraction = )` overrides).
- `validate_compressor(full = TRUE)` now verifies every payload file's byte
  count and SHA-256 and the aggregate `payload_sha256`, as documented; the
  native index is checked against its recorded SHA-256 whenever it is parsed.
- Batched reads that share decoded identity between stores verify each
  member's position/substitution streams against its own manifest first
  (cached per path, size and mtime); a corrupt or half-copied store now
  errors instead of being served another store's keys.
- Every forked read (`pcodec_parallel_lapply()`) checks its results
  centrally: worker errors are re-raised with the store they came from, and a
  worker that died (NULL result) is an error. Region candidate reads could
  previously lose rows silently when a key-decode worker died.
- Native code no longer calls `Rf_error()` from inside C++ catch handlers or
  with C++ objects alive (leaked memory on every failure); errors are raised
  after the C++ state is destroyed.
- The commit-timing manifest update after the atomic commit is written to a
  temporary file and renamed; `overwrite = FALSE` is re-checked at commit.
- Thread counts are no longer written to `native.index.json` and are ignored
  by `canonical_sha256`: payload and canonical hashes no longer depend on
  `threads`. Payload streams are unchanged, but new stores' index bytes and
  `payload_sha256` differ from earlier versions; earlier stores still
  validate.
- Regions without a `chr` prefix (`"1:100-200"`) and `c(chr, start, end)`
  character vectors are accepted (the old pattern required `"ch"`);
  non-integer row IDs are an error instead of being truncated; the unused
  native `compressor_read_pcodec_bridge` reader is removed.
  `read_sumstats_batch()` documents that a store with neither variants nor a
  region is read in full.

- Fix silent candidate loss in `read_candidates_batch()` (0.7.0). When a batch
  mixed variant sets and an identity group (stores with the same variant
  panel) had no candidate rows, stage 2 assigned `shared[[i]] <- NULL`, which
  deletes the list slot and shifted every later store's shared key slice
  onto the wrong store. Those stores came back partial or empty (on 2,940
  UKB-PPP stores, FIS1 returned 1 of 367 flagged rows and EFNA4 0 of 246), or
  the read failed with "subscript out of bounds". Batches of stores sharing
  one variant set were not affected. Every store in every batch composition
  now returns exactly the per-store `read_candidates()` result.
- `read_candidates()` and `read_candidates_batch()` now stop instead of
  returning partial data when the decode step loses selected rows (no
  region), or when a `strategy = "pvalue_flag"` read disagrees with the
  flagged-row count recorded in the store manifest. Batched readers also treat
  a missing worker result (a forked worker that died) as an error.
  `compressor_capabilities()` reports `"candidates_batch_rows_checked"`.

- Faster `compress_sumstats()`; the store payload (every payload file and the
  payload hash) is byte-for-byte unchanged, and for file inputs so is the
  manifest apart from timings and build information:
  - structural QC records each reason as the rows it rejects, evaluating a
    full-length mask only when a cheap aggregate test cannot rule the reason
    out, and finds duplicate keys by an ordered adjacent comparison; the QC
    report (counts and first examples) is unchanged;
  - the variant identity QC already encoded is reused instead of encoding it
    a second time, and the statistic re-validation after compact QC (which
    could only pass) is no longer run;
  - EAF coverage is computed once, p-values are derived from Z only for rows
    without a valid supplied p, within-block position gaps are vectorised,
    and normalisation skips columns that are already canonical;
  - the native manifest is built in memory and written once (plus a final
    write of the commit timing) instead of six write/read round trips.
- The provenance hash of a data-frame input is now a per-column XXH64
  content hash of the columns read, combined with SHA-256
  (`source$hash`, `source$hash_algorithm =
  "xxh64_per_column_sha256_combined_v1"`), computed in parallel and
  independent of thread count and platform. It replaces `source$sha256`
  (the SHA-256 of the whole serialized object) for new data-frame stores, so
  their canonical manifest hash differs from earlier versions; payloads are
  unchanged and stores written earlier open and validate as before. File
  inputs keep `source$sha256` of the file bytes.
- New `allele_columns` argument of `compress_sumstats()`, e.g.
  `allele_columns = c(ref = "other_allele", alt = "effect_allele")` for a
  GWAS-SSF table, declares which input columns hold REF and ALT, so an
  in-memory table can be written without adding REF/ALT columns or writing a
  prepared file. The payload equals that of the prepared table.
- `qc = "none"` is documented as the fast path for tables already validated
  upstream.
- Faster canonical-key and batched reads; outputs are unchanged
  (`identical()` to the previous release on the 10M-row FinnGen store and
  20 simulated stores, including error messages for malformed keys). Times are
  medians of 3 fresh-process runs on BluePebble (Cascade Lake), before -> after:
  - canonical keys are parsed natively (strict `chr:pos:REF:ALT`, with the
    `trimws()` and de-duplication the read path already applied); any other
    spelling falls back to the R parser, so validation and messages are
    unchanged. Within a key block the selection is a galloping merge join
    instead of one binary search per row. `read_sumstats()` with 100k keys:
    1.28 -> 0.53 s (1 thread), 0.94 -> 0.33 s (8 threads); 1M keys:
    13.3 -> 4.8 s, with 25-30% lower peak RSS;
  - `read_sumstats_batch()` normalises and parses a shared key list once,
    hashes each distinct request once, reads stores alone in their panel group
    in one pass each in parallel, and no longer decodes the value pass with
    one thread. 20 shared-panel stores x 100k keys: 4.67 -> 2.35 s
    (8 threads); 5 x 100k: 1.78 -> 0.91 s; 1 store x 100k: 1.20 -> 0.48 s;
    20 stores with distinct panels x 100k keys (one warm process): 16.0 ->
    2.2 s;
  - with `threads = NULL`, key and row-ID reads use one thread per touched key
    block (at most four; one when at most two blocks are touched): 1k keys
    0.36 -> 0.07 s in a warm process. Region reads keep one thread (more
    threads did not help).
- Less aggressive default quantisation. New stores (native format
  `0.4.6-pcodec-native`) use the semantic profile `z10/eaf8/se8+xse`: Z has
  1022 central bins over [-3.5, 3.5) (was 510), the SE log2 residual has 254
  bins (was 62), and rows with |Z| >= 3.5, which already carry a float32
  exception record, take SE from that record as well. The central Z error
  halves (max 0.0034, was 0.0069), relative SE error falls from 1.12% to
  0.27%, and instrument rows are exact to float32 in both Z and SE. Measured
  on BluePebble (FinnGen 10M-row GWAS and two simulated 9M-row traits):
  stores are 9.2-9.5% larger (4.27 bytes per variant, was 3.91); write and
  full-read times are unchanged within run-to-run noise; membership flips at
  reconstructed p <= 0.01 drop from 1.2-1.3% to 0.1% of rows (none at 5e-8
  or 1e-5 either way). In the 10 x 10 fastMR study, MR estimates from `.cpr`
  vs the exact TSV values moved by a median of 0.0015 SE for IVW (was 0.0062;
  max 0.009, was 0.120), 0.0019 for MR-Egger (was 0.0166), 0.0032 for the
  weighted median (was 0.0211) and 0.0040 for the weighted mode (was 0.0156).
  The trade-off grid is in the pull request.
- Z and SE bit widths are recorded in full in the manifest (`z_count`,
  `z_missing`, `z_exception`, `se_count`, ..., `exception_se`) and every
  reader (full, selective, candidates, p-value reconstruction, code-domain
  validation) takes them from there. Stores written by earlier releases are
  Z9/EAF8/SE6 and decode exactly as before; a fixture store written by
  a27b32d is checked in `tests/testthat/test-legacy-store.R`. Releases that
  only know Z9/SE6 refuse 0.4.6 stores rather than misread them. The internal
  option `CompreSSoR.native_profile` (for example `"z9/eaf8/se6"`) writes
  another profile, e.g. for byte-compatible legacy stores; the public API is
  unchanged.
- The SE bin-boundary test no longer hard-codes an SE value whose residual is
  just below the range top on one platform's `log2()`; it searches for one,
  which fixes the failure at `test-quantisation-boundaries.R:68` on Linux.
- Faster store writes with a lower memory peak; the stores are byte-for-byte
  unchanged. On the 10M-row FinnGen table on BluePebble (8 threads, default
  settings including the p-value order domain), `compress_sumstats()` takes
  38.7 s instead of 96.8 s (median of 3 fresh-process runs) and peaks at
  4,670 MiB instead of 5,269 MiB RSS. The changes:
  - every frame of a stream, and every exception frame, is compressed in one
    native call on `threads` threads, instead of forking one R worker per
    block (forking a multi-GB R process cost about 20 s per 10M rows);
  - the `data.frame` input fingerprint (the same SHA-256 as
    `digest::digest()`) is computed while streaming the serialization, without
    materialising it, on a second thread when `threads > 1`;
  - structural QC evaluates allele and text checks once per distinct value,
    skips bookkeeping for rejection reasons no row triggers, and looks up
    chromosome lengths and offsets with `match()`;
  - the flag and order p-value domains share one p-value resolution, and
    already-sorted or fully kept input is not copied.
- Fix: a finite Z one ulp below 3.5 (for example 0.0875 / 0.025, which
  prints as 3.5) was written with Z code 510, the missing sentinel, and no
  exception record, so z, beta and p read back as NA. The SE stream had the
  same collision just below the top of its residual range. Central codes are
  now clamped into the central range. The on-disk format is unchanged, but
  stores written before this fix lost these values at write time and must be
  rewritten from the source to recover them (4 of 10M rows in the FinnGen
  benchmark store).
- Fix: region queries are clamped to the chromosome span, so an end past the
  chromosome end (or a start below 1) no longer returns rows from the next
  (or previous) chromosome.
- Faster full reads: identity columns and stream code-domain checks are built
  natively. On the 10M-row FinnGen store on BluePebble, a fresh-process
  `read_sumstats()` of the eight logical columns drops from 3.4 s to 1.2 s
  (8 threads) and numeric-only reads from 1.1 s to 0.6 s; results are
  identical.

# CompreSSoR 0.6.0

- The exact `pvalue_order` domain is now written by default
  (`pvalue_order = TRUE`, `pvalue_order_threshold = 0.01`; `FALSE` omits it).
  It adds about 1.0% to a 10M-row synthetic store (406 kB on 40.4 MB, 120,904
  candidates) with no measurable compress-time change. It is an optional side
  domain: no native format-version or core-stream change, and stores written
  without it (including `pvalue_order = FALSE`) keep identical payload and
  canonical hashes. `read_candidates(strategy = "exact_order")` adds
  source-p membership (threshold must equal the domain threshold); the
  default `auto` membership is still reconstructed-p. The README gains the
  membership/order contract table and measured quantisation bounds
  (#45, #46).

- One-pass candidates: `read_candidates()` returns any `read_sumstats()` column
  plus a canonical `key` column and (with `order = "exact"`) the exact rank,
  decoding only the key/value/rank blocks that contain candidate rows.
  `read_candidates_batch()` decodes key blocks once per same-panel group. The
  selective (`variants`/`region`) readers now reconstruct `p_value` with the
  decoder's erfc path, so p is bit-identical across full, selective and
  candidate reads (previously selective p was `2 * pnorm(-abs(z))`, which
  could differ in the last ulp). z/beta/se of selective reads are unchanged.

- Performance batch: vectorised key and row matching in the native store
  reader; ingest fast paths for position, allele, chromosome and p-value
  coercion; batched reads reuse variant identity resolution across stores that
  share a panel; a native selective block reader and content-keyed
  (manifest sha256, not mtime) store/index caches. `read_sumstats_batch()`
  gains `region` and accepts `variants = NULL` and zero-based row IDs.
  `read_candidates()` gains a `strategy` override and never auto-selects the
  writer-time p-value flag. Adds a BP-ready read/ingest benchmark harness.

- Adds `read_candidates()` and `read_candidates_batch()` for threshold and
  region candidate extraction from native Pcodec stores. The strategy is
  chosen per store and recorded in `attr(x, "candidate_strategy")`: the
  Z-exception sidecar alone when the threshold is below the p-value of the outermost
  central Z bin (every row outside the central range is an exact float32
  exception); otherwise the Z stream only. Candidate rows are fetched with a
  block-selective reader that opens each payload once. `order = "exact"` uses
  the `pvalue_order` domain; `order = "reconstructed"` is labelled approximate. Membership always equals a full read filtered at
  `p <= threshold`; the writer-time p-value flag is used only with an explicit
  `strategy = "pvalue_flag"` (flag membership follows a supplied p-value and
  may differ).

- Adds an opt-in, explicitly versioned `pvalue_order` side domain for exact
  supplied-p/exact-prepared-Z candidate ordering without changing the locked
  native core streams. `read_pvalue_order()` fails safely when exact ordering
  is unavailable and requires an explicit reconstructed-p fallback.
- Tightens source-package hygiene, native build cleanup, macOS Rust deployment
  targeting, benchmark-path portability, and Rust formatting.

# CompreSSoR 0.5.0

- Refactors the installed package around a strict compression-only contract.
  Inputs must already contain chromosome, position, explicit REF/ALT,
  effect/non-effect alleles, beta, and standard error. The core performs
  structural QC and compression only; it does not look up references, flip
  alleles, resolve rsIDs, or liftover coordinates.
- Adds native GRCh37/hg19 stores alongside GRCh38/hg38. Input and store builds
  must match; cross-build conversion is intentionally outside this package.
- Moves the former reference, harmonisation, liftover, and legacy workflow
  implementation into `archive/harmonisation/` for recovery without loading it
  as part of the installed package.

# CompreSSoR 0.4.0

- Adds the final same-data BP FinnGen chr1 Pareto record: nine formats, five
  access repetitions, explicit self-contained/reference-anchored contracts,
  package round-trip validation, and a README plot generated from the checked
  CSV record.
- Records theoretical validation bounds for the standard lossy Z9/EAF8/SE6
  profile in the manifest, including the SE-relative and beta reconstruction
  tolerances.
- Lowers the Arrow compatibility floor to the Arrow 22 API available on the
  supported BluePebble R environment.
- New stores use the native Pcodec backend when Rust/Cargo is available at
  installation time. The 0.4 block-indexed format keeps the self-contained
  GRCh38 identity key and semantic Z/EAF/SE profile while avoiding the Python
  process entirely for ordinary reads and writes.
- The historical Python Pcodec backend is archived outside the installed
  package; Rust/Cargo is now the sole Pcodec build requirement.
- BP FinnGen chr1 benchmarks selected 131,072-row identity frames and
  131,072-row Pcodec pages: 4,297,931 bytes self-contained versus 15,186,281
  bytes for the source TSV.gz.

# CompreSSoR 0.3.0

- Introduces the wrapped Pcodec v0.3 store: fixed 262,144-row chunks with
  independently readable 4,096-row pages for position, substitution, Z9,
  EAF8, and SE6 streams.
- Adds an integrity-protected binary page index, CRC32 for every Pcodec header,
  chunk metadata record and page, per-file SHA-256, strict four-bit
  substitution validation, and full semantic validation.
- Adds a persistent cached reader, bounded decompressed-page LRU, local binary
  bridges, grouped sparse access, and preallocated whole-stream decoding.
- Keeps v0.2 stores readable while writing only v0.3 stores.
- Makes palindromic alleles fail closed, validates liftover source builds and
  GRCh38 primary targets, and serializes same-store threaded reads.
- Adds real 14,923,434-row FinnGen benchmarks, random sparse/region stress
  tests, and direct compressed-input FastMR integration.

# CompreSSoR 0.2.1

- Moves the default identity-frame geometry to the measured 8,192-row Pareto
  point and value frames to 32,768 rows while retaining compatibility with
  larger 0.2 stores.
- Adds `read_sumstats_batch()` so multi-GWAS analyses pay Python startup once,
  cache store indexes, and coalesce identical canonical-key reads.
- Decouples physical value frames from SE centre blocks and validates every
  self-described frame-size field.
- Adds corrected REF/ALT FinnGen, full-frame, five-repeat access, Tabix, and
  direct 5 x 5 FastMR benchmarks from the Mac mini.

# CompreSSoR 0.2.0

- Makes the self-contained Pcodec Z9/EAF8/SE6 format the default backend.
- Adds public import, GRCh38 liftover, conservative BP-style harmonisation,
  compression, decompression, projection, regional/sparse reads, and full
  integrity validation.
- Adds lossless GRCh38 position/REF/ALT identity keys without per-row rsIDs or
  a decode-time reference dependency.
- Adds indexed sparse lookup by canonical `chromosome:position:REF:ALT` key so
  the same variants can be resolved across stores with different row sets.
- Adds block-aligned exception frames, per-frame and per-file SHA-256 checks,
  atomic store replacement, and a compiled full-scan bridge reader.
- Adds adversarial input/codec tests and five-repeat full-FinnGen API
  benchmarks on the Mac mini.

# CompreSSoR 0.1.0

Initial standalone package release candidate.

- Adds explicit `convert`, `qc`, `all`, `core` and `hm3` workflows.
- Preserves unresolved rows by default and records harmonisation status.
- Supports common delimited, gzipped, VCF/VCF.gz and BIM panel inputs.
- Adds exact/lossy/cache round-trip tests and real FinnGen release-gate
  benchmarks run on the Mac mini.
- Adds optional post-harmonisation `mode = "core_plus"`: the configured core
  panel is stored together with variants within a configurable window of
  canonical `p < 1e-5` signals. `mode = "pvalue_regions"` remains available for
  p-value-only selection.
- Adds deterministic chromosome-sharded core/HM3 panel preparation and fast
  identity-only reads, including optional `pigz` decompression through
  `COMPRESSOR_DECOMPRESSOR`.
