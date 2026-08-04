# CompreSSoR 0.4.0

- New stores use the native Pcodec backend when Rust/Cargo is available at
  installation time. The 0.4 block-indexed format keeps the self-contained
  GRCh38 identity key and semantic Z/EAF/SE profile while avoiding the Python
  process for ordinary reads and writes.
- The Python Pcodec backend remains available for older stores, fallback
  installations, and explicit compatibility testing.

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
