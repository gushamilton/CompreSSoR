# Issue #18 roadmap disposition

Audit base: `origin/main` at `f4e783f` (2026-08-06).

Issue #18 was written for a package with two public workflows: an in-package
harmoniser and a compressor. The final product direction is narrower and more
explicit. The installed package is a prepared-input compressor; the old
reference-backed harmonisation/liftover workflow is external or archived.

## 1. Shipped and verified compressor

- `compress_sumstats()` accepts explicitly prepared, build-labelled input and
  writes either the self-contained native Pcodec `Z9/EAF8/SE6` store or the
  explicit Parquet interoperability backend.
- Native identity is the build-specific position plus directed REF→ALT code
  for supported primary-chromosome biallelic A/C/G/T SNVs. The store does not
  need a shared spine, FASTA, dbSNP table, or chain file when it is read.
- `qc = "compact"` is the validating path. `qc = "none"` is a trusted fast
  path for already-canonical input; it bypasses numeric QC and records that
  fact in the manifest. Supplied finite SE and EAF are handled independently.
- Full, core, HM3, and core-plus selection are implemented. Core-plus is the
  deterministic union of the bundled core panel with inclusive regions around
  `p <= 1e-5` associations, using the default inclusive +/-10,000-bp padding.
  Supplied finite p-values are authoritative; exact prepared Z is the recorded
  fallback when p is absent where the selected contract permits it.
- GRCh37 and GRCh38 are explicit same-build storage profiles. The compressor
  does not silently convert between them.
- Projection, structural/QC reporting, timing provenance, deterministic native
  payloads, empty-input handling, report/error behavior, EAF/SE authority, and
  core-plus selection are covered by the current test and package-check
  evidence. The authoritative performance claim remains the five-run FinnGen
  benchmark in `inst/benchmarks/finngen-10m-bp-20260804/`.

## 2. Intentionally external or archived preparation

The installed package does not export `harmonise_sumstats()` and does not
perform reference lookup, allele flipping/alignment, liftover, or rsID
resolution. A GWASLab/MR-Atlas, tidyGWAS, or other upstream workflow may do
those operations, but that workflow owns its reference version, build mapping,
chain file/hash, QC policy, and timings. It must hand CompreSSoR an explicit
prepared table and provenance.

The former implementation and its observability tests remain under
`archive/harmonisation/` as a migration/recovery path, not as installed API.
Issue #23 is therefore an archive/observability disposition, and #24 is
superseded/not planned for the compressor core. Neither reopens the old
harmonisation scope.

## 3. Remaining #27 evidence blocker

Issue #27 remains open solely for exact-current-base BluePebble evidence for
the large UKB-PPP NTRK3 benchmark. The maintained harness is bounded and
commit/toolchain guarded, but the current BP catalogue exposes Rust/Cargo
1.78 while native Pcodec requires >=1.87.0. No exact-`f4e783f` NTRK3 timing or
RSS claim is made until that prerequisite is available and a successful run
records source/commit, rows, columns, timing, RSS, output bytes, and
validation.

Prior NTRK3 runs and the core-plus harness are supporting evidence only; they
do not replace the exact-current-base #27 run. Core-plus itself is implemented
and tested, but no new real-data benchmark claim is made here.

## 4. Deferred work

The original roadmap's proposed two-public-workflow API, post-harmonisation
selection contract, and a general streaming/parallel compression-writer
architecture are not claims of completion. The latter remains future work
after the concrete #27 memory/evidence question; current thread metadata and
successful benchmarks must not be described as proof of a parallel writer.
Any future writer architecture should be specified and benchmarked as a
separate, bounded-memory work item without changing the locked store contract.

## Disposition

The implementation and documentation work represented by the current core
roadmap is complete enough for #18 to be closed as **superseded by the final
prepared-input compressor direction**, with #27 kept open independently.
Closing #18 must not be interpreted as claiming that the old harmoniser is a
public API, that parallel writer architecture is complete, or that the exact
current-base NTRK3 benchmark has passed.

Current issue state at audit time: #23, #24, #25, #26, #28, #30, and #31 are
closed; #27 is open; #18 remains open pending this roadmap disposition.
