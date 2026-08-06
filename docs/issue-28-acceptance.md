# Issue #28 audit disposition

Audit base: `origin/main` `5f553724a9cebac4609447df68c418c70f202e1f`.

The prepared-input compressor refactor is complete at the implementation
level, and issue #28 can close independently. The exact-current-base
BluePebble large-file benchmark remains a separate issue under #27; it is an
infrastructure/toolchain blocker, not a missing compressor feature.

## Acceptance evidence

- **Prepared contract and resolution:** common aliases and deterministic
  beta/SE, OR/SE, optional Z, EAF, and p-value resolution are documented and
  covered by `test-resolution-qc.R` and `test-input-projection.R`. The strict
  writer does not harmonise, liftover, flip alleles, infer builds, or resolve
  rsIDs.
- **Projection and allocations:** native and default Parquet ingestion expose
  `columns_before`, `columns_read`, and `source_columns_read`; the native
  path projects only recognised core fields. `qc="none"` requires exact
  canonical columns and does not create an absent p-value column or row-level
  QC objects.
- **QC modes:** `qc="compact"` is the safe structural/numeric-QC path;
  `qc="none"` is an explicitly trusted fast path with identity safety and
  native duplicate protection, but numeric validation bypassed. Both policies
  and timings are recorded in manifests.
- **Selection:** full, core, HM3, and core-plus are covered by the selection
  tests. Core-plus uses inclusive `p <= 1e-5` regions with inclusive
  `+/-10,000 bp` padding, unions them with core, treats finite supplied
  p-values as authoritative, and records fallback/provenance when p is absent
  or invalid.
- **Stores and builds:** native Pcodec emits the locked self-contained
  `Z9/EAF8/SE6` identity/value streams with tolerances and checksums. Parquet
  remains the interoperability backend. GRCh37 and GRCh38 same-build stores,
  cross-build rejection, native readers, deterministic payloads, and full
  validation are covered by the test suite.
- **Verification:** focused tests, the full `testthat` suite, the issue #27
  commit-guard smoke test, source build, and `R CMD check --no-manual --as-cran`
  all pass. The package check has 0 errors, 2 known warnings, and 4 known
  notes.

## Residual blocker

The exact-current-base NTRK3 benchmark required to close #27 was not claimed.
Read-only BP inspection on the audit date exposed only `rust/1.78.0-n5bm`,
while the native backend requires Rust/Cargo `>=1.87.0`; the benchmark harness
therefore fails before reading the input. The merged harness verifies the
requested checkout commit and toolchain before creating benchmark outputs, so
the prior external NTRK3 timings remain explicitly prior-commit supporting
evidence rather than exact-current-base acceptance evidence.

## Disposition

Issue #28 is complete and can close independently. Issue #27 remains open as a
separate exact-current-base BluePebble benchmark/toolchain item. Closing #28
does not imply that the exact-current-base NTRK3 timing or RSS evidence exists;
the prior-commit timings remain supporting evidence until #27 obtains a
successful benchmark recording its source commit, columns, timings, RSS, rows,
output size, and validation result.
