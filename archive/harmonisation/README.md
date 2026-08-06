# Archived harmonisation and reference workflow

This directory preserves the pre-0.5 CompreSSoR harmonisation and reference
implementation. It is intentionally outside the package's live `R/`, `man/`,
and `tests/testthat/` trees, so it is not installed or loaded by R.

The archived public API included:

- `harmonise_sumstats()` for reference-backed allele alignment and selection;
- `liftover_sumstats()` for GRCh37/hg19 to GRCh38/hg38 chain liftover;
- `build_canonical_reference()` and `build_ebi_reference()`;
- `resolve_reference()`, `download_reference()`, `grch38_reference()`, and
  `reference_cache_dir()`.

The archived R sources, generated man pages, and legacy tests retain the
former implementation and its coverage. They are not part of the installed
compression-core package. To recover the old workflow in a development
branch, restore the archived files to the live `R/`, `man/`, and test trees,
restore their NAMESPACE exports and optional Bioconductor dependencies, and
use the historical `compress_sumstats()` orchestration from the same commit.

The live package accepts only already prepared rows with explicit GRCh37 or
GRCh38 coordinates and REF/ALT orientation. External harmonisers such as the
GWASLab/MR-Atlas workflow should produce that handoff before compression.

## Archived observability

The reference-backed helper retains a quiet telemetry path for bounded
synthetic or recovery-branch runs. `harmonise_sumstats()` attaches a compact
`timings` object by default; use `observability = FALSE` or
`observability = "off"` to disable it. The object reports monotonic elapsed
and available CPU time, rows in/out/dropped, requested/effective workers, and
the phases that exist in the current path, including import, structural QC,
reference discovery/load/normalisation, liftover, chromosome partitioning,
strand handling, matching, final QC/report construction, and output
preparation.

For batch progress, pass
`observability = list(level = "events", callback = fn)` or
`observability = list(level = "events", jsonl = "run.jsonl")`. Events are
bounded metadata only and contain no alleles, coordinates, identifiers, or
statistics. Phase-boundary RSS sampling is opt-in with `memory = TRUE` and is
reported only where `/proc/self/status` is available. The
`scripts/test-harmonisation-observability.R` runner exercises the archived
API on small in-memory fixtures; it does not address the underlying full
reference materialisation/scaling issue tracked separately in #5/#24.
