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
