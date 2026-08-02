# Handoff

Goal: harden CompreSSoR into a standalone GitHub R package that accepts varied
GWAS summary-statistics files, performs explicit conversion/QC choices, and
has evidence-backed round-trip, stress and speed tests.

Completed: core package, reference anchoring, conversion/QC/core/HM3 modes,
initial edge tests and shipped storage benchmarks.

Next: run the expanded audit and all heavy integration/performance tests on the
Mac mini. Keep raw inputs and temporary stores on
`/Volumes/crucial_x9/CompreSSoR-benchmarks/`; return only compact summaries and
reports to `outputs/`.

Run or check: `R CMD check . --no-manual --as-cran` on the mini after each
release-gate change.

Cautions: do not run jobs on BP login nodes; do not put raw GWAS or large
temporary stores in the synced project.
