# Project instructions

- Read `README.md` before working.
- Preserve original data supplied under `inputs/`.
- Put reusable code in `R/`, `src/` and `scripts/`; put generated artifacts in `outputs/`.
- Main verification command: `R CMD check . --no-manual --as-cran` on a machine with `arrow` installed.
- Run all heavy integration, stress and timing work on the Mac mini (`home`), not on a login node or laptop.
- Keep raw GWAS, reference tables, temporary stores, caches and logs outside the synced project. On the mini use `/Volumes/crucial_x9/CompreSSoR-benchmarks/`; keep only compact CSV/JSON/Markdown summaries in `outputs/`.
- The standard public store is q9/SE10 Parquet-backed; the q8 framed cache is optional.
- Never call a lossy store exact; manifests must record profile, codec, reference identity and tolerances.
