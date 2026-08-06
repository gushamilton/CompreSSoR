#!/bin/bash
#SBATCH --job-name=compressor-core-plus-ntrk3
#SBATCH --partition=test
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --account=sscm013902
#SBATCH --output=/user/work/fh6520/CompreSSoR-core-plus-10kb/logs/%x_%j.out
#SBATCH --error=/user/work/fh6520/CompreSSoR-core-plus-10kb/logs/%x_%j.err

set -euo pipefail
repo=/user/work/fh6520/CompreSSoR-core-plus-10kb
candidate_commit="$(git -C "$repo" rev-parse HEAD 2>/dev/null || true)"
if [[ ! "$candidate_commit" =~ ^[0-9a-fA-F]{40}$ ]]; then
  echo "benchmark setup failed: cannot resolve a 40-character commit for $repo" >&2
  exit 2
fi
bench="/user/work/fh6520/CompreSSoR-protein-e2e-20260805/core-plus-${candidate_commit:0:7}"
native_rlib=/user/work/fh6520/CompreSSoR-protein-e2e-20260805/main-tests/Rlib-927535a
mkdir -p "$bench" "$bench/logs" "$bench/stores" "$bench/results"

module purge
module load languages/R/4.5.1
export R_LIBS="$native_rlib"
export R_LIBS_USER="$native_rlib"

Rscript --vanilla -e 'if (!requireNamespace("CompreSSoR", quietly = TRUE)) stop("CompreSSoR is not available from the pinned native library")' \
  > "$bench/logs/setup_${SLURM_JOB_ID}.log" 2>&1

export COMPRESSOR_CORE_PLUS_INPUT=/user/work/fh6520/CompreSSoR-protein-e2e-20260805/main-tests/raw/NTRK3_explicit_core_columns.tsv
export COMPRESSOR_CORE_PLUS_STORE="$bench/stores/NTRK3_core_plus.cpr"
export COMPRESSOR_CORE_PLUS_RESULT="$bench/results/NTRK3_core_plus.json"
export COMPRESSOR_CORE_PLUS_COMMIT="$candidate_commit"
export COMPRESSOR_CORE_PLUS_NATIVE_LIBRARY_COMMIT=927535a
export COMPRESSOR_CORE_PLUS_SOURCE_ROOT="$repo"
export COMPRESSOR_CORE_PLUS_THREADS=4
export COMPRESSOR_CORE_PLUS_PVALUE_THRESHOLD=1e-5
export COMPRESSOR_CORE_PLUS_REGION_PADDING=10000

/usr/bin/time -v Rscript "$repo/scripts/benchmark_core_plus_ntrk3.R" \
  > "$bench/logs/benchmark_${SLURM_JOB_ID}.log" 2>&1
