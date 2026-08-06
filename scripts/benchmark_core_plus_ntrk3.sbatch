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
bench=/user/work/fh6520/CompreSSoR-protein-e2e-20260805/core-plus-6ed3303
rlib="$bench/Rlib"
mkdir -p "$bench" "$bench/logs" "$bench/stores" "$bench/results" "$rlib"

module purge
module load languages/R/4.5.1
module load rust/1.78.0-n5bm
export R_LIBS_USER="$rlib"
R CMD INSTALL --library="$rlib" "$repo" > "$bench/logs/install_${SLURM_JOB_ID}.log" 2>&1

export COMPRESSOR_CORE_PLUS_INPUT=/user/work/fh6520/CompreSSoR-protein-e2e-20260805/main-tests/raw/NTRK3_explicit_core_columns.tsv
export COMPRESSOR_CORE_PLUS_STORE="$bench/stores/NTRK3_core_plus.cpr"
export COMPRESSOR_CORE_PLUS_RESULT="$bench/results/NTRK3_core_plus.json"
export COMPRESSOR_CORE_PLUS_COMMIT=6ed33033a047820689a7f30d6fa0ea25be6be80b
export COMPRESSOR_CORE_PLUS_THREADS=4
export COMPRESSOR_CORE_PLUS_PVALUE_THRESHOLD=1e-5
export COMPRESSOR_CORE_PLUS_REGION_PADDING=10000

/usr/bin/time -v Rscript "$repo/scripts/benchmark_core_plus_ntrk3.R" \
  > "$bench/logs/benchmark_${SLURM_JOB_ID}.log" 2>&1
