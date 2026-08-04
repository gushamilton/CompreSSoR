#!/usr/bin/env bash
#SBATCH --job-name=compressor-panukb-scale
#SBATCH --account=sscm013902
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --array=1-3
#SBATCH --output=/user/work/fh6520/CompreSSoR-bp-thread-test/logs/%x_%A_%a.out
#SBATCH --error=/user/work/fh6520/CompreSSoR-bp-thread-test/logs/%x_%A_%a.err

set -euo pipefail
project=/user/work/fh6520/CompreSSoR-bp-thread-test
cd "$project"
module purge
module load languages/R/4.5.1
module load bcftools/1.19-openblas-5yp2
module load rust/1.78.0-n5bm
module_cargo=$(command -v cargo)
export PATH="$project/rust-1.87.0/rustc/bin:$(dirname "$module_cargo"):$PATH"
export LD_LIBRARY_PATH="$project/rust-1.87.0/rustc/lib:${LD_LIBRARY_PATH:-}"
lib="$project/r-lib-panukb-scaling-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$lib"
export R_LIBS_USER="$lib"
R CMD INSTALL --library="$lib" --no-multiarch repo
external="$project/external/benchmark-panukb"
export COMPRESSOR_SCALING_SOURCE="$external/panukb_scaling_snps.tsv.gz"
export COMPRESSOR_SCALING_RESULT="$project/results/panukb-scaling-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.csv"
export COMPRESSOR_SCALING_SCRATCH="$external/scratch-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
export COMPRESSOR_SCALING_RUNS="3"
export COMPRESSOR_SCALING_PCODEC_THREADS="4"
case "$SLURM_ARRAY_TASK_ID" in
  1) export COMPRESSOR_SCALING_SIZES="1,2,3,4,5,6,7,8,9,10" ;;
  2) export COMPRESSOR_SCALING_SIZES="11,12,13,14,15,16,17,18,19,20" ;;
  3) export COMPRESSOR_SCALING_SIZES="21,22,23,24,25,26,27,28,29" ;;
  *) echo "unexpected array task" >&2; exit 2 ;;
esac
Rscript --vanilla repo/scripts/benchmark-scaling.R
