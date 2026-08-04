#!/usr/bin/env bash
# Fine-grained 10k-to-1m scaling run for the public Pan-UKB-derived benchmark.
# This is deliberately a separate BP compute job so it does not interfere with
# the 1m-to-21m scaling run.
#SBATCH --job-name=compressor-panukb-fine
#SBATCH --account=sscm013902
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/user/work/fh6520/CompreSSoR-bp-thread-test/logs/%x_%j.out
#SBATCH --error=/user/work/fh6520/CompreSSoR-bp-thread-test/logs/%x_%j.err

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
lib="$project/r-lib-panukb-fine-$SLURM_JOB_ID"
mkdir -p "$lib"
export R_LIBS_USER="$lib"
R CMD INSTALL --library="$lib" --no-multiarch repo

external="$project/external/benchmark-panukb"
export COMPRESSOR_SCALING_SOURCE="$external/panukb_scaling_snps.tsv.gz"
export COMPRESSOR_SCALING_RESULT="$project/results/panukb-scaling-fine-$SLURM_JOB_ID.csv"
export COMPRESSOR_SCALING_SCRATCH="$external/scratch-fine-$SLURM_JOB_ID"
export COMPRESSOR_SCALING_RUNS="3"
export COMPRESSOR_SCALING_PCODEC_THREADS="4"
# benchmark-scaling.R interprets this variable in millions of rows.
export COMPRESSOR_SCALING_SIZES="0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1"
Rscript --vanilla repo/scripts/benchmark-scaling.R
