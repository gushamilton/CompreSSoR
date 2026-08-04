#!/usr/bin/env bash
#SBATCH --job-name=compressor-pcodec-thread-compare
#SBATCH --account=sscm013902
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
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
lib="$project/r-lib-pcodec-thread-compare-$SLURM_JOB_ID"
mkdir -p "$lib"
export R_LIBS_USER="$lib"
R CMD INSTALL --library="$lib" --no-multiarch repo
export COMPRESSOR_FINAL_SOURCE="$project/external/benchmark-10m/finngen_10m_snps.tsv.gz"
export COMPRESSOR_PCODEC_STORE="$project/external/benchmark-10m/pcodec_thread_compare.cpr"
export COMPRESSOR_PCODEC_RESULT="$project/results/pcodec-10m-thread-compare-$SLURM_JOB_ID.csv"
Rscript --vanilla repo/scripts/benchmark-pcodec-thread-compare.R
