#!/usr/bin/env bash
# Five-repeat rerun of the locked production candidate only.
# Submit this through Slurm on BP; it is intentionally not a login-node script.
#SBATCH --job-name=compressor-final-selected
#SBATCH --account=sscm013902
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --array=1-5
#SBATCH --output=/user/work/fh6520/CompreSSoR-bp-thread-test/logs/final-selected-%A_%a.out
#SBATCH --error=/user/work/fh6520/CompreSSoR-bp-thread-test/logs/final-selected-%A_%a.err

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
lib="$project/r-lib-final-selected-$SLURM_ARRAY_TASK_ID"
mkdir -p "$lib"
export R_LIBS_USER="$lib"
R CMD INSTALL --library="$lib" --no-multiarch repo

export COMPRESSOR_FINAL_SOURCE="$project/data/finngen_chr1_core_snvs.tsv.gz"
export COMPRESSOR_FINAL_ONLY=pcodec_native
export COMPRESSOR_FINAL_ROOT="${COMPRESSOR_FINAL_ROOT:-$project/results/final-selected-20260804}"
export COMPRESSOR_FINAL_RUN="${COMPRESSOR_FINAL_RUN:-$SLURM_ARRAY_TASK_ID}"
export COMPRESSOR_FINAL_KEEP_FORMATS=1
Rscript --vanilla repo/scripts/benchmark-final-keyed.R
