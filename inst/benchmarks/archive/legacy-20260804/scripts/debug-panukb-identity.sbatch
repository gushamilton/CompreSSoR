#!/usr/bin/env bash
#SBATCH --job-name=compressor-debug-id
#SBATCH --account=sscm013902
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=/user/work/fh6520/CompreSSoR-bp-thread-test/logs/%x_%j.out
#SBATCH --error=/user/work/fh6520/CompreSSoR-bp-thread-test/logs/%x_%j.err
set -euo pipefail
project=/user/work/fh6520/CompreSSoR-bp-thread-test
cd "$project"
module purge
module load languages/R/4.5.1
export COMPRESSOR_SCALING_SOURCE="$project/external/benchmark-panukb/panukb_scaling_snps.tsv.gz"
export COMPRESSOR_DEBUG_ROWS="${COMPRESSOR_DEBUG_ROWS:-2000000}"
Rscript --vanilla repo/scripts/debug-panukb-identity.R
