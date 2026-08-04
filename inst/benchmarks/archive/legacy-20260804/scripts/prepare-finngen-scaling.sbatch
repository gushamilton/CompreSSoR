#!/usr/bin/env bash
#SBATCH --job-name=compressor-prep-scaling
#SBATCH --account=sscm013902
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --output=/user/work/fh6520/CompreSSoR-bp-thread-test/logs/%x_%j.out
#SBATCH --error=/user/work/fh6520/CompreSSoR-bp-thread-test/logs/%x_%j.err

set -euo pipefail
project=/user/work/fh6520/CompreSSoR-bp-thread-test
cd "$project"
module purge
module load languages/R/4.5.1
external="$project/external/benchmark-scaling"
mkdir -p "$external" "$project/logs"
raw="$project/external/benchmark-10m/finngen_r2_ANTIDEPRESSANTS.gz"
export COMPRESSOR_SCALING_RAW_SOURCE="$raw"
export COMPRESSOR_SCALING_SOURCE="$external/finngen_scaling_snps.tsv.gz"
export COMPRESSOR_SCALING_METADATA="$external/finngen_scaling_metadata.json"
export COMPRESSOR_SCALING_ROWS="${COMPRESSOR_SCALING_ROWS:-15000000}"
Rscript --vanilla repo/scripts/prepare-finngen-scaling.R
