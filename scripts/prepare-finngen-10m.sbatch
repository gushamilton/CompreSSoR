#!/usr/bin/env bash
#SBATCH --job-name=compressor-prep-10m
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
module purge
module load languages/R/4.5.1
external="$project/external/benchmark-10m"
mkdir -p "$external" "$project/logs"
raw="$external/finngen_r2_ANTIDEPRESSANTS.gz"
if [[ ! -s "$raw" ]]; then
  curl --fail --location --retry 3 --output "$raw" \
    https://storage.googleapis.com/finngen-public-data-r2/summary_stats/finngen_r2_ANTIDEPRESSANTS.gz
fi
export COMPRESSOR_RAW_SOURCE="$raw"
export COMPRESSOR_10M_SOURCE="$external/finngen_10m_snps.tsv.gz"
export COMPRESSOR_10M_METADATA="$external/finngen_10m_metadata.json"
Rscript --vanilla repo/scripts/prepare-finngen-10m.R
