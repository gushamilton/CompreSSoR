#!/usr/bin/env bash
#SBATCH --job-name=compressor-panukb-prep
#SBATCH --account=sscm013902
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=/user/work/fh6520/CompreSSoR-bp-thread-test/logs/%x_%j.out
#SBATCH --error=/user/work/fh6520/CompreSSoR-bp-thread-test/logs/%x_%j.err

set -euo pipefail
project=/user/work/fh6520/CompreSSoR-bp-thread-test
cd "$project"
module purge
module load languages/R/4.5.1
external="$project/external/benchmark-panukb"
mkdir -p "$external" "$project/logs"
source="$external/biomarkers-30600-both_sexes-irnt.tsv.bgz"
if [[ ! -s "$source" ]]; then
  curl --fail --location --retry 3 \
    --output "$source" \
    https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/biomarkers-30600-both_sexes-irnt.tsv.bgz
fi
export COMPRESSOR_PANUKB_SOURCE="$source"
export COMPRESSOR_SCALING_SOURCE="$external/panukb_scaling_snps.tsv.gz"
export COMPRESSOR_SCALING_METADATA="$external/panukb_scaling_metadata.json"
export COMPRESSOR_SCALING_ROWS="${COMPRESSOR_SCALING_ROWS:-21000000}"
Rscript --vanilla repo/scripts/prepare-panukb-scaling.R
