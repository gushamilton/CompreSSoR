#!/usr/bin/env bash
# Short compute-node diagnostic for the user-facing R read path.
#SBATCH --job-name=compressor-probe-public
#SBATCH --account=sscm013902
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=00:05:00
#SBATCH --output=/user/work/fh6520/CompreSSoR-bp-thread-test/logs/public-thread-%j.out
#SBATCH --error=/user/work/fh6520/CompreSSoR-bp-thread-test/logs/public-thread-%j.err

set -euo pipefail
project=/user/work/fh6520/CompreSSoR-bp-thread-test
cd "$project"
module purge
module load languages/R/4.5.1
module load bcftools/1.19-openblas-5yp2
module load rust/1.78.0-n5bm
export R_LIBS_USER="$project/r-lib-final-selected-1"
export COMPRESSOR_PROBE_STORE="$project/results/final-selected-20260804/run-1/formats/compressor.cpr"
export COMPRESSOR_PROBE_OUTPUT="$project/results/final-selected-20260804/public-thread-probe.csv"
Rscript --vanilla repo/scripts/probe-public-threads.R
