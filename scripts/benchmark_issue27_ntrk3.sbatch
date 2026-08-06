#!/bin/bash
#SBATCH --job-name=compressor-issue27-ntrk3
#SBATCH --partition=test
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --account=sscm013902

set -euo pipefail

repo="${COMPRESSOR_ISSUE27_REPO:?set COMPRESSOR_ISSUE27_REPO}"
bench="${COMPRESSOR_ISSUE27_BENCH:?set COMPRESSOR_ISSUE27_BENCH}"
selection="${COMPRESSOR_ISSUE27_SELECTION:?set COMPRESSOR_ISSUE27_SELECTION}"
qc="${COMPRESSOR_ISSUE27_QC:?set COMPRESSOR_ISSUE27_QC}"
input="${COMPRESSOR_ISSUE27_INPUT:?set COMPRESSOR_ISSUE27_INPUT}"
commit="${COMPRESSOR_ISSUE27_COMMIT:?set COMPRESSOR_ISSUE27_COMMIT}"
rlib="$bench/Rlib-${commit:0:7}"
output="$bench/stores/${selection}-${qc}.cpr"
summary="$bench/results/${selection}-${qc}.json"

mkdir -p "$bench/logs" "$bench/stores" "$bench/results" "$rlib"
module purge
module load languages/R/4.5.1
if [[ -n "${COMPRESSOR_ISSUE27_RUST_MODULE:-}" ]]; then
  module load "$COMPRESSOR_ISSUE27_RUST_MODULE"
fi

if ! command -v cargo >/dev/null 2>&1; then
  module avail 2>&1 | grep -i rust | head -20 >&2 || true
  echo "benchmark setup failed: cargo is unavailable after loading R" >&2
  exit 2
fi

echo "R=$(R --version | head -1)" >&2
echo "rustc=$(rustc --version)" >&2
echo "cargo=$(cargo --version)" >&2
rust_major=$(rustc --version | awk '{print $2}' | cut -d. -f1)
rust_minor=$(rustc --version | awk '{print $2}' | cut -d. -f2)
if [[ -z "$rust_major" || -z "$rust_minor" ||
      "$rust_major" -lt 1 || ( "$rust_major" -eq 1 && "$rust_minor" -lt 87 ) ]]; then
  echo "benchmark setup failed: native Pcodec requires Rust >= 1.87.0; set COMPRESSOR_ISSUE27_RUST_MODULE" >&2
  exit 2
fi

R_LIBS="$rlib" R_LIBS_USER="$rlib" R CMD INSTALL --library="$rlib" "$repo" \
  > "$bench/logs/install_${SLURM_JOB_ID}.log" 2>&1

export R_LIBS="$rlib"
export R_LIBS_USER="$rlib"
export COMPRESSOR_ISSUE27_COMMIT="$commit"
export COMPRESSOR_ISSUE27_NATIVE_LIBRARY_COMMIT="${COMPRESSOR_ISSUE27_NATIVE_LIBRARY_COMMIT:-$commit}"

/usr/bin/time -v Rscript --vanilla "$repo/scripts/benchmark_issue27_ntrk3.R" \
  "$selection" "$qc" "$input" "$output" "$summary"
