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
rust_root="${COMPRESSOR_ISSUE27_RUST_ROOT:-}"

if ! [[ "$commit" =~ ^[0-9a-fA-F]{40}$ ]]; then
  echo "benchmark setup failed: COMPRESSOR_ISSUE27_COMMIT must be a 40-hex commit" >&2
  exit 2
fi

if command -v git >/dev/null 2>&1; then
  actual_commit="$(git -C "$repo" rev-parse --verify HEAD^{commit} 2>/dev/null || true)"
else
  actual_commit="${COMPRESSOR_ISSUE27_ACTUAL_COMMIT:-}"
  if ! [[ "$actual_commit" =~ ^[0-9a-fA-F]{40}$ ]]; then
    echo "benchmark setup failed: git is unavailable; COMPRESSOR_ISSUE27_ACTUAL_COMMIT must be a 40-hex commit resolved before sbatch" >&2
    exit 2
  fi
fi
if ! [[ "$actual_commit" =~ ^[0-9a-fA-F]{40}$ ]] ||
    [[ "$actual_commit" != "$commit" ]]; then
  echo "benchmark setup failed: requested commit $commit does not match repository HEAD $actual_commit" >&2
  exit 2
fi
if [[ "${COMPRESSOR_ISSUE27_GUARD_ONLY:-0}" == "1" ]]; then
  exit 0
fi

mkdir -p "$bench/logs" "$bench/stores" "$bench/results" "$rlib"
module purge
module load languages/R/4.5.1
if [[ -n "${COMPRESSOR_ISSUE27_RUST_MODULE:-}" ]]; then
  module load "$COMPRESSOR_ISSUE27_RUST_MODULE"
fi
if [[ -n "$rust_root" ]]; then
  if [[ ! -x "$rust_root/bin/rustc" || ! -d "$rust_root/lib" ]]; then
    echo "benchmark setup failed: COMPRESSOR_ISSUE27_RUST_ROOT must contain bin/rustc and lib; Cargo is resolved separately from PATH" >&2
    exit 2
  fi
  export PATH="$rust_root/bin:$PATH"
  export LD_LIBRARY_PATH="$rust_root/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
fi

rustc_bin="$(command -v rustc 2>/dev/null || true)"
if [[ -z "$rustc_bin" ]]; then
  echo "benchmark setup failed: rustc is unavailable after loading the configured Rust root/modules" >&2
  exit 2
fi
cargo_bin="$(command -v cargo 2>/dev/null || true)"
if [[ -z "$cargo_bin" ]]; then
  module avail 2>&1 | grep -i rust | head -20 >&2 || true
  echo "benchmark setup failed: cargo is unavailable after loading R" >&2
  exit 2
fi

echo "R=$(R --version | head -1)" >&2
rustc_version="$($rustc_bin --version)"
cargo_version="$($cargo_bin --version)"
echo "rustc=$rustc_version" >&2
echo "cargo=$cargo_version" >&2
rust_major=$(printf '%s\n' "$rustc_version" | awk '{print $2}' | cut -d. -f1)
rust_minor=$(printf '%s\n' "$rustc_version" | awk '{print $2}' | cut -d. -f2)
if [[ -z "$rust_major" || -z "$rust_minor" ||
      "$rust_major" -lt 1 || ( "$rust_major" -eq 1 && "$rust_minor" -lt 87 ) ]]; then
  echo "benchmark setup failed: native Pcodec requires Rust >= 1.87.0; set COMPRESSOR_ISSUE27_RUST_MODULE" >&2
  exit 2
fi
cargo_major=$(printf '%s\n' "$cargo_version" | awk '{print $2}' | cut -d. -f1)
cargo_minor=$(printf '%s\n' "$cargo_version" | awk '{print $2}' | cut -d. -f2)
cargo_version_re='^cargo[[:space:]]([0-9]+)\.([0-9]+)(\.[0-9]+)?([[:space:]]|[-+]|$)'
if [[ ! "$cargo_version" =~ $cargo_version_re ]]; then
  echo "benchmark setup failed: native Pcodec requires Cargo >= 1.78.0; Cargo is resolved independently from COMPRESSOR_ISSUE27_RUST_ROOT" >&2
  exit 2
fi
cargo_major="${BASH_REMATCH[1]}"
cargo_minor="${BASH_REMATCH[2]}"
if [[ "$cargo_major" -lt 1 ||
      ( "$cargo_major" -eq 1 && "$cargo_minor" -lt 78 ) ]]; then
  echo "benchmark setup failed: native Pcodec requires Cargo >= 1.78.0; Cargo is resolved independently from COMPRESSOR_ISSUE27_RUST_ROOT" >&2
  exit 2
fi

R_LIBS="$rlib" R_LIBS_USER="$rlib" R CMD INSTALL --library="$rlib" "$repo" \
  > "$bench/logs/install_${SLURM_JOB_ID}.log" 2>&1

export R_LIBS="$rlib"
export R_LIBS_USER="$rlib"
export COMPRESSOR_ISSUE27_COMMIT="$commit"
export COMPRESSOR_ISSUE27_ACTUAL_COMMIT="$actual_commit"
export COMPRESSOR_ISSUE27_NATIVE_LIBRARY_COMMIT="${COMPRESSOR_ISSUE27_NATIVE_LIBRARY_COMMIT:-$commit}"

time_prefix=()
time_wrapper="none"
if [[ -x /usr/bin/time ]]; then
  if /usr/bin/time -v true >/dev/null 2>&1; then
    time_prefix=(/usr/bin/time -v)
    time_wrapper="verbose"
  elif /usr/bin/time -p true >/dev/null 2>&1; then
    time_prefix=(/usr/bin/time -p)
    time_wrapper="portable"
  fi
fi
echo "time_wrapper=$time_wrapper" >&2

"${time_prefix[@]}" Rscript --vanilla "$repo/scripts/benchmark_issue27_ntrk3.R" \
  "$selection" "$qc" "$input" "$output" "$summary"
