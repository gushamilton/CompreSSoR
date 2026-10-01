#!/bin/bash
# CompreSSoR performance A/B array job for BluePebble.
# Each array task = one (SHA, replicate) pair; tasks interleave SHAs so
# baseline/optimised land on a mix of nodes. Every task installs its SHA into
# a private library, builds synthetic stores on node-local scratch (timing the
# compress path), runs the op suite one-process-per-op under /usr/bin/time -v,
# and writes one CSV (hostname, CPU, nproc, R, SHA, wall, peak RSS, checksum).
#
# Submit (from a login node; do not run on login):
#   export CSR_SHAS="baseline:<40hex> optimised:<40hex>" CSR_OUT=/user/work/fh6520/CompreSSoR-perf/results
#   n=$(( $(wc -w <<<"$CSR_SHAS") * 5 ))   # 5 replicates per SHA
#   sbatch --array=0-$((n-1))%6 --export=ALL csr_perf_array.sbatch
#SBATCH --job-name=csr-perf
#SBATCH --account=sscm013902
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=/user/work/fh6520/CompreSSoR-perf/logs/%x_%A_%a.out
#SBATCH --error=/user/work/fh6520/CompreSSoR-perf/logs/%x_%A_%a.err

set -euo pipefail
root="${CSR_ROOT:-/user/work/fh6520/CompreSSoR-perf}"
repo="${CSR_REPO_SRC:-$root/repo}"              # clone of gushamilton/CompreSSoR
harness="${CSR_HARNESS:-$root/harness}"          # dir holding csr_bench*.R + run_suite.sh
out_dir="${CSR_OUT:?set CSR_OUT}"
read -r -a pairs <<<"${CSR_SHAS:?set CSR_SHAS='label:sha label:sha'}"
reps_per_op="${CSR_OP_REPS:-3}"
n_rows="${CSR_N_ROWS:-1e7}"
multi_rows="${CSR_MULTI_ROWS:-2e6}"
multi_traits="${CSR_MULTI_TRAITS:-8}"
rust_root="${CSR_RUST_ROOT:-/user/work/fh6520/CompreSSoR-bp-thread-test/rust-1.87.0/rustc}"
rust_module="${CSR_RUST_MODULE:-rust/1.78.0-n5bm}"   # supplies cargo only

task="${SLURM_ARRAY_TASK_ID:?must run as an array task}"
npairs=${#pairs[@]}
pair="${pairs[$((task % npairs))]}"
replicate=$((task / npairs + 1))
label="${pair%%:*}"; sha="${pair#*:}"
[[ "$sha" =~ ^[0-9a-fA-F]{40}$ ]] || { echo "SHA for $label must be 40-hex" >&2; exit 2; }

module purge
module load languages/R/4.5.1
# Cargo >= 1.85 is needed for edition-2024 crates in Cargo.lock. Point
# CSR_CARGO_DIR at a newer cargo/bin (e.g. extracted from the rust 1.87 dist
# tarball); otherwise fall back to the (older) rust module's cargo.
if [[ -n "${CSR_CARGO_DIR:-}" ]]; then
  cargo_dir="$CSR_CARGO_DIR"
else
  module load "$rust_module"
  cargo_dir="$(dirname "$(command -v cargo)")"
fi
# Compute nodes may be offline: crates must be prefetched (cargo fetch --locked)
# on the login node into the shared ~/.cargo registry.
export CARGO_NET_OFFLINE="${CARGO_NET_OFFLINE:-true}"
export PATH="$rust_root/bin:$cargo_dir:$PATH"
export LD_LIBRARY_PATH="$rust_root/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
rustc --version >&2; cargo --version >&2
rv=$(rustc --version | awk '{print $2}')
[[ "$(printf '%s\n1.87.0\n' "$rv" | sort -V | head -1)" == "1.87.0" ]] || { echo "rustc >= 1.87 required (got $rv)" >&2; exit 2; }

work="${TMPDIR:-/tmp}/csr-perf-$SLURM_JOB_ID-$task"
mkdir -p "$work" "$out_dir" "$root/logs"
trap 'rm -rf "$work"' EXIT
# Compute nodes may lack git: prefer a source tree exported on the login node
# (git archive <sha> into $root/src-<sha>); otherwise take a private clone.
if [[ -d "$root/src-$sha" ]]; then
  cp -R "$root/src-$sha" "$work/src"
else
  git clone --quiet --shared --no-checkout "$repo" "$work/src"
  git -C "$work/src" checkout --quiet --detach "$sha"
  actual="$(git -C "$work/src" rev-parse HEAD)"
  [[ "$actual" == "$sha" ]] || { echo "checkout mismatch $actual != $sha" >&2; exit 2; }
fi
lib="$work/lib"; mkdir -p "$lib"
R CMD INSTALL --no-multiarch --no-test-load --library="$lib" "$work/src" >&2

export R_LIBS="$lib" R_LIBS_USER="$lib"
export CSR_LABEL="$label" CSR_GIT_SHA="$sha"
csv="$out_dir/csr-perf_${label}_${sha:0:7}_rep${replicate}_${SLURM_JOB_ID}_${task}.csv"

bash "$harness/run_suite.sh" --build "$work/s10m" "$n_rows" 1 "$csv" 4
bash "$harness/run_suite.sh" "$work/s10m" "$csv" "$reps_per_op" \
  full_mr_t1 full_mr_t4 full_p_t4 rows_1k rows_100k keys_1k keys_10k \
  cis50_p01 cand_flag cand_p01_decode_filter order_exact
bash "$harness/run_suite.sh" "$work/s10m" "$csv" 1 cand_full_rowids
bash "$harness/run_suite.sh" --build "$work/multi" "$multi_rows" "$multi_traits" "$csv" 4
bash "$harness/run_suite.sh" "$work/multi" "$csv" "$reps_per_op" \
  batch_keys1k_t1 batch_keys1k_t4 multi_cand_flag multi_cis50_p01
echo "wrote $csv" >&2
