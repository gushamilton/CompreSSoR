#!/bin/bash
# Login-node helper: stage harness, fetch SHAs, then submit the A/B array.
# Usage: submit_csr_perf.sh <replicates> label:sha [label:sha ...]
# Not run automatically; review before use.
set -euo pipefail
root="${CSR_ROOT:-/user/work/fh6520/CompreSSoR-perf}"
reps="$1"; shift
mkdir -p "$root/logs" "$root/results" "$root/harness"
[[ -d "$root/repo/.git" ]] || git clone --quiet https://github.com/gushamilton/CompreSSoR.git "$root/repo"
git -C "$root/repo" fetch --quiet origin '+refs/heads/*:refs/remotes/origin/*'
for pair in "$@"; do git -C "$root/repo" cat-file -e "${pair#*:}^{commit}"; done
here="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cp "$here"/csr_bench.R "$here"/csr_bench_common.R "$here"/run_suite.sh "$root/harness/"
export CSR_SHAS="$*" CSR_OUT="${CSR_OUT:-$root/results/$(date +%Y%m%d-%H%M%S)}"
n=$(( $# * reps ))
echo "submitting $n tasks -> $CSR_OUT"
sbatch --array=0-$((n-1))%8 --export=ALL "$here/bp/csr_perf_array.sbatch"
