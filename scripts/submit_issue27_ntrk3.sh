#!/bin/bash

# Resolve the checked-out commit on a BP login node, then pass it explicitly to
# the Slurm harness. Compute nodes may not have git, so this is part of the
# benchmark provenance contract.
#
# Usage:
#   COMPRESSOR_ISSUE27_REPO=/path/to/checkout \
#   COMPRESSOR_ISSUE27_COMMIT=$(git -C /path/to/checkout rev-parse HEAD) \
#   scripts/submit_issue27_ntrk3.sh [sbatch options] \
#     scripts/benchmark_issue27_ntrk3.sbatch

set -euo pipefail

repo="${COMPRESSOR_ISSUE27_REPO:?set COMPRESSOR_ISSUE27_REPO to the checkout}"
commit="${COMPRESSOR_ISSUE27_COMMIT:?set COMPRESSOR_ISSUE27_COMMIT to the requested 40-hex commit}"

if ! command -v git >/dev/null 2>&1; then
  echo "issue27 submission failed: git is required on the login/submission node" >&2
  exit 2
fi
if ! [[ "$commit" =~ ^[0-9a-fA-F]{40}$ ]]; then
  echo "issue27 submission failed: COMPRESSOR_ISSUE27_COMMIT must be a 40-hex commit" >&2
  exit 2
fi

actual_commit="$(git -C "$repo" rev-parse --verify HEAD^{commit} 2>/dev/null || true)"
if ! [[ "$actual_commit" =~ ^[0-9a-fA-F]{40}$ ]]; then
  echo "issue27 submission failed: could not resolve HEAD in $repo" >&2
  exit 2
fi
if [[ "$actual_commit" != "$commit" ]]; then
  echo "issue27 submission failed: requested commit $commit does not match checkout HEAD $actual_commit" >&2
  exit 2
fi

if [[ "$#" -eq 0 ]]; then
  echo "usage: $0 [sbatch options] scripts/benchmark_issue27_ntrk3.sbatch" >&2
  exit 2
fi

export COMPRESSOR_ISSUE27_ACTUAL_COMMIT="$actual_commit"
exec sbatch --export="ALL,COMPRESSOR_ISSUE27_ACTUAL_COMMIT=$actual_commit" "$@"
