#!/bin/bash

# Focused shell smoke test for the issue-27 Slurm harness commit pinning.
# It uses guard-only mode, so it never loads modules or installs the package.

set -euo pipefail

script_dir="$(CDPATH= cd -- "$(dirname "$0")/../scripts" && pwd)"
harness="$script_dir/benchmark_issue27_ntrk3.sbatch"
work="$(mktemp -d "${TMPDIR:-/tmp}/compressor-issue27-guard.XXXXXX")"
trap 'rm -rf "$work"' EXIT

repo="$work/repo"
bench="$work/benchmark"
mkdir -p "$repo"
git -C "$repo" init -q
git -C "$repo" config user.email issue27-test@example.invalid
git -C "$repo" config user.name issue27-test
printf '%s\n' marker > "$repo/marker"
git -C "$repo" add marker
git -C "$repo" commit -q -m marker

summary="$bench/results/core-none.json"
set +e
COMPRESSOR_ISSUE27_REPO="$repo" \
COMPRESSOR_ISSUE27_BENCH="$bench" \
COMPRESSOR_ISSUE27_SELECTION=core \
COMPRESSOR_ISSUE27_QC=none \
COMPRESSOR_ISSUE27_INPUT="$work/input.tsv" \
COMPRESSOR_ISSUE27_COMMIT=0000000000000000000000000000000000000000 \
COMPRESSOR_ISSUE27_GUARD_ONLY=1 \
bash "$harness" > "$work/stdout" 2> "$work/stderr"
status=$?
set -e

if [[ "$status" -ne 2 ]]; then
  echo "expected all-zero commit to fail with status 2; got $status" >&2
  exit 1
fi
if [[ -e "$summary" || -d "$bench" ]]; then
  echo "commit guard created benchmark output before rejecting the request" >&2
  exit 1
fi

echo "issue27 harness guard smoke test: PASS"
