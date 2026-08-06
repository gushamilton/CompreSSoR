#!/bin/bash

# Focused shell smoke test for the issue-27 Slurm harness commit pinning.
# It uses guard-only mode, so it never loads modules or installs the package.

set -euo pipefail

script_dir="$(CDPATH= cd -- "$(dirname "$0")/../scripts" && pwd)"
harness="$script_dir/benchmark_issue27_ntrk3.sbatch"
work="$(mktemp -d "${TMPDIR:-/tmp}/compressor-issue27-guard.XXXXXX")"
trap 'rm -rf "$work"' EXIT

repo="$work/repo"
mkdir -p "$repo"
git -C "$repo" init -q
git -C "$repo" config user.email issue27-test@example.invalid
git -C "$repo" config user.name issue27-test
printf '%s\n' marker > "$repo/marker"
git -C "$repo" add marker
git -C "$repo" commit -q -m marker

expected="0123456789abcdef0123456789abcdef01234567"
no_git_bin="$work/no-git-bin"
mkdir -p "$no_git_bin"

run_no_git_guard() {
  local bench="$1"
  local actual_commit="${2:-}"
  local command=(
    /usr/bin/env
    "PATH=$no_git_bin"
    "COMPRESSOR_ISSUE27_REPO=$repo"
    "COMPRESSOR_ISSUE27_BENCH=$bench"
    COMPRESSOR_ISSUE27_SELECTION=core
    COMPRESSOR_ISSUE27_QC=none
    "COMPRESSOR_ISSUE27_INPUT=$work/input.tsv"
    "COMPRESSOR_ISSUE27_COMMIT=$expected"
    COMPRESSOR_ISSUE27_GUARD_ONLY=1
  )
  if [[ -n "$actual_commit" ]]; then
    command+=("COMPRESSOR_ISSUE27_ACTUAL_COMMIT=$actual_commit")
  fi
  command+=(/bin/bash "$harness")
  set +e
  "${command[@]}" > "$work/stdout" 2> "$work/stderr"
  local status=$?
  set -e
  printf '%s\n' "$status"
}

correct_bench="$work/benchmark-correct"
if [[ "$(run_no_git_guard "$correct_bench" "$expected")" -ne 0 ]]; then
  echo "expected no-git guard with a pre-resolved matching commit to pass" >&2
  cat "$work/stderr" >&2
  exit 1
fi
if [[ -e "$correct_bench" ]]; then
  echo "no-git guard created benchmark output before guard-only completion" >&2
  exit 1
fi

missing_bench="$work/benchmark-missing"
if [[ "$(run_no_git_guard "$missing_bench")" -ne 2 ]]; then
  echo "expected missing pre-resolved commit to fail with status 2" >&2
  exit 1
fi
if [[ -e "$missing_bench" ]]; then
  echo "missing-commit guard created benchmark output before rejecting the request" >&2
  exit 1
fi

mismatch_bench="$work/benchmark-mismatch"
if [[ "$(run_no_git_guard "$mismatch_bench" ffffffffffffffffffffffffffffffffffffffff)" -ne 2 ]]; then
  echo "expected mismatched pre-resolved commit to fail with status 2" >&2
  exit 1
fi
if [[ -e "$mismatch_bench" ]]; then
  echo "mismatched-commit guard created benchmark output before rejecting the request" >&2
  exit 1
fi

split_root="$work/rust-root"
split_cargo_path="$work/cargo-path"
fake_tools="$work/fake-tools"
split_bench="$work/benchmark-split-toolchain"
split_cargo_used="$work/cargo-used"
split_rustc_used="$work/rustc-used"
split_rscript_used="$work/rscript-used"
mkdir -p "$split_root/bin" "$split_root/lib" "$split_cargo_path" "$fake_tools"
printf '%s\n' \
  '#!/bin/bash' \
  'printf "%s\\n" "$0" > "${ISSUE27_FAKE_RUSTC_USED:?}"' \
  'printf "%s\\n" "rustc 1.87.0 (issue27-test)"' \
  > "$split_root/bin/rustc"
chmod +x "$split_root/bin/rustc"
printf '%s\n' \
  '#!/bin/bash' \
  'printf "%s\\n" "$0" > "${ISSUE27_FAKE_CARGO_USED:?}"' \
  'printf "%s\\n" "cargo ${ISSUE27_FAKE_CARGO_VERSION:-1.78.0-nightly} (issue27-test)"' \
  > "$split_cargo_path/cargo"
chmod +x "$split_cargo_path/cargo"
printf '%s\n' '#!/bin/bash' 'exit 0' > "$fake_tools/module"
printf '%s\n' \
  '#!/bin/bash' \
  'if [[ "$1" == "--version" ]]; then printf "%s\\n" "R version 4.5.1 (issue27-test)"; fi' \
  'exit 0' \
  > "$fake_tools/R"
printf '%s\n' \
  '#!/bin/bash' \
  'printf "%s\\n" "$0" > "${ISSUE27_FAKE_RSCRIPT_USED:?}"' \
  'exit 0' \
  > "$fake_tools/Rscript"
chmod +x "$fake_tools/module" "$fake_tools/R" "$fake_tools/Rscript"

split_actual_commit="$(git -C "$repo" rev-parse --verify HEAD^{commit})"
set +e
PATH="$split_cargo_path:$fake_tools:$PATH" \
COMPRESSOR_ISSUE27_REPO="$repo" \
COMPRESSOR_ISSUE27_BENCH="$split_bench" \
COMPRESSOR_ISSUE27_SELECTION=core \
COMPRESSOR_ISSUE27_QC=none \
COMPRESSOR_ISSUE27_INPUT="$work/input.tsv" \
COMPRESSOR_ISSUE27_COMMIT="$split_actual_commit" \
COMPRESSOR_ISSUE27_ACTUAL_COMMIT="$split_actual_commit" \
COMPRESSOR_ISSUE27_RUST_ROOT="$split_root" \
SLURM_JOB_ID=issue27-split-test \
ISSUE27_FAKE_CARGO_USED="$split_cargo_used" \
ISSUE27_FAKE_RUSTC_USED="$split_rustc_used" \
ISSUE27_FAKE_RSCRIPT_USED="$split_rscript_used" \
ISSUE27_FAKE_CARGO_VERSION=1.78.0-nightly \
/bin/bash "$harness" > "$work/split-stdout" 2> "$work/split-stderr"
status=$?
set -e
if [[ "$status" -ne 0 ]]; then
  echo "expected split Rust-root/Cargo-PATH setup to pass; got $status" >&2
  cat "$work/split-stderr" >&2
  exit 1
fi
if [[ -e "$split_root/bin/cargo" || ! -s "$split_cargo_used" ||
      ! -s "$split_rustc_used" || ! -s "$split_rscript_used" ]]; then
  echo "split toolchain regression did not use separate Cargo and Rust paths" >&2
  exit 1
fi
if grep -q "illegal option" "$work/split-stderr"; then
  echo "portable timing wrapper used an unsupported Linux-only option" >&2
  cat "$work/split-stderr" >&2
  exit 1
fi

for invalid_cargo_version in 1.77.9 not-a-version; do
  invalid_bench="$work/benchmark-invalid-cargo-${invalid_cargo_version//[^[:alnum:]]/_}"
  set +e
  PATH="$split_cargo_path:$fake_tools:$PATH" \
  COMPRESSOR_ISSUE27_REPO="$repo" \
  COMPRESSOR_ISSUE27_BENCH="$invalid_bench" \
  COMPRESSOR_ISSUE27_SELECTION=core \
  COMPRESSOR_ISSUE27_QC=none \
  COMPRESSOR_ISSUE27_INPUT="$work/input.tsv" \
  COMPRESSOR_ISSUE27_COMMIT="$split_actual_commit" \
  COMPRESSOR_ISSUE27_ACTUAL_COMMIT="$split_actual_commit" \
  COMPRESSOR_ISSUE27_RUST_ROOT="$split_root" \
  SLURM_JOB_ID="issue27-invalid-${invalid_cargo_version//[^[:alnum:]]/_}" \
  ISSUE27_FAKE_CARGO_USED="$split_cargo_used" \
  ISSUE27_FAKE_RUSTC_USED="$split_rustc_used" \
  ISSUE27_FAKE_RSCRIPT_USED="$split_rscript_used" \
  ISSUE27_FAKE_CARGO_VERSION="$invalid_cargo_version" \
  /bin/bash "$harness" > "$work/invalid-cargo-stdout" 2> "$work/invalid-cargo-stderr"
  status=$?
  set -e
  if [[ "$status" -ne 2 ]]; then
    echo "expected Cargo version $invalid_cargo_version to fail with status 2; got $status" >&2
    cat "$work/invalid-cargo-stderr" >&2
    exit 1
  fi
done

summary="$work/benchmark-git/results/core-none.json"
set +e
COMPRESSOR_ISSUE27_REPO="$repo" \
COMPRESSOR_ISSUE27_BENCH="$work/benchmark-git" \
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
if [[ -e "$summary" || -d "$work/benchmark-git" ]]; then
  echo "commit guard created benchmark output before rejecting the request" >&2
  exit 1
fi

echo "issue27 harness guard smoke test: PASS"
