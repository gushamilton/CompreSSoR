#!/usr/bin/env bash
set -euo pipefail

AGENT="$1"
LANE="/user/work/fh6520/CompreSSoR-bp-thread-test/opt-agent-threads/${AGENT}"
ROOT="/user/work/fh6520/CompreSSoR-bp-thread-test/opt-agent-threads"
REPO="/user/work/fh6520/CompreSSoR-bp-thread-test/repo"
BIN="$LANE/bin"
mkdir -p "$BIN" "$LANE/results"
module purge
module load languages/R/4.5.1
module load languages/python/3.11.9
module load rust/1.78.0-n5bm
project="/user/work/fh6520/CompreSSoR-bp-thread-test"
module_cargo="$(command -v cargo)"
export PATH="$project/rust-1.87.0/rustc/bin:$(dirname "$module_cargo"):$PATH"
export LD_LIBRARY_PATH="$project/rust-1.87.0/rustc/lib:${LD_LIBRARY_PATH:-}"
RHOME="$(R RHOME)"

g++ -std=c++17 -O3 -march=native -pthread \
  "$LANE/src/native_parallel_micro.cpp" \
  -I"$LANE/src" -L"$RHOME/lib" -L"$REPO/src" \
  -Wl,-rpath,"$RHOME/lib" -Wl,-rpath,"$REPO/src" \
  -Wl,-l:CompreSSoR.so -lR -ldl -o "$BIN/native_parallel_micro"

python3 "$LANE/src/make_manifest.py" "$ROOT" "$LANE/results/blocks.tsv"
printf 'agent\tmode\tthreads\tselection\tseconds\tvalues\tchecksum\n' > "$LANE/results/micro-results.tsv"
run_one() {
  local mode="$1" threads="$2" selection="$3"
  if "$BIN/native_parallel_micro" "$LANE/results/blocks.tsv" "$mode" "$threads" "$selection" >> "$LANE/results/micro-results.tsv" 2>> "$LANE/results/unsupported.log"; then
    return 0
  fi
  printf '%s\t%s\t%s\t%s\n' "$mode" "$threads" "$selection" "SKIPPED" >> "$LANE/results/unsupported.log"
  return 0
}

case "$AGENT" in
  agent-1-cpu)
    # The shared package is already native; this lane isolates caller overhead and
    # records the compiler/CPU-feature build used for the native scheduling probe.
    run_one serial 1 full; run_one stream 1 full; run_one stream 2 full; run_one stream 4 full; run_one stream 8 full
    ;;
  agent-2-stream)
    for t in 1 2 4 8; do run_one stream "$t" full; run_one compact "$t" full; done
    ;;
  agent-3-block)
    for t in 1 2 4 8; do run_one block "$t" full; run_one block "$t" region; done
    ;;
  agent-4-compact)
    for t in 1 2 4 8; do run_one compact "$t" full; run_one compact "$t" sparse; done
    ;;
  agent-5-access)
    for t in 1 2 4 8; do run_one block "$t" region; run_one block "$t" sparse; run_one stream "$t" sparse; done
    ;;
  *) echo "unknown agent $AGENT" >&2; exit 2 ;;
esac

date -u +%FT%TZ > "$LANE/results/completed.utc"
