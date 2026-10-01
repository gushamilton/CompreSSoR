#!/usr/bin/env bash
# Run CompreSSoR benchmark ops, one process per op, recording peak RSS.
# Usage: run_suite.sh <store_dir> <out.csv> <reps> [op ...]
#   build first: run_suite.sh --build <store_dir> <n_rows> <traits> <out.csv> [threads]
# Env passed through: CSR_LABEL, CSR_GIT_SHA, CSR_REPO, RSCRIPT (default Rscript)
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RSCRIPT="${RSCRIPT:-Rscript}"
if [[ "$(uname)" == "Darwin" ]]; then TIMEFLAG=-l; else TIMEFLAG=-v; fi

peak_kb() {  # parse /usr/bin/time output -> KB
  if [[ "$TIMEFLAG" == "-l" ]]; then
    awk '/maximum resident set size/ {printf "%d", $1/1024}' "$1"
  else
    awk -F: '/Maximum resident set size/ {gsub(/ /,"",$2); printf "%d", $2}' "$1"
  fi
}

run_one() {  # out.csv, then Rscript args...
  local out="$1"; shift
  local tmp; tmp="$(mktemp)"; local tlog; tlog="$(mktemp)"
  rm -f "$tmp"
  /usr/bin/time "$TIMEFLAG" "$RSCRIPT" --vanilla "$here/csr_bench.R" "${@:1:$(($#-1))}" "$tmp" 2> "$tlog" \
    || { cat "$tlog" >&2; echo "FAILED: $*" >&2; rm -f "$tmp" "$tlog"; return 1; }
  local rss; rss="$(peak_kb "$tlog")"
  if [[ ! -s "$out" ]]; then head -1 "$tmp" | sed 's/$/,peak_rss_kb/' > "$out"; fi
  tail -n +2 "$tmp" | sed "s/\$/,${rss}/" >> "$out"
  rm -f "$tmp" "$tlog"
}

if [[ "${1:-}" == "--build" ]]; then
  dir="$2"; n="$3"; traits="$4"; out="$5"; threads="${6:-4}"
  # Rscript args: build dir n traits threads <tmp>  -- csr_bench expects out before threads
  tmp="$(mktemp)"; rm -f "$tmp"; tlog="$(mktemp)"
  /usr/bin/time "$TIMEFLAG" "$RSCRIPT" --vanilla "$here/csr_bench.R" build "$dir" "$n" "$traits" "$tmp" "$threads" 2> "$tlog" \
    || { cat "$tlog" >&2; exit 1; }
  rss="$(peak_kb "$tlog")"
  if [[ ! -s "$out" ]]; then head -1 "$tmp" | sed 's/$/,peak_rss_kb/' > "$out"; fi
  tail -n +2 "$tmp" | sed "s/\$/,${rss}/" >> "$out"
  rm -f "$tmp" "$tlog"
  exit 0
fi

dir="$1"; out="$2"; reps="$3"; shift 3
ops=("$@")
if [[ ${#ops[@]} -eq 0 ]]; then
  ops=()
  while IFS= read -r line; do ops+=("$line"); done < <("$RSCRIPT" --vanilla "$here/csr_bench.R" list-ops)
fi
for op in "${ops[@]}"; do
  echo "[$(date +%T)] $op" >&2
  run_one "$out" op "$dir" "$op" "$reps" placeholder || true
done
