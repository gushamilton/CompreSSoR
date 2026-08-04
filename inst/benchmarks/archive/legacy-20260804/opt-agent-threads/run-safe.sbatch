#!/bin/bash
#SBATCH --job-name=compressor-native-safe
#SBATCH --partition=test
#SBATCH --time=00:45:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=8
#SBATCH --account=sscm013902
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
project=/user/work/fh6520/CompreSSoR-bp-thread-test
lane="$project/opt-agent-threads"
mkdir -p "$lane/logs" "$lane/results"
module purge
module load languages/R/4.5.1
module load bcftools/1.19-openblas-5yp2
module load rust/1.78.0-n5bm
module_cargo="$(command -v cargo)"
export PATH="$project/rust-1.87.0/rustc/bin:$(dirname "$module_cargo"):$PATH"
export LD_LIBRARY_PATH="$project/rust-1.87.0/rustc/lib:${LD_LIBRARY_PATH:-}"
for agent in agent-1-cpu agent-2-stream agent-3-block agent-4-compact agent-5-access; do
  bash "$lane/$agent/run_agent.sh" "$agent"
done
date -u +%FT%TZ > "$lane/results/all-complete.utc"
