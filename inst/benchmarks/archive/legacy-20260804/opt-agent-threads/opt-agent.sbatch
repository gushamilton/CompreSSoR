#!/bin/bash
#SBATCH --partition=test
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=8
#SBATCH --account=sscm013902
#SBATCH --output=logs/%x_%A.out
#SBATCH --error=logs/%x_%A.err
set -euo pipefail
LANE_ROOT=/user/work/fh6520/CompreSSoR-bp-thread-test/opt-agent-threads
AGENT="${AGENT:?AGENT is required}"
mkdir -p "$LANE_ROOT/$AGENT/logs"
bash "$LANE_ROOT/$AGENT/src/../run_agent.sh" "$AGENT"
