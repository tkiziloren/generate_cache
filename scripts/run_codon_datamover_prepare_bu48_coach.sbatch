#!/usr/bin/env bash
#SBATCH --job-name=bu48-coach-stage
#SBATCH --partition=datamover
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --output=/homes/tevfik/PHD/slurm_logs/bu48-coach-stage-%j.out
#SBATCH --error=/homes/tevfik/PHD/slurm_logs/bu48-coach-stage-%j.err

set -euo pipefail

REPO_DIR="${REPO_DIR:-/homes/tevfik/PHD/generate_cache}"
DATA_ROOT="${DATA_ROOT:-/nfs/production/arl/chembl/tevfik/DEEP_APBS_DATASETS}"
PYTHON_BIN="${PYTHON_BIN:-$REPO_DIR/.conda/bin/python}"

mkdir -p /homes/tevfik/PHD/slurm_logs "$DATA_ROOT/datasets/external_benchmarks"
cd "$REPO_DIR"

PATH="$REPO_DIR/.conda/bin:$PATH"

args=(--data-root "$DATA_ROOT" --repo-dir "$REPO_DIR")

if [[ "${FORCE:-0}" == "1" ]]; then
  args+=(--force)
fi
if [[ "${SKIP_PREPARED:-0}" == "1" ]]; then
  args+=(--skip-prepared)
fi

echo "BU48/COACH staging"
echo "Host: $(hostname)"
echo "Repo: $REPO_DIR"
echo "Data root: $DATA_ROOT"
echo "Python: $PYTHON_BIN"
echo "Args: ${args[*]}"
echo

"$PYTHON_BIN" prepare_codon_bu48_coach.py "${args[@]}"

echo
echo "BU48/COACH staging finished."
