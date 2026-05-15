#!/usr/bin/env bash
#SBATCH --job-name=deep-apbs-stage-data
#SBATCH --partition=datamover
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --output=/homes/tevfik/PHD/slurm_logs/deep-apbs-stage-data-%j.out
#SBATCH --error=/homes/tevfik/PHD/slurm_logs/deep-apbs-stage-data-%j.err

set -euo pipefail

REPO_DIR="${REPO_DIR:-/homes/tevfik/PHD/generate_cache}"
DATA_ROOT="${DATA_ROOT:-/nfs/production/arl/chembl/tevfik/DEEP_APBS_DATASETS}"
PYTHON_BIN="${PYTHON_BIN:-$REPO_DIR/.conda/bin/python}"

mkdir -p /homes/tevfik/PHD/slurm_logs "$DATA_ROOT"
cd "$REPO_DIR"

PATH="$REPO_DIR/.conda/bin:$PATH"

args=(--data-root "$DATA_ROOT")

if [[ "${FORCE:-0}" == "1" ]]; then
  args+=(--force)
fi
if [[ "${SKIP_SCPDB:-0}" == "1" ]]; then
  args+=(--skip-scpdb)
fi
if [[ "${SKIP_PDBBIND:-0}" == "1" ]]; then
  args+=(--skip-pdbbind)
fi
if [[ "${SKIP_PURESNET:-0}" == "1" ]]; then
  args+=(--skip-puresnet)
fi
if [[ "${SKIP_P2RANK:-0}" == "1" ]]; then
  args+=(--skip-p2rank)
fi

echo "Deep-APBS dataset staging"
echo "Host: $(hostname)"
echo "Repo: $REPO_DIR"
echo "Data root: $DATA_ROOT"
echo "Python: $PYTHON_BIN"
echo "Args: ${args[*]}"
echo

"$PYTHON_BIN" prepare_codon_datasets.py "${args[@]}"

echo
echo "Dataset staging finished."
echo "Data root: $DATA_ROOT"
