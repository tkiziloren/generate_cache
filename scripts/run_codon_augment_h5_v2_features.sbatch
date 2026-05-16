#!/usr/bin/env bash
#SBATCH --job-name=v2-augment
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/homes/tevfik/PHD/slurm_logs/v2-augment-%j.out
#SBATCH --error=/homes/tevfik/PHD/slurm_logs/v2-augment-%j.err

set -euo pipefail

REPO_DIR="${REPO_DIR:-/homes/tevfik/PHD/generate_cache}"
DATA_ROOT="${DATA_ROOT:-/nfs/production/arl/chembl/tevfik/DEEP_APBS_DATASETS}"
PYTHON_BIN="${PYTHON_BIN:-$REPO_DIR/.conda/bin/python}"
PATH="$REPO_DIR/.conda/bin:$PATH"
export PYTHONPATH="$REPO_DIR:${PYTHONPATH:-}"

TARGET_ROOTS="${TARGET_ROOTS:-$DATA_ROOT/cache/work9_external_benchmark_cache_gridfix_v1 $DATA_ROOT/cache/work11_cache_gridfix_v1}"
GLOB_PATTERN="${GLOB_PATTERN:-*.h5}"
OVERWRITE="${OVERWRITE:-0}"
FULL_SIGNED_SCALE="${FULL_SIGNED_SCALE:-150}"
KEEP_LEAKAGE_FEATURES="${KEEP_LEAKAGE_FEATURES:-0}"
MOVE_LEGACY_ELECTROSTATIC_TO_AUXILIARY="${MOVE_LEGACY_ELECTROSTATIC_TO_AUXILIARY:-0}"

mkdir -p /homes/tevfik/PHD/slurm_logs
cd "$REPO_DIR"

overwrite_args=()
if [[ "$OVERWRITE" == "1" ]]; then
  overwrite_args=(--overwrite)
fi

leakage_args=()
if [[ "$KEEP_LEAKAGE_FEATURES" == "1" ]]; then
  leakage_args=(--keep-leakage-features)
fi

legacy_electrostatic_args=()
if [[ "$MOVE_LEGACY_ELECTROSTATIC_TO_AUXILIARY" == "1" ]]; then
  legacy_electrostatic_args=(--move-legacy-electrostatic-to-auxiliary)
fi

echo "H5 v2 derived feature augmentation"
echo "Host: $(hostname)"
echo "Repo: $REPO_DIR"
echo "Data root: $DATA_ROOT"
echo "Target roots: $TARGET_ROOTS"
echo "Glob: $GLOB_PATTERN"
echo "Overwrite: $OVERWRITE"
echo "Full signed scale: $FULL_SIGNED_SCALE"
echo "Keep leakage features: $KEEP_LEAKAGE_FEATURES"
echo "Move legacy electrostatic to auxiliary: $MOVE_LEGACY_ELECTROSTATIC_TO_AUXILIARY"

for target_root in $TARGET_ROOTS; do
  if [[ ! -d "$target_root" ]]; then
    echo "[SKIP] Missing target root: $target_root"
    continue
  fi

  echo
  echo "============================================================"
  echo "Augmenting: $target_root"
  echo "============================================================"

  "$PYTHON_BIN" augment_h5_v2_features.py \
    --input-dir "$target_root" \
    --glob "$GLOB_PATTERN" \
    --in-place \
    --full-signed-scale "$FULL_SIGNED_SCALE" \
    --summary-csv "$target_root/v2_augment_summary.csv" \
    "${overwrite_args[@]}" \
    "${leakage_args[@]}" \
    "${legacy_electrostatic_args[@]}"
done

echo
echo "Done."
