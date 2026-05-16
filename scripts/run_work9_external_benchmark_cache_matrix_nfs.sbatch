#!/usr/bin/env bash
#SBATCH --job-name=work9-bench-cache
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=2-00:00:00
#SBATCH --output=/homes/tevfik/PHD/slurm_logs/work9-bench-cache-%j.out
#SBATCH --error=/homes/tevfik/PHD/slurm_logs/work9-bench-cache-%j.err

set -euo pipefail

REPO_DIR="${REPO_DIR:-/homes/tevfik/PHD/generate_cache}"
DATA_ROOT="${DATA_ROOT:-/nfs/production/arl/chembl/tevfik/DEEP_APBS_DATASETS}"
PYTHON_BIN="${PYTHON_BIN:-$REPO_DIR/.conda/bin/python}"
PATH="$REPO_DIR/.conda/bin:$PATH"
export PYTHONPATH="$REPO_DIR:${PYTHONPATH:-}"

GRID_SPECS="${GRID_SPECS:-36:70 72:120 161:160}"
NPROC="${NPROC:-8}"
LIMIT="${LIMIT:-all}"
DATASETS="${DATASETS:-coach420_puresnet bu48_puresnet}"
PREPARED_ROOT="${PREPARED_ROOT:-$DATA_ROOT/datasets/external_benchmarks/puresnet_prepared}"
OUTPUT_ROOT="${OUTPUT_ROOT:-$DATA_ROOT/cache/work9_external_benchmark_cache_gridfix_v1}"
OVERWRITE="${OVERWRITE:-0}"

mkdir -p "$OUTPUT_ROOT" /homes/tevfik/PHD/slurm_logs
cd "$REPO_DIR"

limit_args=()
if [[ -n "$LIMIT" && "$LIMIT" != "all" ]]; then
  limit_args=(--limit "$LIMIT")
fi

overwrite_args=()
if [[ "$OVERWRITE" == "1" ]]; then
  overwrite_args=(--overwrite)
fi

echo "Work9 external benchmark cache matrix"
echo "Host: $(hostname)"
echo "Repo: $REPO_DIR"
echo "Data root: $DATA_ROOT"
echo "Prepared root: $PREPARED_ROOT"
echo "Output root: $OUTPUT_ROOT"
echo "Grid specs: $GRID_SPECS"
echo "Datasets: $DATASETS"
echo "NPROC: $NPROC"
echo "LIMIT: $LIMIT"
echo "OVERWRITE: $OVERWRITE"

for grid_spec in $GRID_SPECS; do
  if [[ "$grid_spec" == *":"* ]]; then
    box_size="${grid_spec%%:*}"
    target_span="${grid_spec#*:}"
  else
    box_size="$grid_spec"
    target_span="160"
  fi

  resolution="$(awk "BEGIN { printf \"%.10f\", ${target_span} / (${box_size} - 1) }")"
  log_file="$OUTPUT_ROOT/generation_box${box_size}_span${target_span}.log"

  echo
  echo "============================================================"
  echo "External benchmark cache box${box_size}_span${target_span}"
  echo "Resolution: ${resolution} A/voxel"
  echo "============================================================"

  "$PYTHON_BIN" generate_cache_benchmark_gridfix.py \
    --prepared-root "$PREPARED_ROOT" \
    --output-root "$OUTPUT_ROOT" \
    --datasets $DATASETS \
    --box-size "$box_size" \
    --target-span "$target_span" \
    --resolution "$resolution" \
    "${limit_args[@]}" \
    --nproc "$NPROC" \
    "${overwrite_args[@]}" \
    2>&1 | tee "$log_file"
done

echo
echo "Done."
echo "Output root: $OUTPUT_ROOT"
