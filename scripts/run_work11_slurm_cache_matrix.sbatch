#!/usr/bin/env bash
#SBATCH --job-name=work11-cache
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=5-00:00:00
#SBATCH --output=/homes/tevfik/PHD/slurm_logs/work11-cache-%j.out
#SBATCH --error=/homes/tevfik/PHD/slurm_logs/work11-cache-%j.err

set -euo pipefail

REPO_DIR="${REPO_DIR:-/homes/tevfik/PHD/phd_examples/generate_cache}"
cd "$REPO_DIR"

PYTHON_BIN="${PYTHON_BIN:-$REPO_DIR/.conda/bin/python}"
PATH="$REPO_DIR/.conda/bin:$PATH"

GRID_SPECS="${GRID_SPECS:-}"
BOX_SIZES="${BOX_SIZES:-}"
NPROC="${NPROC:-16}"
LIMIT="${LIMIT:-all}"
TARGET_SPAN="${TARGET_SPAN:-160}"
SCPDB_LABEL_SOURCE="${SCPDB_LABEL_SOURCE:-cavity6}"
SCPDB_ROOT="${SCPDB_ROOT:-/hps/nobackup/arl/chembl/tevfik/deep-apbs-data/scPDB}"
PDBBIND_ROOT="${PDBBIND_ROOT:-/hps/nobackup/arl/chembl/tevfik/deep-apbs-data/pdbbind/refined-set}"
OUTPUT_ROOT="${OUTPUT_ROOT:-/hps/nobackup/arl/chembl/tevfik/deep-apbs-data/work11_cache_gridfix_v1}"
PDBBIND_CASE_LIST="${PDBBIND_CASE_LIST:-pdbbind_gridfix160_full_fits.txt}"

limit_args=()
if [[ -n "$LIMIT" && "$LIMIT" != "all" ]]; then
  limit_args=(--limit "$LIMIT")
fi

if [[ -z "$GRID_SPECS" ]]; then
  if [[ -n "$BOX_SIZES" ]]; then
    for box_size in $BOX_SIZES; do
      GRID_SPECS+="${box_size}:${TARGET_SPAN} "
    done
  else
    GRID_SPECS="36:70 72:120 161:160"
  fi
fi

mkdir -p "$OUTPUT_ROOT" /homes/tevfik/PHD/slurm_logs

echo "Work11 SLURM cache matrix"
echo "Host: $(hostname)"
echo "Repo: $REPO_DIR"
echo "Output root: $OUTPUT_ROOT"
echo "Grid specs: $GRID_SPECS"
echo "NPROC: $NPROC"
echo "LIMIT: $LIMIT"

for grid_spec in $GRID_SPECS; do
  if [[ "$grid_spec" == *":"* ]]; then
    box_size="${grid_spec%%:*}"
    target_span="${grid_spec#*:}"
  else
    box_size="$grid_spec"
    target_span="$TARGET_SPAN"
  fi

  scpdb_out="$OUTPUT_ROOT/scpdb/label_${SCPDB_LABEL_SOURCE}/box${box_size}_span${target_span}"
  pdbbind_out="$OUTPUT_ROOT/pdbbind/refined-set/box${box_size}_span${target_span}"
  mkdir -p "$scpdb_out" "$pdbbind_out"

  echo
  echo "============================================================"
  echo "scPDB box${box_size}_span${target_span}"
  echo "============================================================"
  "$PYTHON_BIN" generate_cache_scpdb_gridfix.py \
    --dataset-root "$SCPDB_ROOT" \
    --output-dir "$scpdb_out" \
    --dataset-label-source "$SCPDB_LABEL_SOURCE" \
    --include-all-labels \
    --box-size "$box_size" \
    --target-span "$target_span" \
    "${limit_args[@]}" \
    --nproc "$NPROC" \
    --log-file "$scpdb_out/generation.log"

  echo
  echo "============================================================"
  echo "PDBBind box${box_size}_span${target_span}"
  echo "============================================================"
  "$PYTHON_BIN" generate_cache_pdbbind_gridfix.py \
    --dataset-root "$PDBBIND_ROOT" \
    --output-dir "$pdbbind_out" \
    --case-list "$PDBBIND_CASE_LIST" \
    --box-size "$box_size" \
    --resolution "$(awk "BEGIN { printf \"%.10f\", ${target_span} / (${box_size} - 1) }")" \
    "${limit_args[@]}" \
    --nproc "$NPROC" \
    2>&1 | tee "$pdbbind_out/generation.log"
done

echo
echo "Done."
echo "Output root: $OUTPUT_ROOT"
