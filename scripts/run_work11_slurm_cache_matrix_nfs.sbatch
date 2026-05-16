#!/usr/bin/env bash
#SBATCH --job-name=work11-cache-nfs
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=5-00:00:00
#SBATCH --output=/homes/tevfik/PHD/slurm_logs/work11-cache-nfs-%j.out
#SBATCH --error=/homes/tevfik/PHD/slurm_logs/work11-cache-nfs-%j.err

set -euo pipefail

REPO_DIR="${REPO_DIR:-/homes/tevfik/PHD/generate_cache}"
DATA_ROOT="${DATA_ROOT:-/nfs/production/arl/chembl/tevfik/DEEP_APBS_DATASETS}"
PYTHON_BIN="${PYTHON_BIN:-$REPO_DIR/.conda/bin/python}"
PATH="$REPO_DIR/.conda/bin:$PATH"
export PYTHONPATH="$REPO_DIR:${PYTHONPATH:-}"

GRID_SPECS="${GRID_SPECS:-}"
BOX_SIZES="${BOX_SIZES:-}"
TARGET_SPAN="${TARGET_SPAN:-160}"
NPROC="${NPROC:-16}"
LIMIT="${LIMIT:-all}"
RUN_SCPDB="${RUN_SCPDB:-1}"
RUN_PDBBIND="${RUN_PDBBIND:-1}"
SCPDB_LABEL_SOURCE="${SCPDB_LABEL_SOURCE:-cavity6}"
SCPDB_ROOT="${SCPDB_ROOT:-$DATA_ROOT/datasets/scPDB}"
PDBBIND_ROOT="${PDBBIND_ROOT:-$DATA_ROOT/datasets/pdbbind/refined-set}"
OUTPUT_ROOT="${OUTPUT_ROOT:-$DATA_ROOT/cache/work11_cache_gridfix_v1}"
PDBBIND_CASE_LIST="${PDBBIND_CASE_LIST:-pdbbind_gridfix160_full_fits.txt}"

mkdir -p "$OUTPUT_ROOT" /homes/tevfik/PHD/slurm_logs
cd "$REPO_DIR"

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

echo "Work11 NFS cache matrix"
echo "Host: $(hostname)"
echo "Repo: $REPO_DIR"
echo "Data root: $DATA_ROOT"
echo "Output root: $OUTPUT_ROOT"
echo "Grid specs: $GRID_SPECS"
echo "NPROC: $NPROC"
echo "LIMIT: $LIMIT"
echo "RUN_SCPDB: $RUN_SCPDB"
echo "RUN_PDBBIND: $RUN_PDBBIND"
echo "scPDB root: $SCPDB_ROOT"
echo "PDBBind root: $PDBBIND_ROOT"

for grid_spec in $GRID_SPECS; do
  if [[ "$grid_spec" == *":"* ]]; then
    box_size="${grid_spec%%:*}"
    target_span="${grid_spec#*:}"
  else
    box_size="$grid_spec"
    target_span="$TARGET_SPAN"
  fi

  if [[ "$RUN_SCPDB" == "1" ]]; then
    scpdb_out="$OUTPUT_ROOT/scpdb/label_${SCPDB_LABEL_SOURCE}/box${box_size}_span${target_span}"
    mkdir -p "$scpdb_out"

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
  fi

  if [[ "$RUN_PDBBIND" == "1" ]]; then
    pdbbind_out="$OUTPUT_ROOT/pdbbind/refined-set/box${box_size}_span${target_span}"
    mkdir -p "$pdbbind_out"

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
  fi
done

echo
echo "Done."
echo "Output root: $OUTPUT_ROOT"
