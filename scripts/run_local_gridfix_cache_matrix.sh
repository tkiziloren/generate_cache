#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$REPO_DIR"

PYTHON_BIN="${PYTHON_BIN:-$REPO_DIR/.conda/bin/python}"
PATH="$REPO_DIR/.conda/bin:$PATH"

BOX_SIZES="${BOX_SIZES:-36 72 161}"
NPROC="${NPROC:-2}"
LIMIT="${LIMIT:-20}"
TARGET_SPAN="${TARGET_SPAN:-160}"
SCPDB_LABEL_SOURCE="${SCPDB_LABEL_SOURCE:-cavity6}"
SCPDB_ROOT="${SCPDB_ROOT:-/Users/tevfik/Sandbox/github/PHD/data/scPDB}"
PDBBIND_ROOT="${PDBBIND_ROOT:-/Users/tevfik/Sandbox/github/PHD/data/pdbbind/refined-set}"
LOCAL_OUTPUT_ROOT="${LOCAL_OUTPUT_ROOT:-/Users/tevfik/Sandbox/github/PHD/data/cache_matrix_local_gridfix_v1}"
PDBBIND_CASE_LIST="${PDBBIND_CASE_LIST:-pdbbind_gridfix160_full_fits.txt}"

mkdir -p "$LOCAL_OUTPUT_ROOT"

echo "Local gridfix cache matrix"
echo "Output root: $LOCAL_OUTPUT_ROOT"
echo "Box sizes: $BOX_SIZES"
echo "Limit per dataset: $LIMIT"

for box_size in $BOX_SIZES; do
  scpdb_out="$LOCAL_OUTPUT_ROOT/scpdb/label_${SCPDB_LABEL_SOURCE}/box${box_size}_span${TARGET_SPAN}"
  pdbbind_out="$LOCAL_OUTPUT_ROOT/pdbbind/refined-set/box${box_size}_span${TARGET_SPAN}"

  echo
  echo "============================================================"
  echo "scPDB box${box_size}"
  echo "============================================================"
  "$PYTHON_BIN" generate_cache_scpdb_gridfix.py \
    --dataset-root "$SCPDB_ROOT" \
    --output-dir "$scpdb_out" \
    --dataset-label-source "$SCPDB_LABEL_SOURCE" \
    --include-all-labels \
    --box-size "$box_size" \
    --target-span "$TARGET_SPAN" \
    --limit "$LIMIT" \
    --nproc "$NPROC" \
    --log-file "$scpdb_out/generation.log"

  echo
  echo "============================================================"
  echo "PDBBind box${box_size}"
  echo "============================================================"
  "$PYTHON_BIN" generate_cache_pdbbind_gridfix.py \
    --dataset-root "$PDBBIND_ROOT" \
    --output-dir "$pdbbind_out" \
    --case-list "$PDBBIND_CASE_LIST" \
    --box-size "$box_size" \
    --limit "$LIMIT" \
    --nproc "$NPROC"
done

echo "Done."
