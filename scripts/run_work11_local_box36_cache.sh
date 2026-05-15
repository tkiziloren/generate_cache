#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$REPO_DIR"

PYTHON_BIN="${PYTHON_BIN:-$REPO_DIR/.conda/bin/python}"
PATH="$REPO_DIR/.conda/bin:$PATH"

BOX_SIZE="${BOX_SIZE:-36}"
TARGET_SPAN="${TARGET_SPAN:-70}"
NPROC="${NPROC:-4}"
LIMIT="${LIMIT:-20}"
SCPDB_LABEL_SOURCE="${SCPDB_LABEL_SOURCE:-cavity6}"
SCPDB_ROOT="${SCPDB_ROOT:-/Users/tevfik/Sandbox/github/PHD/data/scPDB}"
PDBBIND_ROOT="${PDBBIND_ROOT:-/Users/tevfik/Sandbox/github/PHD/data/pdbbind/refined-set}"
LOCAL_OUTPUT_ROOT="${LOCAL_OUTPUT_ROOT:-/Users/tevfik/Sandbox/github/PHD/data/work11_local_cache_gridfix_v1}"
PDBBIND_CASE_LIST="${PDBBIND_CASE_LIST:-pdbbind_gridfix160_full_fits.txt}"

limit_args=()
if [[ -n "$LIMIT" && "$LIMIT" != "all" ]]; then
  limit_args=(--limit "$LIMIT")
fi

scpdb_out="$LOCAL_OUTPUT_ROOT/scpdb/label_${SCPDB_LABEL_SOURCE}/box${BOX_SIZE}_span${TARGET_SPAN}"
pdbbind_out="$LOCAL_OUTPUT_ROOT/pdbbind/refined-set/box${BOX_SIZE}_span${TARGET_SPAN}"

mkdir -p "$scpdb_out" "$pdbbind_out"

echo "Work11 local cache generation"
echo "Repo: $REPO_DIR"
echo "Output root: $LOCAL_OUTPUT_ROOT"
echo "Box size: $BOX_SIZE"
echo "Target span: $TARGET_SPAN"
echo "NPROC: $NPROC"
echo "LIMIT: $LIMIT"
pdbbind_resolution="$(awk "BEGIN { printf \"%.10f\", ${TARGET_SPAN} / (${BOX_SIZE} - 1) }")"

echo
echo "============================================================"
echo "scPDB box${BOX_SIZE}"
echo "============================================================"
"$PYTHON_BIN" generate_cache_scpdb_gridfix.py \
  --dataset-root "$SCPDB_ROOT" \
  --output-dir "$scpdb_out" \
  --dataset-label-source "$SCPDB_LABEL_SOURCE" \
  --include-all-labels \
  --box-size "$BOX_SIZE" \
  --target-span "$TARGET_SPAN" \
  "${limit_args[@]}" \
  --nproc "$NPROC" \
  --log-file "$scpdb_out/generation.log"

echo
echo "============================================================"
echo "PDBBind box${BOX_SIZE}"
echo "============================================================"
"$PYTHON_BIN" generate_cache_pdbbind_gridfix.py \
  --dataset-root "$PDBBIND_ROOT" \
  --output-dir "$pdbbind_out" \
  --case-list "$PDBBIND_CASE_LIST" \
  --box-size "$BOX_SIZE" \
  --resolution "$pdbbind_resolution" \
  "${limit_args[@]}" \
  --nproc "$NPROC" \
  2>&1 | tee "$pdbbind_out/generation.log"

echo
echo "Done."
echo "scPDB log: $scpdb_out/generation.log"
echo "scPDB manifest: $scpdb_out/manifest.csv"
echo "PDBBind log: $pdbbind_out/generation.log"
echo "PDBBind failed cases: $pdbbind_out/failed_cases.txt"
