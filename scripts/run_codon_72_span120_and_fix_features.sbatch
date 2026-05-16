#!/usr/bin/env bash
#SBATCH --job-name=cache72-fix
#SBATCH --cpus-per-task=32
#SBATCH --mem=160G
#SBATCH --time=5-00:00:00
#SBATCH --output=/homes/tevfik/PHD/slurm_logs/cache72-fix-%j.out
#SBATCH --error=/homes/tevfik/PHD/slurm_logs/cache72-fix-%j.err

set -euo pipefail

REPO_DIR="${REPO_DIR:-/homes/tevfik/PHD/generate_cache}"
DATA_ROOT="${DATA_ROOT:-/nfs/production/arl/chembl/tevfik/DEEP_APBS_DATASETS}"
PYTHON_BIN="${PYTHON_BIN:-$REPO_DIR/.conda/bin/python}"
PATH="$REPO_DIR/.conda/bin:$PATH"
export PYTHONPATH="$REPO_DIR:${PYTHONPATH:-}"

BOX_SIZE="${BOX_SIZE:-72}"
TARGET_SPAN="${TARGET_SPAN:-120}"
NPROC="${NPROC:-${SLURM_CPUS_PER_TASK:-32}}"
EXTERNAL_NPROC="${EXTERNAL_NPROC:-$NPROC}"
LIMIT="${LIMIT:-all}"
RUN_SCPDB_72="${RUN_SCPDB_72:-1}"
RUN_PDBBIND_72="${RUN_PDBBIND_72:-1}"
RUN_EXTERNAL_72="${RUN_EXTERNAL_72:-1}"
FIX_FEATURES="${FIX_FEATURES:-1}"
FIX_EXISTING_PDBBIND="${FIX_EXISTING_PDBBIND:-1}"
FIX_EXISTING_EXTERNAL="${FIX_EXISTING_EXTERNAL:-1}"
FIX_EXISTING_SCPDB="${FIX_EXISTING_SCPDB:-0}"
OVERWRITE_EXTERNAL="${OVERWRITE_EXTERNAL:-0}"
OVERWRITE_FEATURE_FIX="${OVERWRITE_FEATURE_FIX:-0}"
DROP_LEGACY_ELECTROSTATIC="${DROP_LEGACY_ELECTROSTATIC:-1}"
APPEND_FULL_PROTEIN_APBS_V2="${APPEND_FULL_PROTEIN_APBS_V2:-1}"

SCPDB_LABEL_SOURCE="${SCPDB_LABEL_SOURCE:-cavity6}"
SCPDB_ROOT="${SCPDB_ROOT:-$DATA_ROOT/datasets/scPDB}"
PDBBIND_ROOT="${PDBBIND_ROOT:-$DATA_ROOT/datasets/pdbbind/refined-set}"
PREPARED_ROOT="${PREPARED_ROOT:-$DATA_ROOT/datasets/external_benchmarks/puresnet_prepared}"
WORK11_ROOT="${WORK11_ROOT:-$DATA_ROOT/cache/work11_cache_gridfix_v1}"
WORK9_ROOT="${WORK9_ROOT:-$DATA_ROOT/cache/work9_external_benchmark_cache_gridfix_v1}"
PDBBIND_CASE_LIST="${PDBBIND_CASE_LIST:-pdbbind_gridfix160_full_fits.txt}"
EXTERNAL_DATASETS="${EXTERNAL_DATASETS:-coach420_puresnet bu48_puresnet}"

mkdir -p /homes/tevfik/PHD/slurm_logs "$WORK11_ROOT" "$WORK9_ROOT"
cd "$REPO_DIR"

span_label="$TARGET_SPAN"
resolution="$(awk "BEGIN { printf \"%.10f\", ${TARGET_SPAN} / (${BOX_SIZE} - 1) }")"

limit_args=()
if [[ -n "$LIMIT" && "$LIMIT" != "all" ]]; then
  limit_args=(--limit "$LIMIT")
fi

external_overwrite_args=()
if [[ "$OVERWRITE_EXTERNAL" == "1" ]]; then
  external_overwrite_args=(--overwrite)
fi

feature_fix_overwrite_args=()
if [[ "$OVERWRITE_FEATURE_FIX" == "1" ]]; then
  feature_fix_overwrite_args=(--overwrite)
fi

drop_legacy_apbs_args=()
if [[ "$DROP_LEGACY_ELECTROSTATIC" == "1" ]]; then
  drop_legacy_apbs_args=(--drop-legacy-electrostatic)
fi

echo "Codon 72:${TARGET_SPAN} cache generation and feature-schema fix"
echo "Host: $(hostname)"
echo "Repo: $REPO_DIR"
echo "Data root: $DATA_ROOT"
echo "Box size: $BOX_SIZE"
echo "Target span: $TARGET_SPAN A"
echo "Resolution: $resolution A/voxel"
echo "RUN_SCPDB_72: $RUN_SCPDB_72"
echo "RUN_PDBBIND_72: $RUN_PDBBIND_72"
echo "RUN_EXTERNAL_72: $RUN_EXTERNAL_72"
echo "FIX_FEATURES: $FIX_FEATURES"
echo "FIX_EXISTING_PDBBIND: $FIX_EXISTING_PDBBIND"
echo "FIX_EXISTING_EXTERNAL: $FIX_EXISTING_EXTERNAL"
echo "FIX_EXISTING_SCPDB: $FIX_EXISTING_SCPDB"
echo "APPEND_FULL_PROTEIN_APBS_V2: $APPEND_FULL_PROTEIN_APBS_V2"
echo "DROP_LEGACY_ELECTROSTATIC: $DROP_LEGACY_ELECTROSTATIC"

TOTAL_STAGES=1
if [[ "$RUN_SCPDB_72" == "1" ]]; then TOTAL_STAGES=$((TOTAL_STAGES + 1)); fi
if [[ "$RUN_PDBBIND_72" == "1" ]]; then TOTAL_STAGES=$((TOTAL_STAGES + 1)); fi
if [[ "$RUN_EXTERNAL_72" == "1" ]]; then TOTAL_STAGES=$((TOTAL_STAGES + 1)); fi
if [[ "$FIX_FEATURES" == "1" ]]; then TOTAL_STAGES=$((TOTAL_STAGES + 1)); fi
if [[ "$APPEND_FULL_PROTEIN_APBS_V2" == "1" ]]; then TOTAL_STAGES=$((TOTAL_STAGES + 1)); fi
CURRENT_STAGE=0

progress_stage() {
  CURRENT_STAGE=$((CURRENT_STAGE + 1))
  echo
  echo "############################################################"
  echo "[STAGE ${CURRENT_STAGE}/${TOTAL_STAGES}] $1"
  echo "############################################################"
}

canonicalize_cache_dir() {
  local old_dir="$1"
  local new_dir="$2"

  if [[ "$old_dir" == "$new_dir" || ! -d "$old_dir" ]]; then
    return
  fi

  if [[ ! -d "$new_dir" ]]; then
    mkdir -p "$(dirname "$new_dir")"
    mv "$old_dir" "$new_dir"
    echo "[RENAMED] $old_dir -> $new_dir"
    return
  fi

  echo "[MERGE] $old_dir -> $new_dir"
  shopt -s nullglob
  for file in "$old_dir"/*.h5 "$old_dir"/*.csv "$old_dir"/*.txt "$old_dir"/*.log; do
    target="$new_dir/$(basename "$file")"
    if [[ -e "$target" ]]; then
      echo "[KEEP] Target exists, leaving source file untouched: $file"
    else
      mv "$file" "$target"
      echo "[MOVED] $file -> $target"
    fi
  done
  shopt -u nullglob

  if [[ -d "$old_dir/_tmp" ]]; then
    echo "[KEEP] Temporary directory left in old path: $old_dir/_tmp"
  fi

  if rmdir "$old_dir" 2>/dev/null; then
    echo "[REMOVED] Empty old directory: $old_dir"
  else
    echo "[KEEP] Old directory still has files: $old_dir"
  fi
}

canonicalize_cache_file() {
  local old_file="$1"
  local new_file="$2"

  if [[ "$old_file" == "$new_file" || ! -f "$old_file" ]]; then
    return
  fi

  if [[ -f "$new_file" ]]; then
    echo "[KEEP] Target file exists, leaving source file untouched: $old_file"
  else
    mv "$old_file" "$new_file"
    echo "[RENAMED] $old_file -> $new_file"
  fi
}

canonicalize_existing_dirs() {
  if [[ "$FIX_EXISTING_PDBBIND" == "1" ]]; then
    canonicalize_cache_dir "$WORK11_ROOT/pdbbind/refined-set/box36" "$WORK11_ROOT/pdbbind/refined-set/box36_span70"
    canonicalize_cache_dir "$WORK11_ROOT/pdbbind/refined-set/box72" "$WORK11_ROOT/pdbbind/refined-set/box72_span160"
    canonicalize_cache_dir "$WORK11_ROOT/pdbbind/refined-set/box161" "$WORK11_ROOT/pdbbind/refined-set/box161_span160"
  fi

  if [[ "$FIX_EXISTING_EXTERNAL" == "1" ]]; then
    for dataset in $EXTERNAL_DATASETS; do
      canonicalize_cache_dir "$WORK9_ROOT/$dataset/box36" "$WORK9_ROOT/$dataset/box36_span70"
      canonicalize_cache_dir "$WORK9_ROOT/$dataset/box72" "$WORK9_ROOT/$dataset/box72_span160"
      canonicalize_cache_dir "$WORK9_ROOT/$dataset/box161" "$WORK9_ROOT/$dataset/box161_span160"
    done
    canonicalize_cache_file "$WORK9_ROOT/manifest_box36.csv" "$WORK9_ROOT/manifest_box36_span70.csv"
    canonicalize_cache_file "$WORK9_ROOT/manifest_box72.csv" "$WORK9_ROOT/manifest_box72_span160.csv"
    canonicalize_cache_file "$WORK9_ROOT/manifest_box161.csv" "$WORK9_ROOT/manifest_box161_span160.csv"
  fi
}

progress_stage "Normalize existing cache directory names"
canonicalize_existing_dirs

if [[ "$RUN_SCPDB_72" == "1" ]]; then
  progress_stage "Generate scPDB box${BOX_SIZE}_span${span_label}"
  scpdb_out="$WORK11_ROOT/scpdb/label_${SCPDB_LABEL_SOURCE}/box${BOX_SIZE}_span${span_label}"
  mkdir -p "$scpdb_out"
  echo
  echo "============================================================"
  echo "Generate scPDB box${BOX_SIZE}_span${span_label}"
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
fi

if [[ "$RUN_PDBBIND_72" == "1" ]]; then
  progress_stage "Generate PDBBind box${BOX_SIZE}_span${span_label}"
  pdbbind_out="$WORK11_ROOT/pdbbind/refined-set/box${BOX_SIZE}_span${span_label}"
  mkdir -p "$pdbbind_out"
  echo
  echo "============================================================"
  echo "Generate PDBBind box${BOX_SIZE}_span${span_label}"
  echo "============================================================"
  "$PYTHON_BIN" generate_cache_pdbbind_gridfix.py \
    --dataset-root "$PDBBIND_ROOT" \
    --output-dir "$pdbbind_out" \
    --case-list "$PDBBIND_CASE_LIST" \
    --box-size "$BOX_SIZE" \
    --resolution "$resolution" \
    "${limit_args[@]}" \
    --nproc "$NPROC" \
    2>&1 | tee "$pdbbind_out/generation.log"
fi

if [[ "$RUN_EXTERNAL_72" == "1" ]]; then
  progress_stage "Generate external benchmarks box${BOX_SIZE}_span${span_label}"
  echo
  echo "============================================================"
  echo "Generate external benchmarks box${BOX_SIZE}_span${span_label}"
  echo "============================================================"
  "$PYTHON_BIN" generate_cache_benchmark_gridfix.py \
    --prepared-root "$PREPARED_ROOT" \
    --output-root "$WORK9_ROOT" \
    --datasets $EXTERNAL_DATASETS \
    --box-size "$BOX_SIZE" \
    --target-span "$TARGET_SPAN" \
    --resolution "$resolution" \
    "${limit_args[@]}" \
    --nproc "$EXTERNAL_NPROC" \
    "${external_overwrite_args[@]}" \
    2>&1 | tee "$WORK9_ROOT/generation_box${BOX_SIZE}_span${span_label}.log"
fi

augment_dir() {
  local dir="$1"
  if [[ ! -d "$dir" ]]; then
    echo "[SKIP] Missing directory: $dir"
    return
  fi

  local count
  count="$(find "$dir" -maxdepth 1 -name "*.h5" | wc -l | awk '{print $1}')"
  if [[ "$count" == "0" ]]; then
    echo "[SKIP] No H5 files: $dir"
    return
  fi

  echo
  echo "============================================================"
  echo "Fix feature schema: $dir"
  echo "H5 count: $count"
  echo "============================================================"
  "$PYTHON_BIN" augment_h5_v2_features.py \
    --input-dir "$dir" \
    --glob "*.h5" \
    --in-place \
    --nproc "$NPROC" \
    --summary-csv "$dir/feature_schema_update.csv" \
    "${feature_fix_overwrite_args[@]}" \
    "${drop_legacy_apbs_args[@]}"
}

append_full_protein_apbs_v2_dir() {
  local dir="$1"
  local source_type="$2"
  if [[ "$APPEND_FULL_PROTEIN_APBS_V2" != "1" ]]; then
    return
  fi
  if [[ ! -d "$dir" ]]; then
    echo "[SKIP] Missing directory for full-protein APBS v2: $dir"
    return
  fi

  local count
  count="$(find "$dir" -maxdepth 1 -name "*.h5" | wc -l | awk '{print $1}')"
  if [[ "$count" == "0" ]]; then
    echo "[SKIP] No H5 files for full-protein APBS v2: $dir"
    return
  fi

  echo
  echo "============================================================"
  echo "Append full-protein APBS v2: $dir"
  echo "Source type: $source_type"
  echo "H5 count: $count"
  echo "============================================================"
  "$PYTHON_BIN" append_full_protein_apbs_v2.py \
    --input-dir "$dir" \
    --glob "*.h5" \
    --source-type "$source_type" \
    --pdbbind-root "$PDBBIND_ROOT" \
    --external-prepared-root "$PREPARED_ROOT" \
    --nproc "$NPROC" \
    --summary-csv "$dir/v2_full_protein_apbs_summary.csv"
}

if [[ "$FIX_FEATURES" == "1" ]]; then
  progress_stage "Fix feature schema in selected cache directories"
  echo
  echo "============================================================"
  echo "Feature schema fix"
  echo "============================================================"

  if [[ "$RUN_SCPDB_72" == "1" || "$FIX_EXISTING_SCPDB" == "1" ]]; then
    augment_dir "$WORK11_ROOT/scpdb/label_${SCPDB_LABEL_SOURCE}/box${BOX_SIZE}_span${span_label}"
  fi

  if [[ "$RUN_PDBBIND_72" == "1" || "$FIX_EXISTING_PDBBIND" == "1" ]]; then
    augment_dir "$WORK11_ROOT/pdbbind/refined-set/box${BOX_SIZE}_span${span_label}"
  fi

  if [[ "$FIX_EXISTING_PDBBIND" == "1" ]]; then
    augment_dir "$WORK11_ROOT/pdbbind/refined-set/box36_span70"
    augment_dir "$WORK11_ROOT/pdbbind/refined-set/box72_span160"
    augment_dir "$WORK11_ROOT/pdbbind/refined-set/box161_span160"
  fi

  if [[ "$FIX_EXISTING_EXTERNAL" == "1" ]]; then
    for dataset in $EXTERNAL_DATASETS; do
      augment_dir "$WORK9_ROOT/$dataset/box36_span70"
      augment_dir "$WORK9_ROOT/$dataset/box72_span160"
      augment_dir "$WORK9_ROOT/$dataset/box${BOX_SIZE}_span${span_label}"
      augment_dir "$WORK9_ROOT/$dataset/box161_span160"
    done
  fi

  if [[ "$FIX_EXISTING_SCPDB" == "1" ]]; then
    augment_dir "$WORK11_ROOT/scpdb/label_${SCPDB_LABEL_SOURCE}/box36_span70"
    augment_dir "$WORK11_ROOT/scpdb/label_${SCPDB_LABEL_SOURCE}/box72_span160"
    augment_dir "$WORK11_ROOT/scpdb/label_${SCPDB_LABEL_SOURCE}/box161_span160"
  fi
fi

if [[ "$APPEND_FULL_PROTEIN_APBS_V2" == "1" ]]; then
  progress_stage "Compute and append full-protein APBS v2 for selected cache directories"

  if [[ "$RUN_SCPDB_72" == "1" || "$FIX_EXISTING_SCPDB" == "1" ]]; then
    echo "[SKIP] scPDB full-protein APBS v2 append is intentionally not wired in this job."
  fi

  if [[ "$RUN_PDBBIND_72" == "1" || "$FIX_EXISTING_PDBBIND" == "1" ]]; then
    append_full_protein_apbs_v2_dir "$WORK11_ROOT/pdbbind/refined-set/box${BOX_SIZE}_span${span_label}" "pdbbind"
  fi

  if [[ "$FIX_EXISTING_PDBBIND" == "1" ]]; then
    append_full_protein_apbs_v2_dir "$WORK11_ROOT/pdbbind/refined-set/box36_span70" "pdbbind"
    append_full_protein_apbs_v2_dir "$WORK11_ROOT/pdbbind/refined-set/box72_span160" "pdbbind"
    append_full_protein_apbs_v2_dir "$WORK11_ROOT/pdbbind/refined-set/box161_span160" "pdbbind"
  fi

  if [[ "$FIX_EXISTING_EXTERNAL" == "1" ]]; then
    for dataset in $EXTERNAL_DATASETS; do
      append_full_protein_apbs_v2_dir "$WORK9_ROOT/$dataset/box36_span70" "external"
      append_full_protein_apbs_v2_dir "$WORK9_ROOT/$dataset/box72_span160" "external"
      append_full_protein_apbs_v2_dir "$WORK9_ROOT/$dataset/box${BOX_SIZE}_span${span_label}" "external"
      append_full_protein_apbs_v2_dir "$WORK9_ROOT/$dataset/box161_span160" "external"
    done
  fi
fi

echo
echo "Done."
echo "Work11 root: $WORK11_ROOT"
echo "Work9 root: $WORK9_ROOT"
