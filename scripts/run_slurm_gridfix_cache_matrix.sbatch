#!/usr/bin/env bash
#SBATCH --job-name=gridfix-cache-matrix
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --output=/homes/tevfik/PHD/slurm_logs/gridfix-cache-matrix-%j.out
#SBATCH --error=/homes/tevfik/PHD/slurm_logs/gridfix-cache-matrix-%j.err

set -euo pipefail

REPO_DIR="${REPO_DIR:-/homes/tevfik/PHD/phd_examples/generate_cache}"
cd "$REPO_DIR"

PYTHON_BIN="${PYTHON_BIN:-$REPO_DIR/.conda/bin/python}"
PATH="$REPO_DIR/.conda/bin:$PATH"

BOX_SIZES="${BOX_SIZES:-36 72 161}"
NPROC="${NPROC:-16}"
TARGET_SPAN="${TARGET_SPAN:-160}"
SCPDB_LABEL_SOURCE="${SCPDB_LABEL_SOURCE:-cavity6}"
SCPDB_ROOT="${SCPDB_ROOT:-/hps/nobackup/arl/chembl/tevfik/deep-apbs-data/scPDB}"
PDBBIND_ROOT="${PDBBIND_ROOT:-/hps/nobackup/arl/chembl/tevfik/deep-apbs-data/pdbbind/refined-set}"
OUTPUT_ROOT="${OUTPUT_ROOT:-/hps/nobackup/arl/chembl/tevfik/deep-apbs-data/gridfix_cache_matrix_v1}"
PDBBIND_CASE_LIST="${PDBBIND_CASE_LIST:-pdbbind_gridfix160_full_fits.txt}"

mkdir -p "$OUTPUT_ROOT" /homes/tevfik/PHD/slurm_logs

echo "SLURM gridfix cache matrix"
echo "Host: $(hostname)"
echo "Repo: $REPO_DIR"
echo "Output root: $OUTPUT_ROOT"
echo "Box sizes: $BOX_SIZES"

for box_size in $BOX_SIZES; do
  scpdb_out="$OUTPUT_ROOT/scpdb/label_${SCPDB_LABEL_SOURCE}/box${box_size}_span${TARGET_SPAN}"
  pdbbind_out="$OUTPUT_ROOT/pdbbind/refined-set/box${box_size}_span${TARGET_SPAN}"

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
    --nproc "$NPROC"
done

echo "Done."
