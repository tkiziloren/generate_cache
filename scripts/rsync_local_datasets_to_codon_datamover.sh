#!/usr/bin/env bash

set -euo pipefail

REMOTE_HOST="${REMOTE_HOST:-codon}"
DATA_ROOT="${DATA_ROOT:-/nfs/production/arl/chembl/tevfik/DEEP_APBS_DATASETS}"
LOCAL_PDBBIND="${LOCAL_PDBBIND:-/Users/tevfik/Sandbox/github/PHD/data/pdbbind/refined-set}"
LOCAL_EXTERNAL_BENCHMARKS="${LOCAL_EXTERNAL_BENCHMARKS:-/Users/tevfik/Sandbox/github/PHD/data/external_benchmarks}"
DATAMOVER_MEM="${DATAMOVER_MEM:-8G}"
DATAMOVER_TIME="${DATAMOVER_TIME:-12:00:00}"
DRY_RUN="${DRY_RUN:-0}"
DELETE_REMOTE="${DELETE_REMOTE:-0}"

if [[ ! -d "$LOCAL_PDBBIND" ]]; then
  echo "Missing local PDBBind directory: $LOCAL_PDBBIND" >&2
  exit 1
fi

if [[ ! -d "$LOCAL_EXTERNAL_BENCHMARKS" ]]; then
  echo "Missing local external benchmarks directory: $LOCAL_EXTERNAL_BENCHMARKS" >&2
  exit 1
fi

remote_prepare=$(cat <<EOF
set -euo pipefail
mkdir -p "$DATA_ROOT/datasets/pdbbind"
mkdir -p "$DATA_ROOT/datasets/pdbbind2020"
mkdir -p "$DATA_ROOT/datasets/external_benchmarks"

if [ -L "$DATA_ROOT/datasets/pdbbind/refined-set" ]; then
  rm -f "$DATA_ROOT/datasets/pdbbind/refined-set"
fi
mkdir -p "$DATA_ROOT/datasets/pdbbind/refined-set"

if [ -L "$DATA_ROOT/datasets/pdbbind2020/refined-set" ]; then
  rm -f "$DATA_ROOT/datasets/pdbbind2020/refined-set"
fi
if [ ! -e "$DATA_ROOT/datasets/pdbbind2020/refined-set" ]; then
  ln -s ../pdbbind/refined-set "$DATA_ROOT/datasets/pdbbind2020/refined-set"
fi
EOF
)

echo "Preparing remote NFS directories via datamover..."
remote_prepare_command="srun --partition=datamover --mem=4G --time=00:30:00 bash -lc $(printf '%q' "$remote_prepare")"
ssh "$REMOTE_HOST" "bash -l -c $(printf '%q' "$remote_prepare_command")"

remote_rsync="bash -l -c \"srun --partition=datamover --mem=$DATAMOVER_MEM --time=$DATAMOVER_TIME rsync\""

rsync_args=(
  -avh
  --partial
  --inplace
  --info=progress2
  --exclude ".DS_Store"
  --exclude "__pycache__"
  --exclude ".git"
)

if [[ "$DRY_RUN" == "1" ]]; then
  rsync_args+=(--dry-run)
fi

if [[ "$DELETE_REMOTE" == "1" ]]; then
  rsync_args+=(--delete-after)
fi

echo
echo "Syncing PDBBind refined-set..."
echo "  local : $LOCAL_PDBBIND"
echo "  remote: $REMOTE_HOST:$DATA_ROOT/datasets/pdbbind/refined-set/"
rsync "${rsync_args[@]}" \
  --rsync-path="$remote_rsync" \
  "$LOCAL_PDBBIND"/ \
  "$REMOTE_HOST:$DATA_ROOT/datasets/pdbbind/refined-set/"

echo
echo "Syncing external benchmarks..."
echo "  local : $LOCAL_EXTERNAL_BENCHMARKS"
echo "  remote: $REMOTE_HOST:$DATA_ROOT/datasets/external_benchmarks/"
rsync "${rsync_args[@]}" \
  --rsync-path="$remote_rsync" \
  "$LOCAL_EXTERNAL_BENCHMARKS"/ \
  "$REMOTE_HOST:$DATA_ROOT/datasets/external_benchmarks/"

echo
echo "Done."
echo "Remote PDBBind root: $DATA_ROOT/datasets/pdbbind/refined-set"
echo "Remote PDBBind2020 alias: $DATA_ROOT/datasets/pdbbind2020/refined-set"
echo "Remote external benchmarks root: $DATA_ROOT/datasets/external_benchmarks"
