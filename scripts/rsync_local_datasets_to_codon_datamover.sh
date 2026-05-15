#!/usr/bin/env bash

set -euo pipefail

REMOTE_HOST="${REMOTE_HOST:-codon}"
DATA_ROOT="${DATA_ROOT:-/nfs/production/arl/chembl/tevfik/DEEP_APBS_DATASETS}"
LOCAL_PDBBIND="${LOCAL_PDBBIND:-/Users/tevfik/Sandbox/github/PHD/data/pdbbind/refined-set}"
LOCAL_EXTERNAL_BENCHMARKS="${LOCAL_EXTERNAL_BENCHMARKS:-/Users/tevfik/Sandbox/github/PHD/data/external_benchmarks}"
DATAMOVER_MEM="${DATAMOVER_MEM:-8G}"
DATAMOVER_TIME="${DATAMOVER_TIME:-12:00:00}"
ONLY="${ONLY:-all}"
DRY_RUN="${DRY_RUN:-0}"
DRY_RUN_SAMPLE="${DRY_RUN_SAMPLE:-1}"
DRY_RUN_RSYNC="${DRY_RUN_RSYNC:-0}"
DELETE_REMOTE="${DELETE_REMOTE:-0}"

case "$ONLY" in
  all|pdbbind|external) ;;
  *)
    echo "ONLY must be one of: all, pdbbind, external" >&2
    exit 1
    ;;
esac

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

remote_rsync=""
if [[ "$DRY_RUN" != "1" || "$DRY_RUN_RSYNC" == "1" ]]; then
  remote_srun="$(ssh "$REMOTE_HOST" "bash -l -c 'command -v srun'")"
  if [[ -z "$remote_srun" ]]; then
    echo "Could not find srun on remote host: $REMOTE_HOST" >&2
    exit 1
  fi
  remote_rsync="$remote_srun --partition=datamover --mem=$DATAMOVER_MEM --time=$DATAMOVER_TIME rsync"
fi

rsync_args=(
  -avh
  --partial
  --inplace
  --progress
  --stats
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

sync_directory() {
  local label="$1"
  local local_root="$2"
  local remote_root="$3"
  local entry_file=""
  local cleanup_entry_file=0
  local extra_args=()
  local sample_count=0

  if [[ "$DRY_RUN" == "1" && "$DRY_RUN_SAMPLE" == "1" ]]; then
    entry_file="$(mktemp)"
    cleanup_entry_file=1
    find "$local_root" -mindepth 1 -maxdepth 1 -print | sort > "$entry_file"
    while IFS= read -r path; do
      name="$(basename "$path")"
      if [[ -d "$path" ]]; then
        extra_args+=(--include="/$name/" --include="/$name/**")
      else
        extra_args+=(--include="/$name")
      fi
      sample_count=$((sample_count + 1))
      if [[ "$sample_count" -ge 5 ]]; then
        break
      fi
    done < "$entry_file"
    extra_args+=(--exclude="*")
  fi

  echo
  if [[ "$DRY_RUN" == "1" && "$DRY_RUN_RSYNC" != "1" ]]; then
    echo "Dry-run check for $label..."
  else
    echo "Syncing $label..."
  fi
  echo "  local : $local_root"
  echo "  remote: $REMOTE_HOST:$remote_root/"
  if [[ "$DRY_RUN" == "1" && "$DRY_RUN_RSYNC" != "1" ]]; then
    echo "  mode  : dry-run metadata only; no remote rsync server is started"
    echo "  sample: first 5 top-level local entries"
    find "$local_root" -mindepth 1 -maxdepth 1 -print | sort | sed -n '1,5p'
    echo "  note  : set DRY_RUN_RSYNC=1 only if you want a full rsync dry-run"
    if [[ "$cleanup_entry_file" == "1" ]]; then
      rm -f "$entry_file"
    fi
    return 0
  elif [[ "$DRY_RUN" == "1" && "$DRY_RUN_SAMPLE" == "1" ]]; then
    echo "  mode  : rsync dry-run sample only, first 5 top-level entries"
  fi

  rsync "${rsync_args[@]}" "${extra_args[@]}" \
    --rsync-path="$remote_rsync" \
    "$local_root"/ \
    "$REMOTE_HOST:$remote_root/"

  if [[ "$cleanup_entry_file" == "1" ]]; then
    rm -f "$entry_file"
  fi
}

if [[ "$ONLY" == "all" || "$ONLY" == "pdbbind" ]]; then
  sync_directory \
    "PDBBind refined-set" \
    "$LOCAL_PDBBIND" \
    "$DATA_ROOT/datasets/pdbbind/refined-set"
fi

if [[ "$ONLY" == "all" || "$ONLY" == "external" ]]; then
  sync_directory \
    "external benchmarks" \
    "$LOCAL_EXTERNAL_BENCHMARKS" \
    "$DATA_ROOT/datasets/external_benchmarks"
fi

echo
echo "Done."
echo "Remote PDBBind root: $DATA_ROOT/datasets/pdbbind/refined-set"
echo "Remote PDBBind2020 alias: $DATA_ROOT/datasets/pdbbind2020/refined-set"
echo "Remote external benchmarks root: $DATA_ROOT/datasets/external_benchmarks"
