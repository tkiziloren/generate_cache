#!/usr/bin/env bash

set -euo pipefail

REMOTE_HOST="${REMOTE_HOST:-codon}"
REMOTE_REPO_DIR="${REMOTE_REPO_DIR:-/homes/tevfik/PHD/generate_cache}"
DATA_ROOT="${DATA_ROOT:-/nfs/production/arl/chembl/tevfik/DEEP_APBS_DATASETS}"
REMOTE_INCOMING="${REMOTE_INCOMING:-$DATA_ROOT/archives/uploads}"
LOCAL_PDBBIND_ZIP="${LOCAL_PDBBIND_ZIP:-/Users/tevfik/Sandbox/github/PHD/data/pdbbind/refined-set.zip}"
LOCAL_EXTERNAL_BENCHMARKS_ZIP="${LOCAL_EXTERNAL_BENCHMARKS_ZIP:-/Users/tevfik/Sandbox/github/PHD/data/external_benchmarks.zip}"
DATAMOVER_MEM="${DATAMOVER_MEM:-16G}"
DATAMOVER_TIME="${DATAMOVER_TIME:-4:00:00}"
ONLY="${ONLY:-all}"
DRY_RUN="${DRY_RUN:-0}"
EXTRACT_ONLY="${EXTRACT_ONLY:-0}"
TRANSFER_ONLY="${TRANSFER_ONLY:-0}"
KEEP_REMOTE_ZIP="${KEEP_REMOTE_ZIP:-0}"
CLEAN_STAGING="${CLEAN_STAGING:-1}"
UPLOAD_MODE="${UPLOAD_MODE:-stream_datamover}"

case "$ONLY" in
  all|pdbbind|external) ;;
  *)
    echo "ONLY must be one of: all, pdbbind, external" >&2
    exit 1
    ;;
esac

validate_zip() {
  local archive="$1"
  local label="$2"

  if [[ ! -f "$archive" ]]; then
    echo "Missing $label archive: $archive" >&2
    exit 1
  fi
  if [[ ! -s "$archive" ]]; then
    echo "$label archive is empty: $archive" >&2
    echo "Wait until the zip command finishes, or recreate the archive." >&2
    exit 1
  fi
  if ! unzip -l "$archive" >/dev/null; then
    echo "$label archive is not a valid zip file: $archive" >&2
    exit 1
  fi
}

remote_login() {
  local command="$1"
  ssh "$REMOTE_HOST" "bash -l -c $(printf '%q' "$command")"
}

remote_datamover() {
  local command="$1"
  local wrapped
  wrapped="srun --partition=datamover --mem=$DATAMOVER_MEM --time=$DATAMOVER_TIME bash -lc $(printf '%q' "$command")"
  remote_login "$wrapped"
}

stream_archive_to_datamover() {
  local archive="$1"
  local label="$2"
  local remote_archive="$3"
  local command

  command=$(cat <<EOF
set -euo pipefail
remote_archive=$(printf '%q' "$remote_archive")
mkdir -p "\$(dirname "\$remote_archive")"
tmp="\$remote_archive.part.\$\$"
echo "[STREAM] stdin -> \$remote_archive"
cat > "\$tmp"
mv "\$tmp" "\$remote_archive"
ls -lh "\$remote_archive"
EOF
)

  echo
  echo "Streaming $label archive directly to NFS via datamover..."
  echo "  local : $archive"
  echo "  remote: $REMOTE_HOST:$remote_archive"
  if command -v pv >/dev/null 2>&1; then
    pv "$archive" | remote_datamover "$command"
  else
    echo "  note  : install pv locally for progress; using cat without progress"
    cat "$archive" | remote_datamover "$command"
  fi
}

transfer_archive() {
  local archive="$1"
  local label="$2"
  local remote_archive="$3"

  if [[ "$UPLOAD_MODE" == "stream_datamover" ]]; then
    stream_archive_to_datamover "$archive" "$label" "$remote_archive"
    return
  fi

  if [[ "$UPLOAD_MODE" != "home_rsync" ]]; then
    echo "UPLOAD_MODE must be one of: stream_datamover, home_rsync" >&2
    exit 1
  fi

  echo
  echo "Uploading $label archive to Codon home..."
  echo "  local : $archive"
  echo "  remote: $REMOTE_HOST:$REMOTE_INCOMING/"
  remote_login "mkdir -p $(printf '%q' "$REMOTE_INCOMING")"
  rsync -avh --progress --partial "$archive" "$REMOTE_HOST:$REMOTE_INCOMING/"
}

extract_command() {
  local label="$1"
  local remote_archive="$2"
  local target="$3"

  cat <<EOF
set -euo pipefail

data_root=$(printf '%q' "$DATA_ROOT")
label=$(printf '%q' "$label")
archive=$(printf '%q' "$remote_archive")
target=$(printf '%q' "$target")
clean_staging=$(printf '%q' "$CLEAN_STAGING")
keep_remote_zip=$(printf '%q' "$KEEP_REMOTE_ZIP")
remote_repo_dir=$(printf '%q' "$REMOTE_REPO_DIR")

if [ ! -s "\$archive" ]; then
  echo "Missing or empty remote archive: \$archive" >&2
  exit 1
fi

tmp="\$data_root/staging/upload_extract_\${label}_\$(date +%Y%m%d_%H%M%S)_\$\$"
mkdir -p "\$tmp" "\$(dirname "\$target")"

echo "[EXTRACT] \$archive -> \$tmp"
if command -v unzip >/dev/null 2>&1; then
  unzip -q "\$archive" -d "\$tmp"
elif [ -x "\$remote_repo_dir/.conda/bin/python" ]; then
  "\$remote_repo_dir/.conda/bin/python" -m zipfile -e "\$archive" "\$tmp"
elif command -v python3 >/dev/null 2>&1; then
  python3 -m zipfile -e "\$archive" "\$tmp"
else
  echo "Could not find unzip or python to extract: \$archive" >&2
  exit 1
fi

payload=""
if [ "\$label" = "pdbbind" ]; then
  for candidate in \
    "\$tmp/refined-set" \
    "\$tmp/pdbbind/refined-set" \
    "\$tmp/PDBBind/refined-set" \
    "\$tmp/data/pdbbind/refined-set"; do
    if [ -d "\$candidate" ]; then
      payload="\$candidate"
      break
    fi
  done
else
  for candidate in \
    "\$tmp/external_benchmarks" \
    "\$tmp/data/external_benchmarks"; do
    if [ -d "\$candidate" ]; then
      payload="\$candidate"
      break
    fi
  done
fi

if [ -z "\$payload" ]; then
  payload="\$tmp"
fi

echo "[PAYLOAD] \$payload"
mkdir -p "\$target"
if command -v rsync >/dev/null 2>&1; then
  rsync -a --delete --exclude ".DS_Store" --exclude "__MACOSX" "\$payload"/ "\$target"/
else
  find "\$target" -mindepth 1 -maxdepth 1 -exec rm -rf {} +
  cp -a "\$payload"/. "\$target"/
  find "\$target" -name ".DS_Store" -delete
  find "\$target" -name "__MACOSX" -type d -prune -exec rm -rf {} +
fi

if [ "\$label" = "pdbbind" ]; then
  alias_dir="\$data_root/datasets/pdbbind2020/refined-set"
  mkdir -p "\$data_root/datasets/pdbbind2020"
  if [ -L "\$alias_dir" ]; then
    rm -f "\$alias_dir"
  elif [ -d "\$alias_dir" ]; then
    if [ -z "\$(find "\$alias_dir" -mindepth 1 -maxdepth 1 -print -quit)" ]; then
      rmdir "\$alias_dir"
    fi
  fi
  if [ ! -e "\$alias_dir" ]; then
    ln -s ../pdbbind/refined-set "\$alias_dir"
  fi
fi

echo "[CHECK] \$target"
du -sh "\$target" || true
find "\$target" -mindepth 1 -maxdepth 1 | wc -l || true

if [ "\$clean_staging" = "1" ]; then
  rm -rf "\$tmp"
fi
if [ "\$keep_remote_zip" != "1" ]; then
  rm -f "\$archive"
fi
EOF
}

process_dataset() {
  local label="$1"
  local local_archive="$2"
  local target="$3"
  local remote_archive="$REMOTE_INCOMING/$(basename "$local_archive")"

  validate_zip "$local_archive" "$label"

  echo
  echo "Dataset: $label"
  echo "  local archive : $local_archive"
  echo "  remote archive: $REMOTE_HOST:$remote_archive"
  echo "  target        : $target"
  du -sh "$local_archive" || true

  if [[ "$DRY_RUN" == "1" ]]; then
    echo "  mode          : dry-run metadata only"
    return
  fi

  if [[ "$EXTRACT_ONLY" != "1" ]]; then
    transfer_archive "$local_archive" "$label" "$remote_archive"
  fi

  if [[ "$TRANSFER_ONLY" != "1" ]]; then
    echo
    echo "Extracting $label archive on datamover..."
    remote_datamover "$(extract_command "$label" "$remote_archive" "$target")"
  fi
}

if [[ "$ONLY" == "all" || "$ONLY" == "pdbbind" ]]; then
  process_dataset \
    "pdbbind" \
    "$LOCAL_PDBBIND_ZIP" \
    "$DATA_ROOT/datasets/pdbbind/refined-set"
fi

if [[ "$ONLY" == "all" || "$ONLY" == "external" ]]; then
  process_dataset \
    "external" \
    "$LOCAL_EXTERNAL_BENCHMARKS_ZIP" \
    "$DATA_ROOT/datasets/external_benchmarks"
fi

echo
echo "Done."
