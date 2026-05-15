#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$REPO_DIR"

PYTHON_BIN="${PYTHON_BIN:-$REPO_DIR/.conda/bin/python}"
PATH="$REPO_DIR/.conda/bin:$PATH"

BOX_SIZES="${BOX_SIZES:-36 72 161}"
NPROC="${NPROC:-2}"
LIMIT="${LIMIT:-}"
PREPARED_ROOT="${PREPARED_ROOT:-/Users/tevfik/Sandbox/github/PHD/data/external_benchmarks/puresnet_prepared}"
OUTPUT_ROOT="${OUTPUT_ROOT:-/Users/tevfik/Sandbox/github/PHD/data/external_benchmarks/puresnet_cache_gridfix_v1}"

LIMIT_ARG=()
if [[ "$LIMIT" != "" ]]; then
  LIMIT_ARG=(--limit "$LIMIT")
fi

for box_size in $BOX_SIZES; do
  echo
  echo "============================================================"
  echo "PUResNet benchmark cache box${box_size}"
  echo "============================================================"
  "$PYTHON_BIN" generate_cache_benchmark_gridfix.py \
    --prepared-root "$PREPARED_ROOT" \
    --output-root "$OUTPUT_ROOT" \
    --box-size "$box_size" \
    --nproc "$NPROC" \
    "${LIMIT_ARG[@]}"
done

echo "Done."
