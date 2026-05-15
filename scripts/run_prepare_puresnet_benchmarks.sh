#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$REPO_DIR"

PYTHON_BIN="${PYTHON_BIN:-$REPO_DIR/.conda/bin/python}"
PATH="$REPO_DIR/.conda/bin:$PATH"

OUTPUT_ROOT="${OUTPUT_ROOT:-/Users/tevfik/Sandbox/github/PHD/data/external_benchmarks/puresnet_prepared}"
LIMIT_ARG=()
if [[ "${LIMIT:-}" != "" ]]; then
  LIMIT_ARG=(--limit "$LIMIT")
fi

"$PYTHON_BIN" prepare_puresnet_benchmarks.py \
  --output-root "$OUTPUT_ROOT" \
  --overwrite \
  "${LIMIT_ARG[@]}"

echo "Prepared benchmark root: $OUTPUT_ROOT"
