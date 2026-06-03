#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

PYTHON_BIN="${PYTHON:-python3}"
"$PYTHON_BIN" scripts/validate_setup.py
"$PYTHON_BIN" scripts/00_make_tree.py
"$PYTHON_BIN" scripts/01_make_mutation_params.py

