#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

PYTHON_BIN="${PYTHON:-python3}"
"$PYTHON_BIN" scripts/validate_setup.py
"$PYTHON_BIN" scripts/00_make_tree.py
"$PYTHON_BIN" scripts/01_make_mutation_params.py
"$PYTHON_BIN" scripts/02_pre_evolve_ancestor.py
"$PYTHON_BIN" scripts/03_run_species_tree.py
"$PYTHON_BIN" scripts/04_collect_final_fastas.py
"$PYTHON_BIN" scripts/05_qc_fastas.py
