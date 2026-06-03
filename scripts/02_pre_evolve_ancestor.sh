#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
python3 scripts/02_pre_evolve_ancestor.py "$@"

