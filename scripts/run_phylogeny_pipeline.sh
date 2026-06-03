#!/usr/bin/env bash
set -Eeuo pipefail
cd "$(dirname "$0")/.."

PYTHON_BIN="${PYTHON:-python3}"

"$PYTHON_BIN" scripts/07_prepare_phylo_inputs.py
bash scripts/08_run_mash_bionj.sh
bash scripts/09_run_mafft_iqtree.sh
bash scripts/10_run_mauve_iqtree.sh
"$PYTHON_BIN" scripts/12_score_trees.py
"$PYTHON_BIN" scripts/13_plot_results.py

cat >&2 <<'MSG'
Implemented phylogeny stages finished:
  - clean FASTA preparation
  - Mash distance tree with R ape::bionj()
  - MAFFT -> IQ-TREE2
  - progressiveMauve -> IQ-TREE2
  - topology-first tree scoring
  - summary SVG plots

The remaining methods are still placeholders:
  - optional MUMmer distance tree
MSG
