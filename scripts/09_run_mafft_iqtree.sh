#!/usr/bin/env bash
set -Eeuo pipefail

cd "$(dirname "$0")/.."

PYTHON_BIN="${PYTHON:-python3}"
THREADS="${THREADS:-AUTO}"

if ! command -v mafft >/dev/null 2>&1; then
  echo "ERROR: mafft is not on PATH. Run scripts/06_check_phylo_tools.sh in WSL." >&2
  exit 1
fi

if command -v iqtree2 >/dev/null 2>&1; then
  IQTREE_BIN="iqtree2"
elif command -v iqtree >/dev/null 2>&1; then
  IQTREE_BIN="iqtree"
else
  echo "ERROR: neither iqtree2 nor iqtree is on PATH. Run scripts/06_check_phylo_tools.sh in WSL." >&2
  exit 1
fi

conditions=("$@")
if ((${#conditions[@]} == 0)); then
  conditions=(very_low low normal high very_high)
fi

for condition in "${conditions[@]}"; do
  replicate="rep_001"
  input_fasta="data/clean_fastas/${condition}/${replicate}/all_tips.clean.fa"
  out_dir="results/phylo/mafft_iqtree/${condition}/${replicate}"
  iqtree_dir="${out_dir}/iqtree"
  alignment="${out_dir}/alignment.fa"
  alignment_stats="${out_dir}/alignment_stats.tsv"
  tree_file="${out_dir}/inferred_tree.nwk"
  log_file="${out_dir}/run.log"
  iqtree_prefix="${iqtree_dir}/alignment"

  if [[ ! -s "$input_fasta" ]]; then
    echo "ERROR: missing clean FASTA for ${condition}: ${input_fasta}" >&2
    echo "Run scripts/07_prepare_phylo_inputs.py first." >&2
    exit 1
  fi

  mkdir -p "$out_dir" "$iqtree_dir"
  {
    echo "condition=${condition}"
    echo "replicate=${replicate}"
    echo "input_fasta=${input_fasta}"
    echo "iqtree_bin=${IQTREE_BIN}"
    echo "threads=${THREADS}"
    echo

    set -x
    mafft --auto "$input_fasta" > "$alignment"
    "$PYTHON_BIN" scripts/09_alignment_stats.py --alignment "$alignment" --out "$alignment_stats"
    "$IQTREE_BIN" -s "$alignment" -m MFP -nt "$THREADS" -pre "$iqtree_prefix"
    cp "${iqtree_prefix}.treefile" "$tree_file"
    set +x
  } > "$log_file" 2>&1

  echo "${condition} ${replicate}: wrote ${tree_file}"
done
