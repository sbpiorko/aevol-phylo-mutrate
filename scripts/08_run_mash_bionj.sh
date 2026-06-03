#!/usr/bin/env bash
set -Eeuo pipefail

cd "$(dirname "$0")/.."

PYTHON_BIN="${PYTHON:-python3}"
MASH_K="${MASH_K:-21}"
MASH_SIZE="${MASH_SIZE:-10000}"

if ! command -v mash >/dev/null 2>&1; then
  echo "ERROR: mash is not on PATH. Run scripts/06_check_phylo_tools.sh in WSL." >&2
  exit 1
fi
if ! command -v Rscript >/dev/null 2>&1; then
  echo "ERROR: Rscript is not on PATH. Run scripts/06_check_phylo_tools.sh in WSL." >&2
  exit 1
fi

conditions=("$@")
if ((${#conditions[@]} == 0)); then
  conditions=(very_low low normal high very_high)
fi

for condition in "${conditions[@]}"; do
  replicate="rep_001"
  input_fasta="data/clean_fastas/${condition}/${replicate}/all_tips.clean.fa"
  out_dir="results/phylo/mash_bionj/${condition}/${replicate}"
  sketch_dir="${out_dir}/sketches"
  log_file="${out_dir}/run.log"
  sketch_prefix="${sketch_dir}/all_tips"
  mash_dist="${out_dir}/mash_distances.tsv"
  distance_matrix="${out_dir}/distance_matrix.tsv"
  distance_long="${out_dir}/distance_long.tsv"
  tree_file="${out_dir}/inferred_tree.nwk"

  if [[ ! -s "$input_fasta" ]]; then
    echo "ERROR: missing clean FASTA for ${condition}: ${input_fasta}" >&2
    echo "Run scripts/07_prepare_phylo_inputs.py first." >&2
    exit 1
  fi

  mkdir -p "$sketch_dir"
  {
    echo "condition=${condition}"
    echo "replicate=${replicate}"
    echo "input_fasta=${input_fasta}"
    echo "mash_k=${MASH_K}"
    echo "mash_size=${MASH_SIZE}"
    echo

    set -x
    mash sketch -i -k "$MASH_K" -s "$MASH_SIZE" -o "$sketch_prefix" "$input_fasta"
    mash dist "${sketch_prefix}.msh" "${sketch_prefix}.msh" > "$mash_dist"
    "$PYTHON_BIN" scripts/08_mash_dist_to_matrix.py \
      --fasta "$input_fasta" \
      --mash-dist "$mash_dist" \
      --out-matrix "$distance_matrix" \
      --out-long "$distance_long"
    Rscript scripts/08_bionj_from_matrix.R "$distance_matrix" "$tree_file"
    set +x
  } > "$log_file" 2>&1

  echo "${condition} ${replicate}: wrote ${tree_file}"
done
