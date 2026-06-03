#!/usr/bin/env bash
set -Eeuo pipefail

cd "$(dirname "$0")/.."

PYTHON_BIN="${PYTHON:-python3}"
THREADS="${THREADS:-AUTO}"

if ! command -v progressiveMauve >/dev/null 2>&1; then
  echo "ERROR: progressiveMauve is not on PATH. Run scripts/06_check_phylo_tools.sh in WSL." >&2
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

expected_tips=(S01 S02 S03 S04 S05 S06 S07 S08 S09 S10 S11 S12 S13 S14 S15 S16)
status_file="results/phylo/mauve_iqtree/mauve_iqtree_status.tsv"
mkdir -p "$(dirname "$status_file")"
printf 'condition\treplicate\tstatus\tmessage\n' > "$status_file"

for condition in "${conditions[@]}"; do
  replicate="rep_001"
  input_fasta="data/clean_fastas/${condition}/${replicate}/all_tips.clean.fa"
  out_dir="results/phylo/mauve_iqtree/${condition}/${replicate}"
  split_dir="${out_dir}/split_fastas"
  mauve_dir="${out_dir}/mauve"
  iqtree_dir="${out_dir}/iqtree"
  xmfa="${mauve_dir}/alignment.xmfa"
  guide_tree="${mauve_dir}/guide_tree.nwk"
  backbone="${mauve_dir}/alignment.backbone"
  core_alignment="${out_dir}/core_alignment.fa"
  block_summary="${out_dir}/xmfa_blocks.tsv"
  alignment_stats="${out_dir}/core_alignment_stats.tsv"
  tree_file="${out_dir}/inferred_tree.nwk"
  log_file="${out_dir}/run.log"
  iqtree_prefix="${iqtree_dir}/core_alignment"

  if [[ ! -s "$input_fasta" ]]; then
    echo "ERROR: missing clean FASTA for ${condition}: ${input_fasta}" >&2
    echo "Run scripts/07_prepare_phylo_inputs.py first." >&2
    exit 1
  fi

  mkdir -p "$split_dir" "$mauve_dir" "$iqtree_dir"
  if [[ -s "$tree_file" && "${FORCE:-0}" != "1" ]]; then
    printf '%s\t%s\tskipped_existing\t%s\n' "$condition" "$replicate" "$tree_file" >> "$status_file"
    echo "${condition} ${replicate}: skipped existing ${tree_file}"
    continue
  fi

  set +e
  (
    echo "condition=${condition}"
    echo "replicate=${replicate}"
    echo "input_fasta=${input_fasta}"
    echo "iqtree_bin=${IQTREE_BIN}"
    echo "threads=${THREADS}"
    echo

    set -x
    "$PYTHON_BIN" scripts/10_split_fasta_records.py --input "$input_fasta" --out-dir "$split_dir" &&
    mapfile -t split_fastas < "${split_dir}/split_fastas.txt" &&
    progressiveMauve \
      --output="$xmfa" \
      --output-guide-tree="$guide_tree" \
      --backbone-output="$backbone" \
      "${split_fastas[@]}" &&
    "$PYTHON_BIN" scripts/10_xmfa_to_core_alignment.py \
      --xmfa "$xmfa" \
      --out-fasta "$core_alignment" \
      --out-blocks "$block_summary" \
      --expected-tips "${expected_tips[@]}" &&
    "$PYTHON_BIN" scripts/09_alignment_stats.py --alignment "$core_alignment" --out "$alignment_stats" &&
    "$IQTREE_BIN" -s "$core_alignment" -m MFP -nt "$THREADS" -pre "$iqtree_prefix" &&
    cp "${iqtree_prefix}.treefile" "$tree_file"
    command_status=$?
    set +x
    exit "$command_status"
  ) > "$log_file" 2>&1
  exit_code=$?
  set -e

  if [[ "$exit_code" -eq 0 ]]; then
    printf '%s\t%s\tcompleted\t%s\n' "$condition" "$replicate" "$tree_file" >> "$status_file"
    echo "${condition} ${replicate}: wrote ${tree_file}"
  else
    printf '%s\t%s\tfailed\tsee %s\n' "$condition" "$replicate" "$log_file" >> "$status_file"
    echo "${condition} ${replicate}: failed; see ${log_file}" >&2
  fi
done
