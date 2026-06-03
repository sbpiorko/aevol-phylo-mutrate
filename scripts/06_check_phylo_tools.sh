#!/usr/bin/env bash
set -Eeuo pipefail

cd "$(dirname "$0")/.."

out_dir="results/software_versions"
out_file="$out_dir/phylo_tools.tsv"
mkdir -p "$out_dir"

tmp_file="$(mktemp)"
missing_required=0

clean_one_line() {
  tr '\t' ' ' | tr '\r' ' ' | sed -E 's/[[:space:]]+/ /g; s/^ //; s/ $//' | head -n 1
}

capture_version() {
  local cmd="$1"
  shift
  if command -v "$cmd" >/dev/null 2>&1; then
    timeout 20s "$cmd" "$@" 2>&1 | clean_one_line || true
  fi
}

record_cmd() {
  local group="$1"
  local label="$2"
  local required="$3"
  local cmd="$4"
  shift 4
  local path=""
  local version=""
  local status="missing"

  if command -v "$cmd" >/dev/null 2>&1; then
    path="$(command -v "$cmd")"
    version="$(capture_version "$cmd" "$@")"
    status="ok"
  elif [[ "$required" == "yes" ]]; then
    missing_required=$((missing_required + 1))
  fi

  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$group" "$label" "$status" "$path" "$version" "$required" >> "$tmp_file"
}

record_iqtree() {
  if command -v iqtree2 >/dev/null 2>&1; then
    record_cmd "bioinformatics" "iqtree" "yes" "iqtree2" "-version"
  elif command -v iqtree >/dev/null 2>&1; then
    record_cmd "bioinformatics" "iqtree" "yes" "iqtree" "-version"
  else
    missing_required=$((missing_required + 1))
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "bioinformatics" "iqtree" "missing" "" "" "yes" >> "$tmp_file"
  fi
}

record_python_package() {
  local package="$1"
  local import_name="$2"
  local required="$3"
  local python_bin=".venv_phylo/bin/python"
  local status="missing"
  local version=""
  local path=""

  if [[ -x "$python_bin" ]]; then
    path="$python_bin"
    if version="$("$python_bin" - "$import_name" <<'PY' 2>/dev/null
import importlib
import sys
name = sys.argv[1]
mod = importlib.import_module(name)
print(getattr(mod, "__version__", "installed"))
PY
)"; then
      status="ok"
      version="$(printf '%s' "$version" | clean_one_line)"
    fi
  fi

  if [[ "$status" == "missing" && "$required" == "yes" ]]; then
    missing_required=$((missing_required + 1))
  fi
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "python_package" "$package" "$status" "$path" "$version" "$required" >> "$tmp_file"
}

record_r_package() {
  local package="$1"
  local required="$2"
  local status="missing"
  local version=""
  local path=""

  if command -v Rscript >/dev/null 2>&1; then
    path="$(command -v Rscript)"
    if version="$(Rscript -e "if (!requireNamespace('$package', quietly=TRUE)) quit(status=2); cat(as.character(utils::packageVersion('$package')))" 2>/dev/null)"; then
      status="ok"
      version="$(printf '%s' "$version" | clean_one_line)"
    fi
  fi

  if [[ "$status" == "missing" && "$required" == "yes" ]]; then
    missing_required=$((missing_required + 1))
  fi
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "r_package" "$package" "$status" "$path" "$version" "$required" >> "$tmp_file"
}

printf '%s\t%s\t%s\t%s\t%s\t%s\n' "group" "tool" "status" "path" "version" "required" > "$out_file"

record_cmd "runtime" "python3" "yes" "python3" "--version"
record_cmd "runtime" "Rscript" "yes" "Rscript" "--version"
record_cmd "bioinformatics" "seqkit" "yes" "seqkit" "version"
record_cmd "bioinformatics" "mash" "yes" "mash" "--version"
record_cmd "bioinformatics" "fastme" "yes" "fastme" "--version"
record_cmd "bioinformatics" "mafft" "yes" "mafft" "--version"
record_iqtree
record_cmd "whole_genome_alignment" "progressiveMauve" "no" "progressiveMauve" "--version"
record_cmd "optional_distance" "nucmer" "no" "nucmer" "--version"
record_cmd "optional_distance" "dnadiff" "no" "dnadiff" "--version"

record_python_package "biopython" "Bio" "yes"
record_python_package "pandas" "pandas" "yes"
record_python_package "ete3" "ete3" "no"

record_r_package "ape" "yes"
record_r_package "phangorn" "yes"
record_r_package "tidyverse" "yes"
record_r_package "dendextend" "no"

cat "$tmp_file" >> "$out_file"
rm -f "$tmp_file"

cat "$out_file"

if ((missing_required)); then
  printf '\nERROR: %d required phylogeny setup checks are missing. See %s\n' "$missing_required" "$out_file" >&2
  exit 1
fi

printf '\nWrote %s\n' "$out_file"
