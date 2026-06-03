Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

Set-Location (Join-Path $PSScriptRoot "..")
. (Join-Path $PSScriptRoot "_resolve_python.ps1")
$Python = Resolve-ProjectPython

& $Python scripts/validate_setup.py
& $Python scripts/00_make_tree.py
& $Python scripts/01_make_mutation_params.py
& $Python scripts/02_pre_evolve_ancestor.py
& $Python scripts/03_run_species_tree.py
& $Python scripts/04_collect_final_fastas.py
& $Python scripts/05_qc_fastas.py
