Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"
Set-Location (Join-Path $PSScriptRoot "..")
. (Join-Path $PSScriptRoot "_resolve_python.ps1")
$Python = Resolve-ProjectPython
& $Python scripts/02_pre_evolve_ancestor.py @args
