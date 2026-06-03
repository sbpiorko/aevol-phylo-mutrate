# Mutation-Rate Aevol 4b Simulation Scaffold

This project scaffold prepares the simulation-only part of the benchmark:

1. pre-evolve one baseline ancestor,
2. run the same 16-tip balanced species tree under five mutation-rate regimes,
3. collect final tip genomes as FASTA files,
4. run basic FASTA quality control.

Phylogeny inference is intentionally out of scope for this version.

## Mutation Regimes

Only these local sequence-level rates are scaled:

```text
POINT_MUTATION_RATE
SMALL_INSERTION_RATE
SMALL_DELETION_RATE
```

The default grid is:

```text
very_low   0.0625x
low        0.25x
normal     1.0x
high       4.0x
very_high  16.0x
```

Large rearrangement rates are left at their baseline values.

## Before You Run

Install Aevol 4b and make sure the command-line executables are on your `PATH`.

The scripts try these command names by default:

```text
create:  aevol_4b_create, aevol_create_4b
run:     aevol_4b_run, aevol_run_4b
extract: aevol_4b_misc_extract, aevol_misc_extract,
         aevol_post_extract_genomes_4b, aevol_4b_post_extract_genomes
```

If your local executable names differ, edit `config/run_config.yaml`.

Also inspect `config/base_params.in`. It is a minimal starter file based on Aevol documentation defaults, not a guarantee of the exact best setup for your local Aevol build. If you already have a known-good Aevol 4b parameter file, replace `config/base_params.in` with that file before running simulations.

## Run Order

From the project root:

```bash
cd aevol_phylo_mutrate
python scripts/validate_setup.py
python scripts/00_make_tree.py
python scripts/01_make_mutation_params.py
python scripts/02_pre_evolve_ancestor.py
python scripts/03_run_species_tree.py
python scripts/04_collect_final_fastas.py
python scripts/05_qc_fastas.py
```

On Windows PowerShell, you can run the setup-only stages with:

```powershell
cd aevol_phylo_mutrate
.\scripts\run_setup.ps1
```

To run the full FASTA-production pipeline after Aevol is installed:

```powershell
.\scripts\run_simulations.ps1
```

If your Python is not globally installed, the PowerShell wrappers will also try the Codex bundled Python runtime.

## Expected Final FASTA Outputs

After the full simulation stage, you should have:

```text
data/final_fastas/very_low/rep_001/all_tips.fa
data/final_fastas/low/rep_001/all_tips.fa
data/final_fastas/normal/rep_001/all_tips.fa
data/final_fastas/high/rep_001/all_tips.fa
data/final_fastas/very_high/rep_001/all_tips.fa
```

Each `all_tips.fa` should contain exactly 16 sequences named:

```text
S01, S02, ..., S16
```

## Useful Debug Commands

Generate only tree/config products:

```bash
scripts/run_setup.sh
```

Run only one condition:

```bash
python scripts/03_run_species_tree.py --conditions normal
```

Run only one condition and one replicate:

```bash
python scripts/03_run_species_tree.py --conditions normal --replicates 1
```

Force rerunning branch directories that already have outputs:

```bash
python scripts/03_run_species_tree.py --conditions normal --force
```

## Notes On Aevol Initialization

Each branch is run in its own directory. Child branches are initialized from the parent FASTA using:

```text
aevol_4b_create params.in --fasta parent.fa
```

The script also supports `aevol_create_4b` if that is the executable name in your installation.

## What This Version Does Not Do

This scaffold does not infer trees, run progressiveMauve, run IQ-TREE, compute k-mer distances, or score RF/KF distances. Those belong in the next phase after the FASTA production pipeline is working.
