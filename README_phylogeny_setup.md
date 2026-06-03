# Phylogeny Tool Setup

This step installs and checks the software environment for the next analysis
phase. It does not run phylogeny inference or create downstream result folders.

Run from Ubuntu/WSL at the project root:

```bash
cd /mnt/c/Users/Fiend/Documents/Projects/Aevol/aevol_phylo_mutrate
bash scripts/06_install_phylo_tools_ubuntu.sh
```

The installer attempts to install:

- FASTA/QC and distance tools: `seqkit`, `mash`
- Tree/distance tools: `fastme`, R `ape`, R `phangorn`
- Alignment tools: `mafft`, `iqtree`
- Whole-genome alignment tool: `progressiveMauve`
- Optional distance tool: `mummer`
- Python packages in `.venv_phylo`: `biopython`, `pandas`, `ete3`
- R packages: `ape`, `phangorn`, `tidyverse`, `dendextend`

After installation, or any time you want to audit the environment, run:

```bash
bash scripts/06_check_phylo_tools.sh
```

The checker writes:

```text
results/software_versions/phylo_tools.tsv
```

If the installer reports that `progressiveMauve` is unavailable from your Ubuntu
repositories, continue with the other installed tools for now. We can add a
targeted fallback for Mauve before implementing that specific method.
