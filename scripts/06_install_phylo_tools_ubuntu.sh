#!/usr/bin/env bash
set -Eeuo pipefail

cd "$(dirname "$0")/.."

log() {
  printf '\n==> %s\n' "$*"
}

warn() {
  printf 'WARNING: %s\n' "$*" >&2
}

have_cmd() {
  command -v "$1" >/dev/null 2>&1
}

apt_has_package() {
  apt-cache show "$1" >/dev/null 2>&1
}

install_apt_packages() {
  local package
  for package in "$@"; do
    if apt_has_package "$package"; then
      log "Installing apt package: $package"
      "${SUDO[@]}" apt-get install -y "$package"
    else
      warn "Apt package not available in configured repositories: $package"
      MISSING_APT+=("$package")
    fi
  done
}

if [[ ! -r /etc/os-release ]]; then
  echo "This installer is intended for Ubuntu or Ubuntu-like WSL environments." >&2
  exit 1
fi

# shellcheck disable=SC1091
. /etc/os-release
if [[ "${ID:-}" != "ubuntu" && "${ID_LIKE:-}" != *"ubuntu"* && "${ID_LIKE:-}" != *"debian"* ]]; then
  warn "Detected ${PRETTY_NAME:-unknown OS}; this script is tested for Ubuntu/WSL."
fi

if [[ "${EUID:-$(id -u)}" -eq 0 ]]; then
  SUDO=()
else
  SUDO=(sudo)
fi

MISSING_APT=()

log "Updating apt metadata"
"${SUDO[@]}" apt-get update

log "Enabling Ubuntu universe repository when available"
"${SUDO[@]}" apt-get install -y software-properties-common
if have_cmd add-apt-repository; then
  "${SUDO[@]}" add-apt-repository -y universe || true
  "${SUDO[@]}" apt-get update
fi

log "Installing core system/runtime packages"
install_apt_packages \
  build-essential \
  ca-certificates \
  curl \
  wget \
  unzip \
  tar \
  gzip \
  bzip2 \
  xz-utils \
  pkg-config \
  python3 \
  python3-pip \
  python3-venv \
  r-base \
  r-base-dev \
  libcurl4-openssl-dev \
  libssl-dev \
  libxml2-dev

log "Installing command-line bioinformatics tools"
install_apt_packages \
  seqkit \
  mash \
  fastme \
  mafft \
  iqtree \
  progressivemauve \
  mummer

log "Installing R phylogenetics/data packages from Ubuntu repositories when available"
install_apt_packages \
  r-cran-ape \
  r-cran-phangorn \
  r-cran-tidyverse \
  r-cran-dendextend

log "Creating project Python environment: .venv_phylo"
python3 -m venv .venv_phylo
.venv_phylo/bin/python -m pip install --upgrade pip
.venv_phylo/bin/python -m pip install biopython pandas ete3

log "Installing any missing R packages from CRAN"
Rscript - <<'RSCRIPT'
repos <- "https://cloud.r-project.org"
pkgs <- c("ape", "phangorn", "tidyverse", "dendextend")
installed <- rownames(installed.packages())
missing <- setdiff(pkgs, installed)
if (length(missing)) {
  install.packages(missing, repos = repos, dependencies = TRUE)
}
RSCRIPT

log "Checking installed tools and writing version manifest"
bash scripts/06_check_phylo_tools.sh

if ((${#MISSING_APT[@]})); then
  warn "Some apt packages were unavailable: ${MISSING_APT[*]}"
  warn "If any are marked missing in the version check, we can add a targeted fallback installer."
fi

log "Phylogeny tool setup finished"
