#!/usr/bin/env bash
# Stage 01. Record the tool versions this run used, and create the conda
# environment from config/env_analysis.yml if asked.
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
  cat <<'USAGE'
Usage: scripts/01_install.sh [--create-env] [--help]

Without flags, checks that the required tools are on PATH and writes their
versions to logs/01_install.versions.tsv. The run recorded in this repository
used the base Python of the machine and no dedicated environment, so that is
what the log describes. --create-env builds a conda environment from
config/env_analysis.yml, which lists the same packages with pinned minimums but
has NOT been solved and locked on the machine that made the results.

Tools needed: python (3.10 or later), rdkit, meeko, gemmi, biopython, scipy,
pandas, matplotlib, psutil, requests, Open Babel (obabel), AutoDock Vina.
Optional: quarto (16_report.qmd), shellcheck. No GPU. No licence keys.
USAGE
  exit 0
fi
if [[ "${1:-}" == "--create-env" ]]; then
  conda env create -f "${REPO_ROOT}/config/env_analysis.yml" -n hiv_stage7 || die "env creation failed"
fi
skip_if_done 01_install && exit 0
OUT="${REPO_ROOT}/logs/01_install.versions.tsv"
printf 'tool\tversion\n' > "${OUT}"
VINA_BIN="${VINA:-$(command -v vina || echo /c/vina/vina.exe.exe)}"
{
  printf 'python\t%s\n' "$(python --version 2>&1)"
  printf 'obabel\t%s\n' "$(obabel -V 2>&1 | head -1)"
  printf 'vina\t%s\n' "$("${VINA_BIN}" --version 2>&1 | head -1)"
  python - <<'PY'
import importlib
for m in ["rdkit", "meeko", "gemmi", "Bio", "scipy", "pandas", "numpy", "matplotlib", "psutil", "requests"]:
    try:
        mod = importlib.import_module(m)
        print(f"{m}\t{getattr(mod, '__version__', 'unknown')}")
    except Exception:
        print(f"{m}\tMISSING")
PY
} >> "${OUT}"
grep -q MISSING "${OUT}" && die "missing packages, see ${OUT}"
stamp_done 01_install
log "wrote ${OUT}"
