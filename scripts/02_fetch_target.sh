#!/usr/bin/env bash
# Stage 02. Structures by RCSB id, ChEMBL Env activity export, UniProt Env sequences.
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
  cat <<'USAGE'
Usage: scripts/02_fetch_target.sh [--help]

Downloads 5 Env structures (about 40 MB compressed in total) into data/structures,
the ChEMBL Env-annotated activity export and UniProt Env sequences into data/.
Writes logs/02_fetch.json with dates, versions, counts and file hashes.
No GPU. Needs network access. FORCE=1 repeats a completed stage.
USAGE
  exit 0
fi
skip_if_done 02_fetch_target && exit 0
run_timed 02_fetch_target python "${REPO_ROOT}/scripts/02_fetch_target.py"
stamp_done 02_fetch_target
