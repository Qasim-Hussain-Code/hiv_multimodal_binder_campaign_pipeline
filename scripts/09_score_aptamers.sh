#!/usr/bin/env bash
# Stage 09. Aptamer arm. Not run.
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
  cat <<'USAGE'
Usage: scripts/09_score_aptamers.sh [--help]

Writes results/modality_aptamer_not_run.md, the documented decision not to run
the aptamer arm. No GPU.
USAGE
  exit 0
fi
skip_if_done 09_score_aptamers && exit 0
mkdir -p "${REPO_ROOT}/results"
cat > "${REPO_ROOT}/results/modality_aptamer_not_run.md" <<'MD'
# Aptamer arm: not run

Stage 6 was never run, so no aptamer candidate exists. The control for this
modality, UCLA1 (derived from B40), was not scored.

The repository stanford-rna-fold predicts RNA 3D structure from sequence for a
folding competition. It does not dock or rank aptamers against a protein and it
was not used.

An aptamer arm would be worth running if there were a sequence-to-affinity
model validated on held-out aptamer-protein pairs and a nucleic-acid null
drawn from the same length and composition. Neither exists in this project. An
empty modality reported as empty is better than a padded one.
MD
stamp_done 09_score_aptamers
log "wrote results/modality_aptamer_not_run.md"
