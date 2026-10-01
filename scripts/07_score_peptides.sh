#!/usr/bin/env bash
# Stage 07. Peptide arm. Not run. Writes the reason instead of a ranked table.
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
  cat <<'USAGE'
Usage: scripts/07_score_peptides.sh [--help]

Writes results/modality_peptide_not_run.md. The peptide arm needs the stage 4
peptide screen, which was never run, so there is nothing to import and no
peptide docking function fitted on peptides to call. No GPU. FORCE=1 repeats.
USAGE
  exit 0
fi
skip_if_done 07_score_peptides && exit 0
mkdir -p "${REPO_ROOT}/results"
cat > "${REPO_ROOT}/results/modality_peptide_not_run.md" <<'MD'
# Peptide arm: not run

No peptide candidates were ranked.

This arm was meant to import stage 4's peptide screen. Stage 4
(vina_litpcba_virtual_screening_pipeline) screens small molecules on 15
LIT-PCBA targets and has no peptide arm. Nothing was carried forward, so the
positive control for this modality, enfuvirtide (T-20), was not scored and no
percentile exists for it.

Docking a peptide with a small-molecule scoring function that was never fitted
on peptides would give a number and no meaning. That is the reason the arm was
left empty and not filled with a table.

What would have to be true for the arm to be worth running: a peptide docking
or design pipeline validated on a held-out set of peptide-protein complexes,
with a null arm of random sequences of matched length, run against a gp41
heptad-repeat model whose conformational state has been checked against the
literature. Enfuvirtide's approval details and structural basis were not
verified for this repository and are left out on purpose.
MD
stamp_done 07_score_peptides
log "wrote results/modality_peptide_not_run.md"
