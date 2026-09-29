#!/usr/bin/env bash
# Stage 08. Designed protein binder arm. Not run.
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
  cat <<'USAGE'
Usage: scripts/08_score_binders.sh [--help]

Writes results/modality_binder_not_run.md. The binder arm needs stage 5
(BindCraft and RFdiffusion, a GPU workload) and stage 2 re-prediction of its
designs. Stage 5 was not run, so there are no designs to re-predict.
No GPU was used by this script.
USAGE
  exit 0
fi
skip_if_done 08_score_binders && exit 0
mkdir -p "${REPO_ROOT}/results"
cat > "${REPO_ROOT}/results/modality_binder_not_run.md" <<'MD'
# Designed protein binder arm: not run

No designs were generated, re-predicted or ranked.

Stage 5 (BindCraft and RFdiffusion, GPU) was never run, so there is no design
set. Stage 2 (colabfold_boltz_structure_confidence_pipeline) exists and could
re-predict designs, but it has nothing to re-predict. The positive controls for
this modality, eCD4-Ig and a CD4-mimetic miniprotein, were not scored.

The CATNAP retrospective test the brief describes belongs to this arm. It was
not run either. CATNAP is distributed through a web form at hiv.lanl.gov and
this machine had no pipeline to feed it.

Evaluating an AlphaFold2-family design with AlphaFold2-family re-prediction is
circular. If this arm is run, the independent arm should be a physics-based or
experimental readout and not a second structure predictor.
MD
stamp_done 08_score_binders
log "wrote results/modality_binder_not_run.md"
