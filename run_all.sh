#!/usr/bin/env bash
# Orchestrator. Runs the stages in order, skipping any that already carry a stamp.
source "$(dirname "${BASH_SOURCE[0]}")/scripts/lib.sh"

usage() {
  cat <<'USAGE'
Usage: bash run_all.sh [--from STAGE] [--modality NAME] [--state ID] [--tier A|B|C] [--jobs N] [--help]

  --from STAGE      first stage number to run, for example 06 (default 00)
  --modality NAME   small_molecule, peptide, binder or aptamer (default: all four;
                    peptide, binder and aptamer only write their not-run notes)
  --state ID        one panel state for the docking stage, for example S2_closed_holo
  --tier A|B|C      bench tier. Tier B and C add stage 17 and need an approval
                    reference in config/bench_tier.conf. Default A.
  --jobs N          concurrent Vina processes, default 1. Serial is the default
                    for reproducible timing, not for memory. Each Vina process
                    peaked near 0.5 GB resident in the run recorded here.

Stage 00 needs its own arguments the first time (threads, RAM, disk):
  bash scripts/00_configure.sh --threads N --ram GB --disk GB --tier A --yes

Docking (stage 06) is the long stage: 168 compounds x 4 states x 2 glycan arms
= 1344 Vina runs, about 30 s each on the machine that produced the results.
No stage needs a GPU. Stages that stamp themselves are skipped on a re-run, and
the skip is printed. FORCE=1 repeats everything. The pre-registration stage
refuses to overwrite an existing preregistration.md.
USAGE
}
FROM=0; MOD=all; STATE=""; TIER=""; JOBS_ARG=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --from) FROM=$((10#$2)); shift 2;;
    --modality) MOD="$2"; shift 2;;
    --state) STATE="$2"; shift 2;;
    --tier) TIER="$2"; shift 2;;
    --jobs) JOBS_ARG="$2"; shift 2;;
    --help|-h) usage; exit 0;;
    *) die "unknown option: $1";;
  esac
done
case "${MOD}" in all|small_molecule|peptide|binder|aptamer) ;; *) die "unknown modality ${MOD}";; esac

S="${REPO_ROOT}/scripts"
runpy() {  # stage-id, script, args...
  local id="$1" script="$2"; shift 2
  if skip_if_done "${id}"; then return 0; fi
  log "running ${id}"
  run_timed "${id}" python "${S}/${script}" "$@"
  stamp_done "${id}"
}
want() { [[ $((10#$1)) -ge ${FROM} ]]; }
mod_is() { [[ "${MOD}" == "all" || "${MOD}" == "$1" ]]; }

want 0  && { [[ -f "${CONF_FILE}" ]] || die "project.conf missing. Run scripts/00_configure.sh first (see --help)."; }
load_conf
[[ -n "${TIER}" && "${TIER}" != "${BENCH_TIER}" ]] && die "--tier ${TIER} differs from project.conf (${BENCH_TIER}). Rerun 00_configure.sh --tier ${TIER} --yes."
want 1  && bash "${S}/01_install.sh"
want 2  && bash "${S}/02_fetch_target.sh"
want 3  && runpy 03_build_conformer_panel 03_build_conformer_panel.py
want 4  && runpy 04_prepare_receptor 04_prepare_receptor.py
want 5  && runpy 05_assemble_candidates 05_assemble_candidates.py
if want 6; then
  mod_is small_molecule && { A=(); [[ -n "${STATE}" ]] && A+=(--states "${STATE}"); [[ -n "${JOBS_ARG}" ]] && A+=(--jobs "${JOBS_ARG}"); bash "${S}/06_score_small_molecules.sh" "${A[@]}"; }
  mod_is peptide && bash "${S}/07_score_peptides.sh"
  mod_is binder  && bash "${S}/08_score_binders.sh"
  mod_is aptamer && bash "${S}/09_score_aptamers.sh"
fi
want 10 && runpy 10_conservation 10_conservation.py
want 11 && runpy 11_rank_within_modality 11_rank_within_modality.py
want 12 && runpy 12_holdout_validation 12_holdout_validation.py
if want 13; then
  if [[ -f "${REPO_ROOT}/results/preregistration.md" ]]; then
    log "skipping 13_prereg: results/preregistration.md exists and is never overwritten"
  else
    runpy 13_prereg 13_prereg.py
  fi
fi
want 14 && runpy 14_analyse 14_analyse.py
want 15 && runpy 15_figures 15_figures.py
if want 16; then
  if command -v quarto >/dev/null 2>&1; then
    (cd "${S}" && quarto render 16_report.qmd --output-dir ../results) || warn "quarto render failed"
  else
    warn "quarto not found, skipping 16_report.qmd (results tables are complete without it)"
  fi
fi
if [[ "${BENCH_TIER}" != "A" ]]; then
  log "tier ${BENCH_TIER}: run scripts/17_assay.py --results <assay.tsv> once assay data exist"
else
  log "tier A: the repository ends at results/preregistration.md and the costed plan in results/cost_model.tsv"
fi
