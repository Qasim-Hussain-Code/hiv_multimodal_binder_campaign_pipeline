#!/usr/bin/env bash
# Stage 00. Record the machine and the budget, project the footprint, and
# refuse to start if the projection does not fit.
#
# Glycan modelling and re-prediction are the expensive steps. The projection is
# only as good as the per-candidate figure it multiplies, so that figure comes
# from --per-candidate-mb once the first ten candidates have been measured.
# Until then a placeholder is used and project.conf says the projection is
# unmeasured. An unmeasured projection is a warning, not a guarantee. The real
# refusal happens again after the ten-candidate measurement.

source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"

usage() {
  cat <<'USAGE'
Usage: scripts/00_configure.sh --threads N --ram GB --disk GB --tier A|B|C [options]

Writes project.conf, which every later script sources.

Threads and RAM are detected where possible and must be passed otherwise.

  --threads N            CPU threads available
  --ram GB               RAM available to the pipeline
  --disk GB              working disk this analysis may use (required)
  --tier A|B|C           bench tier (default from config/bench_tier.conf, A)
  --jobs N               concurrent jobs (default 1)
  --per-candidate-mb X   measured working disk per candidate, per state, per
                         glycan arm, from the first ten candidates
  --n-candidates N       candidates per modality (default 800)
  --n-states N           conformational states in config/panel.tsv (default 4)
  --seed N               master random seed (default 20260929)
  --yes                  overwrite an existing project.conf
  --help                 this message

The default for --jobs is 1. It is serial for reproducible timing and not for
memory: parallel wall-clock distributions are not comparable to serial ones.

GPU: this script does not use one. Stages 08 and any structure prediction step
say in their own --help output whether they need one.
USAGE
}

THREADS_ARG=""; RAM_ARG=""; DISK_ARG=""; TIER_ARG=""; JOBS_ARG=1
PER_CAND_MB=""; N_CAND=800; N_STATES=4; SEED=20260929; YES=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --threads) THREADS_ARG="$2"; shift 2;;
    --ram) RAM_ARG="$2"; shift 2;;
    --disk) DISK_ARG="$2"; shift 2;;
    --tier) TIER_ARG="$2"; shift 2;;
    --jobs) JOBS_ARG="$2"; shift 2;;
    --per-candidate-mb) PER_CAND_MB="$2"; shift 2;;
    --n-candidates) N_CAND="$2"; shift 2;;
    --n-states) N_STATES="$2"; shift 2;;
    --seed) SEED="$2"; shift 2;;
    --yes) YES=1; shift;;
    --help|-h) usage; exit 0;;
    *) die "unknown option: $1 (see --help)";;
  esac
done

if [[ "${YES}" != "1" ]] && skip_if_done 00_configure; then exit 0; fi

# Detection is a convenience. Git Bash on Windows has no /proc/meminfo, so an
# explicit flag always wins and a value that cannot be detected is an error.
[[ -n "${THREADS_ARG}" ]] || THREADS_ARG="$(nproc 2>/dev/null || true)"
[[ -n "${RAM_ARG}" ]] || RAM_ARG="$(awk '/MemTotal/ {printf "%d", $2/1048576}' /proc/meminfo 2>/dev/null || true)"
[[ -n "${THREADS_ARG}" ]] || die "cannot detect threads, pass --threads"
[[ -n "${RAM_ARG}" ]] || die "cannot detect RAM, pass --ram"
[[ -n "${DISK_ARG}" ]] || die "pass --disk (GB of working disk this analysis may use)"

# shellcheck source=/dev/null
source "${REPO_ROOT}/config/bench_tier.conf"
TIER="${TIER_ARG:-${BENCH_TIER}}"
case "${TIER}" in A|B|C) ;; *) die "--tier must be A, B or C";; esac
if [[ "${TIER}" != "A" && -z "${BIOSAFETY_APPROVAL_REF}" ]]; then
  die "tier ${TIER} needs BIOSAFETY_APPROVAL_REF and BIOSAFETY_APPROVAL_DATE in config/bench_tier.conf. Ask the institutional biosafety committee. This script does not supply a classification."
fi

if [[ -f "${CONF_FILE}" && "${YES}" != "1" ]]; then
  die "project.conf exists. Pass --yes to overwrite."
fi

# Projection: candidates x 4 modalities x states x 2 glycan arms x MB each.
PROJ_NOTE="measured"
if [[ -z "${PER_CAND_MB}" ]]; then
  PER_CAND_MB=0.1
  PROJ_NOTE="UNMEASURED placeholder of 0.1 MB (a guess, tables and top-N poses only), replace after the first ten candidates"
fi
PROJ_GB="$(awk -v n="${N_CAND}" -v s="${N_STATES}" -v m="${PER_CAND_MB}" \
  'BEGIN{printf "%.2f", n*4*s*2*m/1024}')"
log "projected working disk: ${PROJ_GB} GB (${PROJ_NOTE}); budget ${DISK_ARG} GB"

if awk -v p="${PROJ_GB}" -v d="${DISK_ARG}" 'BEGIN{exit !(p>d)}'; then
  SHORT="$(awk -v p="${PROJ_GB}" -v d="${DISK_ARG}" 'BEGIN{printf "%.2f", p-d}')"
  die "refusing to start: projection exceeds the disk budget by ${SHORT} GB. Reduce --n-candidates or --n-states, or raise --disk."
fi

cat > "${CONF_FILE}" <<CONF
# Written by scripts/00_configure.sh on $(date -Iseconds). Do not edit by hand.
THREADS=${THREADS_ARG}
RAM_GB=${RAM_ARG}
DISK_GB=${DISK_ARG}
JOBS=${JOBS_ARG}
BENCH_TIER=${TIER}
MASTER_SEED=${SEED}
N_CANDIDATES=${N_CAND}
N_STATES=${N_STATES}
PER_CANDIDATE_MB=${PER_CAND_MB}
PROJECTED_DISK_GB=${PROJ_GB}
PROJECTION_NOTE="${PROJ_NOTE}"
CONF
stamp_done 00_configure
log "wrote project.conf (tier ${TIER}, ${THREADS_ARG} threads, ${RAM_ARG} GB RAM)"
