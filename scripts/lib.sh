#!/usr/bin/env bash
# Shared helpers, sourced by every stage. Conventions follow the stage 1 and
# stage 4 repositories: stamps for idempotence, run_timed for elapsed time and
# peak resident memory, and a scratch registry so a killed run leaves no
# half-modelled glycan tree behind.
#
# Stamps do not hash their inputs. Hashing would invalidate every downstream
# stage whenever a comment changed upstream, and glycan modelling is expensive
# enough that a stale stamp is the smaller failure.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export REPO_ROOT
CONF_FILE="${REPO_ROOT}/project.conf"
STAMP_DIR="${REPO_ROOT}/logs/stamps"
mkdir -p "${STAMP_DIR}"

log()  { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*" >&2; }
warn() { printf '[%s] WARNING: %s\n' "$(date +%H:%M:%S)" "$*" >&2; }
die()  { printf '[%s] ERROR: %s\n' "$(date +%H:%M:%S)" "$*" >&2; exit 1; }

load_conf() {
  [[ -f "${CONF_FILE}" ]] || die "project.conf not found. Run scripts/00_configure.sh first."
  # shellcheck source=/dev/null
  source "${CONF_FILE}"
  : "${THREADS:?THREADS missing from project.conf}"
  : "${BENCH_TIER:?BENCH_TIER missing from project.conf}"
}

stamp_done() {
  { date -Iseconds; git -C "${REPO_ROOT}" rev-parse --short HEAD 2>/dev/null || echo "no-git"; } \
    > "${STAMP_DIR}/$1.done"
}
stamp_exists() { [[ -f "${STAMP_DIR}/$1.done" ]]; }

skip_if_done() {
  if stamp_exists "$1" && [[ "${FORCE:-0}" != "1" ]]; then
    log "skipping $1: already complete (stamp logs/stamps/$1.done). FORCE=1 repeats it."
    return 0
  fi
  return 1
}

# GNU time reports peak RSS; the bash builtin does not. On a machine without
# /usr/bin/time the stage still runs and the timing file says NA, rather than
# the stage failing over a measurement.
run_timed() {
  local stage="$1"; shift
  local tf="${REPO_ROOT}/logs/${stage}.timing"
  mkdir -p "$(dirname "${tf}")"
  if [[ -x /usr/bin/time ]]; then
    /usr/bin/time -f "%e %M" -o "${tf}" -- "$@"
  else
    local t0=${SECONDS}
    "$@"
    echo "$((SECONDS - t0)) NA" > "${tf}"
    warn "${stage}: /usr/bin/time not found, peak RSS not recorded"
  fi
}

dir_mb() { du -sm "$1" 2>/dev/null | awk '{print $1}' || echo 0; }
require_file() { [[ -s "$1" ]] || die "required file missing or empty: $1"; }

# Directories registered here are removed on exit or kill. The trap returns 0
# explicitly: an EXIT trap's last status becomes the script's, and a failed
# final test would turn a finished stage into a failed one.
SCRATCH_DIRS=()
register_scratch() { SCRATCH_DIRS+=("$1"); }
cleanup_scratch() {
  local d
  for d in "${SCRATCH_DIRS[@]:-}"; do
    if [[ -n "${d}" && -d "${d}" ]]; then rm -rf "${d}"; fi
  done
  return 0
}
trap cleanup_scratch EXIT INT TERM

MODALITIES=(small_molecule peptide binder aptamer)
GLYCAN_ARMS=(glycan_free glycosylated)
export MODALITIES GLYCAN_ARMS
