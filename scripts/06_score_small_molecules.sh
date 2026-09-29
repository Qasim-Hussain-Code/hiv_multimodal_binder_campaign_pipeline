#!/usr/bin/env bash
# Stage 06. Dock the small-molecule set against the panel, both glycan arms.
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"

usage() {
  cat <<'USAGE'
Usage: scripts/06_score_small_molecules.sh [--jobs N] [--states a,b] [--help]

Measures the first ten candidates in one state, projects the total wall time
and working disk for all candidates, states and arms, refuses to start if the
disk projection exceeds DISK_GB in project.conf, then docks everything with
AutoDock Vina 1.2.7 (scripts/06_dock.py).

The default is serial (--jobs 1) for reproducible timing and not for memory.
No GPU. Peak resident memory per Vina process is logged per docking.
FORCE=1 repeats a completed stage. Interrupted runs resume: finished dockings
are read back from results/dock_scores.tsv and skipped.
USAGE
}
JOBS_ARG=""; STATES=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --jobs) JOBS_ARG="$2"; shift 2;;
    --states) STATES="$2"; shift 2;;
    --help|-h) usage; exit 0;;
    *) die "unknown option: $1";;
  esac
done
load_conf
JOBS_USE="${JOBS_ARG:-${JOBS}}"
skip_if_done 06_score_small_molecules && exit 0
require_file "${REPO_ROOT}/data/candidates/candidates.tsv"

if [[ ! -s "${REPO_ROOT}/logs/06_projection.json" ]]; then
  log "measuring the first ten candidates (S2, glycan_free) to project the full run"
  M="${REPO_ROOT}/logs/scratch/measure_scores.tsv"; mkdir -p "$(dirname "${M}")"; rm -f "${M}"
  register_scratch "${REPO_ROOT}/logs/scratch"
  python "${REPO_ROOT}/scripts/06_dock.py" --limit 10 --states S2_closed_holo --arms glycan_free --jobs 1 --out "${M}"
  python - "${M}" "${REPO_ROOT}" "${DISK_GB}" "${JOBS_USE}" <<'PY'
import sys, json, csv, pathlib
m, root, disk_gb, jobs = sys.argv[1], pathlib.Path(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
rows = list(csv.DictReader(open(m), delimiter="\t"))
n_c = sum(1 for _ in open(root / "data/candidates/candidates.tsv")) - 1
n_s = len(json.loads((root / "config/boxes.json").read_text()))
mean_s = sum(float(r["wall_s"]) for r in rows) / len(rows)
pose_kb = sum(f.stat().st_size for f in (root / "data/poses").glob("*.pdbqt")) / 1024 / max(1, len(rows))
total = n_c * n_s * 2
proj = dict(measured_dockings=len(rows), mean_wall_s=round(mean_s, 1), candidates=n_c, states=n_s, dockings_total=total,
            projected_wall_h=round(total * mean_s / jobs / 3600, 2), pose_kb_each=round(pose_kb, 1),
            projected_disk_gb=round(total * pose_kb / 1024 / 1024, 3), jobs=jobs)
(root / "logs/06_projection.json").write_text(json.dumps(proj, indent=1))
print(json.dumps(proj))
if proj["projected_disk_gb"] > disk_gb:
    sys.exit(f"refusing to start: projected {proj['projected_disk_gb']} GB exceeds budget {disk_gb} GB by {proj['projected_disk_gb']-disk_gb:.3f} GB")
PY
fi

ARGS=(--jobs "${JOBS_USE}")
[[ -n "${STATES}" ]] && ARGS+=(--states "${STATES}")
run_timed 06_score_small_molecules python "${REPO_ROOT}/scripts/06_dock.py" "${ARGS[@]}"
stamp_done 06_score_small_molecules
