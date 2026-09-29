#!/usr/bin/env python3
"""Stage 06 (driver). Dock every candidate into every state and glycan arm.

Machinery. AutoDock Vina 1.2.7, default 'vina' scoring function, one box per
state from config/boxes.json, exhaustiveness 8, one Vina seed for every run.
Preparation conventions (RDKit ETKDG embedding, MMFF relaxation, Meeko
charges, receptor from Open Babel) follow stage 1, vina_gnina_pose_benchmark_pipeline
at commit 313851a, which this repository does not vendor: the stage 1 code is
built around benchmark complexes and has no entry point for a bare SMILES list.
That is a deviation from 'import, do not rewrite' and the README says so.

The number this stage writes is a Vina score, the output of an empirical
function fitted to a training set of protein-ligand complexes, in the units of
that function (labelled kcal/mol, which it is not: it is not a free energy and
was never a binding energy). It is a ranking key and nothing more.

Controls. The script never reads labels.tsv. Temsavir and BMS-378806 sit in
the candidate list under anonymous ids and go through the same functions as
the null molecules. This is the point of the anchor argument.

Exhaustiveness 8 is a compromise for 16 cores and a 12-hour window. Stage 1
measured the cost of higher exhaustiveness for its own benchmark, not for an
induced pocket, so no claim is made that 8 converges here.
"""
import argparse, json, os, subprocess, sys, time, csv, threading
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import psutil
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
VINA = os.environ.get("VINA") or ("C:/vina/vina.exe.exe" if os.name == "nt" else "vina")
BOX = 22.0
EXH = 8
VINA_SEED = 42
EMBED_SEED = 42
FIELDS = ["id", "state", "arm", "score", "heavy_atoms", "wall_s", "peak_rss_mb", "exhaustiveness", "seed"]
lock = threading.Lock()


def prep_ligand(cid, smi, ld):
    f = ld / f"{cid}.pdbqt"
    if f.exists():
        return f, None
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None, "smiles_unparsable"
    m = Chem.AddHs(m)
    if AllChem.EmbedMolecule(m, randomSeed=EMBED_SEED) != 0:
        if AllChem.EmbedMolecule(m, randomSeed=EMBED_SEED, useRandomCoords=True) != 0:
            return None, "embedding_failed"
    try:
        AllChem.MMFFOptimizeMolecule(m, maxIters=500)
    except Exception:
        pass
    try:
        setups = MoleculePreparation().prepare(m)
        s, ok, err = PDBQTWriterLegacy.write_string(setups[0])
        if not ok:
            return None, "meeko_failed"
    except Exception as e:
        return None, "meeko_" + type(e).__name__
    f.write_text(s)
    return f, None


def run_vina(lig, rec, centre, out):
    cmd = [VINA, "--receptor", str(rec), "--ligand", str(lig), "--center_x", str(centre[0]), "--center_y", str(centre[1]),
           "--center_z", str(centre[2]), "--size_x", str(BOX), "--size_y", str(BOX), "--size_z", str(BOX),
           "--exhaustiveness", str(EXH), "--seed", str(VINA_SEED), "--cpu", "1", "--num_modes", "1", "--out", str(out)]
    t0 = time.time()
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    ps = psutil.Process(p.pid); peak = 0
    while p.poll() is None:
        try:
            peak = max(peak, ps.memory_info().rss)
        except psutil.Error:
            pass
        time.sleep(0.2)
    return time.time() - t0, peak / 1e6, p.returncode


def score_of(pose):
    for l in Path(pose).read_text().splitlines():
        if l.startswith("REMARK VINA RESULT"):
            return float(l.split()[3])
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--jobs", type=int, default=1,
                    help="concurrent Vina processes. Default 1 (serial) for reproducible timing, not for memory.")
    ap.add_argument("--limit", type=int, default=0, help="first N candidates only (the ten-candidate measurement)")
    ap.add_argument("--states", default="", help="comma list of state ids")
    ap.add_argument("--arms", default="glycan_free,glycosylated")
    ap.add_argument("--out", default=str(ROOT / "results" / "dock_scores.tsv"))
    a = ap.parse_args()

    boxes = json.loads((ROOT / "config" / "boxes.json").read_text())
    states = [s for s in a.states.split(",") if s] or list(boxes)
    cands = list(csv.DictReader(open(DATA / "candidates" / "candidates.tsv"), delimiter="\t"))
    if a.limit:
        cands = cands[:a.limit]
    ld = DATA / "ligands"; ld.mkdir(exist_ok=True)
    pd_ = DATA / "poses"; pd_.mkdir(exist_ok=True)
    out = Path(a.out)
    done = set()
    if out.exists():
        done = {(r["id"], r["state"], r["arm"]) for r in csv.DictReader(open(out), delimiter="\t")}
    else:
        out.write_text("\t".join(FIELDS) + "\n")
    fail = ROOT / "results" / "dock_failures.tsv"
    if not fail.exists():
        fail.write_text("id\tstate\tarm\treason\n")

    ligs = {}
    for c in cands:
        f, err = prep_ligand(c["id"], c["smiles"], ld)
        if f is None:
            with fail.open("a") as fh:
                fh.write(f"{c['id']}\tall\tall\t{err}\n")
        else:
            ligs[c["id"]] = (f, Chem.MolFromSmiles(c["smiles"]).GetNumHeavyAtoms())
    print(f"[dock] {len(ligs)} of {len(cands)} ligands prepared", file=sys.stderr)

    tasks = [(cid, s, arm) for s in states for arm in a.arms.split(",") for cid in ligs if (cid, s, arm) not in done]
    print(f"[dock] {len(tasks)} dockings to run ({len(done)} already done), jobs={a.jobs}", file=sys.stderr)

    def work(t):
        cid, s, arm = t
        rec = DATA / "receptors" / f"{s}__{arm}.pdbqt"
        pose = pd_ / f"{cid}__{s}__{arm}.pdbqt"
        wall, rss, rc = run_vina(ligs[cid][0], rec, boxes[s]["centre"], pose)
        sc = score_of(pose) if rc == 0 and pose.exists() else None
        with lock:
            if sc is None:
                with fail.open("a") as fh:
                    fh.write(f"{cid}\t{s}\t{arm}\tvina_rc_{rc}\n")
            else:
                with out.open("a") as fh:
                    fh.write("\t".join(map(str, [cid, s, arm, sc, ligs[cid][1], round(wall, 2), round(rss, 1), EXH, VINA_SEED])) + "\n")

    t0 = time.time()
    with ThreadPoolExecutor(a.jobs) as ex:
        list(ex.map(work, tasks))
    print(f"[dock] finished in {time.time()-t0:.0f} s", file=sys.stderr)


if __name__ == "__main__":
    main()
